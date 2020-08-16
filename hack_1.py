import os
import time
from argparse import ArgumentParser
import numpy as np


TRACING = False
KMAD = 1.4826
NSAMPS_TOTAL = 55000
NCHANS = 2048
NSAMPS = 1024


def os_error(path, function, ose):
    msg = "{} failed [{}], errno={}, strerror:{}".format(function, path, ose.errno, ose.strerror)
    raise OSError(msg)


def correct_drops(block, mask_array, mask_offset, thresh):
    '''
    Correct Drops in block

    Parameters
    ----------
    block : float32 array of shape (nsamps, NCHANS)
    thresh :

    Returns
    -------
    mask : All zeros except the to_fill position which = 1

    '''
    ts = block.sum(axis=1)
    med = np.median(ts)
    mad = np.median(np.abs(ts - med))
    ts -= med
    ts /= (KMAD * mad)

    bp_mad = np.median(np.abs(block - np.median(block,axis=0)), axis=0)
    #UNUSED bp_m = block.mean(axis=0)
    bp_m_med = np.median(block, axis=0)

    to_fill = np.where((ts < - thresh) | (ts > thresh))[0]
    block[to_fill] = np.random.normal(bp_m_med,
                                      bp_mad * KMAD,
                                      size=(len(to_fill), len(bp_mad)))
    #UNUSED ndrops = len(to_fill)
    mask_array[mask_offset + to_fill] = 1


#======================= main ========================
def main():
    """
    Usage:
        rfi_mitigation  [<FULL_PATH_TO_BIN_FILE>]  [options]

    Options:
        -h, --help            Show this help message and exit
        -o OUT_DIR, --out_dir=OUT_DIR
                              Location for output files. Default: current dir.
    """
    p = ArgumentParser(description='RFI Mitigation (Wael)')
    p.add_argument('-i', '--in_path', dest='in_path', type=str, default='./block.bin',
                   help='Full/relative path of input file')
    p.add_argument('-o', '--out_dir', dest='out_dir', type=str, default='./',
                   help='Location for output files. Default: local dir. ')
    p.add_argument('-T', '--threshold', dest='threshold', type=int, default=5,
                   help='Threshold (Wael?)')
    args = p.parse_args()
    in_block_path = args.in_path
    out_corr_path = args.out_dir + 'corrected_' + os.path.basename(in_block_path)
    out_mask_path = args.out_dir + 'mask.bin'
    thresh = args.threshold

    print("main: BEGIN .....")
    time_start = time.time()
    mask_array = np.zeros(NSAMPS_TOTAL, dtype=np.int8)
    block_shape = (NSAMPS_TOTAL, NCHANS)

    try:
        nbytes = os.path.getsize(in_block_path)
    except OSError as ose:
        os_error(in_block_path, 'os.path.getsize', ose)

    # Load input block data.
    block = np.fromfile(in_block_path, dtype='float32').reshape(block_shape)

    # Process all of the input blocks.
    offset = 0
    ibulk = 0
    nsamps_to_correct = 0
    while True:
        nsamps_to_correct = min(NSAMPS, 
                                NSAMPS_TOTAL - (ibulk * nsamps_to_correct))
        correct_drops(block[offset:offset + nsamps_to_correct],
                      mask_array,
                      offset,
                      thresh)
        if TRACING:
            print("main:TRACING: block {} completed".format(ibulk))
        if (ibulk * NSAMPS) + nsamps_to_correct == NSAMPS_TOTAL:
            break
        ibulk += 1
        offset += nsamps_to_correct

    # Write out corrected blocks.
    try:
        with open(out_corr_path, "wb") as fd_corr:
            fd_corr.write(block)
    except OSError as ose:
        os_error(out_corr_path, 'write(corr)', ose)
    if TRACING:
        print("main:TRACING: corrected block array written")

    # Write out mask array.
    try:
        with open(out_mask_path, "wb") as fd_mask:
            fd_mask.write(mask_array)
    except OSError as ose:
        os_error(out_mask_path, 'write(mask)', ose)
    if TRACING:
        print("main:TRACING: mask array written")

    time_stop = time.time()
    et = time_stop - time_start
    print("main: END, Total blocks =", ibulk)
    print("main: END, Total elapsed time = {:.2f} seconds".format(et))
    print("main: END, Blocks/second = {:.2f}".format(ibulk / et))
    print("main: END, MB/second = {:.2f}".format((nbytes / et) / 1e6))


if __name__ == "__main__":
    main()