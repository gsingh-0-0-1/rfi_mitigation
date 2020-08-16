import numpy as np
import matplotlib.pyplot as plt
import time


def correct_drops(block, thresh):
    # Block shape is: (nsamps, nchans)
    ts = block.sum(axis=1)
    med = np.median(ts)
    mad = np.median(np.abs(ts-med))
    ts -= med
    ts /= (1.4826*mad)

    bp_mad = np.median(np.abs(block - np.median(block,axis=0)), axis=0)
    bp_m = block.mean(axis=0)
    bp_m_med = np.median(block, axis=0)

    to_fill = np.where((ts < - thresh) | (ts > thresh))[0]
    block[to_fill] = np.random.normal(bp_m_med, bp_mad*1.4826,
            size=(len(to_fill), len(bp_mad)))
    ndrops = len(to_fill)

    mask = np.zeros(block.shape[0], dtype=np.int8)
    mask[to_fill] = 1

    return mask, bp_m_med, bp_mad



block_shape = (55000, 2048)
thresh = 5

# load data
block = np.fromfile('./block.bin', dtype='float32').reshape(block_shape)
block_corrected_c = np.fromfile('corrected_block.bin', 
        dtype='float32').reshape(block_shape)

# plot timeseries pre-clean
plt.figure(1)
plt.title("Time series")
plt.plot(np.arange(len(block.sum(axis=1)))*728e-6,
        block.sum(axis=1), label="pre-clean")

i = 0
ibulk = 0
nsamps_to_correct = 0
mask = []
nsamps = 512
nsamps_total = block_shape[0]
while True:
    nsamps_to_correct = min(nsamps, nsamps_total - ibulk*nsamps_to_correct)
    m, bp_m_med, bp_mad = correct_drops(block[i:i+nsamps_to_correct], thresh)

    if (ibulk*nsamps + nsamps_to_correct == nsamps_total):
        break
    ibulk += 1
    i += nsamps_to_correct
    print (i)

"""
for i in np.linspace(0, block_shape[0], 100).astype(int):
    if i == 0:
        continue
    m, bp_m_med, bp_mad = correct_drops(block[ii:i], thresh)
    mask.append(m)
    ii = i
    print (ii)
print (time.time() - t)
"""


# plot data post-clean + the C version
plt.plot(np.arange(len(block.sum(axis=1)))*728e-6, 
        block.sum(axis=1), 
        label="python clean")
plt.plot(np.arange(len(block.sum(axis=1)))*728e-6,
        block_corrected_c.sum(axis=1), label="C clean")
plt.xlabel("Time (seconds)")

# Plot waterfall
plt.legend()
plt.figure()
plt.imshow(10*np.log(block_corrected_c), interpolation='nearest', aspect='auto')
plt.xlabel("Frequency chanel")
plt.ylabel("Time in samples")
plt.title("Clean using C")

# Plot waterfall
plt.figure()
plt.imshow(10*np.log(block), interpolation='nearest', aspect='auto')
plt.title("Clean using python")
plt.xlabel("Frequency chanel")
plt.ylabel("Time in samples")

plt.show()