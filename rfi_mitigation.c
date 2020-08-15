/*
 * rfi_mitigation.c
 *
*/
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <time.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define ELEM_SWAP(a,b) { register float t=(a);(a)=(b);(b)=t; }

//#define DETAILED_TRACING

static int float_compare(const void * a, const void * b) {
    return (int) ( *(float *)a - *(float*)b );;
}


static double timer_drops = 0.0;
static double timer_median = 0.0;
static double timer_norm = 0.0;


/*
 * Algorithm:
 * http://www.stat.cmu.edu/~ryantibs/median/
 *
 */
float quick_select_median(float arr[], uint32_t n) {
    uint32_t low, high ;
    uint32_t median;
    uint32_t middle, ll, hh;
    clock_t clock_start = clock();
    low = 0 ; high = n-1 ; median = (low + high) / 2;

    for (;;) {
        if (high <= low) /* One element only */ {
            timer_median += (double) (clock() - clock_start) / CLOCKS_PER_SEC;
            return arr[median] ;
        }
        if (high == low + 1) /* Two elements only */ {
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]) ;
            timer_median += (double) (clock() - clock_start) / CLOCKS_PER_SEC;
            return arr[median] ;
        }
        /* Find median of low, middle and high items; swap into position low */
        middle = (low + high) / 2;
        if (arr[middle] > arr[high])
        ELEM_SWAP(arr[middle], arr[high]) ;
        if (arr[low] > arr[high])
        ELEM_SWAP(arr[low], arr[high]) ;
        if (arr[middle] > arr[low])
        ELEM_SWAP(arr[middle], arr[low]) ;
        /* Swap low item (now in position middle) into position (low+1) */
        ELEM_SWAP(arr[middle], arr[low+1]) ;
        /* Nibble from each end towards middle, swapping items when stuck */
        ll = low + 1;
        hh = high;
        for (;;) {
            do ll++; while (arr[low] > arr[ll]) ;
            do hh--; while (arr[hh] > arr[low]) ;
            if (hh < ll)
                break;
            ELEM_SWAP(arr[ll], arr[hh]) ;
        }
        /* Swap middle item (in position low) back into correct position */
        ELEM_SWAP(arr[low], arr[hh]) ;
        /* Re-set active partition */
        if (hh <= median)
            low = ll;
        if (hh >= median)
            high = hh - 1;
    }
//    timer_median += (double) (clock() - clock_start) / CLOCKS_PER_SEC;
//    return arr[median] ;
}


/*
 * Algorithm:
 * https://www.csee.usf.edu/~kchriste/tools/gennorm.c
 *
 */
double norm(double mean, double std_dev) {
    double   u, r, theta;           // Variables for Box-Muller method
    double   x;                     // Normal(0, 1) rv
    double   norm_rv;               // The adjusted normal rv
    clock_t clock_start = clock();

    // Generate u
    u = 0.0;
    while (u == 0.0)
        u = (rand() % 1000) / 1000.0;

    // Compute r
    r = sqrt(-2.0 * log(u));

    // Generate theta
    theta = 0.0;
    while (theta == 0.0)
        theta = 2.0 * 3.1415 * (rand() % 1000) / 1000.0;

    // Generate x value
    x = r * cos(theta);

    // Adjust x value for specified mean and variance
    norm_rv = (x * std_dev) + mean;

    timer_norm += (double) (clock() - clock_start) / CLOCKS_PER_SEC;

    // Return the normally distributed RV value
    return(norm_rv);
}


static inline void correct_drops(float* block, uint64_t nchans, uint64_t nsamps, double thresh, char* mask) {
    //block shape is (nchans, nsamps)
    
    size_t ichan,isamp;

    //ts is the sum of block across axis 0
    float ts[nsamps];
    //ensure that we set everything to start at 0, since we're iteratively adding to the value
    //in the next loop, rather than just setting it
    memset(ts, 0, nsamps * sizeof(*ts));
    
    //axiswise sum of block array to get the timeseries
    for (isamp = 0; isamp < nsamps; isamp++){
        for (ichan = 0; ichan < nchans; ichan++)
            ts[isamp] = ts[isamp] + *(block + (isamp*nchans) + ichan);
    }

    #ifdef DETAILED_TRACING
    fprintf(stderr, "correct_drops: Finished calculating time-series.\n");
    #endif
        
    //sort ts to get median - create a separate array and copy ts into it, so we can preserve
    //the original ts for later
    float ts_sorted[nsamps];
    memcpy(ts_sorted, ts, nsamps * sizeof(*ts));
    qsort(ts_sorted, nsamps, sizeof(*ts), float_compare);
    
    //now get the median
    float median = quick_select_median(ts_sorted, nsamps);
    
    #ifdef DETAILED_TRACING
    fprintf(stderr, "correct_drops: Time-series median: %0.1f \n", median);
    #endif
    
    //calculate mad
    //start by subtracting median and getting absolute value
    float ts_processed[nsamps];
    for (isamp = 0; isamp < nsamps; isamp++)
        ts_processed[isamp] = fabs(ts[isamp] - median);
    
    //now we have the processed array
    //now we actually calculate the mad
    qsort(ts_processed, nsamps, sizeof(*ts_processed), float_compare);
    float mad = quick_select_median(ts_processed, nsamps);

    #ifdef DETAILED_TRACING
    fprintf(stderr, "correct_drops: Time-series MAD: %0.1f \n", mad);
    #endif
        
    //modify the ts using the mad and median
    float factor = 1.0 / (1.4826 * mad);
    for (isamp = 0; isamp < nsamps; isamp++)
        ts[isamp] = (ts[isamp] - median) * factor;
    
    //output the ts
    #ifdef DETAILED_TRACING
    printf("correct_drops: Finished modifying time-series.\n");
    #endif    
    
    //axiswise median of block array to get the bandpass
    
    //create a list to keep the medians in. this will be separate from the one in which
    //we create the mad values
    float bp_med[nchans];
    //mad value list
    float bp_mad[nchans];
    
    //create the channel variable to calculate stats across each channel
    float channel[nsamps];

    
    for (ichan = 0; ichan < nchans; ichan++) {
        for (isamp = 0; isamp < nsamps; isamp++)
            channel[isamp] = *(block + ichan + (isamp * nchans));

        //get the median of the channel
        qsort(channel, nsamps, sizeof(*channel), float_compare);
        bp_med[ichan] = quick_select_median(channel, nsamps);//getMedian(channel, nsamps);
        
        //now we calculate the mad of the channel
        for (isamp = 0; isamp < nsamps; isamp++)
            channel[isamp] = fabs(channel[isamp] - bp_med[ichan]);

        qsort(channel, nsamps, sizeof(*channel), float_compare);
        bp_mad[ichan] = quick_select_median(channel, nsamps);//getMedian(channel, nsamps);
    }

    #ifdef DETAILED_TRACING
    printf("correct_drops: Finished calculating bandpass med/mad.\n");
    #endif    
    
    //now we replace the points in the ts that are outside (-thresh, thresh)
    //iterate through the ts to find any indices where it is outside that interval
    for (isamp = 0; isamp < nsamps; isamp++) {
        //check for points
        if (ts[isamp] < -thresh){// || ts[i] > thresh) {
            //set mask to 1 here
            *(mask + isamp) = 1;
            //if the ts exceeds the threshold, start looping through the block at that index
            //set all the sample values at those channels
            //from a normal distribution with a mean at the bp_med of that sample
            //and a standard deviation of the bp_mad of that sample
            for (ichan = 0; ichan < nchans; ichan++)
                *(block + (isamp * nchans) + ichan) = 
                  norm(bp_med[ichan], 1.4826*bp_mad[ichan]);
        }
    }
}

int main() {
    printf("main: BEGIN .....\n");
    int main_start_time = (unsigned) time(NULL);
    srand(time(NULL));
    
    const int nchans = 2048;
    const int nsamps_total = 55000;
    const int nsamps = 1024;
    
    double thresh = 5;

    char in_file_name[64] = "block.bin";
    char out_file_name[64] = "corrected_block.bin";
    char mask_file_name[64] = "mask.bin";

    // ### NOT USED ### int chan_display_max = 5;
    // ### NOT USED ### int samp_display_max = 200;

    FILE *ptr_in;
    FILE *ptr_out;
    FILE *ptr_mask;
    
    float *block;
    block = calloc(nchans*nsamps, sizeof *block);

    char *mask;
    mask = calloc(nsamps, sizeof *mask);

    ptr_in = fopen(in_file_name, "rb");
    if (ptr_in == NULL) {
        fprintf(stderr, "Error: Input file not open...\n");
        return EXIT_FAILURE;
    }
    
    ptr_out = fopen(out_file_name, "wb");
    if (ptr_out == NULL) {
        fprintf(stderr, "Error: output file not open...\n");
        return EXIT_FAILURE;
    }
    
    ptr_mask = fopen(mask_file_name, "wb");
    if (ptr_mask == NULL) {
        fprintf(stderr, "Error: mask file not open...\n");
        return EXIT_FAILURE;
    }

    int nsamps_to_read = 0;
    int ibulk = 0;
    size_t nread = 0; // ### fread SET BUT NOT USED ###
    long nbytes = 0L;
    clock_t clock_start;

    #ifdef DETAILED_TRACING
    int block_start_time, block_end_time;
    #endif

    while (1) {
        #ifdef DETAILED_TRACING
        fprintf(stderr, "block: ibulk number: %i\n", ibulk);
        #endif
        
        nsamps_to_read = MIN(nsamps, nsamps_total - ibulk*nsamps);
        nread = fread(block, sizeof(*block), nchans*nsamps_to_read, ptr_in);

        #ifdef DETAILED_TRACING
        printf("block: Block loaded, block[0]: %.1f.\n", block[0]);
        block_start_time = (unsigned) time(NULL);
        #endif
                
        memset(mask, 0, nsamps * sizeof(*mask));
        
        clock_start = clock();
        correct_drops(block, nchans, nsamps_to_read, thresh, mask);
        timer_drops += (double) (clock() - clock_start) / CLOCKS_PER_SEC;
        
        #ifdef DETAILED_TRACING
        block_end_time = (unsigned) time(NULL);
        printf("block: Total processing time: %d seconds\n", block_end_time - block_start_time);
        printf("block: Finished correcting block, created mask.\n");
        #endif
                
        fwrite(block, sizeof(*block), nchans*nsamps_to_read, ptr_out);
        fwrite(mask, sizeof(*mask), nsamps_to_read, ptr_mask);

        #ifdef DETAILED_TRACING
        printf("block: Finished writing to corrected_block.bin\n");
        #endif

        if (ibulk*nsamps + nsamps_to_read == nsamps_total)
            break;

        ++ibulk;
        nbytes += nread;
    }
    
    fclose(ptr_in);
    fclose(ptr_out);
    fclose(ptr_mask);
    free(block);
    free(mask);

    int main_stop_time = (unsigned) time(NULL);
    int et = main_stop_time - main_start_time;
    float float_et = (float) et;
    
    float pct_median = timer_median * 100.0 / float_et;
    float pct_norm = timer_norm * 100.0 / float_et;
    float pct_drops = (timer_drops - timer_median - timer_norm) * 100.0 / float_et;
    
    printf("main: END, Total blocks = %d\n", ibulk);
    printf("main: END, Total elapsed time = %d seconds\n", et);
    printf("main: END, Blocks/second = %.2f\n", (float) ibulk / float_et);
    printf("main: END, MB/second = %.2f\n", ((float) nbytes / float_et) / 1e6);
    printf("main: END, Pct time in quick_select_median() = %.2f %%\n", pct_median);
    printf("main: END, Pct time in norm() = %.2f %%\n", pct_norm);
    printf("main: END, Pct time in correct_drops() = %.2f %%\n", pct_drops);

}
