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
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

int compare(const void * a, const void * b){
    return (int) ( *(float *)a - *(float*)b );;
}

float quick_select_median(float arr[], uint16_t n){
    uint16_t low, high ;
    uint16_t median;
    uint16_t middle, ll, hh;
    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
    if (high <= low) /* One element only */
    return arr[median] ;
    if (high == low + 1) { /* Two elements only */
    if (arr[low] > arr[high])
    ELEM_SWAP(arr[low], arr[high]) ;
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
    return arr[median] ;
}

/*float getMedian(float *arr, int len){
    float median;
    if (len % 2 == 0){
        int center = len / 2;
        median = 0.5 * (arr[center] + arr[center - 1]);
    }
    else{
        int center = 1 + (len / 2);
        median = 1.0 * arr[center];
    }
    
    return median;
}*/

double norm(double mean, double std_dev)
{
  double   u, r, theta;           // Variables for Box-Muller method
  double   x;                     // Normal(0, 1) rv
  double   norm_rv;               // The adjusted normal rv

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

  // Return the normally distributed RV value
  return(norm_rv);
}

void correct_drops(float* block, int64_t nchans, int64_t nsamps, double thresh, char* mask){
    //block shape is (nchans, nsamps)
    
    //ts is the sum of block across axis 0
    float ts[nsamps];
    //ensure that we set everything to start at 0, since we're iteratively adding to the value
    //in the next loop, rather than just setting it
    memset(ts, 0, nsamps * sizeof(*ts));
    
    //axiswise sum of block array to get the timeseries
    for (int i = 0; i < nsamps; i++){
        for (int j = 0; j < nchans; j++){
            ts[i] = ts[i] + *(block + (i*nchans) + j);
        }
    }

    fprintf(stderr, "Finished calculating time-series.\n");
    
    //sort ts to get median - create a separate array and copy ts into it, so we can preserve
    //the original ts for later
    float ts_sorted[nsamps];
    memcpy(ts_sorted, ts, nsamps * sizeof(*ts));
    qsort(ts_sorted, nsamps, sizeof(*ts), compare);
    
    //now get the median
    float median;
    median = quick_select_median(ts_sorted, nsamps);
    fprintf(stderr, "Time-series median: %0.1f \n", median);
    
    //calculate mad
    //start by subtracting median and getting absolute value
    float mad;
    float ts_processed[nsamps];
    for (int i = 0; i < nsamps; i++){
        float temp = ts[i] - median;
        if (temp < 0){ temp = temp * -1; }
        ts_processed[i] = temp;
    }
    
    //now we have the processed array
    //now we actually calculate the mad
    qsort(ts_processed, nsamps, sizeof(*ts_processed), compare);
    mad = quick_select_median(ts_processed, nsamps);
    fprintf(stderr, "Time-series MAD: %0.1f \n", mad);
    
    //modify the ts using the mad and median
    for (int i = 0; i < nsamps; i++){
        ts[i] = ts[i] - median;
        ts[i] = ts[i] / (1.4826 * mad);
    }
    
    //output the ts
    printf("Finished modifying time-series.\n");
    
    
    //start bandpass-related calculations
    float bp[nchans];
    //initialize to 0
    for (int i = 0; i < nchans; i++){
        bp[i] = 0;
    }
    
    //axiswise median of block array to get the bandpass
    
    //create a list to keep the medians in. this will be separate from the one in which
    //we create the mad values
    float bp_med[nchans];
    //mad value list
    float bp_mad[nchans];
    
    //create the channel variable to calculate stats across each channel
    float channel[nsamps];
    
    for (int i = 0; i < nchans; i++){
        for (int j = 0; j < nsamps; j++){
            channel[j] = *(block + i + (j * nchans));
        }
        //get the median of the channel
        qsort(channel, nsamps, sizeof(*channel), compare);
        float channel_median = quick_select_median(channel, nsamps);//getMedian(channel, nsamps);
        bp_med[i] = channel_median;
        
        //now we calculate the mad of the channel
        for (int k = 0; k < nsamps; k++){
            channel[k] = channel[k] - channel_median;
            if (channel[k] < 0) { channel[k] = channel[k] * -1; }
        }
        qsort(channel, nsamps, sizeof(*channel), compare);
        float channel_mad = quick_select_median(channel, nsamps);//getMedian(channel, nsamps);
        bp_mad[i] = channel_mad;
    }

    printf("Finished calculating bandpass med/mad.\n");
    
    //now we replace the points in the ts that are outside (-thresh, thresh)
    //iterate through the ts to find any indices where it is outside that interval
    for (int i = 0; i < nsamps; i++){
        //check for points
        if (ts[i] < -thresh){// || ts[i] > thresh){
            //set mask to 1 here
            *(mask + i) = 1;
            //if the ts exceeds the threshold, start looping through the block at that index
            //set all the sample values at those channels
            //from a normal distribution with a mean at the bp_med of that sample
            //and a standard deviation of the bp_mad of that sample
            for (int ch = 0; ch < nchans; ch++){
                *(block + (i * nchans) + ch) = norm(bp_med[ch], 1.4826*bp_mad[ch]);
            }
        }
    }
    
}

int main(){
    srand(time(NULL));

    int nchans = 2048;
    int nsamps_total = 55000;
    int nsamps = 1024;
    
    double thresh = 5;

    char in_file_name[64] = "block.bin";
    char out_file_name[64] = "corrected_block.bin";
    char mask_file_name[64] = "mask.bin";

    int chan_display_max = 5;
    int samp_display_max = 200;

    FILE *ptr_in;
    FILE *ptr_out;
    FILE *mask_ptr;
    
    float *block;
    block = calloc(nchans*nsamps, sizeof *block);

    char *mask;
    mask = calloc(nsamps, sizeof *mask);

    ptr_in = fopen(in_file_name, "rb");
    if (ptr_in == NULL){
      fprintf(stderr, "Error: Input file not open...\n");
      return EXIT_FAILURE;
    }
    
    ptr_out = fopen(out_file_name, "wb");
    if (ptr_out == NULL){
      fprintf(stderr, "Error: output file not open...\n");
      return EXIT_FAILURE;
    }
    
    mask_ptr = fopen(mask_file_name, "wb");
    if (mask_ptr == NULL){
      fprintf(stderr, "Error: mask file not open...\n");
      return EXIT_FAILURE;
    }

    int nsamps_to_read = 0;
    int ibulk = 0;
    size_t nread = 0;
    while (1){
        fprintf(stderr, "ibulk number: %i\n", ibulk);
        nsamps_to_read = MIN(nsamps, nsamps_total - ibulk*nsamps);
        nread = fread(block, sizeof(*block), nchans*nsamps_to_read, ptr_in);

        printf("Block loaded, block[0]: %.1f.\n", block[0]);

        int start_time = (unsigned)time(NULL);
        printf("Processing start time: %d\n", start_time);
        
        memset(mask, 0, nsamps * sizeof(*mask));
        
        correct_drops(block, nchans, nsamps_to_read, thresh, mask);
        
        int end_time = (unsigned)time(NULL);
        printf("Processing end time: %d\n", end_time);
        printf("\n");

        printf("Total processing time: %d\n", end_time - start_time);

        printf("Finished correcting block, created mask.\n");

        fwrite(block, sizeof(*block), nchans*nsamps_to_read, ptr_out);
        fwrite(mask, sizeof(*mask), nsamps_to_read, mask_ptr);

        printf("Finished writing to corrected_block.bin\n");
        if (ibulk*nsamps + nsamps_to_read == nsamps_total){
            break;
        }
        ibulk += 1;
    }
    
    fclose(ptr_in);
    fclose(ptr_out);
    fclose(mask_ptr);
    free(block);
    free(mask);
}
