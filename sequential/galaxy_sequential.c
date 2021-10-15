#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "../common/galaxy_utils.h"

float *real_rasc, *real_decl, *rand_rasc, *rand_decl;
long int MemoryAllocatedCPU = 0L;
int numberOfGalaxies = 100000;

int main(int argc, char* argv[]) 
    {

    struct timeval _ttime;
    struct timezone _tzone;
    float  pif = acosf(-1.0f);


    gettimeofday(&_ttime, &_tzone);
    double time_start = (double)_ttime.tv_sec + (double)_ttime.tv_usec/1000000.;

    // store right ascension and declination for real galaxies here
    // Note: indices run from 0 to 99999 = 100000-1: realrasc[0] -> realrasc[99999] 
    // realrasc[100000] is out of bounds for allocated memory!
    real_rasc        = (float *)calloc(numberOfGalaxies, sizeof(float));
    real_decl        = (float *)calloc(numberOfGalaxies, sizeof(float));

    // store right ascension and declination for synthetic random galaxies here
    rand_rasc        = (float *)calloc(numberOfGalaxies, sizeof(float));
    rand_decl        = (float *)calloc(numberOfGalaxies, sizeof(float));

    MemoryAllocatedCPU += 10L*numberOfGalaxies*sizeof(float);

    // read input data from files given on the command line
    if ( parseargs_readinput(argc, argv, numberOfGalaxies, real_rasc, real_decl, rand_rasc, rand_decl) != 0 ) {printf("   Program stopped.\n");return(0);}
    printf("   Input data read, now calculating histograms\n");

    long int histogram_DD[360] = {0L};
    long int histogram_DR[360] = {0L};
    long int histogram_RR[360] = {0L};
    MemoryAllocatedCPU += 3L*360L*sizeof(long int);
    
    float quantum = .25;

      // compute DR
    for(int i=0; i<numberOfGalaxies; ++i)
        for(int j=0; j<numberOfGalaxies; ++j){
            float d = angle_between(real_rasc[i], real_decl[i], rand_rasc[j], rand_decl[j]) * 180 / pif;
            ++histogram_DR[(int)floor(d/quantum)];

        }


    // compute DD
    histogram_DD[0] = numberOfGalaxies;
    for(int i=0; i<numberOfGalaxies; ++i){
        float d = angle_between(real_rasc[i], real_decl[i], real_rasc[i], real_decl[i]);
        for(int j=i+1; j<numberOfGalaxies; ++j){
            float d = angle_between(real_rasc[i], real_decl[i], real_rasc[j], real_decl[j]) * 180 / pif;
            histogram_DD[(int)floor(d/quantum)]+=2;

        }
    }
    // compute RR
    histogram_RR[0] = numberOfGalaxies;
    for(int i=0; i<numberOfGalaxies; ++i){
        for(int j=i+1; j<numberOfGalaxies; ++j){
            float d = angle_between(rand_rasc[i], rand_decl[i], rand_rasc[j], rand_decl[j]) * 180 / pif;    
            histogram_RR[(int)floor(d/quantum)]+=2;
        }
    }


    // check point: the sum of all historgram entries should be 10 000 000 000
    long int histsum = 0L;
    int      correct_value=1;
    for ( int i = 0; i < 360; ++i ) histsum += histogram_DD[i];
    printf("   Histogram DD : sum = %ld\n",histsum);
    if ( histsum != (numberOfGalaxies*numberOfGalaxies)) correct_value = 0;

    histsum = 0L;
    for ( int i = 0; i < 360; ++i ) histsum += histogram_DR[i];
    printf("   Histogram DR : sum = %ld\n",histsum);
    if ( histsum != numberOfGalaxies*numberOfGalaxies ) correct_value = 0;

    histsum = 0L;
    for ( int i = 0; i < 360; ++i ) histsum += histogram_RR[i];
    printf("   Histogram RR : sum = %ld\n",histsum);
    if ( histsum != (numberOfGalaxies*numberOfGalaxies)) correct_value = 0;

    if ( correct_value != 1 ) 
       {printf("   Histogram sums should be %d. Ending program prematurely\n", numberOfGalaxies*numberOfGalaxies);return(0);}


    printf("   Omega values for the histograms:\n");
    float omega[360];
    for ( int i = 0; i < 10; ++i ) 
        if ( histogram_RR[i] != 0L )
           {
           omega[i] = (histogram_DD[i] - 2L*histogram_DR[i] + histogram_RR[i])/((float)(histogram_RR[i]));
           if ( i < 10 ) printf("      angle %.2f deg. -> %.2f deg. : %.3f\n", i*0.25, (i+1)*0.25, omega[i]);
           }

    FILE *out_file = fopen(argv[3],"w");
    if ( out_file == NULL ) printf("   ERROR: Cannot open output file %s\n",argv[3]);
    else
       {
       for ( int i = 0; i < 360; ++i ) 
           if ( histogram_RR[i] != 0L )
              fprintf(out_file,"%.2f  : %.3f\n", i*0.25, omega[i] ); 
       fclose(out_file);
       printf("   Omega values written to file %s\n",argv[3]);
       }
       

    free(real_rasc); free(real_decl);
    free(rand_rasc); free(rand_decl);

    printf("   Total memory allocated = %.1lf MB\n",MemoryAllocatedCPU/1000000.0);
    gettimeofday(&_ttime, &_tzone);
    double time_end = (double)_ttime.tv_sec + (double)_ttime.tv_usec/1000000.;

    printf("   Wall clock run time    = %.1lf secs\n",time_end - time_start);

    return(0);
}
