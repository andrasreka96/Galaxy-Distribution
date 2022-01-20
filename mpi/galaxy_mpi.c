#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <mpi.h>

#include "../common/galaxy_utils.h"


int main(int argc, char* argv[]) 
{

  int np, myid;
  double starttime, endtime;

  MPI_Init(&argc, &argv);                /* Initialize MPI */
  MPI_Comm_size(MPI_COMM_WORLD, &np);    /* Get nr of processes */
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);    /* Get own identifier */

  if (np < 3) {
    if (myid == 0) printf("You have to run this program with at least 3 processes\n");
    MPI_Finalize();
    exit(0);
  }

 
  float *real_rasc, *real_decl, *rand_rasc, *rand_decl, *real_rasc_partial;
  long int numberOfGalaxies = 1000;
  int num_of_angle_intervals = 360;
  float quantum_reciprocal = 4;
  float radian_to_angle = 180 / acosf(-1.0f);

  // store right ascension and declination for real galaxies here
  real_rasc = (float *)calloc(numberOfGalaxies, sizeof(float));
  real_decl = (float *)calloc(numberOfGalaxies, sizeof(float));

  // store right ascension and declination for synthetic random galaxies here
  rand_rasc = (float *)calloc(numberOfGalaxies, sizeof(float));
  rand_decl = (float *)calloc(numberOfGalaxies, sizeof(float));

  if ( parseargs_readinput(argc, argv, numberOfGalaxies, real_rasc, real_decl, rand_rasc, rand_decl) != 0 ) {printf("   Program stopped.\n");return(0);}
    printf("   Input data read, now calculating histograms\n");

  //init global histograms, these will store the result of reduce
  long int histogram_DD_g[360] = {0L};
  long int histogram_DR_g[360] = {0L};
  long int histogram_RR_g[360] = {0L};

 
  //init local histograms, these will store the partial calculations
  long int histogram_DD_l[360] = {0L};
  long int histogram_DR_l[360] = {0L};
  long int histogram_RR_l[360] = {0L};

  printf("--------------%f", real_rasc[1]);

  /* Start measuring time */
  if (myid == 0) {
    starttime = MPI_Wtime();
  }

  int p0_n = (np-1)/3 + 1;
  int p1_n = (np-2)/3 + 1;
  int p2_n = (np-3)/3 + 1;


  if(myid%3 == 0){
    //calculates DR

    int starti = (myid/3) * (numberOfGalaxies/p0_n);
    for(int i=starti; i<starti+numberOfGalaxies/p0_n && i<numberOfGalaxies; ++i)
        for(int j=0; j<numberOfGalaxies; ++j){
          float angle = angle_between(real_rasc[i], real_decl[i], rand_rasc[j], rand_decl[j]) * radian_to_angle;
          ++histogram_DR_l[ (int)(angle*quantum_reciprocal)];
        }
  }
  else if (myid%3 == 1)
  {
    //calculates DD


    int starti = (myid/3) * (numberOfGalaxies/p1_n);
    for(int i=starti; i<starti+numberOfGalaxies/p1_n && i<numberOfGalaxies; ++i)       
      for(int j=i+1; j<numberOfGalaxies; ++j){
        float angle = angle_between(real_rasc[i], real_decl[i], real_rasc[j], real_decl[j]) * radian_to_angle;
        histogram_DD_l[ (int)(angle*quantum_reciprocal)] += 2;
      }
  }

  else{
      
    //calculates RR

    int starti = (myid/3) * (numberOfGalaxies/p2_n);
    for(int i=starti; i<starti+numberOfGalaxies/p2_n && i<numberOfGalaxies; ++i)
        for(int j=i+1; j<numberOfGalaxies; ++j){
          float angle = angle_between(rand_rasc[i], rand_decl[i], rand_rasc[j], rand_decl[j]) * radian_to_angle;
          histogram_RR_l[ (int)(angle*quantum_reciprocal)] += 2;
        }
    
  }

  MPI_Reduce(histogram_DD_l, histogram_DD_g, 360 ,MPI_LONG ,MPI_SUM ,0, MPI_COMM_WORLD);
  MPI_Reduce(histogram_DR_l, histogram_DR_g, 360 ,MPI_LONG ,MPI_SUM ,0, MPI_COMM_WORLD);
  MPI_Reduce(histogram_RR_l, histogram_RR_g, 360 ,MPI_LONG ,MPI_SUM ,0, MPI_COMM_WORLD);


  if(myid==0)

    {

    // galaxie distance of itself will be 0    
    histogram_DD_g[0] += numberOfGalaxies;
    histogram_RR_g[0] += numberOfGalaxies;

         // check point: the sum of all historgram entries should be galaxy_num**2
    long int histsum = 0L;
    int      correct_value=1;
    for ( int i = 0; i < 360; ++i ) histsum += histogram_DD_g[i];
    printf("   Histogram DD : sum = %ld\n",histsum);
    if ( histsum != (numberOfGalaxies*numberOfGalaxies)) {printf("   Histogram sums should be %ld. Ending program prematurely\n", numberOfGalaxies*numberOfGalaxies);return(0);}


    histsum = 0L;
    for ( int i = 0; i < 360; ++i ) histsum += histogram_DR_g[i];
    printf("   Histogram DR : sum = %ld\n",histsum);
    if ( histsum != numberOfGalaxies*numberOfGalaxies ) {printf("   Histogram sums should be %ld. Ending program prematurely\n", numberOfGalaxies*numberOfGalaxies);return(0);}


    histsum = 0L;
    for ( int i = 0; i < 360; ++i ) histsum += histogram_RR_g[i];
    printf("   Histogram RR : sum = %ld\n",histsum);
    if ( histsum != (numberOfGalaxies*numberOfGalaxies)) {printf("   Histogram sums should be %ld. Ending program prematurely\n", numberOfGalaxies*numberOfGalaxies);return(0);}


    printf("   Omega values for the histograms:\n");
    float omega[360];
    for ( int i = 0; i < 10; ++i ) 
        if ( histogram_RR_g[i] != 0L )
           {
           omega[i] = (histogram_DD_g[i] - 2L*histogram_DR_g[i] + histogram_RR_g[i])/((float)(histogram_RR_g[i]));
           if ( i < 10 ) printf("      angle %.2f deg. -> %.2f deg. : %.3f\n", i*0.25, (i+1)*0.25, omega[i]);
           }

    FILE *out_file = fopen(argv[3],"w");
    if ( out_file == NULL ) printf("   ERROR: Cannot open output file %s\n",argv[3]);
    else
       {
       for ( int i = 0; i < 360; ++i ) 
           if ( histogram_RR_g[i] != 0L )
              fprintf(out_file,"%.2f  : %.3f\n", i*0.25, omega[i] ); 
       fclose(out_file);
       printf("   Omega values written to file %s\n",argv[3]);
       }

    free(real_rasc); free(real_decl);
    free(rand_rasc); free(rand_decl);


    }
    MPI_Finalize();

    

  return (0);
}
