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
 
  float *real_rasc, *real_decl, *rand_rasc, *rand_decl, *real_rasc_partial;
  long int numberOfGalaxies = 100000;
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


  /* Start measuring time */
  if (myid == 0) {
    starttime = MPI_Wtime();
  }

  int p0_n =  np%3 == 0 ? np/3 : np/3 + np%3;
  int p1_n = np/3;
  int p2_n = p1_n;

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
            float angle = angle_between(real_rasc[i], real_decl[i], rand_rasc[j], rand_decl[j]) * radian_to_angle;
            histogram_DD_l[ (int)(angle*quantum_reciprocal)] += 2;
         }

  }
  else{
      
    //calculates RR

    int starti = (myid/3) * (numberOfGalaxies/p2_n);
    for(int i=starti; i<starti+numberOfGalaxies/p2_n && i<numberOfGalaxies; ++i)
        for(int j=i+1; j<numberOfGalaxies; ++j){
          float angle = angle_between(real_rasc[i], real_decl[i], rand_rasc[j], rand_decl[j]) * radian_to_angle;
          histogram_RR_l[ (int)(angle*quantum_reciprocal)] += 2;
        }


  }

  MPI_Reduce(histogram_DD_l, histogram_DD_g, numberOfGalaxies ,MPI_INT ,MPI_SUM ,0, MPI_COMM_WORLD);
  MPI_Reduce(histogram_DR_l, histogram_DR_g, numberOfGalaxies ,MPI_INT ,MPI_SUM ,0, MPI_COMM_WORLD);
  MPI_Reduce(histogram_RR_l, histogram_RR_g, numberOfGalaxies ,MPI_INT ,MPI_SUM ,0, MPI_COMM_WORLD);


  if(myid==0)
    {
        printf("globalsum1 = %ld \n",histogram_DD_g[0]);
        printf("globalsum2 = %ld \n",histogram_DD_g[1]);
        printf("globalsum3 = %ld \n",histogram_DD_g[2]);

    }

    MPI_Finalize();

  return (0);
}
