#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <mpi.h>

#include "../common/galaxy_utils.h"


int main(int argc, char* argv[]) 
{
 
  float *real_rasc, *real_decl, *rand_rasc, *rand_decl;
  long int MemoryAllocatedCPU = 0L;
  long int numberOfGalaxies = 100000;
  int num_of_angle_intervals = 360;

  float  pif = acosf(-1.0f);

  /* this will store partial histogram calculations
  by summing up the columns groupping the indixes in mod 3, we will get the final DR,DD and RR histograms */
  long int ** histograms;

  int np, me;
  double starttime, endtime;

  MPI_Init(&argc, &argv);                /* Initialize MPI */
  MPI_Comm_size(MPI_COMM_WORLD, &np);    /* Get nr of processes */
  MPI_Comm_rank(MPI_COMM_WORLD, &me);    /* Get own identifier */

    
  /* Check that we run on at least four processes */
  if (np<4 ) {

    if (me == 0) {
      printf("You have to run the program on at least four processes\n");
    }
    MPI_Finalize();
    exit(0);
  }

  // store right ascension and declination for real galaxies here
  real_rasc        = (float *)calloc(numberOfGalaxies, sizeof(float));
  real_decl        = (float *)calloc(numberOfGalaxies, sizeof(float));

  // store right ascension and declination for synthetic random galaxies here
  rand_rasc        = (float *)calloc(numberOfGalaxies, sizeof(float));
  rand_decl        = (float *)calloc(numberOfGalaxies, sizeof(float));
  
  if (me == 0) {

      // read input data from files given on the command line
      if ( parseargs_readinput(argc, argv, numberOfGalaxies, real_rasc, real_decl, rand_rasc, rand_decl) != 0 ) {printf("   Program stopped.\n");return(0);}
      printf("   Input data read, now calculating histograms\n");

      starttime = MPI_Wtime();
  }
  
  //broadcast to all other processes from root
  MPI_Bcast(real_rasc, numberOfGalaxies, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(real_decl, numberOfGalaxies, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(rand_rasc, numberOfGalaxies, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(rand_decl, numberOfGalaxies, MPI_FLOAT, 0, MPI_COMM_WORLD);

  if (me != 0)
    printf("   process %d running rand_decl %f\n", me, rand_decl[3]);

    return(0);
}
