#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>


//histogram size
size_t HIST_SIZE = 360 * sizeof(unsigned long long int);
//galaxy number
long int N = 100000;
//input data
int    NoofReal;
int    NoofRand;
float *real_rasc, *real_decl;
float *rand_rasc, *rand_decl;
long int CPUMemory = 0L;
long int GPUMemory = 0L;
int  readdata(char *, char *);


__global__  void calcHist(unsigned long long int *DR,unsigned long long int *DD, unsigned long long int *RR, float *real_rasc, float *real_decl, float *rand_rasc, float *rand_decl, long int N)
{
   long int blockId = blockIdx.x
      + blockIdx.y * gridDim.x
      + gridDim.x * gridDim.y * blockIdx.z;

   long int threadId = blockId*(blockDim.x*blockDim.y*blockDim.z)
      + (threadIdx.z *(blockDim.x * blockDim.y))
      + (threadIdx.y * blockDim.x)
      + threadIdx.x;

   // printf("Thread id: %d\n", threadId);

   if (threadId >= 2*N*N) return;

   //calculate DR
   if (threadId < N * N){
      // printf("Calculate DR (%d, %d) pair from total %d\n", threadId/N, threadId%N, N);

      float acos_in = sin(real_decl[threadId/N])*sin(rand_decl[threadId%N])+cos(real_decl[threadId/N])*cos(rand_decl[threadId%N])*cos(real_rasc[threadId/N]-rand_rasc[threadId%N]);
      if (acos_in >= 0)
         acos_in = fmin(acos_in,1);
      else
         acos_in = fmax(acos_in,-1);
      
      atomicAdd(&DR[(int)(acos(acos_in) * 180 / acosf(-1.0f)*4)], 1);
   }else{

      long int i = (threadId - N*N)/N;
      long int j = threadId % N;

      if (i == j) return;

      if (i < j){
         // compute DD

         // printf("Calculate DD (%d, %d) pair from total %d\n", i, j, N);

         float acos_in = sin(real_decl[i])*sin(real_decl[j])+cos(real_decl[i])*cos(real_decl[j])*cos(real_rasc[i]-real_rasc[j]);
         if (acos_in >= 0)
            acos_in = fmin(acos_in,1);
         else
            acos_in = fmax(acos_in,-1);
         
         atomicAdd(&DD[(int)(acos(acos_in)* 180 / acosf(-1.0f) * 4)], 2);

      }
      else{
         // compute RR

         // printf("Calculate RR (%d, %d) pair from total %d\n", j, i, N);

         float acos_in = sin(rand_decl[j])*sin(rand_decl[i])+cos(rand_decl[j])*cos(rand_decl[i])*cos(rand_rasc[j]-rand_rasc[i]);
         if (acos_in >= 0)
            acos_in = fmin(acos_in,1);
         else
            acos_in = fmax(acos_in,-1);
         
         atomicAdd(&RR[(int)(acos(acos_in)* 180 / acosf(-1.0f) *  4)], 2);

      }

   }

}
void calcHistWrep(unsigned long long int *DR,unsigned long long int *DD, unsigned long long int *RR)
{
   unsigned long long int *histogramDR_gpu, *histogramDD_gpu, *histogramRR_gpu;
   float *real_rasc_gpu, *real_decl_gpu;
   float *rand_rasc_gpu, *rand_decl_gpu;
   long int *n;
  
   // input data is available in the arrays float real_rasc[], real_decl[], rand_rasc[], rand_decl[];
   // allocate memory on the GPU for input data
   cudaMalloc(&real_rasc_gpu, N * sizeof(float));
   cudaMalloc(&real_decl_gpu, N * sizeof(float));
   cudaMalloc(&rand_rasc_gpu, N * sizeof(float));
   cudaMalloc(&rand_decl_gpu, N * sizeof(float));
   
   // the number of galaxies is needed because of the work distribution logic
   cudaMalloc((void**)&n, sizeof(long int));

    // allocate memory on the GPU for histograms
   cudaMalloc((void **)&histogramDR_gpu, HIST_SIZE);
   cudaMalloc((void **)&histogramDD_gpu, HIST_SIZE);
   cudaMalloc((void **)&histogramRR_gpu, HIST_SIZE);
   cudaMemset(histogramDR_gpu, 0, HIST_SIZE);
   cudaMemset(histogramDD_gpu, 0, HIST_SIZE);
   cudaMemset(histogramRR_gpu, 0, HIST_SIZE);

   // initialize the data on GPU by copying the real and rand data to the GPU
   cudaMemcpy(real_rasc_gpu, real_rasc, N * sizeof(float), cudaMemcpyHostToDevice);
   cudaMemcpy(real_decl_gpu, real_decl, N * sizeof(float), cudaMemcpyHostToDevice);
   cudaMemcpy(rand_rasc_gpu, rand_rasc, N * sizeof(float), cudaMemcpyHostToDevice);
   cudaMemcpy(rand_decl_gpu, rand_decl, N * sizeof(float), cudaMemcpyHostToDevice);
   cudaMemcpy(n, &N, sizeof(long int), cudaMemcpyHostToDevice);
 
   long int threadsInBlock = 32*32;
   long int blocksInGrid = ( 2*N*N + threadsInBlock - 1 )/threadsInBlock;
   printf("Size of the blocks in grid: %ld\nThread number in each block:%ld\n",  blocksInGrid, threadsInBlock );
   calcHist<<<  blocksInGrid, threadsInBlock >>>(histogramDR_gpu,histogramDD_gpu,histogramRR_gpu, real_rasc_gpu, real_decl_gpu, rand_rasc_gpu, rand_decl_gpu, N);

   //move back the results from the gpu memory
   cudaMemcpy(DR, histogramDR_gpu, HIST_SIZE, cudaMemcpyDeviceToHost);
   cudaMemcpy(DD, histogramDD_gpu, HIST_SIZE, cudaMemcpyDeviceToHost);
   cudaMemcpy(RR, histogramRR_gpu, HIST_SIZE, cudaMemcpyDeviceToHost);
  
}

int main(int argc, char** argv){

   struct timeval _ttime;
   struct timezone _tzone;

   gettimeofday(&_ttime, &_tzone);
   double time_start = (double)_ttime.tv_sec + (double)_ttime.tv_usec/1000000.;

   //read real and random galaxies
   if ( readdata(argv[1], argv[2]) != 0 ) return(-1);

   // init global histogram arrays
   unsigned long long int *histogramDR, *histogramDD, *histogramRR;

   //allocate memory for the global histograms
   histogramDR=(unsigned long long int *) malloc(HIST_SIZE);
   histogramDD=(unsigned long long int *) malloc(HIST_SIZE);
   histogramRR=(unsigned long long int *) malloc(HIST_SIZE);

   //run the GPU kernel
   calcHistWrep(histogramDR, histogramDD, histogramRR);

   // galaxie distance of itself will be 0    
   histogramDD[0] += N;
   histogramRR[0] += N;    

    
   // check point: the sum of all historgram entries should be galaxy_num**2 
   long int histsum = 0L;

    for ( int i = 0; i < 360; ++i ) histsum += histogramDD[i];
    printf("   Histogram DD : sum = %ld\n",histsum);
    if ( histsum != (N*N)) {printf("   Histogram sums should be %ld. Ending program prematurely\n", N*N);return(0);}

    histsum = 0L;
    for ( int i = 0; i < 360; ++i ) histsum += histogramRR[i];
    printf("   Histogram RR : sum = %ld\n",histsum);
    if ( histsum != (N*N)) {printf("   Histogram sums should be %ld. Ending program prematurely\n", N*N);return(0);}

    histsum = 0L;
    for ( int i = 0; i < 360; ++i ) histsum += histogramDR[i];
    printf("   Histogram DR : sum = %ld\n",histsum);
    if ( histsum != N*N ) {printf("   Histogram sums should be %ld. Ending program prematurely\n", N*N);return(0);}

    printf("   Omega values for the histograms:\n");
    float omega[360];
    for ( int i = 0; i < 10; ++i ) 
        if ( (long long int)histogramRR[i] != 0LL )
           {
           omega[i] = ((long long int)histogramDD[i] - 2L*(long long int)histogramDR[i] + (long long int)histogramRR[i])/((float)(histogramRR[i]));
           if ( i < 10 ) printf("      angle %.2f deg. -> %.2f deg. : %.3f\n", i*0.25, (i+1)*0.25, omega[i]);
           }

    FILE *out_file = fopen(argv[3],"w");
    if ( out_file == NULL ) printf("   ERROR: Cannot open output file %s\n",argv[3]);
    else
       {
       for ( int i = 0; i < 360; ++i ) 
           if ( (long long int)histogramRR[i] != 0LL )
              fprintf(out_file,"%.2f  : %.3f\n", i*0.25, omega[i] ); 
       fclose(out_file);
       printf("   Omega values written to file %s\n",argv[3]);
       }

    free(real_rasc); free(real_decl);
    free(rand_rasc); free(rand_decl);

    gettimeofday(&_ttime, &_tzone);
    double time_end = (double)_ttime.tv_sec + (double)_ttime.tv_usec/1000000.;

    printf("   Wall clock run time    = %.1lf secs\n",time_end - time_start);


   return 0;

}

int readdata(char *argv1, char *argv2)
{
  int    i,linecount;
  char   inbuf[80];
  double ra, dec, dpi;
  FILE  *infil;
                                         
  printf("   Assuming data is in arc minutes!\n");
                          // phi   = ra/60.0 * dpi/180.0;
                          // theta = (90.0-dec/60.0)*dpi/180.0;
                          // otherwise use 
                          // phi   = ra * dpi/180.0;
                          // theta = (90.0-dec)*dpi/180.0;

  dpi = acos(-1.0);
  infil = fopen(argv1,"r");
  if ( infil == NULL ) {printf("Cannot open input file %s\n",argv1);return(-1);}

  linecount =0;
  while ( fgets(inbuf,80,infil) != NULL ) ++linecount;
  rewind(infil);

  printf("   %s contains %d galaxies\n",argv1, linecount-1);

  NoofReal = linecount-1;

  if ( NoofReal != 100000 ) {printf("Incorrect number of galaxies\n");return(1);}

  real_rasc = (float *)calloc(NoofReal,sizeof(float));
  real_decl = (float *)calloc(NoofReal,sizeof(float));
  CPUMemory += 2L*NoofReal*sizeof(float);

  fgets(inbuf,80,infil);
  sscanf(inbuf,"%d",&linecount);
  if ( linecount != 100000 ) {printf("Incorrect number of galaxies\n");return(1);}

  i = 0;
  while ( fgets(inbuf,80,infil) != NULL )
      {
      if ( sscanf(inbuf,"%lf %lf",&ra,&dec) != 2 ) 
         {
         printf("   Cannot read line %d in %s\n",i+1,argv1);
         fclose(infil);
         return(-1);
         }
      real_rasc[i] = (float)( ra/60.0*dpi/180.0);
      real_decl[i] = (float)(dec/60.0*dpi/180.0);
      ++i;
      }

  fclose(infil);

  if ( i != NoofReal ) 
      {
      printf("   Cannot read %s correctly\n",argv1);
      return(-1);
      }

  infil = fopen(argv2,"r");
  if ( infil == NULL ) {printf("Cannot open input file %s\n",argv2);return(-1);}

  linecount =0;
  while ( fgets(inbuf,80,infil) != NULL ) ++linecount;
  rewind(infil);

  printf("   %s contains %d galaxies\n",argv2, linecount-1);

  NoofRand = linecount-1;
  if ( NoofRand != 100000 ) {printf("Incorrect number of random galaxies\n");return(1);}

  rand_rasc = (float *)calloc(NoofRand,sizeof(float));
  rand_decl = (float *)calloc(NoofRand,sizeof(float));
  CPUMemory += 2L*NoofRand*sizeof(float);

  fgets(inbuf,80,infil);
  sscanf(inbuf,"%d",&linecount);
  if ( linecount != 100000 ) {printf("Incorrect number of random galaxies\n");return(1);}

  i =0;
  while ( fgets(inbuf,80,infil) != NULL )
      {
      if ( sscanf(inbuf,"%lf %lf",&ra,&dec) != 2 ) 
         {
         printf("   Cannot read line %d in %s\n",i+1,argv2);
         fclose(infil);
         return(-1);
         }
      rand_rasc[i] = (float)( ra/60.0*dpi/180.0);
      rand_decl[i] = (float)(dec/60.0*dpi/180.0);
      ++i;
      }

  fclose(infil);

  if ( i != NoofReal ) 
      {
      printf("   Cannot read %s correctly\n",argv2);
      return(-1);
      }

  return(0);
}