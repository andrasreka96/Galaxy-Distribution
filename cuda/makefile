galaxy_cuda: galaxy_cuda.cu
	nvcc -O3 -arch=sm_70 --ptxas-options=-v -o galaxy galaxy_cuda.cu -lm
