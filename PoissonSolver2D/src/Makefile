# --- Compilers & Tools ---
CC = gcc
#MPICC = mpicc
MPICC = mpixlc
MPIXL = mpixlc
NVCC = nvcc

# --- Paths & Libraries ---
MATH_LIB = -lm
CUDA_PATH ?= /usr/local/cuda-11.2
CUDA_LIBS = -L$(CUDA_PATH)/lib64 -lcudadevrt -lcudart -lstdc++  


# --- Flags ---
OPTFLAGS ?= -O3
ARCH_FLAGS = -arch=sm_70

ALL_CFLAGS = $(OPTFLAGS) $(CFLAGS)

# --- Default target ---
.PHONY: all clean serial mpi cuda  mpi-cuda
all: serial mpi cuda mpi-cuda

# --- Targets ---
# 1. Serial Version
serial: poisson2d_serial
poisson2d_serial: poisson2d_serial.c
	$(CC) $(ALL_CFLAGS) $< -o $@ $(MATH_LIB)

# 2. Only MPI Version
mpi: poisson2d_mpi
poisson2d_mpi: poisson2d_mpi.c
	$(MPICC) $(ALL_CFLAGS) $< -o $@ $(MATH_LIB)

# 3. Only CUDA Version
cuda: poisson2d_cuda
poisson2d_cuda: poisson2d_cuda.cu
	$(NVCC) $(ALL_CFLAGS) $(ARCH_FLAGS) $< -o $@ $(MATH_LIB)

# 4. MPI + CUDA (Power9 / XL Compiler)
mpi-cuda: poisson2d_mpi_cuda
poisson2d_mpi_cuda: poisson2d_mpi_cuda.c poisson2d_mpi_cuda.cu
	$(MPIXL) $(ALL_CFLAGS) poisson2d_mpi_cuda.c -c -o poisson2d_mpi_cuda_xlc.o
	$(NVCC)  $(ALL_CFLAGS) $(ARCH_FLAGS) poisson2d_mpi_cuda.cu -c -o poisson2d_mpi_cuda_nvcc.o
	$(MPIXL) $(ALL_CFLAGS) poisson2d_mpi_cuda_xlc.o poisson2d_mpi_cuda_nvcc.o -o $@ $(CUDA_LIBS) $(MATH_LIB)

# --- Utility ---
clean:
	rm -f poisson2d_serial poisson2d_mpi poisson2d_cuda poisson2d_mpi_cuda *.o





# # CFLAGS =-O3

# mpi-cuda: poisson2d_mpi_cuda.c poisson2d_mpi_cuda.cu
# 	mpixlc -O3 poisson2d_mpi_cuda.c -c -o poisson2d_mpi_cuda_xlc.o
# 	nvcc -O3 -arch=sm_70 poisson2d_mpi_cuda.cu -c -o poisson2d_mpi_cuda_nvcc.o 
# 	mpixlc -O3 poisson2d_mpi_cuda_xlc.o poisson2d_mpi_cuda_nvcc.o -o poisson2d_mpi_cuda -L/usr/local/cuda-11.2/lib64/ -lcudadevrt -lcudart -lstdc++ 

# serial: poisson2d_serial.c
# 	$(CC) $(OPTFLAGS) $(CFLAGS) poisson2d_serial.c -o poisson2d_serial

# serial_compute: poisson2d_serial.c
# 	$(CC) $(CFLAGS) -DTEST=COMPUTE_SOLUTION $< -o poisson2d_serial

# cuda: poisson2d_cuda.cu
# 	nvcc $(OPTFLAGS) $(CFLAGS) poisson2d_cuda.cu -o poisson2d_cuda

# mpi: poisson2d_mpi.c
# 	mpicc $(OPTFLAGS) $(CFLAGS) $< -o poisson2d_mpi



# serial_exact_poly: poisson2d_serial.c
# 	$(CC) $(CFLAGS) -DTEST=EXACT_SOLUTION -DSOLUTION=POLY $< -o poisson2d_serial

# serial_exact_gaussian: poisson2d_serial.c
# 	$(CC) $(CFLAGS) -DTEST=EXACT_SOLUTION -DSOLUTION=GAUSSIAN $< -o poisson2d_serial

# 	