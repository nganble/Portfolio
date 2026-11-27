# Makefile for the 2D Poisson solvers (SOR): serial, MPI, CUDA, MPI+CUDA
# Usage:
#   make                 serial and MPI (any machine with gcc and an MPI library)
#   make cuda            CUDA (needs nvcc)
#   make mpi-cuda        MPI+CUDA (needs nvcc)
#   make POWER9=1 all    IBM Power9 cluster (XL compiler, Power9 timer)
#   make CFLAGS="-DTEST=EXACT_SOLUTION -DSOLUTION=GAUSSIAN"   other test case
#   make clean

POWER9 ?= 0

# --- Compilers ---
CC = gcc
ifeq ($(POWER9),1)
  MPICC     = mpixlc
  CUDA_PATH ?= /usr/local/cuda-11.2
else
  MPICC     = mpicc
endif
NVCC = nvcc

# --- Flags ---
OPTFLAGS ?= -O3
ALL_CFLAGS = $(OPTFLAGS) -DUSE_POWER9=$(POWER9) $(CFLAGS)
ARCH_FLAGS ?= -arch=sm_70   # V100 GPU
CUDA_PATH  ?= /usr/local/cuda
CUDA_LIBS  = -L$(CUDA_PATH)/lib64 -lcudadevrt -lcudart -lstdc++
MATH_LIB   = -lm

HEADERS = poisson2d_common.h poisson2d_cuda.h poisson2d_mpi.h timing.h

.PHONY: default all serial mpi cuda mpi-cuda clean
default: serial mpi
all: serial mpi cuda mpi-cuda

# --- Serial ---
serial: poisson2d_serial
poisson2d_serial: poisson2d_serial.c
	$(CC) $(ALL_CFLAGS) $< -o $@ $(MATH_LIB)

# --- MPI ---
mpi: poisson2d_mpi
poisson2d_mpi: poisson2d_mpi.c
	$(MPICC) $(ALL_CFLAGS) $< -o $@ $(MATH_LIB)

# --- CUDA ---
cuda: poisson2d_cuda
poisson2d_cuda: poisson2d_cuda.cu
	$(NVCC) $(ALL_CFLAGS) $(ARCH_FLAGS) $< -o $@ $(MATH_LIB)

# --- MPI + CUDA: host code with MPICC, kernels with NVCC, linked with MPICC ---
mpi-cuda: poisson2d_mpi_cuda
poisson2d_mpi_cuda: poisson2d_mpi_cuda.c poisson2d_mpi_cuda.cu $(HEADERS)
	$(MPICC) $(ALL_CFLAGS) -c poisson2d_mpi_cuda.c -o poisson2d_mpi_cuda_host.o
	$(NVCC)  $(ALL_CFLAGS) $(ARCH_FLAGS) -c poisson2d_mpi_cuda.cu -o poisson2d_mpi_cuda_dev.o
	$(MPICC) $(ALL_CFLAGS) poisson2d_mpi_cuda_host.o poisson2d_mpi_cuda_dev.o -o $@ $(CUDA_LIBS) $(MATH_LIB)

clean:
	rm -f poisson2d_serial poisson2d_mpi poisson2d_cuda poisson2d_mpi_cuda *.o
