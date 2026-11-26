# Makefile for the MPI anisotropic wave solver
# Usage: make | make clean | make CXX=<compiler>   (override the MPI compiler wrapper)

CXX      = mpicxx
CXXFLAGS = -O3 -Wall

all: wave2dMPI

wave2dMPI: wave2dMPI.o
	$(CXX) $(CXXFLAGS) -o $@ $^

wave2dMPI.o: wave2dMPI.C ../common/parseCommand.h ../common/getCPU.h ../common/testCases.h getLocalIndexBounds.h writeMatlabArray.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f *.o wave2dMPI

.PHONY: all clean
