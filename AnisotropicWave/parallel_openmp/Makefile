# Makefile for the OpenMP anisotropic wave solver
# Usage: make | make clean | make CXX=<compiler>   (override the compiler)

# macOS: Homebrew GCC (Apple clang has no OpenMP); Linux: system GCC
ifeq ($(shell uname -s),Darwin)
  CXX = g++-16
else
  CXX = g++
endif
CXXFLAGS = -O3 -Wall -fopenmp

all: wave2dOpenMP

wave2dOpenMP: wave2dOpenMP.o
	$(CXX) $(CXXFLAGS) -o $@ $^

wave2dOpenMP.o: wave2dOpenMP.C ../common/parseCommand.h ../common/getCPU.h ../common/testCases.h ../common/grid.h ../common/stencil.h ../common/output.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f *.o wave2dOpenMP

.PHONY: all clean
