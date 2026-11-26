# Makefile for the serial anisotropic wave solver
# Usage: make | make clean | make CXX=clang++ (override the compiler)

CXX= g++
CXXFLAGS = -O3 -Wall


# First target is made by default when using "make", traditionally named "all"
all: wave2dSerial

wave2dSerial: wave2dSerial.o
	$(CXX) $(CXXFLAGS) -o $@ $^

wave2dSerial.o: wave2dSerial.C ../common/parseCommand.h ../common/getCPU.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<


clean:
	rm -f *.o wave2dSerial

.PHONY: all clean
