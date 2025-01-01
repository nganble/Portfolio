#
# THIS IS THE MAKE FILE FOR OPENMP CODE:
#
#

# First target is made by default when using "make", traditionally named "all"
all = wave2dOpenMP


# Set names of compilers once in case we need to change them
CC = gcc
#CXX = g++

# For parallel computing
CXX = g++ -fopenmp

opt = -O3

# implicit (generic) rule to compile .C files
# $@ = file name of the target
# $< = name of the first prerequisite
%.o : %.C
	$(CXX) $(opt) -o $@ -c $<

#using open MP
wave2dOpenMPFiles = wave2dOpenMP.o
wave2dOpenMP: $(wave2dOpenMPFiles)
	$(CXX) $(opt) -o $@ $(wave2dOpenMPFiles) $(LIBS)

clean:; rm *.o