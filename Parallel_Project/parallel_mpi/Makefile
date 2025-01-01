#
# THIS MAKE FILE USES THE FOLLOWING ENVIROMENTAL VARIABLE:
#

# First target is made by default when using "make", traditionally named "all"
# all = heat2d
all = wave2dMPI

# Set names of compilers once in case we need to change them
CC = gcc
CXX = mpic++

opt = -O3


# include -DUSE_PPP with compile flags:
CCFLAGS= -DUSE_PPP -fPIC $(opt)

# link to libPpp instead of libApp
#AppLibraries = -Wl,-rpath,$(APlusPlus)/lib -L$(APlusPlus)/lib -lPpp -lPpp_static

# List of all libraries
LIBS = $(AppLibraries)

# implicit (generic) rule to compile .C files
# $@ = file name of the target
# $< = name of the first prerequisite
%.o : %.C
	$(CXX) $(CCFLAGS) -o $@ -c $<

# 2D heat equation, explicit time-stepping, A++ arrays
wave2dMPIFiles = wave2dMPI.o
wave2dMPI: $(wave2dMPIFiles)
	$(CXX) $(CCFLAGS) -o $@ $(wave2dMPIFiles) $(LIBS)
	
clean:; rm *.o