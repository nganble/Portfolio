#
# THIS IS THE MAKE FILE FOR SERIAL CODE:
#
#

# First target is made by default when using "make", traditionally named "all"
all = wave2dSerial


# Set names of compilers once in case we need to change them
CC = gcc
#CXX = g++

# For parallel computing
CXX = g++

opt = -O3

# implicit (generic) rule to compile .C files
# $@ = file name of the target
# $< = name of the first prerequisite
%.o : %.C
	$(CXX) $(opt) -o $@ -c $<

#using open MP
wave2dSerialFiles = wave2dSerial.o
wave2dSerial: $(wave2dSerialFiles)
	$(CXX) $(opt) -o $@ $(wave2dSerialFiles) $(LIBS)

clean:; rm *.o