# This is a makefile template


# Compiler
# --------
#CC=g++
# to compile for ICE uncomment
CC = g++ -DICE -DUSE_PAPI -D_OPENMP
#CC=icpc -DICE

include files.mk


# Compiler flags
# -------------------------------------------------------------------------
## intel flags
#CFLAGS=-I/opt/include -g -O3 -openmp
## gcc flags
CFLAGS=-I/opt/include -g -O3 
## gcc flags for gprof
#CFLAGS=-I/opt/include -g -O3 -pg

# Linker flags
# ------------
#g++
LDFLAGS= -L/opt/lib -L$$HOME/lib -L/lrz/sys/libraries/xerces/3.1/lib -lpapi -lXrandr -lxerces-c -lrt -llog4cxx -lcppunit -lglfw -lGL -lX11 -lpthread -lGLU -fopenmp 

#icc
#LDFLAGS= -static-intel -L/opt/lib -L$$HOME/lib -L/lrz/sys/libraries/xerces/3.1/lib -lirc -lpapi -lXrandr -lxerces-c -lrt -llog4cxx -lcppunit -lglfw -lGL -lX11 -lpthread -lGLU -openmp 
## linker flags for gprof
#LDFLAGS=-L/opt/lib -lxerces-c -lrt -llog4cxx -lcppunit -pg

INCLUDES= -I./src -I$$HOME/include

OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=molsim

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@  
clean:
	rm $(OBJECTS)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

