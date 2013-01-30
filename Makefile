# default make file for linux/mac platform

# Compiler
# --------
CC=g++
#CC=icpc

include files.mk


# Compiler flags
# -------------------------------------------------------------------------
## intel flags
#CFLAGS=-I/opt/include -g -O3 -openmp
## g++ flags
CFLAGS=-I/opt/include -g -O3 -fopenmp
## g++ flags for gprof
#CFLAGS=-I/opt/include -g -O3 -pg

# Linker flags
# ------------
LDFLAGS= -L/opt/lib -lxerces-c -lrt -llog4cxx -lcppunit -lglfw -lGL -lX11 -lpthread -lGLU -fopenmp 

## linker flags for gprof
#LDFLAGS=-L/opt/lib -lxerces-c -lrt -llog4cxx -lcppunit -pg

INCLUDES= -I./src 

OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=molsim

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@ 


clean:
	rm $(OBJECTS)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@


