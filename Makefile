# This is a makefile template


# Compiler
# --------
#CC=g++
CC=icpc

include files.mk


# Compiler flags
# -------------------------------------------------------------------------
## intel flags
CFLAGS=-I/opt/include -g -O3 -ip -ipo -fast
## gcc flags
#CFLAGS=-I/opt/include -g -O3
## gcc flags for gprof
#CFLAGS=-I/opt/include -g -O3 -pg

# Linker flags
# ------------
LDFLAGS= -L/opt/lib -lxerces-c -lrt -llog4cxx -lcppunit 

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

