# This is a makefile template


# Compiler
# --------
CC=g++
#CC=icpc

include files.mk


# Compiler flags
# -------------------------------------------------------------------------
CFLAGS=-I/opt/include -g -O3

# Linker flags
# ------------
LDFLAGS= -L/opt/lib -lxerces-c -lrt -llog4cxx -lcppunit 



INCLUDES= -I./src -I./libxsd

OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=molsim

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@ 

clean:
	rm $(OBJECTS)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

