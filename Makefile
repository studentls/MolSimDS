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
LDFLAGS= -L/opt/lib -lxerces-c 

SOURCES=\
src/MolSim.cpp\
src/outputWriter/XYZWriter.cpp\
src/outputWriter/VTKWriter.cpp\
src/outputWriter/vtk-unstructured.cpp\
src/FileReader.cpp\
src/Particle.cpp\

INCLUDES= -I./src -I./libxsd

OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=MolSim

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@ 

clean:
	rm $(OBJECTS)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

