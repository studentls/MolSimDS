//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File MolSim.cpp
// contains main class MolSim
//------------------------------------------------------------------------------------------------

#ifndef MOLSIM_HEADER_
#define MOLSIM_HEADER_

#include "Simulation.h"

#include "outputWriter/XYZWriter.h"
#include "outputWriter/VTKWriter.h"
#include "FileReader.h"
#include "utils/utils.h"

#include <list>
#include <cstring>
#include <cstdlib>
#include <iostream>

///
/// application's main class, handling all logic stuff
///
class MolSim
{
private:

	Simulation	*sim;

	err_type	parseLine(int argc, char *argsv[]);
	void		printHelloMessage();

public:
	MolSim():sim(NULL)			{}

	/// sets up simulator and parses any given arguments
	err_type Init(int argc, char *argsv[]);

	/// runs a simulation
	err_type Run();

	/// cleans up
	err_type Release();
	
};


#endif
