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

#include "Logging.h"
#include "Simulation.h"

#include "outputWriter/XYZWriter.h"
#include "outputWriter/VTKWriter.h"
#include "FileReader.h"
#include "utils/utils.h"

#include <list>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <string>

//include tests
#include "ParticleContainerTest.h"

/// specifies in which state the application is currently
enum ApplicationState
{
	AS_NONE,		/// no state 
	AS_SIMULATION,	/// a simulation will be run
	AS_TESTS,		/// run all tests
	AS_SINGLETEST,	/// run single test
	AS_HELP,		/// show help
	AS_SHOWTESTS	/// show all avaliable tests
};

/// maximal level to display test structure
#define MAX_LEVEL	5

/// application's main class, handling all logic stuff
class MolSim
{
private:

	Simulation			*sim;

	ApplicationState	state;

	std::string			strTestCase;

	/// parse Argument line, check for variable ranges
	err_type			parseLine(int argc, char *argsv[]);

	/// print sime nice logo
	void				printHelloMessage();

	/// print Usage
	void				printUsage();

	/// run Tests
	void				RunTests();

	/// recursive function to display test tree
	void				showTestTree(CppUnit::Test *node, int level, int maxlevel);

	/// print all avaliable tests
	void				showTests();

	/// run single test case
	void				runSingleTest(std::string s);

	/// show help, list all avaliable commands
	void				showHelp();

public:

	MolSim():sim(NULL), state(AS_NONE)			{}

	/// sets up simulator and parses any given arguments
	err_type			Init(int argc, char *argsv[]);

	/// runs a simulation
	err_type			Run();

	/// cleans up
	err_type			Release();	
};

#endif
