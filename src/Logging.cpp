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
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------

#include "Logging.h"

using namespace log4cxx;

void configureLoggers() {
	// Load configuration file
	// this currently just sets the root logger to output to the console
	// and also to a file called "molsim.log"
	// this makes the differentiation between loggers meaningless for now
	// but if we decide to make several types of output go to different files
	// this will be very easy and fast to change
	PropertyConfigurator::configure("config.property");
}