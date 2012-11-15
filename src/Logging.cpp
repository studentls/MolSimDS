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

LoggerPtr initializationLogger(Logger::getLogger("initializationLogger"));
LoggerPtr particleGenerationLogger(Logger::getLogger("initializationLogger.particleGenerationLogger"));
LoggerPtr fileReaderLogger(Logger::getLogger("initializationLogger.particleGenerationLogger.fileReaderLogger"));
LoggerPtr simulationInitializationLogger(Logger::getLogger("initializationLogger.simulationInitializationLogger"));
LoggerPtr generalOutputLogger(Logger::getLogger("generalOutputLogger"));
LoggerPtr testLogger(Logger::getLogger("generalOutputLogger.testLogger"));
LoggerPtr simulationLogger(Logger::getLogger("simulationLogger"));

void configureLoggers() {
	// this currently just sets the root logger to output to the console
	BasicConfigurator::configure();
	// this makes the differentiation between loggers meaningless for now
	// but if we decide to make several types of output go to different files
	// this will be very easy and fast to change
}
