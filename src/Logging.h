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
#ifndef LOGGING_H_
#define LOGGING_H_

#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/helpers/exception.h>

using namespace log4cxx;
using namespace log4cxx::helpers;

/// define all loggers used in this application
static LoggerPtr initializationLogger(Logger::getLogger("initializationLogger"));
static LoggerPtr particleGenerationLogger(Logger::getLogger("initializationLogger.particleGenerationLogger"));
static LoggerPtr fileReaderLogger(Logger::getLogger("initializationLogger.particleGenerationLogger.fileReaderLogger"));
static LoggerPtr simulationInitializationLogger(Logger::getLogger("initializationLogger.simulationInitializationLogger"));
static LoggerPtr generalOutputLogger(Logger::getLogger("generalOutputLogger"));
static LoggerPtr testLogger(Logger::getLogger("generalOutputLogger.testLogger"));
static LoggerPtr simulationLogger(Logger::getLogger("simulationLogger"));

/// configure the loggers defined above
/// a configuration file may be used later, but as long as there is so little
/// it is more convenient to define it all right in the project itself
void configureLoggers();

#endif