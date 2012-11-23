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
#include "utils/Helper.h"

using namespace log4cxx;

static const char g_szDefaultConfigFilename[] = "config.property";

/// generates a simple config file
void generateDefaultConfigFile()
{		
	FILE *pFile = fopen(g_szDefaultConfigFilename, "w");

	// sucessful?
	if(!pFile)return; // maybe do some error handling here...

	// simple config layout...
	fprintf(pFile,	"log4j.rootLogger=trace, stdout, R\n" \
					"\n" \
					"log4j.appender.stdout=org.apache.log4j.ConsoleAppender\n" \
					"log4j.appender.stdout.layout=org.apache.log4j.PatternLayout\n" \
					"\n" \
					"# Pattern to output the caller's file name and line number.\n" \
					"log4j.appender.stdout.layout.ConversionPattern=%%m%%n\n" \
					"\n" \
					"log4j.appender.R=org.apache.log4j.RollingFileAppender\n" \
					"log4j.appender.R.File=molsim.log\n" \
					"\n" \
					"log4j.appender.R.MaxFileSize=100KB\n" \
					"# Keep one backup file" \
					"log4j.appender.R.MaxBackupIndex=1\n" \
					"\n" \
					"log4j.appender.R.layout=org.apache.log4j.PatternLayout\n" \
					"log4j.appender.R.layout.ConversionPattern=%%5p [%%t] (%%F:%%L) - %%m%%n\n");

	fclose(pFile);	
}

/// configures Loggers according to the config.property file
void configureLoggers() {

	//it can happen, that maybe a user deletes the config file, in order to prevent program from crashing, write a simple config file
	if(!utils::fileExists(g_szDefaultConfigFilename))generateDefaultConfigFile();


	// Load configuration file
	// this currently just sets the root logger to output to the console
	// and also to a file called "molsim.log"
	// this makes the differentiation between loggers meaningless for now
	// but if we decide to make several types of output go to different files
	// this will be very easy and fast to change
	PropertyConfigurator::configure(g_szDefaultConfigFilename);
}