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

#include "MolSim.h"

using namespace std;
using namespace utils;


//CPPUNIT's ugly static variable hiding concept...
CPPUNIT_TEST_SUITE_REGISTRATION(ParticleContainerTest);
CPPUNIT_TEST_SUITE_REGISTRATION(XMLFileReaderTest);

err_type MolSim::Init(int argc, char *argsv[])
{
	Timer timer;

	// parse command line arguments
	if(FAILED(parseLine(argc, argsv)))return E_INVALIDPARAM;

	// print hello message
	printHelloMessage();

	// test mode?
	if(state != AS_SIMULATION)return S_OK;

	// Init Simulation Data
	sim = new Simulation();
	if(!sim)
	{
		LOG4CXX_ERROR(initializationLogger, " >> failed to allocate mem");
		return E_OUTOFMEMORY;
	}

	// init random generator, use defined timer for that and offset it a bit,
	// just to make it look cool!
	srand((unsigned int)timer.getElapsedTime() + 0x362f0824);

	// the value desc should be later removed from Init!

	// set desc
	SimulationDesc desc;
	
	if(FAILED(sim->Init(desc)))return E_UNKNOWN;

	// use xml file
	if(FAILED(sim->CreateSimulationFromXMLFile(simFileName.c_str())))return E_INVALIDPARAM;

	// if viewer is active, set interior particle array size
	int size = sim->getParticleCount();
	
	// not for ICE
#ifndef ICE
	if(Viewer::Instance().IsRunning())Viewer::Instance().setParticleArraySize(sim->getParticleCount());
#endif
	// set grid
	ParticleContainer *pc = sim->getParticleContainer();
	
#ifndef ICE
	// valid?
	// beware, everything currently only in 2D!
	if(pc && Viewer::Instance().IsRunning())
	{
		utils::BoundingBox bb;
		bb = pc->getBoundingBox();
		int xcount = 1;
		int ycount = 1;

		// linked cell?
		if(pc->getType() == PCT_LINKEDCELL)
		{
			xcount = ((LinkedCellParticleContainer*)pc)->getCellExtent()[0];
			ycount = ((LinkedCellParticleContainer*)pc)->getCellExtent()[1];
		}

		Viewer::Instance().setGrid(bb.center[0] - bb.extent[0] * 0.5, bb.center[1] - bb.extent[1] * 0.5, bb.extent[0], bb.extent[1], xcount, ycount);
	}
#endif
	return S_OK;
}

// worker Thread, runs simulation
void MolSim::workerMain(MolSim *molsim)
{
	// simply run it
	molsim->sim->Run();
}

err_type MolSim::Run()
{
	switch(state)
	{
	case AS_NONE:
		return E_UNKNOWN;
		break;
	case AS_HELP:
		showHelp();
		break;
	case AS_SHOWTESTS:
		showTests();
		break;
	case AS_TESTS:
		RunTests();	
		break;
	case AS_SINGLETEST:
		runSingleTest(strTestCase);
		break;
	case AS_PTEST:
		{
			PerformanceTest PTest;

			// run performance tests...
			PTest.Run(simFileName.c_str());

			break;
		}
	case AS_SIMULATION:
		{
			// valid pointer?
			if(!sim)return E_NOTINITIALIZED;
#ifndef ICE	
			// if viewer is active(--viewer option enabled)
			// run multithreaded
			// else simply run simulation in one thread
			if(Viewer::Instance().IsRunning())
			{	
#ifdef USE_BOOST
				// create boost thread
				 boost::thread workerThread(workerMain, this);

				// GLFW Thread has to be the main thread, because Cocoa seems to have problems otherwise
				// run viewer's message loop
				Viewer::Instance().MessageLoop();

				// join threads
				workerThread.join();
#else
				GLFWthread worker;
			
				// create glfw based thread
				worker = glfwCreateThread((GLFWthreadfun)workerMain, this);

				// GLFW Thread has to be the main thread, because Cocoa seems to have problems otherwise
				// run viewer's message loop
				Viewer::Instance().MessageLoop();

				// join threads
				glfwWaitThread(worker, GLFW_WAIT);
#endif
			}
			else
			{
				// run simulation...
				err_type e = sim->Run();
				if(FAILED(e))return e;
			}	
#else
			err_type e = sim->Run();
			if(FAILED(e))return e;
#endif
			// show statistics...
			SimulationStatistics &stats = sim->getStatistics();
			showStatistics(stats);
			break;
		}
	}
	return S_OK;
}

err_type MolSim::Release()
{
	SAFE_DELETE(sim);

	return S_OK;
}

/// parses and verifies command line arguments
/// @param argc count of arguments
/// @param argsv string arguments
/// @return returns E_INVALIDPARAM if argument count is wrong(not four) or second
///					or third argument are of non number format
///			returns E_FILENOTFOUND if no file exists
err_type MolSim::parseLine(int argc, char *argsv[])
{
	state = AS_NONE;
	
	// Syntax is molsim scene.xml (--viewer)
	// where t_end and delta_t denote a floating point value
	if(argc != 2 && argc != 3)
	{
		LOG4CXX_ERROR(simulationInitializationLogger, ">> error: invalid count of arguments");
		printUsage();
		return E_INVALIDPARAM;
	}

	// check if called with -test
	if(argc == 2)
	{
		if(strcmp(argsv[1], "-test") == 0)
			state = AS_TESTS;
		else if(strcmp(argsv[1], "-help") == 0)
			state = AS_HELP;
		else if(strcmp(argsv[1], "-showtests") == 0)
			state = AS_SHOWTESTS;
		else if(strcmp(argsv[1], "-ptest") == 0)
		{
			LOG4CXX_ERROR(generalOutputLogger, ">> error: please specify a test file!");
			return E_INVALIDPARAM;
		}
	}
	else if(argc == 3)
	{
		if(strcmp(argsv[1], "-test") == 0)
		{
			state = AS_SINGLETEST;
			strTestCase = argsv[2];
		}
		// can be displayed before or after the filename
		else if(strcmp(argsv[1], "--viewer") == 0 || strcmp(argsv[2], "--viewer") == 0)
		{
#ifndef ICE	
			// initialize viewer
			Viewer::Instance().InitAndDisplay();
#else
			LOG4CXX_ERROR(generalOutputLogger, ">> sorry, but currently X11 output does not work for ICE!");
			return E_INVALIDPARAM;
#endif		
		}
		else if(strcmp(argsv[1], "-ptest") == 0)
		{
			state = AS_PTEST;
			simFileName = argsv[2];
			if(strcmp(utils::getFileExtension(argsv[2]), "xml") != 0)
			{
				LOG4CXX_ERROR(generalOutputLogger, ">> error: invalid file format!");
				return E_INVALIDPARAM;
			}
		}
		else
		{
			LOG4CXX_ERROR(simulationInitializationLogger, ">> error: invalid argument");
			printUsage();
			return E_INVALIDPARAM;
		}
	}

	if(state == AS_NONE)
		if(argc == 2 || argc == 3 )
		{
			int index = 1;
		
			//set index according where a possible -- option may be set
			if(argsv[1][0] == '-')index = 2;
			else index = 1;

			// parse file

			// check if file exists
			if(!fileExists(argsv[index]))
			{
				LOG4CXX_ERROR(simulationInitializationLogger, "error: file doesn't exist!");
				printUsage();
				return E_FILENOTFOUND;
			}

			// has file correct ending?
			if(strcmp(utils::getFileExtension(argsv[index]), "xml") != 0)
			{
				LOG4CXX_ERROR(simulationInitializationLogger, "error: file has no .xml extension!");
				printUsage();
				return E_FILEERROR;
			}

			// everything ok...
			// set name
			simFileName = argsv[index];

			//run Simulation
			state = AS_SIMULATION;
		}

	return S_OK;
}

void MolSim::showHelp()
{
	LOG4CXX_INFO(generalOutputLogger, " >> "<<"options\t\t\tdescription");
	LOG4CXX_INFO(generalOutputLogger, "    "<<" <file> "<<"\t"<<"run simulation according to <file>");
	LOG4CXX_INFO(generalOutputLogger, "    "<<"-help"<<"\t\t\t"<<"show help");
	LOG4CXX_INFO(generalOutputLogger, "    "<<"-test <name>"<<"\t\t"<<"run single test case or leave\n\t\t\t\t<name> blank to run all tests");
	LOG4CXX_INFO(generalOutputLogger, "    "<<"-showtests"<<"\t\t\t"<<"list all avaliable tests by name");
	LOG4CXX_INFO(generalOutputLogger, "    "<<"-ptest <file>"<<"\t\t\t"<<"run performance tests for the given file");
	LOG4CXX_INFO(generalOutputLogger, "    "<<"--viewer"<<"\t\t\t"<<"start with builtin OpenGL viewer");
	
}

void MolSim::showTestTree(CppUnit::Test *node, int level, int maxlevel)
{
	// no children?
	if(node->getChildTestCount() == 0)
	{
		cout<<"    ";

		// space indent
		for(int i = 0; i < level; i++)cout<<"  ";

		LOG4CXX_INFO(generalOutputLogger, "- "<<node->getName());
	}
	else
	{
		cout<<"    ";

		// space indent
		for(int i = 0; i < level; i++)cout<<"  ";
		
		// max level reached?
		if(level >= maxlevel)cout<<"..."<<endl;
		else
		{
			int count = node->getChildTestCount();
			LOG4CXX_INFO(generalOutputLogger, "+ "<<node->getName());

			//print children
			for(int i = 0; i < count; i++)
				showTestTree(node->getChildTestAt(i), level + 1, maxlevel);
		}
			
	}
}

void MolSim::showTests()
{
	CppUnit::TestSuite *suite = NULL;
	
	// show some info
	LOG4CXX_INFO(generalOutputLogger, " >> displaying test structure, to run a single test\n    or a bevy of tests call molsym -test <name>");
	
	// recursive tree display
	CppUnit::Test* root = CppUnit::TestFactoryRegistry::getRegistry().makeTest();

	showTestTree(root, 0, MAX_LEVEL);
}

void MolSim::printHelloMessage()
{
	line();
	LOG4CXX_INFO(generalOutputLogger, "MolSim for PSE");
#ifdef DEBUG
	cout << "compiled: "<<__DATE__<<"  "<<__TIME__<<endl;
#endif

	line();
	cout << endl;
	LOG4CXX_INFO(generalOutputLogger, "(c) 2012 by F.Dietz & L.Spiegelberg"); 
	cout << endl;
	LOG4CXX_INFO(generalOutputLogger, "Molecular Simulator handling *.txt files");
	line();
	cout << endl;
}

void MolSim::printUsage()
{
	LOG4CXX_INFO(generalOutputLogger, "   usage: molsim <file>");
	LOG4CXX_INFO(generalOutputLogger, "          molsim -help");
}

void MolSim::RunTests()
{
	LOG4CXX_INFO(testLogger, " >> start running test suite");

	CppUnit::TextUi::TestRunner runner;
	CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
	runner.addTest( registry.makeTest() );
	bool wasSuccessful = runner.run( "", false );
	
	if(wasSuccessful) {
		LOG4CXX_INFO(testLogger, " >> all tests succeeded! ");
	}
	else LOG4CXX_INFO(testLogger, " >> test suite failed! ");

}

void MolSim::runSingleTest(string s)
{
	//find test according to name
	CppUnit::Test* root = CppUnit::TestFactoryRegistry::getRegistry().makeTest();

	CppUnit::Test* test = NULL;

	try
	{
		// find test
		// maybe later some custom string matching routine would be desirable,
		// as the testcase has to be specified very, very exactly
		test = root->findTest(s);
	}
	catch(...)
	{
		// no test found
		LOG4CXX_ERROR(testLogger, " >> "<<"no test "<<s<<" could be found! ");
	}

	if(test)
	{
		// run single test
		CppUnit::TextUi::TestRunner runner;
		runner.addTest(test);
		bool wasSuccessful = runner.run( "", false);
	
		if(wasSuccessful) {
			LOG4CXX_INFO(testLogger, " >> test succeeded! ");
		}
		else LOG4CXX_ERROR(testLogger, " >> test failed! ");
	}
}

void MolSim::showStatistics(SimulationStatistics& s)
{
	LOG4CXX_INFO(generalOutputLogger, "Statistics");
	utils::line();
	LOG4CXX_INFO(generalOutputLogger, "total particle count:\t"<<s.particle_count);
	LOG4CXX_INFO(generalOutputLogger, "iteration count: \t"<<s.step_count);
	
	// for less than 10 seconds print detailed time out
	if(s.time < 10.0)
	{
		LOG4CXX_INFO(generalOutputLogger, "total time:\t\t"<<s.time * 1000<<" ms");
	}
	else
	{
		LOG4CXX_INFO(generalOutputLogger, "total time:\t\t"<<utils::secondsToHMS((int)s.time).c_str());
	}
	
	LOG4CXX_INFO(generalOutputLogger, "timer per step:\t\t"<<1000 * s.timeperstep<<" ms");
}
