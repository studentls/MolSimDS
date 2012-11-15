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

	// set desc
	SimulationDesc desc;

	desc.output_fmt = SOF_VTK;
	desc.start_time = 0.0;
	desc.end_time = atof(argsv[2]);
	desc.delta_t = atof(argsv[3]);

	// test if input values make sense...
	int nSteps = (desc.end_time - desc.start_time) / desc.delta_t;

	//if more than 1000000 steps warn user
	if(nSteps > 1000000)
	{
		LOG4CXX_INFO(generalOutputLogger, " >> The Simulation will need approximately "<<nSteps<<" steps to finish");
		LOG4CXX_INFO(generalOutputLogger, "    finishing calculation may take a very long time, proceed ? (y/n)");
		std::string s;
		cin>>s;
		LOG4CXX_INFO(generalOutputLogger, "input: " << s);

		if(s[0] == 'y') {
			//...
		}
		else if(s[0] == 'n')return E_INVALIDPARAM;
		else 
		{
			LOG4CXX_INFO(generalOutputLogger, " >> please enter next time (y/n)");
			return E_INVALIDPARAM;
		}
		
	}

	
	if(FAILED(sim->Init(desc)))return E_UNKNOWN;

	if(FAILED(sim->AddParticlesFromFile(argsv[1])))return E_INVALIDPARAM;

	return S_OK;
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
	case AS_SIMULATION:
		// valid pointer?
		if(!sim)return E_NOTINITIALIZED;
	
		// run simulation...
		return sim->Run();
		break;
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
	// Syntax is molsim scene.txt t_end delta_t
	// where t_end and delta_t denote a floating point value
	if(argc != 4 && argc != 2 && argc != 3)
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
		else
		{
			LOG4CXX_ERROR(simulationInitializationLogger, ">> error: invalid argument");
			printUsage();
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
		else
		{
			LOG4CXX_ERROR(simulationInitializationLogger, ">> error: invalid argument");
			printUsage();
			return E_INVALIDPARAM;
		}
	}
	else
	{
		// check if endtime, delta are numbers and file exists
		if(!fileExists(argsv[1]))
		{
			LOG4CXX_ERROR(simulationInitializationLogger, "error: file doesn't exist!");
			printUsage();
			return E_FILENOTFOUND;
		}

		if(!strIsNumber(argsv[2]))
		{
			LOG4CXX_ERROR(simulationInitializationLogger, "error: endtime not a valid number");
			printUsage();
			return E_INVALIDPARAM;
		}

		if(!strIsNumber(argsv[3]))
		{
			LOG4CXX_ERROR(simulationInitializationLogger, "error: delta_t not a valid number");
			printUsage();
			return E_INVALIDPARAM;
		}

		//run Simulation
		state = AS_SIMULATION;
	}
	// parse args...


	return S_OK;
}

void MolSim::showHelp()
{
	LOG4CXX_INFO(generalOutputLogger, " >> "<<"options\t\t\tdescription");
	LOG4CXX_INFO(generalOutputLogger, "    "<<" <file> <endtime> <delta_t>"<<"\t"<<"run simulation according to file");
	LOG4CXX_INFO(generalOutputLogger, "    "<<"-help"<<"\t\t\t"<<"show help");
	LOG4CXX_INFO(generalOutputLogger, "    "<<"-test <name>"<<"\t\t"<<"run single test case or leave\n\t\t\t\t<name> blank to run all tests");
	LOG4CXX_INFO(generalOutputLogger, "    "<<"-showtests"<<"\t\t\t"<<"list all avaliable tests by name");
	
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
	LOG4CXX_INFO(generalOutputLogger, "   usage: molsim file endtime delta_t");
	LOG4CXX_INFO(generalOutputLogger, "   usage: molsim -test");
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
	
		if(wasSuccessful)LOG4CXX_INFO(testLogger, " >> test succeeded! ");
		else LOG4CXX_ERROR(testLogger, " >> test failed! ");
	}
}