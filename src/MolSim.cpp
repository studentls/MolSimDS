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


err_type MolSim::Init(int argc, char *argsv[])
{
	// parse command line arguments
	if(FAILED(parseLine(argc, argsv)))return E_INVALIDPARAM;

	// print hello message
	printHelloMessage();

	// Init Simulation Data
	sim = new Simulation();
	if(!sim)
	{
		cout<<"failed to allocate mem"<<endl;
		return E_OUTOFMEMORY;
	}

	// set desc
	SimulationDesc desc;

	desc.output_fmt = SOF_VTK;
	desc.start_time = 0.0;
	desc.end_time = atof(argsv[2]);
	desc.delta_t = atof(argsv[3]);
	
	if(FAILED(sim->Init(desc)))return E_UNKNOWN;

	if(FAILED(sim->AddParticlesFromFile(argsv[1])))return E_INVALIDPARAM;

	return S_OK;
}

err_type MolSim::Run()
{
	if(!sim)return E_NOTINITIALIZED;

	// run simulation...
	return sim->Run();
}

err_type MolSim::Release()
{
	DELETE(sim);

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
	if(argc != 4)
	{
		cout<<"error: invalid count of arguments"<<endl;
		cout<<"usage: molsim file endtime delta_t"<<endl;
		return E_INVALIDPARAM;
	}

	// check if endtime, delta are numbers and file exists
	if(!fileExists(argsv[1]))
	{
		cout<<"error: file doesn't exist!"<<endl;
		cout<<"usage: molsim file endtime delta_t"<<endl;
		return E_FILENOTFOUND;
	}

	if(!strIsNumber(argsv[2]))
	{
		cout<<"error: endtime not a valid number"<<endl;
		cout<<"usage: molsim file endtime delta_t"<<endl;
		return E_INVALIDPARAM;
	}

	if(!strIsNumber(argsv[3]))
	{
		cout<<"error: delta_t not a valid number"<<endl;
		cout<<"usage: molsim file endtime delta_t"<<endl;
		return E_INVALIDPARAM;
	}

	// parse args...


	return S_OK;
}


void MolSim::printHelloMessage()
{
	line();
	cout << "MolSim for PSE" << endl;
#ifdef DEBUG
	cout << "compiled: "<<__DATE__<<"  "<<__TIME__<<endl;
#endif

	line();
	cout << endl;
	cout << "(c) 2012 by F.Dietz & L.Spiegelberg" << endl; 
	cout << endl;
	cout << "Molecular Simulator handling *.txt files" << endl;
	line();
	cout << endl;
}

/// console tool
void line()
{
	for(int i = 0; i < 40; i++)cout<<"-";
	cout<<endl;
}

