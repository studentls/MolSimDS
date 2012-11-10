//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File main.cpp
// contains main function
//------------------------------------------------------------------------------------------------


#include "MolSim.h"

// entry point
int main(int argc, char* argsv[]) {

	MolSim *molsim = new MolSim();

	if(FAILED(molsim->Init(argc, argsv)))
	{
		std::cout<<">> Initialization of Molecular Simulator failed"<<std::endl;
		std::cout<<">> quitting program..."<<std::endl;
		return 0;
	}

	molsim->Run();

	molsim->Release();

	SAFE_DELETE(molsim);

	return 0;
}
