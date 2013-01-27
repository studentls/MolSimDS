//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File PerformanceTest.cpp
// contains class PerformanceTest
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------
#include "PerformanceTest.h"


err_type PerformanceTest::Run(const char *xmlFile)
{
	err_type res = S_OK;
	int perfstepcount = 1000; // use 1000 steps to analyse Performance
	int maximumthreads = 16; // perf up to 16 threads

	LOG4CXX_INFO(generalOutputLogger, ">> init performance tests for file "<<xmlFile);

	// create and init Simulation
	Simulation sim;

	res = sim.CreateSimulationFromXMLFile(xmlFile);
	if(FAILED(res))return res;

	// modify simulation desc
	sim.desc.output_fmt = SOF_NONE;			// no output for performance analysis
	sim.desc.iterationsperoutput = 9999;	// no output

	sim.desc.end_time = sim.desc.start_time +  sim.desc.delta_t * perfstepcount; // use perfsteps
	
	double speedup = 1.0;
	double basespeed = 1.0;

	LOG4CXX_INFO(generalOutputLogger, ">> starting performance tests for file "<<xmlFile);
	// performe analysis
	for(int i = 1; i <= maximumthreads; i++)
	{
		sim.particles->setThreadCount(i);
		LOG4CXX_INFO(generalOutputLogger, ">> running performance test "<<i<<" out of "<<maximumthreads);
		utils::Timer timer;
		sim.Run();
		double time = timer.getElapsedTime();
		if(i == 1)basespeed = time; // store for one thread

		// calc speed up
		speedup = basespeed / time;

		LOG4CXX_INFO(generalOutputLogger, ">> "<<time << " s @ "<<i<<" threads , speedup factor: "<<speedup);

		// push back data
		PerformanceTestData entry;
		entry.time = time;
		entry.iterations = perfstepcount;
		entry.num_threads = i;
		entry.speedup = speedup;

		data.push_back(entry);
	}


	// save data
	std::string filename = xmlFile;
	size_t pos = filename.rfind(".");
	if(pos != std::string::npos)filename = filename.substr(0, pos) + ".csv";
	res = saveCSV(filename.c_str());

	return res;
}

err_type PerformanceTest::saveCSV(const char *filename)
{
	FILE *pFile = NULL;

	// open file (ASCII, replace)
	pFile = fopen(filename, "w"); 
	if(!pFile)return E_UNKNOWN;

	// print header
	fprintf(pFile, "number of threads, iterations, time, speedup\n");

	// print data
	for(int i = 0; i < data.size(); i++)
	{
		fprintf(pFile, "%d, %d, %lf, %lf\n", data[i].num_threads, data[i].iterations, data[i].time, data[i].speedup);
	}


	fclose(pFile);
}