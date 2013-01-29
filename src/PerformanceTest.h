//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File PerformanceTest.h
// contains class PerformanceTest
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------
#ifndef PERFORMANCE_TEST_HEADER_
#define PERFORMANCE_TEST_HEADER_

#include "Simulation.h"
#include "utils/utils.h"

/// struct to hold measured data for one performance test
struct PerformanceTestData
{
	int		num_threads;	/// how many threads made this data?
	int		iterations;		/// how many iterations have been used?
	double	time;			/// how much time was needed?
	double  speedup;		/// speed up factor
};

/// class to measure automatically Multithread Performance
class PerformanceTest
{
private:
	utils::TFastArray<PerformanceTestData> data;

	/// helper function to save accumulated data to a Comma Separated Value file
	err_type	saveCSV(const char *filename);	
public:

	/// Run Performance Tests, internal the maximum steps will be cut at 1000 steps to ensure a good runtime
	/// @param inputFile for which performance shall be evaluated
	err_type Run(const char *xmlFile);
};

#endif