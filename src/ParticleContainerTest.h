//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File ParticleContainerTest.h
// tests a Particle Container
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------


#ifndef PARTICLE_CONTAINER_TEST_H_
#define PARTICLE_CONTAINER_TEST_H_

//CPP Unit
#include <cppunit/TestCase.h>
#include <cppunit/TestFixture.h>
#include <cppunit/ui/text/TextTestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>
#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/XmlOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>


#include "ParticleContainer.h"
#include "ListParticleContainer.h"

// see http://cppunit.sourceforge.net/doc/lastest/cppunit_cookbook.html for details

/// Particle Container test class
class ParticleContainerTest : public CppUnit::TestFixture
{	
	//fast setup macros
	CPPUNIT_TEST_SUITE(ParticleContainerTest);
	CPPUNIT_TEST(testListParticleContainerIteration);
	CPPUNIT_TEST(testListParticleContainerIterationPairwise);
	CPPUNIT_TEST_SUITE_END();
private:

	/// test function, summing up data
	static void funcIterationSum(void *data, Particle& p)
	{
		double *d = ((double*)data);

		//add masses...
		*d += p.m;
	}

	/// test function, summing up data pairwise
	static void funcIterationPWSum(void *data, Particle& p1, Particle& p2)
	{
		double *d = ((double*)data);

		//add masses...
		*d += 1.0;
	}

public:
	void setUp()	{}

	void tearDown()	{}

	void testListParticleContainerIteration()
	{
		ListParticleContainer PC;
		
		//add 10 Particles with position, velocity and mass
		int iParticleCount = 10;
		for(int i = 0; i < iParticleCount; i++)
		{
			Particle p1(Vec3(1.0), Vec3(1.0), 1.0);
			PC.AddParticle(p1);
		}

		//now simple iterate over the Particles and add their mass
		double sum = 0.0;
		PC.Iterate(funcIterationSum, &sum);

		//the final value in sum should be something around iParticleCount 
		//because of floating point arithmetics
		
		double epsilon = 0.00001;

		CPPUNIT_ASSERT( -epsilon < (sum - (double)iParticleCount)  && (sum - (double)iParticleCount) < epsilon);
	}

	void testListParticleContainerIterationPairwise()
	{
		ListParticleContainer PC;

		//add 10 Particles with position, velocity and mass
		int iParticleCount = 10;
		for(int i = 0; i < iParticleCount; i++)
		{
			Particle p1(Vec3(1.0), Vec3(1.0), 1.0);
			PC.AddParticle(p1);
		}

		//now simple iterate over the Particles and add their mass
		double sum = 0.0;
		PC.IteratePairwise(funcIterationPWSum, &sum);

		// Iterate Pairwise shall be called for every unsorted pair once!
		// also pairs like (1,1) are not allowed!
		// so sum should return n * (n - 1 ) / 2
		//because of floating point arithmetics
		
		double epsilon = 0.00001;
		double result = (double)(iParticleCount * (iParticleCount - 1) / 2);
		CPPUNIT_ASSERT( -epsilon < (sum - result)  && (sum - result) < epsilon);
	
	}
};


#endif