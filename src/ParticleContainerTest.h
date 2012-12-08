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
#include "LinkedCellParticleContainer.h"
#include "ParticleGenerator.h"

// see http://cppunit.sourceforge.net/doc/lastest/cppunit_cookbook.html for details

/// Particle Container test class
class ParticleContainerTest : public CppUnit::TestFixture
{	
	//fast setup macros
	CPPUNIT_TEST_SUITE(ParticleContainerTest);
	CPPUNIT_TEST(testListParticleContainerIteration);
	CPPUNIT_TEST(testListParticleContainerIterationPairwise);/*
	CPPUNIT_TEST(testLinkedCellParticleContainerGetHalo2D);	
	CPPUNIT_TEST(testLinkedCellParticleContainerGetHalo3D);*/
	CPPUNIT_TEST(testLinkedCellParticleContainerGetHaloBig);	
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

	/// add to particle container from array
	/// positions has to be array of 3*count elements!
	void addSimpleParticlesToParticleContainerFromArray(std::vector<Particle> &particles, double *positions, const int count, const double mass, const int type)
	{
		Particle p;
		p.m = mass;
		p.type = type;
		for(int i = 0; i < count; i++)
		{
			p.x[0] = positions[3 * i + 0];
			p.x[1] = positions[3 * i + 1];
			p.x[2] = positions[3 * i + 2];
			
			particles.push_back(p);
		}
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

	void testLinkedCellParticleContainerGetHalo2D()
	{
		// 2D test

		// construct particles
		std::vector<Particle> particles;
		utils::Vector<double, 3> corner;
		utils::Vector<double, 3> extent;

		//set 1x1 simulation area
		for(int i = 0; i < 2; i++)
		{
			corner[i] = 0;
			extent[i] = 1;
		}

		// place 8 halo particles, in each of the 8 possible halo regions!
		
		double positions[] = {-1, -1, 0,
							  -0.5, -4, 0,
							  5, -2, 0,
							  3, 0.8, 0,
							  25, 13, 0,
							  0.3, 15, 0,
							  -1, 6, 0,
							  -2, 0.6, 0};

		
		// add to array
		addSimpleParticlesToParticleContainerFromArray(particles, positions, 8, 1.0, 1);

		LinkedCellParticleContainer *pc = new LinkedCellParticleContainer(2, particles, 1, corner, extent, BC_NONE, 1.0);

		// get halo particles
		std::vector<Particle> halo = pc->getHaloParticles();

		CPPUNIT_ASSERT(!halo.empty());

		int sum = 0;
		// perform sum check
		for(std::vector<Particle>::iterator it = halo.begin(); it != halo.end(); it++)
		{
			sum += it->type;
		}

		SAFE_DELETE(pc);

		// there shall be 8 particles contained
		CPPUNIT_ASSERT(sum == 8);
	}

	void testLinkedCellParticleContainerGetHalo3D()
	{
		// 3D test

		// construct particles
		std::vector<Particle> particles;
		utils::Vector<double, 3> corner;
		utils::Vector<double, 3> extent;

		//set 1x1x1 simulation area
		for(int i = 0; i < 3; i++)
		{
			corner[i] = 0;
			extent[i] = 1;
		}

		// place 8 halo particles, in each of the 8 possible halo regions!
		
		double positions[] = {-1, -1, 0.3,
							  -0.5, -4, 0.3,
							  5, -2, 0.4,
							  3, 0.8, 0.2,
							  25, 13, 0.5,
							  0.3, 15, 0.6,
							  -1, 6, 0.1,
							  -2, 0.6, 0.3,
							   0.5, 0.3, 0.2};

		// add to array
		addSimpleParticlesToParticleContainerFromArray(particles, positions, 8, 1.0, 1);

		// lower plane
		for(int i = 0; i < 9; i++)positions[3*i + 2] -= 12;	// z positions shift
		
		addSimpleParticlesToParticleContainerFromArray(particles, positions, 9, 1.0, 1);

		// upper plane
		for(int i = 0; i < 9; i++)positions[3*i + 2] += 38; // z positions shift
		
		addSimpleParticlesToParticleContainerFromArray(particles, positions, 9, 1.0, 1);


		//construct container...
		LinkedCellParticleContainer *pc = new LinkedCellParticleContainer(2, particles, 1, corner, extent, BC_NONE, 1.0);

		// get halo particles
		std::vector<Particle> halo = pc->getHaloParticles();

		CPPUNIT_ASSERT(!halo.empty());

		int sum = 0;
		// perform sum check
		for(std::vector<Particle>::iterator it = halo.begin(); it != halo.end(); it++)
		{
			sum += it->type;
		}

		SAFE_DELETE(pc);

		// there shall be 26 ( = 9 + 8 + 9 particles contained
		CPPUNIT_ASSERT(sum == 26);		
	}

	void testLinkedCellParticleContainerGetHaloBig()
	{
		
		utils::Vector<unsigned int, 3> N;

		int n = 3;

		

		// construct particles
		std::vector<Particle> particles;
		utils::Vector<double, 3> corner;
		utils::Vector<double, 3> lowercorner;
		utils::Vector<double, 3> extent;
		utils::Vector<double, 3> vel;

		// for 2D and 3D (i = 0 => 2D, i = 1 => 3D)
		for(int i = 0; i < 2; i++)
		{
			int dim = 2 + i;

			lowercorner[0] = -0.5;
			lowercorner[1] = -0.5;
			lowercorner[2] = -0.5 * i; // set according to used func
			
			N[0] = n;
			N[1] = n;
			N[2] = 1 + (n-1) * i;
			
			// set simulation area
			// we will place a cuboid over all cells
			for(int i = 0; i < dim; i++)
			{
				corner[i] = 0;
				extent[i] = N[i] - 2;
			}

			ListParticleContainer lpc;
			ParticleGenerator::makeCuboid(lpc, lowercorner, N, 1.0, 1.0, vel, 0.0, 1);

			lowercorner[0] = 0.0;
			lowercorner[1] = 0.0;
			lowercorner[2] = 0.0; 

			//construct container...
			LinkedCellParticleContainer *pc = new LinkedCellParticleContainer(dim, lpc.getParticles(), 1.0, lowercorner, extent, BC_NONE, 1.0);

			// get halo particles
			std::vector<Particle> halo = pc->getHaloParticles();

			CPPUNIT_ASSERT(!halo.empty());

			int sum = 0;
			// perform sum check
			for(std::vector<Particle>::iterator it = halo.begin(); it != halo.end(); it++)
			{
				sum += it->type;
			}

			SAFE_DELETE(pc);

			// number of halo particles should be N_0 * N_1 - (N_0 - 2) * (N_1 - 2)
			// analog in 3D
			if( i == 0)
			{
				
				int calcsum = N[0] * N[1] - (N[0] - 2)*(N[1] - 2);			
				CPPUNIT_ASSERT(sum == calcsum);		
			}
			else
			{
				int calcsum = N[0] * N[1] * N[2] - (N[0] - 2)*(N[1] - 2)*(N[2] - 2);			
				CPPUNIT_ASSERT(sum == calcsum);		
			}
		}
	}
};


#endif