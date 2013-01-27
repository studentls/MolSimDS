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
	CPPUNIT_TEST(testListParticleContainerIterationPairwise);
	CPPUNIT_TEST(testLinkedCellParticleContainerGetHalo);	
	CPPUNIT_TEST(testLinkedCellParticleContainerGetBoundary);	
	CPPUNIT_TEST(testLinkedCellParticleContainerIndices);
	CPPUNIT_TEST(testLinkedCellOpenMPIndices);
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

	/// create simple LinkedCellContainer with a cuboid in it
	LinkedCellParticleContainer* createSimpleLC(const unsigned int dimx, const unsigned int dimy, const unsigned int dimz)
	{
		ListParticleContainer lpc;
		utils::Vector<double, 3> vel;
		utils::Vector<double, 3> extent;
		utils::Vector<unsigned int, 3> N;
		utils::Vector<double, 3> lowercorner;
		ParticleGenerator::makeCuboid(lpc, lowercorner, N, 1.0, 1.0, vel, 0.0, 1);

		int dim = dimz == 0 ? 2 : 3;
		lowercorner[0] = 0.0;
		lowercorner[1] = 0.0;
		lowercorner[2] = 0.0; 

		N[0] = dimx;
		N[1] = dimy;
		N[2] = dimz;

		for(int i = 0; i < dim; i++)extent[i] = N[i];

		//construct container...
		LinkedCellParticleContainer *pc = new LinkedCellParticleContainer(dim, lpc.getParticles(), 1.0, lowercorner, extent, BC_NONE, 1.0);

		return pc;
	}

public:
	void setUp()	{}

	void tearDown()	{}

	/// test iteration function by a simple sum function
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

	/// test pairwise iteration function by a simple sum function
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

	/// test if halo particles are retunred properly for constructed linked cell particle container
	void testLinkedCellParticleContainerGetHalo()
	{
		
		utils::Vector<unsigned int, 3> N;

		// set as you want, it should work fine :)
		int n = 10;

		

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

	/// test if boundary particles are retunred properly for constructed linked cell particle container
	void testLinkedCellParticleContainerGetBoundary()
	{
		
		utils::Vector<unsigned int, 3> N;

		// set as you want, it should work fine :)
		int n = 10;

		

		// construct particles
		std::vector<Particle> particles;
		utils::Vector<double, 3> corner;
		utils::Vector<double, 3> lowercorner;
		utils::Vector<double, 3> extent;
		utils::Vector<double, 3> vel;

		// for 2D and 3D (i = 0 => 2D, i = 1 => 3D)
		for(int i = 1; i < 2; i++)
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
			std::vector<Particle> boundary = pc->getBoundaryParticles();

			CPPUNIT_ASSERT(!boundary.empty());

			int sum = 0;
			// perform sum check
			for(std::vector<Particle>::iterator it = boundary.begin(); it != boundary.end(); it++)
			{
				sum += it->type;
			}

			SAFE_DELETE(pc);

			// number of boundary particles should be (N_0 - 2) * (N_1 - 2) - (N_0 - 4) * (N_1 - 4)
			// note that it is only ensured, that N_i >= 3
			// so special case has to be regarded
			// analog in 3D
			if( i == 0)
			{
				
				int calcsum = (N[0] - 2)*(N[1] - 2);			
				if(N[0] > 3 && N[1] > 3)calcsum -= (N[0] - 4)*(N[1] - 4);

				CPPUNIT_ASSERT(sum == calcsum);		
			}
			else
			{
				int calcsum = (N[0] - 2)*(N[1] - 2)*(N[2] - 2);	
				if(N[0] > 3 && N[1] > 3 && N[2] > 3)calcsum -= (N[0] - 4)*(N[1] - 4)*(N[2] - 4);
				CPPUNIT_ASSERT(sum == calcsum);		
			}
		}

	}

	/// test if indices are constructed correctly for linked cell particle container
	void testLinkedCellParticleContainerIndices()
	{
		using namespace utils;
		int index = 0;
		int dimx = 20;
		int dimy = 20;
		int dimz = 20;

		// create a linkedcell particle container
		LinkedCellParticleContainer *pc2D = this->createSimpleLC(dimx, dimy, 0);
		LinkedCellParticleContainer *pc3D = this->createSimpleLC(dimx, dimy, dimz);

		// 2D
		for(int x = 0; x < dimx; x++)
			for(int y = 0; y < dimy; y++)
			{
				index = pc2D->Index2DTo1D(x, y);
				CPPUNIT_ASSERT(pc2D->Index1DTo2D(index)[0] == x);
				CPPUNIT_ASSERT(pc2D->Index1DTo2D(index)[1] == y);
			}

		// 3D
		for(int x = 0; x < dimx; x++)
			for(int y = 0; y < dimy; y++)
				for(int z = 0; z < dimz; z++)
				{
					index = pc3D->Index3DTo1D(x, y, z);
					CPPUNIT_ASSERT(pc3D->Index1DTo3D(index)[0] == x);
					CPPUNIT_ASSERT(pc3D->Index1DTo3D(index)[1] == y);
					CPPUNIT_ASSERT(pc3D->Index1DTo3D(index)[2] == z);
				}		

		SAFE_DELETE(pc2D);
		SAFE_DELETE(pc3D);
	}

	/// tests if the construction of stripped indices works properly
	void testLinkedCellOpenMPIndices()
	{
		// create a indexstrip
		IndexStrip strip;

		// construct for a domain of 10x7x1 indices for the first strip
		strip.constructVerticalStripIndices(0, 7, utils::Vector<unsigned int, 3>(10, 7, 1));

		// there should be 4(n-1) + n pairs, where n denotes cellCount[1] for one strip
		CPPUNIT_ASSERT(strip.size() == 4*(7 - 1) + 7);

		// now test for several domain sizes in 2D
		for(int nx = 1; nx < 5; nx++)
			for(int ny = 1; ny < 5; ny++)
			{
				// construct every time a linked cell particle container
				LinkedCellParticleContainer *pc2D = this->createSimpleLC(nx, ny, 0);

				// halo cells shall also be noticed!
				int dimx = nx + 2;
				int dimy = ny + 2;

				// compare sizes
				// every strip should have 4(dimy - 1) +  dimy pairs contained
				// and oddstrips + evenstrips should equal dimx - 1
				CPPUNIT_ASSERT(pc2D->evenStrips.size() + pc2D->oddStrips.size() == dimx - 1);

				for(int i = 0; i < pc2D->evenStrips.size(); i++)
					CPPUNIT_ASSERT(pc2D->evenStrips[i].size() == 4 * (dimy - 1) + dimy);

				for(int i = 0; i < pc2D->oddStrips.size(); i++)
					CPPUNIT_ASSERT(pc2D->oddStrips[i].size() == 4 * (dimy - 1) + dimy);

				SAFE_DELETE(pc2D);
			}

		// maybe include here a correctness check...
		// ...
		LOG4CXX_INFO(generalOutputLogger, "add correctness check for strip Indices");
	}
};


#endif