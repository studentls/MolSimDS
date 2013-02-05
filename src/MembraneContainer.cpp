//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File LinkedCellParticleContainer.cpp
// contains class LinkedCellParticleContainer
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------

#include <cstdio>
#include <cstring>

#include "MembraneContainer.h"

MembraneContainer::MembraneContainer(unsigned int pullIterations) {
	particles = NULL;
	particleCount = 0;
	dimensions = 0;
	this->pullIterations = pullIterations;
}

void MembraneContainer::AddParticle(const Particle& particle) {
	throw "AddParticle() does not work for MembraneContainers, since they work differently from other containers";
}

void MembraneContainer::Iterate(void(*func)(void*, Particle&), void *data) {
	// iterate over all Particles and call the function on it
	for (int i = 0; i < particleCount; i++)
	{
		Particle& p = particles[i];
		(*func)(data, p);
	}
}

void MembraneContainer::IteratePairwise(void(*func)(void *data, Particle&, Particle&), void *data) {
	// iterate over all Particles
	for (int i = 0; i < particleCount; i++)
		for (int j = i + 1; j < particleCount; j++)
		{
			Particle& p1 = particles[i];
			Particle& p2 = particles[j];
			// check that the particles are close enough
			double dist = p1.x.distance(p2.x);
			if (abs(dist) < 1.225)
				// call the function on the pair of Particles
				(*func)(data, p1, p2);
		}
}

void MembraneContainer::ApplyMembraneForces()
{
	// iterate over all Particles
	for (int i = 0; i < particleCount; i++)
	{
		Particle& p1 = particles[i];
		// do the direct neighbors
		std::vector<int> neighbors = p1.getDirectNeighbors();
		for (std::vector<int>::iterator it = neighbors.begin(); it < neighbors.end(); it++)
		{
			int j = *it;
			Particle& p2 = particles[j];
			// call the function on the pair of Particles
			ApplyMembraneForce(1.0, p1, p2);
		}
		// do the diagonal neighbors
		neighbors = p1.getDiagonalNeighbors();
		double distanceFactor = sqrt(2.0);
		for (std::vector<int>::iterator it = neighbors.begin(); it < neighbors.end(); it++)
		{
			int j = *it;
			Particle& p2 = particles[j];
			// call the function on the pair of Particles
			ApplyMembraneForce(distanceFactor, p1, p2);
		}
	}
}

void MembraneContainer::ApplyMembraneForce(double distanceFactor, Particle& p1, Particle& p2)
{
	double k = 300.0;
	double r0 = 2.2;

	double r = r0 * distanceFactor;
	double dst = p1.x.distance(p2.x);
	
	double factor = k * (dst - r);

	// the division by the distance is for normalization of the force vector
	utils::Vector<double, 3> force = (p2.x - p1.x) * (factor / dst);
	
	// add individual particle to particle force to sum
	p1.addForce(force);
	
	// new function to avoid unnecessary object construction
	p2.substractForce(force);
}

void MembraneContainer::ApplyReflectiveBoundaryAtBottom(void(*func)(void*, Particle&, Particle&), void *data)
{
	for (int i = 0; i < particleCount; i++)
	{
		Particle& p = particles[i];
		if (p.x[2] < 1.1225)
		{
			Particle p2(p);
			p2.x[2] = 0.0;
			(*func)(data, p, p2);
		}
	}
}

void MembraneContainer::SetMembrane(const Vec3& vLowerLeftFrontCorner,
					   const utils::Vector<unsigned int, 2>& dimensions,
					   const double meshWidth)
{
	// initialize particles
	particleCount = dimensions[0] * dimensions[1];
	this->dimensions = dimensions;
	particles = new Particle[particleCount];

	// acknowledge that default constructor sets forces and velocities to zero...
	Particle p;
	// go through mesh grid...
	for(unsigned int x = 0; x < dimensions[0]; x++)
		for(unsigned int y = 0; y < dimensions[1]; y++)
		{
			p.x[0]	= vLowerLeftFrontCorner[0] + x * meshWidth;
			p.x[1]	= vLowerLeftFrontCorner[1] + y * meshWidth;
			p.x[2]	= vLowerLeftFrontCorner[2];
			// the type is used to determine if the particle is to be pulled up in the beginning
			if ((x == 16 || x == 17) &&
				(y == 23 || y == 24))
				p.type = 1;
			/*if (x == 5 && y == 5)
				p.type = 1;*/
			else
				p.type	= 0;
			
			int thisIndex = x + y * dimensions[0];
			// store all neighbors with a lower index
			std::vector<int> directN = std::vector<int>();
			std::vector<int> diagonalN = std::vector<int>();
			for (int xp = -1; xp < 2; xp++)
				for (int yp = -1; yp < 2; yp++)
				{
					if (xp == 0 && yp == 0)
						continue;
					int xt = x + xp;
					int yt = y + yp;
					if (xt < 0 || xt >= dimensions[0])
						continue;
					if (yt < 0 || yt >= dimensions[1])
						continue;
					int i = xt + yt * dimensions[0];
					if (i < thisIndex)
						if (xp == 0 || yp == 0)
							directN.push_back(i);
						else
							diagonalN.push_back(i);
				}
			p.setDirectNeighbors(directN);
			p.setDiagonalNeighbors(diagonalN);

			// add it
			particles[thisIndex] = p;
		}
}

void MembraneContainer::AddParticlesFromFile(const char *filename)
{
	throw "AddParticlesFromFile() does not work for MembraneContainers, since they work differently from other containers";
}

/// Quick 'n' Dirty Method
bool MembraneContainer::AddParticlesFromFileNew(const char *filename)
{
	throw "AddParticlesFromFileNew() does not work for MembraneContainers, since they work differently from other containers";
}


utils::BoundingBox MembraneContainer::getBoundingBox()
{
	using namespace utils;

	// some work has to be done here,
	// calc bounding box
	BoundingBox bb;

	Vector<double, 3> vmin(999999.9f);	//+infty
	Vector<double, 3> vmax(-9999999.9f); //-infty

	for(int index = 0; index < particleCount; index++)
	{
		Particle& it = particles[index];
		for(int i = 0; i < 3; i++)
		{
			vmin[i] = min(vmin[i], it.x[i]);
			vmax[i] = max(vmax[i], it.x[i]);
		}			
	}

	// calc area(bounding box)
	bb.extent = vmax - vmin;
	bb.center = (vmax + vmin) * 0.5;

	return bb;
}
