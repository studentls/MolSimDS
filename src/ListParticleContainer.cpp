//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File ListParticleContainer.cpp
// contains class ListParticleContainer
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------

#include <cstdio>
#include <cstring>

#include "ListParticleContainer.h"
#include "FileReader.h"
#include "ParticleGenerator.h"

// gcc has currently problems, to compile sources with indirect includes for STL!
// so to fix this, use direct include 
// and include old c libraries in C++ style, as internet sources state
// in the current gcc version seems to reorganize files

ListParticleContainer::ListParticleContainer() {
	particles = std::vector<Particle>();
}

ListParticleContainer::ListParticleContainer(const std::vector<Particle>& particles_args) {
	// copy data
	particles = particles_args;
}

void ListParticleContainer::AddParticle(const Particle& particle) {
	particles.push_back(particle);
}

void ListParticleContainer::Iterate(void(*func)(void*, Particle&), void *data) {
	// iterate over all Particles and call the function on it
	for (std::vector<Particle>::iterator it = particles.begin() ; it < particles.end(); it++)
	{
		Particle& p = *it;
		(*func)(data, p);
	}
}

void ListParticleContainer::IteratePairwise(void(*func)(void *data, Particle&, Particle&), void *data) {
	// iterate over all Particles
	for (std::vector<Particle>::iterator it1 = particles.begin() ; it1 < particles.end(); it1++)
		for (std::vector<Particle>::iterator it2 = it1 + 1; it2 < particles.end(); it2++)
			// make sure that a Particle is not paired with itself
			if (it1 != it2)
			{
				// call the function on the pair of Particles
				Particle& p1 = *it1;
				Particle& p2 = *it2;
				(*func)(data, p1, p2);
			}
}

utils::BoundingBox ListParticleContainer::getBoundingBox()
{
	using namespace utils;

	// some work has to be done here,
	// calc bounding box
	BoundingBox bb;

	Vector<double, 3> vmin(999999.9f);	//+infty
	Vector<double, 3> vmax(-9999999.9f); //-infty

	if(!particles.empty())
	{
		for(std::vector<Particle>::iterator it = particles.begin(); it != particles.end(); it++)
		{
			for(int i = 0; i < 3; i++)
			{
				vmin[i] = min(vmin[i], it->x[i]);
				vmax[i] = max(vmax[i], it->x[i]);
			}			
		}
	}

	// calc area(bounding box)
	bb.extent = vmax - vmin;
	bb.center = (vmax + vmin) * 0.5;

	return bb;
}
