//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File ParticleContainer.cpp
// contains class ParticleContainer
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------

#include "ParticleContainer.h"

ParticleContainer::ParticleContainer() {
	particles = std::vector<Particle>();
}

ParticleContainer::ParticleContainer(std::vector<Particle> particles_args) {
	particles = particles_args;
}

void ParticleContainer::AddParticle(Particle& particle) {
	particles.push_back(particle);
}

void ParticleContainer::Iterate(void(*func)(Particle&)) {
	/// iterate over all Particles and call the function on it
	for (std::vector<Particle>::iterator it = particles.begin() ; it < particles.end(); it++)
	{
		Particle& p = *it;
		(*func)(p);
	}
}

void ParticleContainer::IteratePairwise(void(*func)(Particle&, Particle&)) {
	/// iterate over all Particles
	for (std::vector<Particle>::iterator it1 = particles.begin() ; it1 < particles.end(); it1++)
		for (std::vector<Particle>::iterator it2 = it1 + 1; it2 < particles.end(); it2++)
			/// make sure that a Particle is not paired with itself
			if (it1 != it2)
			{
				/// call the function on the pair of Particles
				Particle& p1 = *it1;
				Particle& p2 = *it2;
				(*func)(p1, p2);
			}
}