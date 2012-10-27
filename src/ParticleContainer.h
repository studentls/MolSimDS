//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File ParticleContainer.h
// contains class ParticleContainer
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------


#ifndef PARTICLE_CONTAINER_H_
#define PARTICLE_CONTAINER_H_

#include <vector>
#include "Particle.h"

/// a class that is used to store Particles and iterate over them
class ParticleContainer {

private:
	/// the vector<Particles> that is used by the ParticleContainer to store the Particles
	std::vector<Particle> particles;

public:
	/// a constructor that creates the value particles from scratch
	ParticleContainer();

	/// a constructor that takes one argument, which is the vector<Particles>
	ParticleContainer(std::vector<Particle> particles);

	/// a function to add a Particle to the ParticleContainer
	void AddParticle(Particle particle);

	/// a function that takes a void(*func)(Particle) and uses it to iterate over all Particles
	void Iterate(void(*func)(Particle));
	
	/// a function that takes a void(*func)(Particle, Particle) and uses them to iterate over all pairs of Particles
	/// each symmetrical pair is only taken once to reduce redundancy
	void IteratePairwise(void(*func)(Particle, Particle));

};

#endif /* PARTICLE_CONTAINER_H_ */