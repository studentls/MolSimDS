//------------------------------------------------------------------------------------------------
// File ListParticleContainer.h
// contains class ListParticleContainer
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------


#ifndef PARTICLE_CONTAINER_H_
#define PARTICLE_CONTAINER_H_

#include <vector>
#include "Logging.h"
#include "Particle.h"

/// an abstract class that is used to store Particles and iterate over them
class ParticleContainer {
public:
	/// a method to add a Particle to the ListParticleContainer
	virtual void AddParticle(const Particle& particle) = 0;

	/// a method that takes a void(*func)(void*, Particle) and uses it to iterate over all Particles
	/// @param data additional data given to func
	virtual void Iterate(void(*func)(void*, Particle&), void *data) = 0;
	
	/// a method that takes a void(*func)(void*, Particle, Particle) and uses it to iterate over all pairs of Particles (each symmetrical pair is only taken once to reduce redundancy)
	/// @param data additional data given to func
	virtual void IteratePairwise(void(*func)(void*, Particle&, Particle&), void *data) = 0;

	/// add particles from *.txt file
	virtual void AddParticlesFromFile(const char *filename) = 0;

	/// our new fileformat, replace later AddParticlesFromFile
	/// @return return true if file could be read
	virtual bool AddParticlesFromFileNew(const char *filename) = 0;

	/// removes all particles
	virtual void Clear() = 0;

	/// are any particles contained?
	virtual bool IsEmpty() = 0;
	// this method shall be later removed...
	///returns ListParticleContainer's internal container
	virtual const std::vector<Particle>& getParticles() = 0;

};

#endif /* PARTICLE_CONTAINER_H_ */

