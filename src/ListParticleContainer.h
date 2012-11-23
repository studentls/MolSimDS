//------------------------------------------------------------------------------------------------
// File ListParticleContainer.h
// contains class ListParticleContainer
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------


#ifndef LIST_PARTICLE_CONTAINER_H_
#define LIST_PARTICLE_CONTAINER_H_

#include <vector>
#include "Logging.h"
#include "Particle.h"
#include "ParticleContainer.h"

/// a class that is used to store Particles and iterate over them
class ListParticleContainer : public ParticleContainer {

private:
	/// the vector<Particles> that is used by the ParticleContainer to store the Particles
	std::vector<Particle> particles;

public:
	/// a constructor that creates the value particles from scratch
	ListParticleContainer();

	/// a constructor that takes one argument, which is the vector<Particles>
	ListParticleContainer(const std::vector<Particle>& particles);

	/// a method to add a Particle to the ListParticleContainer
	void AddParticle(const Particle& particle);

	/// a method that takes a void(*func)(void*, Particle) and uses it to iterate over all Particles
	/// @param data additional data given to func
	void Iterate(void(*func)(void*, Particle&), void *data);
	
	/// a method that takes a void(*func)(void*, Particle, Particle) and uses it to iterate over all pairs of Particles (each symmetrical pair is only taken once to reduce redundancy)
	/// @param data additional data given to func
	void IteratePairwise(void(*func)(void*, Particle&, Particle&), void *data);

	/// add particles from *.txt file
	void AddParticlesFromFile(const char *filename);

	/// our new fileformat, replace later AddParticlesFromFile
	/// @return return true if file could be read
	bool AddParticlesFromFileNew(const char *filename);

	/// removes all particles
	void Clear()	{if(!particles.empty())particles.clear();}

	/// are any particles contained?
	bool IsEmpty()				{return particles.empty();}

	/// returns how many particles are contained in this container
	/// @return returns count of particles in this container
	unsigned int					getParticleCount()	{return particles.size();}

	// this method shall be later removed...
	///returns ListParticleContainer's internal container
	const std::vector<Particle>& getParticles()	{return particles;}

};

#endif 

