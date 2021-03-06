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
#include "utils/utils.h"

/// used to identify ParticleContainer
enum ParticleContainerType
{
	PCT_UNKOWN,
	PCT_LIST,
	PCT_LINKEDCELL,
	PCT_MEMBRANE
};

/// an abstract class that is used to store Particles and iterate over them
class ParticleContainer {
protected:
	int										threadCount;
public:

	/// constructor
	ParticleContainer()
	{
#ifdef OPENMP
		threadCount = omp_get_max_threads(); // use per default maximum avaliable threads
#else
		threadCount = 1;
#endif
	}

	/// a method to add a Particle to the ListParticleContainer
	virtual void							AddParticle(const Particle& particle) = 0;

	/// a method that takes a void(*func)(void*, Particle) and uses it to iterate over all Particles
	/// @param data additional data given to func
	virtual void							Iterate(void(*func)(void*, Particle&), void *data) = 0;
	
	/// a method that takes a void(*func)(void*, Particle, Particle) and uses it to iterate over all pairs of Particles (each symmetrical pair is only taken once to reduce redundancy)
	/// @param data additional data given to func
	virtual void							IteratePairwise(void(*func)(void*, Particle&, Particle&), void *data) = 0;

	/// removes all particles
	virtual void							Clear() = 0;

	/// are any particles contained?
	/// return true if at least one particle is contained in this container, false otherwise
	virtual bool							IsEmpty() = 0;

	/// returns how many particles are contained in this container
	/// @return returns count of particles in this container
	virtual unsigned int					getParticleCount() = 0;

	/// this method shall be later removed...
	/// @return returns ListParticleContainer's internal container
	virtual const std::vector<Particle>&	getParticles() = 0;

	/// method to identify container
	/// @return returns member of ParticleContainerType
	virtual ParticleContainerType			getType() = 0;

	/// method to retrieve a Bounding Box, which surrounds all particles
	/// @return returns a BoundingBox, which defines extent and center of all particles in the container(bounding box)
	virtual utils::BoundingBox				getBoundingBox() = 0;

	/// set how many threads shall be used
	/// @param threadCount minimum 1
	void									setThreadCount(const int threadCount)	{this->threadCount = threadCount < 1 ? 1 : threadCount;}
};

#endif /* PARTICLE_CONTAINER_H_ */

