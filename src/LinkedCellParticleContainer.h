//------------------------------------------------------------------------------------------------
// File LinkedCellParticleContainer.h
// contains class LinkedCelltParticleContainer
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------


#ifndef LINKEDCELL_PARTICLE_CONTAINER_H_
#define LINKEDCELL_PARTICLE_CONTAINER_H_

#include <vector>
#include "Logging.h"
#include "Particle.h"
#include "ParticleContainer.h"

/// a class that is used to store Particles and iterate over them
template<typename T> class LinkedCellParticleContainer : public ParticleContainer {

private:

	// ...

public:
	
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
	void Clear();

	/// are any particles contained?
	bool IsEmpty();

	/// returns how many particles are contained in this container
	/// @return returns count of particles in this container
	unsigned int					getParticleCount();

	// this method shall be later removed...
	///returns ListParticleContainer's internal container
	const std::vector<Particle>& getParticles();

	/// included for use in XMLFileReader...
	LinkedCellParticleContainer(const unsigned int dim,
		const std::vector<Particle>& particles,
		const double cutoffDistance,
		utils::Vector<double, 3> frontLowerLeftCorner,
		utils::Vector<double, 3> simulationAreaExtent,
		int iterationsPerParticleToCellReassignment,
		bool leftReflectiveBoundary, bool rightReflectiveBoundary,
		bool frontReflectiveBoundary, bool backReflectiveBoundary, // these two will be ignored in the two-dimensional case		
		bool bottomReflectiveBoundary, bool topReflectiveBoundary,			
		double sigma)	
	{
	}

};


// implement here...

template<typename T> void LinkedCellParticleContainer<T>::AddParticle(const Particle& particle)
{
	// functions have to look like this...
	// first template<typename T> ... LinkedCellParticleContainer<T>::...
}

#endif 

