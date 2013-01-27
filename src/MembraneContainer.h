//------------------------------------------------------------------------------------------------
// File LinkedCellParticleContainer.h
// contains class LinkedCelltParticleContainer
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------


#ifndef MEMBRANE_CONTAINER_H_
#define MEMBRANE_CONTAINER_H_

#include <vector>
#include "Logging.h"
#include "Particle.h"
#include "ParticleContainer.h"

/// a class that is used to simulate a membrane
/// this class works differently from the other containers
/// because particles have to be linked to each other and are therefore created within the container
/// instead of by a ParticleGenerator
class MembraneContainer : public ParticleContainer {

private:
	/// the array that is used by the ParticleContainer to store the Particles
	Particle *particles;

	/// the number of particles in the mesh
	int particleCount;

	/// the number of particles in each dimension
	utils::Vector<unsigned int, 2> dimensions;

public:

	unsigned int pullIterations;

	/// a constructor that takes the number of iterations for which the pull is active
	MembraneContainer(unsigned int pullIterations);

	/// a method to add a Particle to the ListParticleContainer
	void AddParticle(const Particle& particle);

	/// a method that takes a void(*func)(void*, Particle) and uses it to iterate over all Particles
	/// @param data additional data given to func
	void Iterate(void(*func)(void*, Particle&), void *data);
	
	/// a method that takes a void(*func)(void*, Particle, Particle) and uses it to iterate over all pairs of Particles (each symmetrical pair is only taken once to reduce redundancy)
	/// @param data additional data given to func
	void IteratePairwise(void(*func)(void*, Particle&, Particle&), void *data);
	
	/// calculate and apply the force between a pair of neighboring particles in a membrane. Used in calculateF()
	void ApplyMembraneForces();
	
	/// calculate and apply the force between a pair of neighboring particles in a membrane. Used in ApplyMembraneForces()
	void ApplyMembraneForce(double distanceFactor, Particle& p1, Particle& p2);

	/// applies a reflective boundary at the bottom
	void ApplyReflectiveBoundaryAtBottom(void(*func)(void*, Particle&, Particle&), void *data);

	/// a method that initializes the mesh the MembraneContainer holds
	void SetMembrane(const Vec3& vLowerLeftFrontCorner,
					   const utils::Vector<unsigned int, 2>& dimensions,
					   const double meshWidth,
					   const double mass);

	/// add particles from *.txt file
	void AddParticlesFromFile(const char *filename);

	/// our new fileformat, replace later AddParticlesFromFile
	/// @return return true if file could be read
	bool AddParticlesFromFileNew(const char *filename);

	/// removes all particles
	void Clear()	{ particles = NULL; particleCount = 0; dimensions = 0;}

	/// are any particles contained?
	bool IsEmpty()				{return particleCount == 0;}

	/// returns how many particles are contained in this container
	/// @return returns count of particles in this container
	unsigned int					getParticleCount()	{return particleCount;}
	
	// Quick 'n' dirty hack
	std::vector<Particle> p;

	// this method shall be later removed...
	///returns MembraneContainer's internal container
	const std::vector<Particle>& getParticles()
	{
		// clear if necessary
		if(!p.empty())p.clear();

		assert(particles);
		
		for (int i = 0; i < particleCount; i++)
			p.push_back(particles[i]);
		return p;
	}

	/// method to identify container
	/// @return returns PCT_LIST
	ParticleContainerType			getType() {return PCT_MEMBRANE;}

	/// method to retrieve a Bounding Box, which surrounds all particles
	/// @return returns a BoundingBox, which defines extent and center of all particles in the container(bounding box)
	utils::BoundingBox				getBoundingBox();
};

#endif 