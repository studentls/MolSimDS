/*
 * ParticleContainer.h
 *
 *  Created on: 26.10.2012
 *      Authors: Florian Dietz, Leonhard Spiegelberg
 */


#ifndef PARTICLE_CONTAINER_H_
#define PARTICLE_CONTAINER_H_

#include "utils/Vector.h"
#include "Particle.h"


class ParticleContainer {

private:
	utils::Vector<Particle> particles;

public:
	ParticleContainer();

	ParticleContainer(util::vector<Particle> particles);

	void AddParticle(Particle particle);

	void Iterate(void(*func)(Particle));
	
	void IteratePairwise(void(*func)(Particle, Particle));

};

#endif /* PARTICLE_CONTAINER_H_ */