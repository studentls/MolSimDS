/*
 * ParticleContainer.h
 *
 *  Created on: 26.10.2012
 *      Authors: Florian Dietz, Leonhard Spiegelberg
 */

#include "ParticleContainer.h"

ParticleContainer() {
	particles = new util::vector<Particle>();
};

ParticleContainer(utils::Vector<Particle> particles) {
	this.particles = particles;
};

void AddParticle(Particle particle) {
	particles.push_back(particle);
};

void Iterate(void(*func)(Particle)) {
	for (utils::Vector<Particle>::iterator it=particles.begin() ; it < particles.end(); it++)
	{
		Particle p = *it;
		(*func)(p);
	}
};

void IteratePairwise(void(*func)(Particle, Particle)) {
	for (utils::Vector<Particle>::iterator it1=particles.begin() ; it1 < particles.end(); it1++)
		for (utils::Vector<Particle>::iterator it2=particles.begin() ; it2 < particles.end(); it2++)
			if (it1 != it2)
			{
				Particle p1 = *it1;
				Particle p2 = *it2;
				(*func)(p1, p2);
			}
};