//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File Simulation.cpp
// contains class Simulation
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------

#ifndef PARTICLECUBOID_H_
#define PARTICLECUBOID_H_

#include "Particle.h"
#include "ParticleContainer.h"
#include "MaxwellBoltzmannDistribution.h"

/// adds a cuboid of particles to a ParticleContainer
/// @param cont the ParticleContainer to add to
/// @param lowerLeftFrontCorner the lower left front corner of the cuboid
/// @param n1 the number of particles on the x-axis
/// @param n2 the number of particles on the y-axis
/// @param n3 the number of particles on the z-axis
/// @param particleDistance the distance between neighbouring particles
/// @param mass the mass of each particle
/// @param initialParticleVelocity the initial velocity of particles
/// @param desc the simulation description containing relevant parameters, in particular the brownianMotionFactor
/// @param type the type of Particle to create, which is currently an unused value
void AddParticleCuboid(ParticleContainer& cont,
	utils::Vector<double,3>& lowerLeftFrontCorner,
	int n1,
	int n2,
	int n3,
	double particleDistance,
	double mass,
	utils::Vector<double, 3>& initialParticleVelocity,
	const SimulationDesc& desc,
	int particleType);


#endif /* PARTICLECUBOID_H_ */