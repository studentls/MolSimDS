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

#include "Simulation.h"

void AddParticleCuboid(ParticleContainer& cont,
	utils::Vector<double,3>& lowerLeftFrontCorner,
	int particlesPerDimension,
	double particleDistance,
	double mass,
	utils::Vector<double, 3>& initialParticleVelocity,
	const SimulationDesc& desc,
	int particleType) {
		// for every combination of valid x, y and z coordinates...
		for (int x = 0; x < particlesPerDimension; x++)
			for (int y = 0; y < particlesPerDimension; y++)
				for (int z = 0; z < particlesPerDimension; z++)
				{
					// find the position
					double posArr [] = {
						lowerLeftFrontCorner[0] + x * particleDistance,
						lowerLeftFrontCorner[1] + y * particleDistance,
						lowerLeftFrontCorner[2] + z * particleDistance
					};
					utils::Vector<double,3> pos(posArr);
					// uses a copy of initialParticleVelocity as the base velocity, without Brownian motion
					utils::Vector<double,3> vel(initialParticleVelocity);
					Particle particle(pos, vel, mass, particleType);
					// apply Brownian motion via the Boltzmann distribution
					MaxwellBoltzmannDistribution(particle, desc.brownianMotionFactor, 3);
					cont.AddParticle(particle);
				}
}