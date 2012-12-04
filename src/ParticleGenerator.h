//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File MolSim.cpp
// contains main class MolSim
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------


#ifndef PARTICLE_GENERATOR_H_
#define PARTICLE_GENERATOR_H_

#include "ParticleContainer.h"
#include "MaxwellBoltzmannDistribution.h"

/// a class to generate particles
class ParticleGenerator
{
public:
	/// adds a cuboid of particles to a ParticleContainer
	/// @param pc the ParticleContainer where particles will be added to
	/// @param vLowerLeftFrontCorner the lower left front corner of the cuboid
	/// @param dimensions specify the dimensions according to SizeX x SizeY x SizeZ
	/// @param meshWidth distance between neighbouring particles
	/// @param mass the mass of each particle
	/// @param vInitialParticleVelocity the initial velocity of each particle
	/// @param dBrownianMotion brownianMotionFactor
	/// @param type ParticleType
	static void	makeCuboid(ParticleContainer& pc,
					   const Vec3& vLowerLeftFrontCorner,
					   const utils::Vector<unsigned int, 3>& dimensions,
					   const double meshWidth,
					   const double mass,
					   const Vec3& vInitialVelocity,
					   const double dBrownianMotion = 0.1,// set this later individual, use default value at the moment
					   const int type = 0) 
	{
		// acknowledge that default constructor sets forces to zero...
		Particle p;

		// go through mesh grid...
		for(unsigned int x = 0; x < dimensions[0]; x++)
			for(unsigned int y = 0; y < dimensions[1]; y++)
				for(unsigned int z = 0; z < dimensions[2]; z++)
				{
					p.x[0]	= vLowerLeftFrontCorner[0] + x * meshWidth;
					p.x[1]	= vLowerLeftFrontCorner[1] + y * meshWidth;
					p.x[2]	= vLowerLeftFrontCorner[2] + z * meshWidth;

					p.v		= vInitialVelocity;
					p.type	= type;
					p.m		= mass;

					// apply Brownian motion via Boltzmann distribution

					// at the moment only 2D !!!
					MaxwellBoltzmannDistribution(p, dBrownianMotion, 2);
					
					// add it
					pc.AddParticle(p);
				}
	}

};

#endif

