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
	/// @param vInitialVelocity the initial velocity of each particle
	/// @param epsilon the epsilon value to use for the particles
	/// @param sigma the sigma value to use for the particles
	/// @param dBrownianMotion brownianMotionFactor
	/// @param type ParticleType
	static void	makeCuboid(ParticleContainer& pc,
					   const Vec3& vLowerLeftFrontCorner,
					   const utils::Vector<unsigned int, 3>& dimensions,
					   const double meshWidth,
					   const double mass,
					   const Vec3& vInitialVelocity,
					   const double epsilon,
					   const double sigma,
					   const double dBrownianMotion = 0.1,// set this later individual, use default value at the moment
					   const int type = 0)
	{
		// acknowledge that default constructor sets forces to zero...
		Particle p;

		int dim = dimensions[2] <= 1 ? 2 : 3;

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

					p.epsilon = epsilon;
					p.sigma = sigma;

					// apply Brownian motion via Boltzmann distribution

					// at the moment only 2D !!!
					MaxwellBoltzmannDistribution(p, dBrownianMotion, dim);
					
					// add it
					pc.AddParticle(p);
				}
	}

	/// adds a grid forming a sphere of particles to a ParticleContainer
	/// @param pc the ParticleContainer where particles will be added to
	/// @param vCenter the center of the sphere
	/// @param radius radius of the sphere measured in particles
	/// @param meshWidth distance between neighbouring particles
	/// @param mass the mass of each particle
	/// @param vInitialVelocity the initial velocity of each particle
	/// @param epsilon the epsilon value to use for the particles
	/// @param sigma the sigma value to use for the particles
	/// @param dimensions valid values are 2 or 3. For dimensions = 2, this will generate a flat disc, for 3 a real sphere
	/// @param dBrownianMotion brownianMotionFactor
	/// @param type ParticleType
	static void makeSphere(ParticleContainer& pc,
							const Vec3& vCenter,
							const Vec3& vInitialVelocity,
							const double mass,
							const unsigned int radius,
							const double meshWidth,
						    const double epsilon,
						    const double sigma,
							const unsigned int dimensions = 2,
							const double dBrownianMotion = 0.1,// set this later individual, use default value at the moment
							const int type = 0)
	{
		assert(dimensions >=2 && dimensions <= 3);

		int r = radius;

		if(dimensions == 2)
		{
			for(int x = -r; x <= r; x++)
				for(int y = -r; y <= r; y++)
				{
					// is coordinate in sphere?
					// x^2 + y^2 <= r^2!
					if(x*x + y*y <= r * r)
					{
						Particle p;
						p.x[0]		= (double)x * meshWidth + vCenter[0];
						p.x[1]		= (double)y * meshWidth + vCenter[1];
						p.x[2]		= 0.0;

						p.v			= vInitialVelocity;
						p.type		= type;
						p.m			= mass;

						p.epsilon = epsilon;
						p.sigma = sigma;

						// apply Brownian motion via Boltzmann distribution
						MaxwellBoltzmannDistribution(p, dBrownianMotion, dimensions);
						pc.AddParticle(p);
					}
				}
		}
		else if(dimensions == 3)
		{
			for(int z = -r; z <= r; z++)
				for(int x = -r; x <= r; x++)
					for(int y = -r; y <= r; y++)
					{
						// is coordinate in sphere?
						// x^2 + y^2 + z^2 <= r^2!
						if(x*x + y*y + z*z <= radius * radius)
						{
							Particle p;
							p.x[0]		= (double)x * meshWidth + vCenter[0];
							p.x[1]		= (double)y * meshWidth + vCenter[1];							
							p.x[2]		= (double)z * meshWidth + vCenter[2];

							p.v			= vInitialVelocity;
							p.type		= type;
							p.m			= mass;

							// apply Brownian motion via Boltzmann distribution
							MaxwellBoltzmannDistribution(p, dBrownianMotion, dimensions);
							pc.AddParticle(p);
						}
					}
		}
	}



};

#endif

