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

#ifndef SIMULATION_HEADER_
#define SIMULATION_HEADER_

#include "utils\utils.h"
#include "Particle.h"
#include "ParticleContainer.h"
#include <vector>



///
/// Output Format
///
enum SimulationOutputFormat
{
	SOF_NONE,
	SOF_VTK,
	SOF_XYZ
};

/// a struct that describes common Simulation params
/// @param delta_t step size
/// @param start_time start time of simulation
/// @param end_time end time of simulation
/// @param gravitational_constant the gravitational constant used in the Simulation
/// @param output_fmt specify OutputFormat, use SOF_NONE for no output
struct SimulationDesc
{
	double	delta_t;
	double	start_time;
	double	end_time;
	/// this value may be renamed appropriately in a later version
	/// if the particles are affected by magnetic instead of gravitational forces
	double	gravitational_constant;

	SimulationOutputFormat output_fmt;

	/// constructor that creates a default SimulationDesc object
	SimulationDesc()
	{
		delta_t = 0.01;
		start_time = 0.0;
		end_time = 1.0;
		output_fmt = SOF_NONE;
		gravitational_constant = 1.0;
	}
};


/// a class that is used to represent a simulation
class Simulation
{
private:
	/// stores values that are used for calculations for this Simulation
	SimulationDesc desc;

	/// stores the particles used in this Simulation
	ParticleContainer particles;

	/// performs one time step based on delta_t
	void performStep();

	/// calculate the force for all particles
	void calculateF();

	/// resets the force on a particle for a new iteration. Used in calculateF()
	static void forceResetter(void*, Particle& p);

	/// calculate and apply the force between a pair of particles. Used in calculateF()
	static void forceCalculator(void*, Particle&, Particle&);

	/// calculate the position for all particles
	void calculateX();

	/// calculate the new position of a particle for a new iteration. Used in calculateX()
	static void posCalculator(void*, Particle&);

	/// calculate the velocity for all particles
	void calculateV();

	/// calculate the new velocity of a particle for a new iteration. Used in calculateV()
	static void velCalculator(void*, Particle& p);

	/// plot the particles to a xyz-file
	/// @param iteration the number of this iteration
	void plotParticles(int iteration);

	static void test(Particle& p);


public:
	Simulation()	{}

	~Simulation()	{ Release(); }
	
	/// set up simulation using description desc
	/// @param desc simulation description
	err_type Init(const SimulationDesc& desc);

	/// adds Particles to particle list from file
	/// @param filename particle file, has to use .txt format
	err_type AddParticlesFromFile(const char *filename);

	/// Run Simulation
	/// @param showStatistics print out statistic Data at end
	err_type Run(bool showStatistics = false);

	/// cleans memory and flushes Simulation Data
	err_type Release();
};


#endif