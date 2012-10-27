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

///
/// describes common Simulation params
/// @param delta_t step size
/// @param start_time start time of simulation
/// @param end_time end time of simulation
/// @param output_fmt specify OutputFormat, use SOF_NONE for no output
///
struct SimulationDesc
{
	double	delta_t;
	double	start_time;
	double	end_time;
	double	gravitational_constant;

	SimulationOutputFormat output_fmt;

	//also structs can contain a constructor...
	SimulationDesc()
	{
		delta_t = 0.01;
		start_time = 0.0;
		end_time = 1.0;
		output_fmt = SOF_NONE;
		gravitational_constant = 1.0;
	}
};


///
/// Molecule Simulation class
///
class Simulation
{
private:
	SimulationDesc	desc;

	ParticleContainer particles;

	/// performs one time step based on delta_t
	void performStep();


	//prototypes...
	/**
	 * calculate the force for all particles
	 */
	void calculateF();

	/**
	 * calculate the position for all particles
	 */
	void calculateX();

	/**
	 * calculate the velocity for all particles
	 */
	void calculateV();

	/**
	 * plot the particles to a xyz-file
	 */
	void plotParticles(int iteration);


public:
	Simulation()	{}

	~Simulation()	{ Release(); }
	///
	/// set up simulation using description desc
	// @param desc simulation description
	///
	err_type Init(const SimulationDesc& desc);

	///
	/// adds Particles to particle list from file
	/// @param filename particle file, has to use .txt format
	///
	err_type AddParticlesFromFile(const char *filename);

	///
	/// Run Simulation
	/// @param showStatistics print out statistic Data at end
	///
	err_type Run(bool showStatistics = false);

	///
	/// cleans memory and flushes Simulation Data
	/// 
	err_type Release();
};


#endif