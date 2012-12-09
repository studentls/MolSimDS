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

#include "Logging.h"
#include "utils/utils.h"
#include "Particle.h"
#include "LinkedCellParticleContainer.h"
#include "ListParticleContainer.h"
#include "ParticleGenerator.h"
#include <vector>

/// Output Format
enum SimulationOutputFormat
{
	SOF_NONE,
	SOF_VTK,
	SOF_XYZ
};

/// struct to hold information about a material, that is assigned to a particle with type
struct Material
{
	double	sigma;
	double	epsilon;
	std::string name;

	/// constructor setting values to default values
	Material()
	{
		epsilon = 5.0;
		sigma	= 1.0; 
	}
};

/// a struct that describes common Simulation params
/// @param delta_t step size
/// @param start_time start time of simulation
/// @param end_time end time of simulation
/// @param gravitational_constant the gravitational constant used in the Simulation
/// @param output_fmt specify OutputFormat, use SOF_NONE for no output
/// @param iterationsperoutput after iterationsperoutput iterations, the simulation will output data
/// @param outname filename for output
struct SimulationDesc
{
	double				delta_t;
	double				start_time;
	double				end_time;
	
	double				brownianMotionFactor;

	unsigned int		iterationsperoutput;

	std::string			outname;

	/// stores all materials in the simulation
	std::vector<Material>	materials;

	// DEPRECATED
	// this value may be renamed appropriately in a later version
	// if the particles are affected by magnetic instead of gravitational forces
	//double	gravitational_constant;

	SimulationOutputFormat output_fmt;

	/// constructor that creates a default SimulationDesc object
	SimulationDesc()
	{
		delta_t = 0.0002;
		start_time = 0.0;
		end_time = 5.0;
		output_fmt = SOF_NONE;
		
	    brownianMotionFactor = 0.1;

		iterationsperoutput = 10;

		outname = "out";
		// DEPRECATED
		//gravitational_constant = 1.0;
	}
};


/// a struct to hold statistical data of a simulation
struct SimulationStatistics
{
	unsigned int	particle_count; /// total particles in simulation
	unsigned int	step_count;		/// how many steps?
	double			time;			/// contains how many seconds the simulation lasted
	double			timeperstep;	/// how long did a step last?

	/// constructor to set everything to zero
	SimulationStatistics()
	{
		particle_count	= 0;
		step_count		= 0;
		time			= 0.0;
		timeperstep		= 0.0;
	}
};

/// a class that is used to represent a simulation
class Simulation
{
private:
	/// stores values that are used for calculations for this Simulation
	SimulationDesc			desc;

	/// stores a simulation's statistical data
	SimulationStatistics	statistics;

	/// stores the particles used in this Simulation
	ParticleContainer		*particles;

	/// performs one time step based on delta_t
	void					performStep();

	/// calculate the force for all particles
	void					calculateF();

	/// resets the force on a particle for a new iteration. Used in calculateF()
	static void				forceResetter(void*, Particle& p);

	/// calculate and apply the force between a pair of particles. Used in calculateF()
	static void				forceCalculator(void*, Particle&, Particle&);

	/// calculate the position for all particles
	void					calculateX();

	/// calculate the new position of a particle for a new iteration. Used in calculateX()
	static void				posCalculator(void*, Particle&);

	/// calculate the velocity for all particles
	void					calculateV();

	/// calculate the new velocity of a particle for a new iteration. Used in calculateV()
	static void				velCalculator(void*, Particle& p);

	/// plot the particles to a xyz-file
	/// @param iteration the number of this iteration
	void					plotParticles(int iteration);

public:
	Simulation():particles(NULL)			{}

	~Simulation()							{ Release(); }
	
	/// sets up simulation using description desc, any old Simulation will be lost when this function is called
	/// @param desc simulation description
	/// @return returns always S_OK, 
	err_type				Init(const SimulationDesc& desc);

	/// parses an .xml file according to simulationfile.xsd and creates based on this file a simulation with settings & data
	/// @param filename xml file, must have .xml ending
	/// @return returns S_OK if file could be loaded successfully, otherwise E_FILEERROR
	err_type				CreateSimulationFromXMLFile(const char *filename);

	/// runs Simulation according to settings given by Init
	/// @return returns S_OK or E_NOTINITIALIZED if particle data is empty
	err_type				Run();

	/// cleans memory and flushes Simulation Data
	/// @return returns always S_OK
	err_type				Release();

	/// get statistics
	/// @return statistical data about the simulation
	SimulationStatistics&	getStatistics()	{return this->statistics;}
};


#endif
