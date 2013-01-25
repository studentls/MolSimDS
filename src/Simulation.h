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
#include "MembraneContainer.h"
#include "ParticleGenerator.h"
#include <vector>
#include "Viewer.h"

/// Output Format
enum SimulationOutputFormat
{
	SOF_NONE,
	SOF_VTK,
	SOF_XYZ,
	SOF_TXT
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

/// struct, used to cache materials fast
struct CachedMaterial
{
	double sigmaSq; // sigma * sigma
	double epsilon24; // 24.0 * epsilon
};

/// a class that describes common Simulation params
/// @param delta_t step size
/// @param start_time start time of simulation
/// @param end_time end time of simulation
/// @param brownianMotionFactor the factor for the Brownian motion
/// @param epsilon the epsilon value
/// @param sigma the sigma value
/// @param dimensions the number of dimensions
/// @param timestepsPerThermostatApplication the number of iterations/timesteps before the heat is normalized again
/// @param targetTemperature the temperature at which the particles are kept by the heat normalizer
/// @param temperatureDifferenceStepSize the maximum difference in heat that may be made per heat normalization
/// @param gravitational_constant the gravitational constant used in the Simulation
/// @param output_fmt specify OutputFormat, use SOF_NONE for no output
/// @param iterationsperoutput after iterationsperoutput iterations, the simulation will output data
/// @param outname filename for output
class SimulationDesc
{
private:
	CachedMaterial	*cached_mat;
	unsigned int	cachesize;
	unsigned int	cacheline;
public:
	double					delta_t;
	double					start_time;
	double					end_time;
	
	double					brownianMotionFactor;

	unsigned int			iterationsperoutput;

	unsigned int			iterationsTillThermostatApplication;		/// steps, till thermostat is applied, 0 if no thermostat exists
	double					temperature;							/// current temperature in Kelvin
	double					targetTemperature;						/// target temperature
	double					temperatureStepSize;					/// after timestepsTillThermostatApplication the temperature will be increased till it reaches targetTemperature

	unsigned int			dimensions;								/// Dimensions of the simulation, can be 2 or 3

	std::string				outname;
		
	/// stores all materials in the simulation
	std::vector<Material>	materials;

	double	gravitational_constant;

	SimulationOutputFormat output_fmt;

	/// constructor that creates a default SimulationDesc object
	SimulationDesc()
	{
		delta_t = 0.0002;
		start_time = 0.0;
		end_time = 5.0;
		output_fmt = SOF_NONE;
		
	    brownianMotionFactor = 0.1;
		
		iterationsTillThermostatApplication = 0;
		temperature = 0;
		targetTemperature = 0;
		temperatureStepSize = 0;

		iterationsperoutput = 10;
		dimensions = 0;

		outname = "out";
		
		// no gravitational force per default
		gravitational_constant = 0.0;

		cached_mat = NULL;
		cachesize = 0;
		cacheline = 0;
	}

	/// copy constructor
	SimulationDesc(const SimulationDesc& desc)
	{
		delta_t = desc.delta_t;
		start_time = desc.start_time;
		end_time = desc.end_time;
		output_fmt = desc.output_fmt;
		
		brownianMotionFactor = desc.brownianMotionFactor;
		
		iterationsTillThermostatApplication = desc.iterationsTillThermostatApplication;
		temperature = desc.temperature;
		targetTemperature = desc.targetTemperature;
		temperatureStepSize = desc.temperatureStepSize;


		dimensions = desc.dimensions;

		iterationsperoutput = desc.iterationsperoutput;

		outname =  desc.outname;
		
		gravitational_constant = desc.gravitational_constant;

		// copy materials
		materials = desc.materials;

		// do not copy cache!
		cached_mat = NULL;
	}

	~SimulationDesc()
	{
		SAFE_DELETE_A(cached_mat);
	}

	SimulationDesc& operator=(const SimulationDesc& desc)
	{
		delta_t = desc.delta_t;
		start_time = desc.start_time;
		end_time = desc.end_time;
		output_fmt = desc.output_fmt;
		
		brownianMotionFactor = desc.brownianMotionFactor;
		
		iterationsTillThermostatApplication = desc.iterationsTillThermostatApplication;
		temperature = desc.temperature;
		targetTemperature = desc.targetTemperature;
		temperatureStepSize = desc.temperatureStepSize;


		dimensions = desc.dimensions;

		iterationsperoutput = desc.iterationsperoutput;

		outname =  desc.outname;
		
		gravitational_constant = desc.gravitational_constant;

		// copy materials
		materials = desc.materials;

		// do not copy cache!
		cached_mat = NULL;

		return *this;
	}

	/// function to test if a thermostat is present and should be applied
	/// @return true if thermostat is present and should be applied
	bool		applyThermostat()	{return iterationsTillThermostatApplication != 0;}

	/// generate Cache values
	void		generateCache()	
	{
		// delete if necessary
		SAFE_DELETE_A(cached_mat);

		// use a squared array to store values
		cachesize = materials.size() * materials.size(); // sacve size
		cacheline = materials.size();
		cached_mat = new CachedMaterial[cachesize];

		// go through mats and calculate pairs(actually one could save half the space, as only an upper/lower triangle matrix is needed, but for performance we
		// don't care for this tiny space saving
		int i = 0, j = 0;
		for(int i = 0; i < cacheline; i++)
		{
			for(int j = 0; j < cacheline; j++)
			{
				// mix materials
				double epsilon = sqrt(materials[i].epsilon * materials[j].epsilon);
				double sigma = (materials[i].sigma  + materials[j].sigma) * 0.5;

				// calc factors
				cached_mat[i + j * cacheline].epsilon24 = 24.0 * epsilon;
				cached_mat[i + j * cacheline].sigmaSq = sigma * sigma;
			}
		}
	}

	inline double getSigmaSq(const unsigned int& type1, const unsigned int& type2)	{assert(type1 < materials.size()); assert(type2 < materials.size()); return cached_mat[type1 + cacheline * type2].sigmaSq;}
	inline double getEpsilon24(const unsigned int& type1, const unsigned int& type2)	{assert(type1 < materials.size()); assert(type2 < materials.size()); return cached_mat[type1 + cacheline * type2].epsilon24;}
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

	/// apply a pulling force to a membrane. Used in Run()
	static void				forceCalculatorMembranePull(void*, Particle&);

	/// calculate Gravity force for each particle
	static void				gravityCalculator(void*, Particle&);

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

	/// calculates kinetic energy of the system and outputs into a double
	static void				kineticEnergyCalculator(void* data, Particle& p);

	/// calculates summed mass of all particles, data is a double pointer
	static void				totalMassCalculator(void *data, Particle& p);

	/// set Brownianmotion hard
	static void				setBrownianMotionCalculator(void *data, Particle& p);

	/// initialize particles to given temperature
	void					initializeThermostat();

	/// adjust temperature
	void					adjustThermostat();

	/// apply beta(Velocity Scale Factor) to particles
	static void				applyTemperatureScalingFactor(void *data, Particle& p);

	/// send particles to viewer
	void					notifyViewer();

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

	/// @return returns count of particles in the container
	int						getParticleCount()	{return (particles) ? particles->getParticleCount() : 0;}

	/// @return returns simulation area extent, including center
	utils::BoundingBox		getSimulationAreaExtent()	{return particles ? particles->getBoundingBox() : utils::BoundingBox();}

	/// @return returns the particle container pointer
	ParticleContainer*		getParticleContainer()	{return particles;}
};


#endif
