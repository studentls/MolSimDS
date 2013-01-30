//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File SimulationDesc.h
// contains helper classes for Simulation
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------

#ifndef SIMULATIONDESCRIPTION_HEADER_
#define SIMULATIONDESCRIPTION_HEADER_

#include "Logging.h"
#include "utils/utils.h"
#include "Particle.h"
#include <vector>


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
	double	mass; /// mass of each particle that has this material assigned
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

/// enum to hold allowed force methods
enum EPotentialForce
{
	PF_LJ,		/// lennard jones potential
	PF_SLJ,		/// smoothed lennard jones potential
	PF_GRAVITY /// gravity force
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
	double			*cached_halfinvmass; // 0.5 * 1.0 / mass
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
		
	double					totalMass;								/// totalMass
	double					kineticEnergy;							/// calculated kinetic Energy

	double					SLJfactor;								/// smoothed Lennard Jones potential factor
	double					cutoffRadius;							/// cutoffRadius

	EPotentialForce			potentialForce;

	unsigned int			iterationsTillTStatisticsCalculation;			/// steps, till thermostatistical data is generated
	unsigned int			iterationsPerTStatisticsCalculation;			/// how many steps shall be averaged?
	double					rdfdelta_t;

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
		cached_halfinvmass = NULL;
		cachesize = 0;
		cacheline = 0;

		totalMass = 0.0;
		kineticEnergy = 0.0;

		SLJfactor = 1.0;
		cutoffRadius = 1.1225;

		// no statistics per default
		iterationsTillTStatisticsCalculation = 0;	
		iterationsPerTStatisticsCalculation = 0;
		rdfdelta_t = 0.1;

		// use LJ as a standard
		potentialForce = PF_LJ;
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
		cached_halfinvmass = NULL;

		totalMass = 0.0;
		kineticEnergy = 0.0;

		SLJfactor = desc.SLJfactor;
		cutoffRadius = desc.cutoffRadius;
		potentialForce = desc.potentialForce;

		iterationsTillTStatisticsCalculation = desc.iterationsTillTStatisticsCalculation;	
		iterationsPerTStatisticsCalculation  = desc.iterationsPerTStatisticsCalculation;	
		rdfdelta_t = desc.rdfdelta_t;
	}											

	~SimulationDesc()
	{
		SAFE_DELETE_A(cached_mat);
		SAFE_DELETE_A(cached_halfinvmass);
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

		totalMass = desc.totalMass;
		kineticEnergy = desc.kineticEnergy;

		SLJfactor = desc.SLJfactor;
		cutoffRadius = desc.cutoffRadius;
		potentialForce = desc.potentialForce;

		
		iterationsTillTStatisticsCalculation = desc.iterationsTillTStatisticsCalculation;	
		iterationsPerTStatisticsCalculation  = desc.iterationsPerTStatisticsCalculation;	
		rdfdelta_t = desc.rdfdelta_t;

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
		cached_halfinvmass = new double[cacheline];

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

			// inv mass
			cached_halfinvmass[i] = 0.5 * 1.0 / materials[i].mass;
		}
	}

	inline double getSigmaSq(const int type1, const int type2) const	{assert(type1 < materials.size()); assert(type2 < materials.size()); return cached_mat[type1 + cacheline * type2].sigmaSq;}
	inline double getEpsilon24(const int type1, const int type2) const 	{assert(type1 < materials.size()); assert(type2 < materials.size()); return cached_mat[type1 + cacheline * type2].epsilon24;}
	inline double getHalfInvMass(const int type) const {assert(type < materials.size()); return cached_halfinvmass[type];}
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


#endif