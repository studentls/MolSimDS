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
#include "PerformanceTest.h"

#include "SimulationDesc.h"

/// helper struct, stores thermodynamical statistics
struct ThermostatisticalData
{
	/// initial positions(necessary for diffusion computation
	utils::TFastArray<utils::Vector<double, 3> > initialPositions;

	/// diffusion (first coordinate is time t_i, second the diffusion value for this time
	/// a writeoutput everytime the statistical data is generated is suggested
	utils::TFastArray<utils::Vector<double, 2> > diffusionValues;
};

/// a class that is used to represent a simulation
class Simulation
{
private:
	/// stores values that are used for calculations for this Simulation
	SimulationDesc				desc;

	/// stores a simulation's statistical data
	SimulationStatistics		statistics;

	/// Thermodynamical statistics
	ThermostatisticalData		*thermostatistics;

	/// stores the particles used in this Simulation
	ParticleContainer			*particles;

	/// performs one time step based on delta_t
	void						performStep();

	/// calculate the force for all particles
	void						calculateF();

	/// resets the force on a particle for a new iteration. Used in calculateF()
	static void					forceResetter(void*, Particle& p);

	/// calculate and apply the force between a pair of particles. Used in calculateF()
	/// use Lennard Jones (12, 6) potential for inter molecular forces
	static void					forceLJCalculator(void*, Particle&, Particle&);

	/// calculate and apply the force between a pair of particles. Used in calculateF()
	/// use Smoothed Lennard Jones (12, 6) potential for inter molecular forces
	static void					forceSLJCalculator(void*, Particle&, Particle&);

	/// calculate and apply the force between a pair of particles. Used in calculateF()
	/// use gravitational forces for inter molecular forces
	static void					forceGravitationCalculator(void*, Particle&, Particle&);

	/// member variable holding a possible inter molecular force calculation function
	void						(*forceCalculator)(void*, Particle&, Particle&);

	/// apply a pulling force to a membrane. Used in Run()
	static void					forceCalculatorMembranePull(void*, Particle&);

	/// calculate Gravity force for each particle
	static void					gravityCalculator(void*, Particle&);

	/// helper function to calculate Diffusion (takes as argument a ThermodynamicalStatistics data structure)
	static void					diffusionCalculator(void*, Particle&);

	/// helper to give each particle an id
	static void					idAssigner(void*, Particle&);

	/// calculate the position for all particles
	void						calculateX();

	/// calculate the new position of a particle for a new iteration. Used in calculateX()
	static void					posCalculator(void*, Particle&);

	/// calculate the velocity for all particles
	void						calculateV();

	/// calculate the new velocity of a particle for a new iteration. Used in calculateV()
	static void					velCalculator(void*, Particle& p);
	
	/// plot the particles to a xyz-file
	/// @param iteration the number of this iteration
	void						plotParticles(int iteration);

	/// calculates kinetic energy of the system and outputs into a double
	static void					kineticEnergyCalculator(void* data, Particle& p);

	/// calculates summed mass of all particles, data is a double pointer
	static void					totalMassCalculator(void *data, Particle& p);

	/// set Brownianmotion hard
	static void					setBrownianMotionCalculator(void *data, Particle& p);

	/// initialize particles to given temperature
	void						initializeThermostat();

	/// adjust temperature
	void						adjustThermostat();

	/// calculate thermostatistical data
	void						calculateThermostatisticalData(const double t);

	/// writes CSV files for thermostatistical data
	err_type						writeThermostatisticalDataToFile();

	/// apply beta(Velocity Scale Factor) to particles
	static void					applyTemperatureScalingFactor(void *data, Particle& p);

	/// send particles to viewer
	void						notifyViewer();

public:
	Simulation():particles(NULL), thermostatistics(NULL)			{}

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

	// befriended with performance test class
	friend class PerformanceTest;
};


#endif
