//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File Simulation.h
// contains class Simulation
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------

#include "Simulation.h"
#include "outputWriter/VTKWriter.h"
#include "outputWriter/XYZWriter.h"
#include "XMLFileReader.h"


using namespace std;


err_type Simulation::Init(const SimulationDesc& desc)
{
	this->desc = desc;

	// clear particles if it is not already empty
	Release();


	return S_OK;
}

err_type Simulation::CreateSimulationFromXMLFile(const char *filename)
{

	XMLFileReader fr;

	if(FAILED(fr.readFile(filename)))return E_FILEERROR;

	this->desc = fr.getDescription();

	// get container
	Release(); // to free mem

	fr.makeParticleContainer(&particles);
	
	// is particles a valid pointer? - If not file has not been parsed correctly...
	if(!particles)return E_FILEERROR;

	//calc initial forces...
	calculateF();

	return S_OK;
}


void Simulation::performStep()
{
	// calculate new x
	calculateX();
	
	// calculate new f
	calculateF();

	// calculate new v
	calculateV();
}


err_type Simulation::Run()
{
	assert(particles);

	// check if the particles are valid
	if(particles->IsEmpty())return E_NOTINITIALIZED;

	// initialize some values
	double current_time = desc.start_time;
	int iteration = 0;
	updatingHeat = false;

	// set common statistical values...
	statistics.particle_count = particles->getParticleCount();
	statistics.step_count = (desc.end_time - desc.start_time) / desc.delta_t;

	// output that calculation have started ("starting calculation...")
	LOG4CXX_TRACE(simulationLogger, "starting calculation...");

	//start timer
	utils::Timer timer;

	// initiate heat for all particles
	initiateHeat();

	// iterate until the end time is reached...
	while (current_time < desc.end_time) {

		// perform one iteration step
		performStep();

	    // normalize v to get the heat you want
	    // do this only every thousandth iteration
		// unless the updating takes several iterations
		if (desc.timestepsPerThermostatApplication > 0)
		{
	        if (iteration % desc.timestepsPerThermostatApplication == 0)
		    	updatingHeat = true;
		    if (updatingHeat)
	            normalizeHeat();
		}

		// TODO: reset to 100
		// plot the particles on every hundredth iteration, beginning with the first
		if (iteration % desc.iterationsperoutput == 0) {
			plotParticles(iteration);
			
#ifndef DEBUG
			// output that an iteration has finished
			if(iteration % 25 == 0)
			LOG4CXX_TRACE(simulationLogger, "Iteration " << iteration << " finished.");
#else
			// output that an iteration has finished
			if(iteration % 5 == 0)
			LOG4CXX_TRACE(simulationLogger, "Iteration " << iteration << " finished.");
#endif
		}
		
		// increment loop values
		iteration++;
		current_time += desc.delta_t;
	}

	//generate statistical data
	statistics.time = timer.getElapsedTime();
	statistics.timeperstep = statistics.time / statistics.step_count;

	// output that the output has finished
	cout << endl;
	LOG4CXX_TRACE(simulationLogger, "output written. Terminating...");

	return S_OK;
}

err_type Simulation::Release()
{
	// delete Particle data
	if(particles)particles->Clear();

	SAFE_DELETE(particles);

	return S_OK;
}

void Simulation::calculateF() {

	// reassign particles, if LinkedCell Algorithm is used...
	if(particles->getType() == PCT_LINKEDCELL)
		((LinkedCellParticleContainer*)particles)->ReassignParticles();

	// call particles.Iterate() on forceResetter
	particles->Iterate(forceResetter, (void*)&desc);

	// call particles.IteratePairwise() on forceCalculator
	particles->IteratePairwise(forceCalculator, (void*)&desc);

	// apply Boundary conditions, if LinkedCell Algorithm is used...
	if(particles->getType() == PCT_LINKEDCELL)
		((LinkedCellParticleContainer*)particles)->ApplyReflectiveBoundaryConditions(forceCalculator, (void*)&desc);
}

void Simulation::forceResetter(void* data, Particle& p) {
	// call particlesresetForce(), which also ensures that old_f gets changed as well as f
	p.resetForce();
}

void Simulation::forceCalculator(void* data, Particle& p1, Particle& p2)
{
#ifdef DEBUG
	// assert data is a valid pointer!
	if(!data)
	{
		LOG4CXX_ERROR(simulationLogger, "error: data is not a valid pointer!");
		return;
	}
#endif

	SimulationDesc *desc = (SimulationDesc*)data;

	// calculate force via Lennard-Jones Potential
	// these calculations are rather straightforward, looking at the formula
	double dist = p1.x.distance(p2.x);
	double temp1 = desc->sigma / dist;
	// this is faster than power functions or manual multiplication
	double pow2 = temp1 * temp1;
	double pow3 = pow2 * temp1;
	double pow6 = pow3 * pow3;
	double pow12 = pow6 * pow6;

	double temp2 = pow6 - 2.0 * pow12;
	double prefactor = 24.0 * desc->epsilon / dist / dist;
	double factor = prefactor * temp2;
	
	// DEPRECATED
	// using simple force calculation model
	//double invdist = 1.0 / p1.x.distance(p2.x);
	// use for speed up
	//double denominator = invdist * invdist * invdist; 
	//double factor = p1.m * p2.m * denominator * desc->gravitational_constant;
	

	utils::Vector<double, 3> force = (p2.x - p1.x) * factor;
	
	// add individual particle to particle force to sum
	p1.addForce(force);
	
	
	p2.addForce(-1.0 * force);
}

void Simulation::calculateX() {
	
	assert(particles);

	// call particles.Iterate() on posCalculator
	particles->Iterate(posCalculator, (void*)&desc);
}

void Simulation::posCalculator(void* data, Particle& p) {
	#ifdef DEBUG
	// assert data is a valid pointer!
	if(!data)
	{
		LOG4CXX_ERROR(simulationLogger, "error: data is not a valid pointer!");
		return;
	}
#endif

	SimulationDesc *desc = (SimulationDesc*)data;
	
	// calc new pos
	// base calculation on Velocity-Störmer-Verlet-Algorithm
	// x_i ( t^{n+1} ) = x_i (t^n) + delta_t * p_i (t^n) + delta_t^2 * f_i (t^n) * 0.5 / m_i
	p.x = (p.x + desc->delta_t * p.v + desc->delta_t * desc->delta_t * p.getF() * 0.5 * (1.0 / p.m));
}

void Simulation::calculateV() {
	// call particles.Iterate() on velCalculator
	particles->Iterate(velCalculator, (void*)&desc);
}

void Simulation::velCalculator(void* data, Particle& p) {
	#ifdef DEBUG
	// assert data is a valid pointer!
	if(!data)
	{
		LOG4CXX_ERROR(simulationLogger, "error: data is not a valid pointer!");
		return;
	}
#endif

	SimulationDesc *desc = (SimulationDesc*)data;

	// calc new vel
	// base calculation on Velocity-Störmer-Verlet-Algorithm
	// v_i ( t^{n+1} ) = v_i(t^n) + dt * (F_i(t^n) + F_i(T^{n+1} ) ) / 2m_i
	p.v = (p.v +  desc->delta_t * (p.getF() + p.getOldF() ) * 0.5 * (1.0 / p.m ));
}

bool Simulation::updatingHeat = false;
double Simulation::sumOfKineticEnergies = 0.0;
double Simulation::heatNormalizationFactor = 0.0;
double Simulation::dBrownianMotionMassless = 0.0;

void Simulation::initiateHeat() {
	dBrownianMotionMassless = sqrt(desc.initialTemperature * desc.dimensions);
	particles->Iterate(heatInitializer, (void*)&desc);
}

void Simulation::heatInitializer(void* data, Particle& p) {
	#ifdef DEBUG
	// assert data is a valid pointer!
	if(!data)
	{
		LOG4CXX_ERROR(simulationLogger, "error: data is not a valid pointer!");
		return;
	}
#endif

	SimulationDesc *desc = (SimulationDesc*)data;

	// apply Brownian motion via Boltzmann distribution
	MaxwellBoltzmannDistribution(p, dBrownianMotionMassless / sqrt(p.m), desc->dimensions);
}

void Simulation::normalizeHeat() {
	// determine the current sum of kinetic energies
	sumOfKineticEnergies = 0.0;
	particles->Iterate(velocitySumAcquirer, (void*)&desc);
	// determine the current heat
	double currentHeat = sumOfKineticEnergies / particles->getParticleCount() / desc.dimensions * 2.0;
	double currentTargetTemperature = desc.targetTemperature;
	// ensure that the heat isn't adapted by too much at once
	if (currentTargetTemperature > currentHeat + desc.temperatureDifferenceStepSize)
		currentTargetTemperature = currentHeat + desc.temperatureDifferenceStepSize;
	else if (currentTargetTemperature < currentHeat - desc.temperatureDifferenceStepSize)
		currentTargetTemperature = currentHeat - desc.temperatureDifferenceStepSize;
	else
		updatingHeat = false;
	// determine the factor and update the velocities
	heatNormalizationFactor = sqrt(currentTargetTemperature / currentHeat);
	std::cout << "Heat: " << currentHeat << std::endl;
	std::cout << "Normalizing heat. Factor: " << heatNormalizationFactor << std::endl;
	particles->Iterate(heatNormalizer, (void*)&desc);
}

void Simulation::velocitySumAcquirer(void* data, Particle& p) {
	sumOfKineticEnergies += p.v.L2NormSq() * p.m * 0.5;
}

void Simulation::heatNormalizer(void* data, Particle& p) {
	p.v = p.v * heatNormalizationFactor;
}

void Simulation::plotParticles(int iteration) {

	assert(particles);
	
	// switch between VTK and XYZ output
	// depending on the value of desc.output_fmt
	switch(desc.output_fmt)
	{
	case SOF_VTK:
		{
			// VTK Output
			outputWriter::VTKWriter writer;
			writer.plotParticles(particles->getParticles(), desc.outname, iteration);
			break;
		}
	case SOF_XYZ:
		{
			// XYZ Output
			outputWriter::XYZWriter writer;
			writer.plotParticles(particles->getParticles(), desc.outname, iteration);
			break;
		}
	default:
		{
		}
	}
}
