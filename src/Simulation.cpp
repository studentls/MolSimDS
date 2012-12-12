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

	desc = fr.getDescription();

	// get container
	Release(); // to free mem

	fr.makeParticleContainer(&particles);
	
	// is particles a valid pointer? - If not file has not been parsed correctly...
	if(!particles)return E_FILEERROR;

	// calc initial forces...
	calculateF();

	// is a thermostat present? if yes, initialize it
	if(desc.applyThermostat())initializeThermostat();

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

	// set common statistical values...
	statistics.particle_count = particles->getParticleCount();
	statistics.step_count = (desc.end_time - desc.start_time) / desc.delta_t;

	// output that calculation have started ("starting calculation...")
	LOG4CXX_TRACE(simulationLogger, "starting calculation...");

	//start timer
	utils::Timer timer;

	// iterate until the end time is reached...
	while (current_time < desc.end_time) {

		// thermostat present? if yes apply(for iteration 0 thermostat is already set)
		if(desc.applyThermostat())
			if( iteration % desc.iterationsTillThermostatApplication == 0 && iteration != 0)
				adjustThermostat();
		
		// perform one iteration step
		performStep();

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

	// call gravityCalculator to add gravity force for each particle
	if(desc.gravitational_constant != 0.0) // bad floating point variable check, but in this case o.k.
		particles->Iterate(gravityCalculator, (void*)&desc);

	// apply Boundary conditions, if LinkedCell Algorithm is used...
	if(particles->getType() == PCT_LINKEDCELL)
	{
		((LinkedCellParticleContainer*)particles)->ApplyReflectiveBoundaryConditions(forceCalculator, (void*)&desc);
		((LinkedCellParticleContainer*)particles)->ApplyPeriodicBoundaryConditionsForce(forceCalculator, (void*)&desc);
	}
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

	// assert indices
	assert(p1.type >= 0);
	assert(p2.type >= 0);
	assert(p1.type < desc->materials.size());
	assert(p2.type < desc->materials.size());

	// cache this later for better performance...
	double epsilon = sqrt(desc->materials[p1.type].epsilon * desc->materials[p2.type].epsilon);
	double sigma = (desc->materials[p1.type].sigma  + desc->materials[p2.type].sigma) * 0.5;

	// calculate force via Lennard-Jones Potential
	// these calculations are rather straightforward, looking at the formula
	double dist = p1.x.distance(p2.x);
	double temp1 = sigma / dist;
	// this is faster than power functions or manual multiplication
	double pow2 = temp1 * temp1;
	double pow3 = pow2 * temp1;
	double pow6 = pow3 * pow3;
	double pow12 = pow6 * pow6;

	double temp2 = pow6 - 2.0 * pow12;
	double prefactor = 24.0 * epsilon / dist / dist;
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

void Simulation::gravityCalculator(void *data, Particle& p)
{
	SimulationDesc *desc = (SimulationDesc*)data;

	// add gravitational force, based on G = m * g
	utils::Vector<double, 3> grav_force;

	// only y - component is affected...
	grav_force[1] = p.m * desc->gravitational_constant;

	p.addForce(grav_force);
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

void Simulation::kineticEnergyCalculator(void *data, Particle& p)
{
	double *eout = (double*)data;

	// calc according to
	// Task 1, Sheet 4
	*eout += p.m * p.v.L2NormSq() * 0.5;
}

void Simulation::totalMassCalculator(void *data, Particle& p)
{
	double *sum = (double*)data;

	*sum += p.m;
}

void Simulation::setBrownianMotionCalculator(void *data, Particle& p)
{
	double brownianMotionFactor = ((SimulationDesc*)data)->brownianMotionFactor;
	unsigned int dim = ((SimulationDesc*)data)->dimensions;

	// apply brownian Motion
	MaxwellBoltzmannDistribution(p, brownianMotionFactor, dim);
}

void Simulation::applyTemperatureScalingFactor(void *data, Particle& p)
{
	double beta = *((double*)data);
	p.v = p.v * beta;
}

void Simulation::initializeThermostat()
{
	assert(particles);

	// calc kinetic Energy of the whole system
	// temperature is start temperature
	// currently we use a dimensionless method
	double EnergyProposed = 0.5 * desc.dimensions * particles->getParticleCount() *desc.temperature;

	// we need the total mass of all particles...
	double totalMass = 0.0;
	particles->Iterate(totalMassCalculator, (void*)&totalMass);

	// calculate inital velocity that is set for each particle
	desc.brownianMotionFactor = sqrt(2.0 * EnergyProposed / (desc.dimensions * totalMass));

	// set brownianMotion to every particle
	particles->Iterate(setBrownianMotionCalculator, (void*)&desc);

}

void Simulation::adjustThermostat()
{
	assert(particles);
	
	// inc temperature
	desc.temperature += desc.temperatureStepSize;

	// secure that temperature is not above target temperature
	if(desc.temperature > desc.targetTemperature)desc.temperature = desc.targetTemperature;
	// calc beta

	//get energy of the system
	double EnergyAtTheMoment = 0.0;
	particles->Iterate(kineticEnergyCalculator, (double*)&EnergyAtTheMoment);

	// proposed energy
	double EnergyProposed = 0.5 * desc.dimensions * particles->getParticleCount() *desc.temperature;

	// beta
	double beta = sqrt(EnergyProposed / EnergyAtTheMoment);

	particles->Iterate(applyTemperatureScalingFactor, (double*)&beta);

	
}
