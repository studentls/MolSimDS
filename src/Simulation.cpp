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

using namespace std;

err_type Simulation::Init(const SimulationDesc& desc)
{
	this->desc = desc;

	// clear particles if it is not already empty
	if(particles)if(!particles->IsEmpty())particles->Clear();


	return S_OK;
}

err_type Simulation::AddParticlesFromFile(const char *filename)
{
	//assert(particles);

	// read the file New !
	//if(!particles->AddParticlesFromFileNew(filename))return E_FILEERROR;
	
	// set up particles with some data...
	utils::Vector<double, 2> extent;
	extent[0] = 180;
	extent[1] = 90;

	ListParticleContainer PC;
	PC.AddParticlesFromFileNew(filename);

	//method
	
	// use LinkedCell
	particles = new LinkedCellParticleContainer<2>(PC.getParticles(), 3.0, utils::Vector<double,2>(0.0), extent, 1,
		true, true, true, true, true, true, 1.0); 

	// use standard
	//particles = new ListParticleContainer(PC.getParticles());


	// call calculateF() because the forces are needed to calculate x, but are not given in the input file.
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

	// set common statistical values...
	statistics.particle_count = particles->getParticleCount();
	statistics.step_count = (desc.end_time - desc.start_time) / desc.delta_t;

	// output that calculation have started ("starting calculation...")
	LOG4CXX_TRACE(simulationLogger, "starting calculation...");

	//start timer
	utils::Timer timer;

	// iterate until the end time is reached...
	while (current_time < desc.end_time) {

		// perform one iteration step
		performStep();

		// TODO: reset to 100
		// plot the particles on every hundredth iteration, beginning with the first
		if (iteration % 10 == 0) {
			plotParticles(iteration);
			
			// output that an iteration has finished
			if(iteration % 25 == 0)
			LOG4CXX_TRACE(simulationLogger, "Iteration " << iteration << " finished.");
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
	particles->Clear();

	SAFE_DELETE(particles);

	return S_OK;
}

void Simulation::calculateF() {
	// call particles.Iterate() on forceResetter
	particles->Iterate(forceResetter, (void*)&desc);
	// call particles.IteratePairwise() on forceCalculator
	particles->IteratePairwise(forceCalculator, (void*)&desc);

	// also add forces from boundary conditions if necessary

	// TODO: typechecking in c++? only do this if it's a linkedCell
	if ()
		particles->ApplyReflectiveBoundaryConditions();
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
	// base calculation on Velocity-St�rmer-Verlet-Algorithm
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
	// base calculation on Velocity-St�rmer-Verlet-Algorithm
	// v_i ( t^{n+1} ) = v_i(t^n) + dt * (F_i(t^n) + F_i(T^{n+1} ) ) / 2m_i
	p.v = (p.v +  desc->delta_t * (p.getF() + p.getOldF() ) * 0.5 * (1.0 / p.m ));
}

void Simulation::plotParticles(int iteration) {

	assert(particles);
	
	string out_name("MD_vtk");
	// switch between VTK and XYZ output
	// depending on the value of desc.output_fmt
	switch(desc.output_fmt)
	{
	case SOF_VTK:
		{
			// VTK Output
			outputWriter::VTKWriter writer;
			
			writer.plotParticles(particles->getParticles(), out_name, iteration);
			break;
		}
	case SOF_XYZ:
		{
			// XYZ Output
			outputWriter::XYZWriter writer;
			writer.plotParticles(particles->getParticles(), out_name, iteration);
			break;
		}
	default:
		{
		}
	}
}
