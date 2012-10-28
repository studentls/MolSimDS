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
#include "FileReader.h"
#include "outputWriter\VTKWriter.h"
#include "outputWriter\XYZWriter.h"

using namespace std;

err_type Simulation::Init(const SimulationDesc& desc)
{
	this->desc = desc;

	//clear
	if(!particles.getParticles().empty())particles.getParticles().clear();

	return S_OK;
}

err_type Simulation::AddParticlesFromFile(const char *filename)
{
	FileReader fileReader;
	fileReader.readFile(particles.getParticles(), filename);
	
	// the forces are needed to calculate x, but are not given in the input file.
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

err_type Simulation::Run(bool showStatistics)
{	
	active_desc = desc;

	//particles valid?
	if(particles.getParticles().empty())return E_NOTINITIALIZED;

	double current_time = desc.start_time;

	int iteration = 0;

	cout<<"starting calculation..."<<endl;

	 // for this loop, we assume: current x, current f and current v are known
	while (current_time < desc.end_time) {

		//step Simulation forward
		performStep();
		
		//we want to start plotting from initial configuration onwards!!!
		if (iteration % 10 == 0) {
			plotParticles(iteration);
		}
		
		iteration++;
		
		cout << "Iteration " << iteration << " finished." << endl;

		current_time += desc.delta_t;
	}
	cout << endl;
	cout << "output written. Terminating..." << endl;

	return S_OK;
}

err_type Simulation::Release()
{
	//delete Particle data...
	if(!particles.getParticles().empty())particles.getParticles().clear();

	return S_OK;
}

void Simulation::calculateF() {
	particles.Iterate(forceResetter);
	particles.IteratePairwise(forceCalculator);
}

void Simulation::forceResetter(Particle& p) {
	p.resetForce();
}

void Simulation::forceCalculator(Particle& p1, Particle& p2)
	//using simple force calculation model
	double invdist = 1.0 / p1.x.distance(p2.x);
	//use for speed up
	double denominator = invdist * invdist * invdist; 
	double factor = p1.m * p2.m * denominator * active_desc.gravitational_constant;
	utils::Vector<double, 3> force = (*p2.x - *p1.x) * factor;
	//add individual particle to particle force to sum
	p1.addForce(force);
	p2.addForce(-1 * force);
}

void Simulation::calculateX() {
	particles.Iterate(posCalculator);
}

void Simulation::posCalculator(Particle& p) {
	// TODO: accidentally deleted the description here. put it back in
	// when merging with main branch
	p.x = (p.x + active_desc.delta_t * p.v + active_desc.delta_t * active_desc.delta_t * p.getF() * 0.5 * (1.0 / p.m));
}

void Simulation::calculateV() {
	particles.Iterate(velCalculator);
}

void Simulation::velCalculator(Particle& p) {
	//base calculation on Velocity-Störmer-Verlet-Algorithm
	//v_i ( t^{n+1} ) = v_i(t^n) + dt * (F_i(t^n) + F_i(T^{n+1} ) ) / 2m_i
	p.v = (p.v +  active_desc.delta_t * (p.getF() + p.getOldF() ) * 0.5 * (1.0 / p.m ));
}


void Simulation::plotParticles(int iteration) {

	string out_name("MD_vtk");

	switch(desc.output_fmt)
	{
	case SOF_VTK:
		{
			//VTK Output
			outputWriter::VTKWriter writer;
			writer.plotParticles(particles.getParticles(), out_name, iteration);
			break;
		}
	case SOF_XYZ:
		{
			//XYZ Output
			outputWriter::XYZWriter writer;
			writer.plotParticles(particles.getParticles(), out_name, iteration);
			break;
		}

	default:
		{
		}
	}
}
