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
	if(!particles.empty())particles.clear();

	return S_OK;
}

err_type Simulation::AddParticlesFromFile(const char *filename)
{
	FileReader fileReader;
	fileReader.readFile(particles, filename);
	
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
	//particles valid?
	if(particles.empty())return E_NOTINITIALIZED;

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
	if(!particles.empty())particles.clear();

	return S_OK;
}

void Simulation::calculateF() {
	particles.Iterate(forceResetter);
	particles.IteratePairwise(forceCalculator);
}

void Simulation::forceResetter(Particle p) {
	p.resetForce();
}

void Simulation::forceCalculator(Particle p1, Particle p2)
	vector<Particle>::iterator iterator;
	//using simple force calculation model
	double invdist = 1.0 / p1.getX().distance(p2.getX());
	//use for speed up
	double denominator = invdist * invdist * invdist; 
	double factor = p1.getM() * p2.getM() *denominator * desc.gravitational_constant;
	utils::Vector<double, 3> force = (p2.getX() - p1.getX()) * factor;
	//add individual particle to particle force to sum
	p1.addForce(force);
	p2.addForce(-force)
}

void Simulation::calculateX() {
	particles.Iterate(posCalculator);
}

void Simulation::posCalculator(Particle p) {
	//base calculation on Velocity-Störmer-Verlet-Algorithm
	//v_i ( t^{n+1} ) = v_i(t^n) + dt * (F_i(t^n) + F_i(T^{n+1} ) ) / 2m_i
	p.setX(p.getX() + desc.delta_t * p.getV() + desc.delta_t * desc.delta_t * p.getF() * 0.5 * (1.0 / p.getM()));
}

void Simulation::calculateV() {
	particles.Iterate(velCalculator);
}

void Simulation::velCalculator(Particle p) {
	//base calculation on Velocity-Störmer-Verlet-Algorithm
	//v_i ( t^{n+1} ) = v_i(t^n) + dt * (F_i(t^n) + F_i(T^{n+1} ) ) / 2m_i
	p.setV(p.getV() +  desc.delta_t * (p.getF() + p.getOldF() ) * 0.5 * (1.0 / p.getM() ));
}


void Simulation::plotParticles(int iteration) {

	string out_name("MD_vtk");

	switch(desc.output_fmt)
	{
	case SOF_VTK:
		{
			//VTK Output
			outputWriter::VTKWriter writer;
			writer.plotParticles(particles, out_name, iteration);
			break;
		}
	case SOF_XYZ:
		{
			//XYZ Output
			outputWriter::XYZWriter writer;
			writer.plotParticles(particles, out_name, iteration);
			break;
		}

	default:
		{
		}
	}
}
