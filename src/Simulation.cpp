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
#include "outputWriter/VTKWriter.h"
#include "outputWriter/XYZWriter.h"

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
	vector<Particle>::iterator iterator;
	iterator = particles.begin();

	//go through particles
	while (iterator != particles.end()) {
		vector<Particle>::iterator innerIterator = particles.begin();
		utils::Vector<double, 3> forceAcc = 0.0;

		Particle& p1 = *iterator;
		while (innerIterator != particles.end()) {
			if (innerIterator != iterator) {

				Particle& p2 = *innerIterator;

				//using simple force calculation model
				double invdist = 1.0 / p1.getX().distance(p2.getX());

				//use for speed up
				double denominator = invdist * invdist * invdist; 

				double factor = p1.getM() * p2.getM() *denominator * desc.gravitational_constant;

				utils::Vector<double, 3> force = (p2.getX() - p1.getX()) * factor;

				//add individual particle to particle force to sum
				forceAcc = forceAcc + force;
			}
			++innerIterator;
		}

		//set new force (sets internally old force to the now old value)
		p1.setForce(forceAcc);
		++iterator;
	}
}

void Simulation::calculateX() {
	vector<Particle>::iterator iterator = particles.begin();
	while (iterator != particles.end()) {

		Particle& p = *iterator;

		//base calculation on Velocity-Störmer-Verlet-Algorithm

		//x_i ( t^{n+1} ) = x_i(T^n) + dt * v_i(t^n) + (dt)^2 * F_i(t^n) / 2m_i

		p.setX(p.getX() + desc.delta_t * p.getV() + desc.delta_t * desc.delta_t * p.getF() * 0.5 * (1.0 / p.getM()));

		++iterator;
	}
}


void Simulation::calculateV() {
	vector<Particle>::iterator iterator = particles.begin();
	while (iterator != particles.end()) {

		Particle& p = *iterator;

		//base calculation on Velocity-Störmer-Verlet-Algorithm

		//v_i ( t^{n+1} ) = v_i(t^n) + dt * (F_i(t^n) + F_i(T^{n+1} ) ) / 2m_i

		p.setV(p.getV() +  desc.delta_t * (p.getF() + p.getOldF() ) * 0.5 * (1.0 / p.getM() ));

		++iterator;
	}
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
