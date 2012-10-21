
#include "outputWriter/XYZWriter.h"
#include "outputWriter/VTKWriter.h"
#include "FileReader.h"

#include <list>
#include <cstring>
#include <cstdlib>
#include <iostream>

using namespace std;

/**** forward declaration of the calculation functions ****/

/**
 * calculate the force for all particles
 */
void calculateF();

/**
 * calculate the position for all particles
 */
void calculateX();

/**
 * calculate the velocity for all particles
 */
void calculateV();

/**
 * plot the particles to a xyz-file
 */
void plotParticles(int iteration);



//use better default values!!! tooo long!
double start_time = 0;
double end_time = 10;
double delta_t = 0.014;

std::list<Particle> particles;

//console tool
void line()
{
	for(int i = 0; i < 40; i++)cout<<"-";
	cout<<endl;
}


int main(int argc, char* argsv[]) {

	//added header
	line();
	cout << "MolSim for PSE" << endl;
	line();
	cout << endl;
	cout << "(c) 2012 by F.Dietz & L.Spiegelberg" << endl; 
	cout << endl;
	cout << "Molecular Simulator handling *.txt files" << endl;
	line();
	cout << endl;


	//test for correct argument count
	if (argc != 2) {
		cout << "Errounous programme call! " << endl;
		cout << "usage ./molsym filename" << endl;

		//quit program to prevent bad memory access!
		return 0;
	}
	
	FileReader fileReader;
	fileReader.readFile(particles, argsv[1]);
	// the forces are needed to calculate x, but are not given in the input file.
	calculateF();

	double current_time = start_time;

	int iteration = 0;

	 // for this loop, we assume: current x, current f and current v are known
	while (current_time < end_time) {
		// calculate new x
		calculateX();

		// calculate new f
		calculateF();
		// calculate new v
		calculateV();
		
		//we want to start plotting from initial configuration onwards!!!
		if (iteration % 10 == 0) {
			plotParticles(iteration);
		}
		
		iteration++;
		
		cout << "Iteration " << iteration << " finished." << endl;

		current_time += delta_t;
	}

	cout << "output written. Terminating..." << endl;
	return 0;
}


// gets the distance between two Particles
// use double instead of long! long is an integer,here a floating point type is needed
//revisited code
double GetDistance(Particle p1, Particle p2) {
	
	//distance between p1 and p2
	utils::Vector<double, 3> distance = p1.getX() - p2.getX();

	return distance.L2Norm();
}

// the strength of gravity
const double gravitationalConstant = 1.0;

void calculateF() {
	list<Particle>::iterator iterator;
	iterator = particles.begin();

	//go through particles
	while (iterator != particles.end()) {
		list<Particle>::iterator innerIterator = particles.begin();
		utils::Vector<double, 3> forceAcc = 0.0;

		Particle& p1 = *iterator;
		while (innerIterator != particles.end()) {
			if (innerIterator != iterator) {

				Particle& p2 = *innerIterator;

				//using simple force calculation model
				double invdist = 1.0 / GetDistance(p1, p2);

				//use for speed up
				double denominator = invdist * invdist * invdist; 

				double factor = p1.getM() * p2.getM() *denominator * gravitationalConstant;

				utils::Vector<double, 3> force = (p1.getX() - p2.getX()) * factor;

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


void calculateX() {
	list<Particle>::iterator iterator = particles.begin();
	while (iterator != particles.end()) {

		Particle& p = *iterator;

		//base calculation on Velocity-Störmer-Verlet-Algorithm

		//x_i ( t^{n+1} ) = x_i(T^n) + dt * v_i(t^n) + (dt)^2 * F_i(t^n) / 2m_i

		p.setX(p.getX() + delta_t * p.getV() + delta_t * delta_t * p.getF() * 0.5 *(1.0 / p.getM()));

		++iterator;
	}
}


void calculateV() {
	list<Particle>::iterator iterator = particles.begin();
	while (iterator != particles.end()) {

		Particle& p = *iterator;

		//base calculation on Velocity-Störmer-Verlet-Algorithm

		//v_i ( t^{n+1} ) = v_i(t^n) + dt * (F_i(t^n) + F_i(T^{n+1} ) ) / 2m_i

		p.setV(p.getV() +  delta_t * (p.getF() + p.getOldF() ) * 0.5 * (1.0 / p.getM() ));

		++iterator;
	}
}


void plotParticles(int iteration) {

	string out_name("MD_vtk");

	//xyz output...
	//outputWriter::XYZWriter writer;
	//writer.plotParticles(particles, out_name, iteration);

	//VTK Output
	outputWriter::VTKWriter writer;
	writer.plotParticles(particles, out_name, iteration);
}
