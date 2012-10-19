
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


double start_time = 0;
double end_time = 1000;
double delta_t = 0.014;

std::list<Particle> particles;


int main(int argc, char* argsv[]) {
	
	cout << "Hello from MolSim for PSE!" << endl;
	if (argc != 2) {
		cout << "Errounous programme call! " << endl;
		cout << "./molsym filename" << endl;

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

		iteration++;
		if (iteration % 10 == 0) {
			plotParticles(iteration);
		}
		cout << "Iteration " << iteration << " finished." << endl;

		current_time += delta_t;
	}

	cout << "output written. Terminating..." << endl;
	return 0;
}


void calculateF() {
	list<Particle>::iterator iterator;
	iterator = particles.begin();

	while (iterator != particles.end()) {
		list<Particle>::iterator innerIterator = particles.begin();

		while (innerIterator != particles.end()) {
			if (innerIterator != iterator) {

				Particle& p1 = *iterator;
				Particle& p2 = *innerIterator;

				// insert calculation of force here!

			}
			++innerIterator;
		}
		++iterator;
	}
}


void calculateX() {
	list<Particle>::iterator iterator = particles.begin();
	while (iterator != particles.end()) {

		Particle& p = *iterator;

		//base calculation on Velocity-St�rmer-Verlet-Algorithm

		//x_i ( t^{n+1} ) = x_i(T^n) + dt * v_i(t^n) + (dt)^2 * F_i(t^n) / 2m_i

		p.setX(p.getX() + delta_t * p.getV() + delta_t * delta_t * p.getF() * 0.5 *(1.0 / p.getM()));

		++iterator;
	}
}


void calculateV() {
	list<Particle>::iterator iterator = particles.begin();
	while (iterator != particles.end()) {

		Particle& p = *iterator;

		//base calculation on Velocity-St�rmer-Verlet-Algorithm

		//v_i ( t^{n+1} ) = v_i(t^n) + dt * (F_i(t^n) + F_i(T^{n+1} ) ) / 2m_i

		p.setV(p.getV() +  delta_t * (p.getF() + p.getOldF() ) * 0.5 * (1.0 / p.getM() ));

		++iterator;
	}
}


void plotParticles(int iteration) {

	string out_name("MD_vtk");

	outputWriter::XYZWriter writer;
	writer.plotParticles(particles, out_name, iteration);
}
