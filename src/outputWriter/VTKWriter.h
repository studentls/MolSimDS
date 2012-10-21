/*
 * VTKWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#ifndef VTKWRITER_H_
#define VTKWRITER_H_

#include "../Particle.h"
#include "vtk-unstructured.h"

#include <list>

namespace outputWriter {

/**
 * This class implements the functionality to generate vtk output from particles.
 */
class VTKWriter {

public:
	VTKWriter();

	virtual ~VTKWriter();

	/**
	 * set up internal data structures and prepare to plot a particle.
	 */
	void initializeOutput(int numParticles);

	/**
	 * plot type, mass, position, velocity and force of a particle.
	 *
	 * @note: initializeOutput() must have been called before.
	 */
	void plotParticle(Particle& p);

	/**
	 * writes the final output file.
	 *
	 * @param filename the base name of the file to be written.
	 * @param iteration the number of the current iteration,
	 *        which is used to generate an unique filename
	 */
	void writeFile(const std::string& filename, int iteration);

   /**
	* plots particle list to file
	*
	* @param particles list of particles
	* @param filename the output file
	* @param iteration the number of the current iteration used to 
	*        create a unique filename
	*/
	void plotParticles(std::list<Particle> particles, const std::string& filename, int iteration);

private:
	VTKFile_t* vtkFile;
};

}

#endif /* VTKWRITER_H_ */
