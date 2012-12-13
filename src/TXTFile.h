//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File TXTFile.h
// contains class to handle input/output of txt files
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------
#ifndef TXTFILE_H_
#define TXTFILE_H_

#include "utils/utils.h"
#include "ParticleContainer.h"
#include "Simulation.h"

/// class to handle input/output of .txt files in a very simple format

class TXTFile
{
private:
	// a vector of particles
	std::vector<Particle> particles;

	// a vector of materials
	std::vector<Material> materials;
public:

	// read file
	err_type readFile(const char *filename);

	// write file
	err_type writeFile(const char *filename, const std::vector<Particle>& particles, const std::vector<Material>& materials);

	// get Particles
	std::vector<Particle>& getParticles()	{return particles;}

	// get Materials
	std::vector<Material>& getMaterials()	{return materials;}

};

#endif