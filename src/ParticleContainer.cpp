//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File ParticleContainer.cpp
// contains class ParticleContainer
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------

#include "ParticleContainer.h"
#include "FileReader.h"
#include "ParticleGenerator.h"

ParticleContainer::ParticleContainer() {
	particles = std::vector<Particle>();
}

ParticleContainer::ParticleContainer(const std::vector<Particle>& particles_args) {
	// copy data
	particles = particles_args;
}

void ParticleContainer::AddParticle(const Particle& particle) {
	particles.push_back(particle);
}

void ParticleContainer::Iterate(void(*func)(void*, Particle&), void *data) {
	// iterate over all Particles and call the function on it
	for (std::vector<Particle>::iterator it = particles.begin() ; it < particles.end(); it++)
	{
		Particle& p = *it;
		(*func)(data, p);
	}
}

void ParticleContainer::IteratePairwise(void(*func)(void *data, Particle&, Particle&), void *data) {
	// iterate over all Particles
	for (std::vector<Particle>::iterator it1 = particles.begin() ; it1 < particles.end(); it1++)
		for (std::vector<Particle>::iterator it2 = it1 + 1; it2 < particles.end(); it2++)
			// make sure that a Particle is not paired with itself
			if (it1 != it2)
			{
				// call the function on the pair of Particles
				Particle& p1 = *it1;
				Particle& p2 = *it2;
				(*func)(data, p1, p2);
			}
}

void ParticleContainer::AddParticlesFromFile(const char *filename)
{
	// read the file
	FileReader fileReader;
	fileReader.readFile(particles, filename);
}

/// Quick 'n' Dirty Method
bool ParticleContainer::AddParticlesFromFileNew(const char *filename)
{
	// implement later a separate class

	FILE *pFile = NULL;
	char buffer[512]; // only 512 characters per line allowed

	pFile = fopen(filename, "r");

	// valid handle?
	if(!pFile)return false;

	while(!feof(pFile))
	{
		memset(buffer, 0, 512 * sizeof(char));
		fgets(buffer, 512, pFile);

		// if line starts with # or is an empty line ignore
		if(buffer[0] == '#' || buffer[0] == '\n')continue;
		
		// include in doxygen
		/// the only allowed syntax will be 
		/// cuboid x1 x2 x3 v1 v2 v3 mass n1 n2 n3 meshwidth
		/// where x1, x2, x3 are the componets of the position vector
		///  of the lower left front corner of the cuboid
		/// v1 v2 v3 are the initial velocity parameters
		/// n1 n2 n3 are the dimensions of the cuboid
		/// meshwidth is the distance between two particles
		/// mass the initial mass
		
		utils::Vector<unsigned int, 3> dim;
		char cmd[256];
		Vec3 v;
		Vec3 x;
		double mass;
		double meshwidth;
		
		if(sscanf(buffer, "%s %lf %lf %lf %lf %lf %lf %lf %d %d %d  %lf", 
			&cmd, &x[0], &x[1], &x[2], &v[0], &v[1], &v[2], &mass, &dim[0], &dim[1], &dim[2], &meshwidth))
		{
			if(strcmp(cmd, "cuboid") == 0)
			{
				// assert values
				// max range for dim should be 1000x1000x1000
				for(int i = 0; i < 3; i++)assert(dim[i] < 1000);

				ParticleGenerator::makeCuboid(*this,
					x, dim, meshwidth, mass, v);
			}
		}
		else LOG4CXX_ERROR(particleGenerationLogger, " >> matching error, file corrupted?");
		
	}


	fclose(pFile);

	return true;
}