//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File TXTFile.cpp
// contains class to handle input/output of txt files
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------

#include "TXTFile.h"

#define MODE_DATA		0x1
#define MODE_MATERIAL	0x2

// small helper to delete comment line starting with #
void deletecomment(char *str, int length)
{
	int pos = 0;
	bool found = false;
	while(pos < length && str[pos] != '\0')
	{
		// comment?
		if(str[pos] == '#')found = true;

		if(found)str[pos] = '\0';

		pos++;
	}
}

err_type TXTFile::readFile(const char *filename)
{
	FILE		*pFile = NULL;
	char		buffer[512];
	int			mode = 0;
	int			line = 0;

	// does file exist?
	if(!utils::fileExists(filename))
	{
		LOG4CXX_ERROR(generalOutputLogger, "file "<<filename <<" not found");
		return E_FILENOTFOUND;
	}

	// open file (in ASCII mode)
	pFile = fopen(filename, "r");
	if(!pFile)return E_UNKNOWN;

	// clear internal containers...
	if(!particles.empty())particles.clear();
	if(!materials.empty())materials.clear();

	// go through file...
	while(!feof(pFile))
	{
		memset(buffer, 0, 512 * sizeof(char));
		fgets(buffer, 512, pFile);

		// empty line? skip it if so
		if(buffer[0] == '\n')continue;

		// search for first #, because then it is a comment line
		int pos = 0;
		while(buffer[pos] == ' ' && pos < 512)pos++; // skip spaces

		// skip comment line
		if(buffer[pos] == '#')continue;

		// now it is a real line
		// determine if the line contains
		// material:
		// or
		// data:
		// first set # and everything behind to zero
		deletecomment(buffer, 512);

		//use strcmp
		if(0 == strcmp(buffer, "material:"))
		{
			mode = MODE_MATERIAL;
		}

		if(0 == strcmp(buffer, "data:"))
		{
			mode = MODE_DATA;
		}

		// scan file
		if(MODE_MATERIAL == mode)
		{
			Material mat;
			char strbuf[256];
			memset(buffer, 0, 256 * sizeof(char));

			// syntax is
			//	   name epsilon sigma
			// e.g mat 1.0 2.0
			if(sscanf(buffer, "%s %lf %lf", strbuf, &mat.epsilon, &mat.sigma) <= 0)
			{
				LOG4CXX_ERROR(generalOutputLogger, "failed to parse file "<<filename<<" at line "<<line<<" correctly");
				return E_PARSEERROR;
			}

			mat.name = strbuf;

			materials.push_back(mat);
		}

		if(MODE_DATA == mode)
		{
			// only particles currently supported for simplicity reasons
			Particle p;

			// syntax is
			//	   xyz, velocity, mass, type
			//	   where xyz, velocity are 3D vectors and type is an index to the i-th material
			// e.g 1.0 1.0 1.0, 0.0 -10.0 0.0, 1.0, 0
			if(sscanf(buffer, "%lf %lf %lf, %lf %lf %lf, %lf, %d",
				&p.x[0], &p.x[1], &p.x[2], &p.v[0], &p.v[1], &p.v[2],
				&p.m, &p.type) <= 0)
			{
				LOG4CXX_ERROR(generalOutputLogger, "failed to parse file "<<filename<<" at line "<<line<<" correctly");
				return E_PARSEERROR;
			}

			// push back
			particles.push_back(p);

		}

		// inc line pointer
		line++;

	}

	// close file
	fclose(pFile);


	// assert types
	for(std::vector<Particle>::iterator it = particles.begin(); it != particles.end(); it++)
	{
		// any materials?
		if(materials.empty())it->type = 0;

		// secure index of type is valid, otherwise set it to zero
		if(it->type < 0)
		{
			it->type = 0;
			LOG4CXX_INFO(generalOutputLogger, "invalid index in "<<filename<<" found");
		}
		if(it->type >= materials.size())
		{
			it->type = 0;
			LOG4CXX_INFO(generalOutputLogger, "invalid index in "<<filename<<" found");
		}
	}

	// evrything parsed, o.k.
	return S_OK;
}

err_type TXTFile::writeFile(const char *filename, const std::vector<Particle>& particles, const std::vector<Material>& materials)
{
	FILE *pFile = NULL;
	
	// correct containers?
	if(particles.empty() || materials.empty())
	{
		LOG4CXX_ERROR(generalOutputLogger, "tried to write empty containers to file "<<filename);
		return E_INVALIDPARAM;
	}
	
#ifdef DEBUG
	// assert types
	for(std::vector<Particle>::const_iterator it = particles.begin(); it != particles.end(); it++)
	{

		// secure index of type is valid, otherwise set it to zero
		if(it->type < 0)
		{
			LOG4CXX_INFO(generalOutputLogger, "invalid index for writing "<<filename<<" found");
			return E_INVALIDPARAM;
		}
		if(it->type >= materials.size())
		{
			LOG4CXX_INFO(generalOutputLogger, "invalid index for writing "<<filename<<" found");
			return E_INVALIDPARAM;
		}
	}
#endif

	// open file (ASCII, replace)
	pFile = fopen(filename, "w"); 
	if(!pFile)return E_UNKNOWN;

	// print description
	fprintf(pFile,	"#\n" \
					"# file format:\n"\
					"# anything after a '#' character is treated as a comment\n"\
					"# empty lines are skipped\n"\
					"# there are two possible chunks that have to exist\n"
					"# a line starting with \"material:\" identifies the material chunk\n"\
					"# where materials have the form\n"\
					"#      materialname epsilon sigma\n"\
					"# e.g. material     1.0     4.0\n"\
					"# next is data chunk, starting with \"data:\"\n"\
					"# to encode single particles in the form\n"\
					"#      xyz-coord,    velocity,     mass, type\n"\
					"# e.g. 1.0 1.0 0.0,  0.0 1.0 0.0,  1.0,  0\n"\
					"# where type is indexing the i-th material described in material chunk\n"\
					"#\n");

	// print materials
	fprintf(pFile, "material:\n");
	for(std::vector<Material>::const_iterator it = materials.begin(); it != materials.end(); it++)
	{
		fprintf(pFile, "%s %lf %lf\n", it->name.c_str(), it->epsilon, it->sigma);
	}

	// print particles...
	fprintf(pFile, "\ndata:\n");
	for(std::vector<Particle>::const_iterator it = particles.begin(); it != particles.end(); it++)
	{
		fprintf(pFile, "%lf %lf %lf, %lf %lf %lf, %lf, %d\n", it->x[0], it->x[1], it->x[2],
			it->v[0], it->v[1], it->v[2], it->m, it->type);
	}
	
	// close, everything o.k.
	fclose(pFile);

	return S_OK;
}