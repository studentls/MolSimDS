//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File XMLFileReader.h
// contains class XMLFileReader
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------

#ifndef XMLFILEREADER_HEADER_
#define XMLFILEREADER_HEADER_

#include "utils/utils.h"

#include "Simulation.h"
#include "ParticleContainer.h"

// generated definitions
#include "XMLFile/simulationfile.h"
#include <memory>
#include <iostream>


class XMLFileReader
{
private:
	SimulationDesc					desc;
	
	std::auto_ptr<simulationfile_t> file;
	bool							fileParsed;

	/// creates a simple default material
	/// @return a default material
	inline Material						createDefaultMaterial()
	{
		Material mat;

		mat.name = "default";
		mat.epsilon = 1.0;
		mat.sigma = 1.0;
		return mat;
	}

	/// function to search for the index in the materials' array of a simulation description
	/// @param strMatName the materialname to search for
	/// @param desc a simulation description
	/// @return returns the index or -1 if a material with name str is not found
	int  getMaterialIndex(const SimulationDesc& desc, const std::string & strMatName)
	{
		int i = 0;

		// empty?
		if(desc.materials.empty())return -1;

		for(std::vector<Material>::const_iterator it = desc.materials.begin();
			it != desc.materials.end(); it++)
		{
			// strcmp and return index if match
			if(strMatName.compare(it->name) == 0)return i;

			i++;
		}

		// return -1 if not found
		return -1;
	}
public:
	XMLFileReader() : fileParsed(false)	{}

	/// reads an XML file
	/// @param filename filename of .xml file
	/// @param validate if true, the file will be validated against simulationfile.xsd
	/// @return returns S_OK on success or
	/// E_FILENOTFOUND if the file is not avaliable
	/// or E_FILEERROR if validation fails
	err_type readFile(const char *filename, bool validate = false);

	/// getter
	/// @return returns description of simulation(all simulation params)
	SimulationDesc getDescription()	{return desc;}


	/// outputs a particle container with all information that could be retrieved from the file
	/// @param out a pointer reference where XMLFileReader will allocate the new container to
	/// @return returns S_OK or E_OUTOFMEMORY if allocating mem failed
	/// or E_UNKOWN if out is not a NULL-Pointer
	/// don't forget to delete what is given to out!
	err_type makeParticleContainer(ParticleContainer **out);

};

#endif
