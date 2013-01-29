//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File XMLFileReader.cpp
// contains class XMLFileReader
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------

#include "XMLFileReader.h"

#include "LinkedCellParticleContainer.h"
#include "ListParticleContainer.h"
#include "ParticleGenerator.h"
#include "TXTFile.h"

//delete later
static const double maxsigma = 1.0; 


err_type XMLFileReader::readFile(const char *filename, bool validate)
{
	// exists file?
	if(!utils::fileExists(filename))return E_FILENOTFOUND;

	// parse like in the Hello example from XSD
	
	// try - catch block to avoid mistakes...
	try
	{
		::xml_schema::flags flags;
		flags = !validate ? ::xml_schema::flags::dont_validate : 0;

		file = std::auto_ptr<simulationfile_t>(simulationfile(filename, flags));
	
	}
	catch(const xml_schema::exception& e)
	{
	  LOG4CXX_ERROR(generalOutputLogger, e.what());
	  return E_FILEERROR;
	}
	
	// is there a brownian Motion factor encoded? if yes set to the value, otherwise to zero
	if(file->params().brownianMotionFactor().present())
		desc.brownianMotionFactor	= file->params().brownianMotionFactor().get();
	else desc.brownianMotionFactor = 0.0;

	desc.delta_t				= file->params().delta_t();
	desc.end_time				= file->params().t_end();

	
	//compare strings
	if(strcmp(file->params().outputfmt().c_str(), "VTK") == 0)desc.output_fmt = SOF_VTK;
	else if(strcmp(file->params().outputfmt().c_str(), "XYZ") == 0)desc.output_fmt = SOF_XYZ;
	else if(strcmp(file->params().outputfmt().c_str(), "TXT") == 0)desc.output_fmt = SOF_TXT;
	else desc.output_fmt = SOF_NONE;
	
	desc.start_time				= file->params().t_start();

	desc.iterationsperoutput	= file->params().iterationsperoutput();

	desc.outname				= file->params().output();
	
	// add a default material
	desc.materials.push_back(createDefaultMaterial());

	// go through materials, if no material exists, create default material
	for( ::data_t::material_const_iterator it = file->data().material().begin();
		it != file->data().material().end(); it++)
	{
		Material mat;
		// get values
		mat.epsilon = it->epsilon().get();
		mat.sigma = it->sigma().get();
		mat.name = it->name().get();
		mat.mass = it->mass().get();

		// add
		desc.materials.push_back(mat);
	}
	
	// thermostat
	// is timeStepsTillThermostatApplication and the rest of factors present?
	if(file->params().iterationsTillThermostatApplication().present() && file->params().initialTemperature().present()
		&& file->params().targetTemperature().present() && file->params().temperatureStepSize().present())
	{
		desc.iterationsTillThermostatApplication = file->params().iterationsTillThermostatApplication().get();
		// set current temperature to initial temperature
		// currently dimensionless!
		desc.temperature = file->params().initialTemperature().get();
		desc.targetTemperature = file->params().targetTemperature().get();
		desc.temperatureStepSize = file->params().temperatureStepSize().get();
	}
	else
	{
		// no thermostat application...
		desc.iterationsTillThermostatApplication = 0;
		desc.temperature = 0;
		desc.targetTemperature = 0;
		desc.temperatureStepSize = 0;
	}
	
	// dimension
	desc.dimensions = file->params().dimension();

	// assert dimension is 2 or 3
	if(desc.dimensions != 2 && desc.dimensions != 3)
		LOG4CXX_ERROR(generalOutputLogger, "unsupported dimension, please use only 2D or 3D");

	// set gravity factor
	if(file->params().gravity().present())
		desc.gravitational_constant = file->params().gravity().get();
	else desc.gravitational_constant = 0.0;

	//file parsed...
	fileParsed = true;

	// rest will be done in makeParticleContainer...

	return S_OK;
}

err_type XMLFileReader::makeParticleContainer(ParticleContainer **out)
{
	using namespace std;

	// Checks...
	if(!fileParsed)return E_UNKNOWN;

	if(*out)return E_INVALIDPARAM;

	
	// collect Data...

	// big particle container to store everything
	vector<Particle> particles;

	
	// at the moment inputfile is a simple particle list, materials have to be defined in the xml file separately!
	// go through input files...
	for(::data_t::inputfile_const_iterator it = file->data().inputfile().begin();
		it != file->data().inputfile().end(); it++)
	{
		// use new txtfile
		// at the moment we simply ignore the material chunk
		TXTFile txt;
		if(FAILED(txt.readFile(it->c_str())))
		{
			LOG4CXX_ERROR(generalOutputLogger, "failed to read inputfile");
			return E_FILEERROR;
		}

		// test for correct particle index
		for(std::vector<Particle>::const_iterator it = txt.getParticles().begin(); it != txt.getParticles().end(); it++)
		{
			if(it->type < 0){LOG4CXX_ERROR(generalOutputLogger, "invalid index!");return E_FILEERROR;}
			if(it->type >= desc.materials.size()){LOG4CXX_ERROR(generalOutputLogger, "material not known for type = "<<it->type);return E_FILEERROR;}
		}
		
		// add particles...
		particles.insert(particles.begin(), txt.getParticles().begin(), txt.getParticles().end());
		
	}

	// go through particles...
	for(::data_t::particle_const_iterator it = file->data().particle().begin();
		it != file->data().particle().end(); it++)
	{
		Particle p;
		p.x[0] = it->X().at(0);
		p.x[1] = it->X().at(1);
		p.x[2] = (it->X().size() > 2) ? it->X().at(2) : 0;

		p.v[0] = it->V().at(0);
		p.v[1] = it->V().at(1);
		p.v[2] = (it->V().size() > 2) ? it->V().at(2) : 0;
		
		// set type to corresponding material
		if(it->material().present())
		{
			int index = getMaterialIndex(desc, it->material().get().c_str());
			p.type = index >= 0 ? index : 0;
		}
		else p.type = 0; // set to default material

		// add to vector
		particles.push_back(p);
	}

	// go through cuboids...
	for(::data_t::cuboid_const_iterator it = file->data().cuboid().begin();
		it != file->data().cuboid().end(); it++)
	{
		ListParticleContainer pc;

		// construct variables
		utils::Vector<double, 3> corner;
		utils::Vector<double, 3> v;
		double h;
		utils::Vector<unsigned int, 3> dim;

		corner[0] = it->X().at(0);
		corner[1] = it->X().at(1);
		corner[2] = (it->X().size() > 2) ? it->X().at(2) : 0;

		v[0] = it->V().at(0);
		v[1] = it->V().at(1);
		v[2] = (it->V().size() > 2) ? it->V().at(2) : 0;
		
		h = it->h();

		dim[0] = it->N().at(0);
		dim[1] = it->N().at(1);
		dim[2] = (it->N().size() > 2) ? it->N().at(2) : 1;


		int type = 0;

		// set type to corresponding material if avaliable
		if(it->material().present())
		{
			int index = getMaterialIndex(desc, it->material().get().c_str());
			type = index >= 0 ? index : 0;
		}

		// make Cuboid and add to particle
		ParticleGenerator::makeCuboid(pc, corner, dim, h, v, type, desc.brownianMotionFactor);
		vector<Particle> temp = pc.getParticles();
		particles.insert(particles.begin(), temp.begin(), temp.end());
	
	}

	// go through spheres...
	for(::data_t::sphere_const_iterator it = file->data().sphere().begin();
		it != file->data().sphere().end(); it++)
	{
		ListParticleContainer pc;

		// construct variables
		utils::Vector<double, 3> center;
		utils::Vector<double, 3> v;
		double m;
		double h;
		unsigned int dim;
		unsigned int radius;

		center[0] = it->X().at(0);
		center[1] = it->X().at(1);
		center[2] = (it->X().size() > 2) ? it->X().at(2) : 0;

		v[0] = it->V().at(0);
		v[1] = it->V().at(1);
		v[2] = (it->V().size() > 2) ? it->V().at(2) : 0;
		
		h = it->h();

		dim = it->dimensions();

		radius = it->r();

		int type = 0;

		// set type to corresponding material if avaliable
		if(it->material().present())
		{
			int index = getMaterialIndex(desc, it->material().get().c_str());
			type = index >= 0 ? index : 0;
		}

		// make Cuboid and add to particle
		ParticleGenerator::makeSphere(pc, center, v,radius, h, dim, desc.brownianMotionFactor, type);
		vector<Particle> temp = pc.getParticles();
		particles.insert(particles.begin(), temp.begin(), temp.end());
	
	}


	// now construct Particle Container...
	ParticleContainer *container = NULL;

	// is a LinkedCell Container present?
	if(file->params().algorithm().LinkedCell().present())
	{
		LinkedCell_t& lc = file->params().algorithm().LinkedCell().get();

		// set dimensions to 2D if z is set to zero
		unsigned int dim = lc.sizeofdomainZ() == 0 ? 2 : 3;

		double cutoffDistance = lc.cutoff_radius();

		Vec3 frontLowerLeftCorner;
		Vec3 simulationAreaExtent;
		
		frontLowerLeftCorner[0] = lc.offset().at(0);
		frontLowerLeftCorner[1] = lc.offset().at(1);
		frontLowerLeftCorner[2] = lc.offset().at(2);

		simulationAreaExtent[0] = lc.sizeofdomainX();
		simulationAreaExtent[1] = lc.sizeofdomainY();
		simulationAreaExtent[2] = lc.sizeofdomainZ();


		// for LinkedCell, optional boundary conditions can be specified
		unsigned int bc = BC_NONE;

		// are any conditions contained?
		if(lc.conditions().present())
		{
			// go through conditions...
			for(::conditions_t::condition_const_iterator it = lc.conditions().get().condition().begin();
				it != lc.conditions().get().condition().end(); it++)
			{

				std::string val = it->value();

				// is a periodic boundary present?
				if(it->type().present())
				{
					if(strcmp(it->type().get().c_str(), "periodic") == 0)
					{
						// set flags
						if(strcmp(val.c_str(), "left")		== 0)bc	|= BC_PERIODIC_XAXIS;
						if(strcmp(val.c_str(), "right")		== 0)bc	|= BC_PERIODIC_XAXIS;
						if(strcmp(val.c_str(), "top")		== 0)bc	|= BC_PERIODIC_ZAXIS;
						if(strcmp(val.c_str(), "bottom")	== 0)bc	|= BC_PERIODIC_ZAXIS;
						if(strcmp(val.c_str(), "front")		== 0)bc	|= BC_PERIODIC_YAXIS;				
						if(strcmp(val.c_str(), "back")		== 0)bc	|= BC_PERIODIC_YAXIS;
					}
					else if(strcmp(it->type().get().c_str(), "reflective") == 0)
					{
						// set flags
						if(strcmp(val.c_str(), "left")		== 0)bc	|= BC_LEFT;
						if(strcmp(val.c_str(), "right")		== 0)bc	|= BC_RIGHT;
						if(strcmp(val.c_str(), "top")		== 0)bc	|= BC_TOP;
						if(strcmp(val.c_str(), "bottom")	== 0)bc	|= BC_BOTTOM;
						if(strcmp(val.c_str(), "front")		== 0)bc	|= BC_FRONT;				
						if(strcmp(val.c_str(), "back")		== 0)bc	|= BC_BACK;
					}
					else
						LOG4CXX_ERROR(generalOutputLogger, "unkown boundary type");
				}
				else
				{			
					// otherwise use reflective boundaries as default

					// set flags
					if(strcmp(val.c_str(), "left")		== 0)bc	|= BC_LEFT;
					if(strcmp(val.c_str(), "right")		== 0)bc	|= BC_RIGHT;
					if(strcmp(val.c_str(), "top")		== 0)bc	|= BC_TOP;
					if(strcmp(val.c_str(), "bottom")	== 0)bc	|= BC_BOTTOM;
					if(strcmp(val.c_str(), "front")		== 0)bc	|= BC_FRONT;				
					if(strcmp(val.c_str(), "back")		== 0)bc	|= BC_BACK;							
				}
				
				// special cases
				if(strcmp(val.c_str(), "all")		== 0)bc	= BC_ALL;
				if(strcmp(val.c_str(), "none")		== 0)bc	= BC_NONE;


				if(strcmp(val.c_str(), "outflow")		== 0)bc	|= BC_OUTFLOW;
			}

		}

		

		container = new LinkedCellParticleContainer(dim, particles, cutoffDistance, frontLowerLeftCorner,
			simulationAreaExtent, bc, maxsigma);
		
		(*out) = container;

	}
	else if(file->params().algorithm().Membrane().present())
	{
		// use the membrane list
		Membrane_t& lc = file->params().algorithm().Membrane().get();
		unsigned int pullIterations = lc.pull_iterations();;
		container = new MembraneContainer(pullIterations);

		// hardcoded contents for the container
		utils::Vector<unsigned int, 2> dimensions(50, 50);
		Vec3 lowerLeftFrontCorner(15.0, 15.0, 1.5);
		((MembraneContainer*)container)->SetMembrane(lowerLeftFrontCorner, dimensions, 2.2);

		(*out) = container;
	}
	else
	{
		// use in this case, the simple ListParticleContainer
		container = new ListParticleContainer(particles);

		(*out) = container;
	}


	// secure check, see if types are all valid!
	std::vector<Particle> tmp = (*out)->getParticles();
	if(tmp.empty())return S_OK;
	else
	{
		// go through particles and check type
		for(std::vector<Particle>::const_iterator it = tmp.begin(); it != tmp.end(); ++it)
		{
			if(it->type < 0 || it->type >= desc.materials.size())
			{
				LOG4CXX_ERROR(generalOutputLogger, ">> error: invalid type found! total failure!");
				return E_UNKNOWN;
			}
		}
	}

	return S_OK;
}
