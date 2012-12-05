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

	
	//set desc
	desc.brownianMotionFactor	= file->params().brownianMotionFactor();
	desc.delta_t				= file->params().delta_t();
	desc.end_time				= file->params().t_end();
	desc.epsilon				= file->params().epsilon();
	
	//compare strings
	if(strcmp(file->params().outputfmt().c_str(), "VTK") == 0)desc.output_fmt = SOF_VTK;
	else if(strcmp(file->params().outputfmt().c_str(), "XYZ") == 0)desc.output_fmt = SOF_XYZ;
	else desc.output_fmt = SOF_NONE;

	desc.sigma					= file->params().sigma();
	desc.start_time				= file->params().t_start();

	desc.iterationsperoutput	= file->params().iterationsperoutput();

	desc.outname				= file->params().output();
	
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

	// go through input files...
	for(::data_t::inputfile_const_iterator it = file->data().inputfile().begin();
		it != file->data().inputfile().end(); it++)
	{
		// quick read per ListParticleContainer
		// better make a new class for reading .txt files
		ListParticleContainer pc;
		pc.AddParticlesFromFileNew(it->c_str());

		vector<Particle> temp = pc.getParticles();
		particles.insert(particles.begin(), temp.begin(), temp.end());
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

		p.m = it->m();

		// type is not supported yet...

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
		double m;
		double h;
		utils::Vector<unsigned int, 3> dim;

		corner[0] = it->X().at(0);
		corner[1] = it->X().at(1);
		corner[2] = (it->X().size() > 2) ? it->X().at(2) : 0;

		v[0] = it->V().at(0);
		v[1] = it->V().at(1);
		v[2] = (it->V().size() > 2) ? it->V().at(2) : 0;

		m = it->m();

		h = it->h();

		dim[0] = it->N().at(0);
		dim[1] = it->N().at(1);
		dim[2] = (it->N().size() > 2) ? it->N().at(2) : 1;


		// make Cuboid and add to particle
		ParticleGenerator::makeCuboid(pc, corner, dim, h, m, v, desc.brownianMotionFactor);
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

		m = it->m();

		h = it->h();

		dim = it->dimensions();

		radius = it->r();

		// make Cuboid and add to particle
		ParticleGenerator::makeSphere(pc, center, v, m, radius, h, dim, desc.brownianMotionFactor);
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


		container = new LinkedCellParticleContainer(dim, particles, cutoffDistance, frontLowerLeftCorner,
			simulationAreaExtent, false, false, false, false, false, false, desc.sigma);
		
		(*out) = container;

	}
	else
	{
		// use in this case, the simple ListParticleContainer
		container = new ListParticleContainer(particles);

		(*out) = container;
	}

	return S_OK;
}