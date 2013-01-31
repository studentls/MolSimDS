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
#include "outputWriter/VTKWriter.h"
#include "outputWriter/XYZWriter.h"
#include "XMLFileReader.h"
#include "TXTFile.h"

using namespace std;


err_type Simulation::Init(const SimulationDesc& desc)
{
	this->desc = desc;

	// clear particles if it is not already empty
	Release();


	return S_OK;
}

err_type Simulation::CreateSimulationFromXMLFile(const char *filename)
{
    // does file exist?
    if(!utils::fileExists(filename))return E_FILENOTFOUND;
    
	XMLFileReader fr;

	if(FAILED(fr.readFile(filename)))return E_FILEERROR;

	desc = fr.getDescription();

	// generate Cache
	desc.generateCache();

	// get container
	Release(); // to free mem

	fr.makeParticleContainer(&particles);

	// set a force calculation method...
	if(desc.potentialForce == PF_SLJ)forceCalculator = Simulation::forceSLJCalculator;
	else if(desc.potentialForce == PF_GRAVITY)forceCalculator = Simulation::forceGravitationCalculator;
	else forceCalculator = Simulation::forceLJCalculator;
	
	// initialize thermostatisticaldata if specified in desc
	if(desc.iterationsTillTStatisticsCalculation > 0)
	{
		this->thermostatistics = new ThermostatisticalData;

		// set id
		int id = 0;
		particles->Iterate(Simulation::idAssigner, (void*)&id);

		// set capacity
		thermostatistics->initialPositions.setsize(particles->getParticleCount());

		// set rdf intervals		
		int countofrdfintervals = desc.cutoffRadius / desc.rdfdelta_t;
		//adjust delta
		thermostatistics->delta_t = desc.cutoffRadius / (double)countofrdfintervals;
		thermostatistics->distanceCounter.setsize(countofrdfintervals);
		for(int i = 0; i < thermostatistics->distanceCounter.size(); i++)thermostatistics->distanceCounter[i] = 0;

		// copy initial positions...
		std::vector<Particle> tmp = particles->getParticles();
		for(std::vector<Particle>::const_iterator it = tmp.begin(); it != tmp.end(); it++)
		{
			assert(it->id >= 0 && it->id < thermostatistics->initialPositions.capacity());
			thermostatistics->initialPositions[it->id] = (it->x);
		}
	}

	// is particles a valid pointer? - If not file has not been parsed correctly...
	if(!particles)return E_FILEERROR;

	// calc initial forces...
	calculateF();

	// is a thermostat present? if yes, initialize it
	if(desc.applyThermostat())if(desc.initThermostat)initializeThermostat();

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


err_type Simulation::Run()
{
	assert(particles);

	// check if the particles are valid
	if(particles->IsEmpty())return E_NOTINITIALIZED;

	// initialize some values
	double current_time = desc.start_time;
	int iteration = 0;

	// set common statistical values...
	statistics.particle_count = particles->getParticleCount();
	statistics.step_count = (desc.end_time - desc.start_time) / desc.delta_t;

	// output that calculation have started ("starting calculation...")
	LOG4CXX_TRACE(simulationLogger, ">> starting calculation...");
	LOG4CXX_TRACE(simulationLogger, "\n   progress\t\t\telapsed time / remaining time\n");

	//start timer
	utils::Timer timer;
	utils::Timer progress_timer; // update progress view

	// additional running variable to finish if events occur
	bool running = true;

#ifndef ICE
	// has the user requested a viewer to run?
	bool viewer = Viewer::Instance().IsRunning();
#endif
	// calculate how many steps will be 1%
	int stepsperpercent = statistics.step_count / 100;

	// iterate until the end time is reached...
	while (current_time < desc.end_time && running) {

		// thermostat present? if yes apply(for iteration 0 thermostat is already set)
		if(desc.applyThermostat())
			if( iteration % desc.iterationsTillThermostatApplication == 0 && iteration != 0)
				adjustThermostat();
		
		if (particles->getType() == PCT_MEMBRANE && iteration < ((MembraneContainer*)particles)->pullIterations)
			particles->Iterate(forceCalculatorMembranePull, (void*)&desc);
		
		// perform one iteration step
		performStep();

		// plot the particles on every desiredth iteration, beginning with the first
		if (iteration % desc.iterationsperoutput == 0) {
			plotParticles(iteration);
		}
		
		// calculate statistical data if wished
		if(desc.iterationsTillTStatisticsCalculation > 0)
		{
			if(iteration % desc.iterationsTillTStatisticsCalculation < desc.iterationsPerTStatisticsCalculation)
				calculateThermostatisticalData((double)iteration * desc.delta_t + desc.start_time);

			// output data
			if(iteration != 0 && iteration % desc.iterationsTillTStatisticsCalculation == 0)writeThermostatisticalDataToFile();
		}

		// output info every 5s
		if(progress_timer.getElapsedTime() > 5.0f)
		{
			stringstream str;

			str<<"   [";

			int max = 15;
			int count = max * iteration / statistics.step_count;

			for(int i = 0; i < count; i++)
				str<<".";
			for(int i = count; i < max; i++)
				str<<" ";
			// calc how many seconds iterations will approximately last from now onwards
			int seconds = (int)(timer.getElapsedTime() / (double)iteration * statistics.step_count);

			str<<"]  \t"<<(int)(100 * iteration / statistics.step_count)<<"%  \t"<<utils::secondsToHMS((int)timer.getElapsedTime())<<" / "<<utils::secondsToHMS(seconds);

			LOG4CXX_TRACE(simulationLogger, str.str().c_str());

			progress_timer.reset();
		}
		
#ifndef ICE	
		// is a viewer active?
		if(viewer)
		{
			// update running dependend to the status of the viewer
			running = Viewer::Instance().IsRunning();

			// if running, notify(send particle data to viewer)
			if(running)notifyViewer();
		}
#endif
		// increment loop values
		iteration++;
		current_time += desc.delta_t;
	}

	// has the user aborted? is runnign set to false?
	if(!running)LOG4CXX_INFO(simulationLogger, ">> simulation aborted by user");

	// generate statistical data
	statistics.time = timer.getElapsedTime();
	statistics.timeperstep = statistics.time / statistics.step_count;

	// write statistics to file
	if(desc.iterationsTillTStatisticsCalculation > 0)	
			// output data
			if(iteration != 0 && iteration % desc.iterationsTillTStatisticsCalculation == 0)writeThermostatisticalDataToFile();


	// output that the output has finished
	cout << endl;
	LOG4CXX_TRACE(simulationLogger, ">> output written. Terminating...\n");

	return S_OK;
}

err_type Simulation::Release()
{	
	// delete thermostatistics
	SAFE_DELETE(thermostatistics);

	// delete Particle data
	if(particles)particles->Clear();

	SAFE_DELETE(particles);

	return S_OK;
}

void Simulation::calculateF() {

	// reassign particles, if LinkedCell Algorithm is used...
	if(particles->getType() == PCT_LINKEDCELL)
		((LinkedCellParticleContainer*)particles)->ReassignParticles();

	// call particles.Iterate() on forceResetter
	particles->Iterate(forceResetter, (void*)&desc);

	// if it's a membrane, calculate the membrane forces
	if (particles->getType() == PCT_MEMBRANE)
		((MembraneContainer*)particles)->ApplyMembraneForces();
	
	// call particles.IteratePairwise() on forceCalculator
	particles->IteratePairwise(forceCalculator, (void*)&desc);

	// call gravityCalculator to add gravity force for each particle
	if(desc.gravitational_constant != 0.0) // bad floating point variable check, but in this case o.k.
		particles->Iterate(gravityCalculator, (void*)&desc);

	// apply Boundary conditions, if LinkedCell Algorithm is used...
	if(particles->getType() == PCT_LINKEDCELL)
		((LinkedCellParticleContainer*)particles)->ApplyBoundaryConditions(forceCalculator, (void*)&desc);

	// apply Boundary conditions, if Membrane Algorithm is used...
	if(particles->getType() == PCT_MEMBRANE)
		((MembraneContainer*)particles)->ApplyReflectiveBoundaryAtBottom(forceCalculator, (void*)&desc);
}

void Simulation::forceResetter(void* data, Particle& p) {
	// call particlesresetForce(), which also ensures that old_f gets changed as well as f
	p.resetForce();
}

void Simulation::forceSLJCalculator(void* data, Particle& p1, Particle& p2)
{
	using namespace utils;

	SimulationDesc *desc = (SimulationDesc*)data;

	// assert indices
	assert(p1.type >= 0);
	assert(p2.type >= 0);
	assert(p1.type < desc->materials.size());
	assert(p2.type < desc->materials.size());


	/// the force will be calculated as
	// F = U_LJ(r) * S'(r) + U_LJ'(r) * S(r)
	// where r = ||x_i - x_j||_2	
	
	double eps24 = desc->getEpsilon24(p1.type, p2.type);
	double sigmaSq = desc->getSigmaSq(p1.type, p2.type);		

	double ULennardJones = 0.0;
	double ULennardJonesGradient = 0.0;
	double Smoothing = 0.0;
	double SmoothingGradient = 0.0;

	double r_c = desc->cutoffRadius;
	double r_l = desc->SLJfactor;
	double eps = eps24 / 24.0;	

	// calc vector	
	Vector<double, 3> dir = p2.x - p1.x;
	double norm = dir.L2Norm();
	double r = norm;
	double invnorm = 1.0 / norm; // inverted norm
	dir = dir * (invnorm); // normalize vector


	// calc U_LJ, S, U_LJ', S'

	// U_LJ
	double invdistSq  = 1.0 / p1.x.distanceSq(p2.x);
	double temp1 = sigmaSq * invdistSq; 
	double pow6 = temp1 * temp1 * temp1;
	double pow12 = pow6 * pow6;
	ULennardJones = 4.0 * eps * (pow12 - pow6);

	// U_LJ'
	ULennardJonesGradient = 24.0 * eps * invnorm * (pow6 - 2.0 * pow12);

	double rdiff = r_c-r_l;

	// S
	if(r < r_l)Smoothing = 1.0;
	else if(r > r_c)Smoothing = 0.0;
	else
	{
		Smoothing = 1.0 - (r - r_l) * (r - r_l) * (3.0 * r_c - r_l - 2.0 * r)/(rdiff*rdiff*rdiff);
	}

	// S'
	//SmoothingGradient = - 2.0 / (rdiff) * (-rdiff + rdiff * r_l / norm + 3.0 * norm - 4.0 * r_l - r_l * r_l / norm) * xDif;
	double rtmp1 = r - r_c;
	double rtmp2 = r_c-r_l;
	SmoothingGradient = 6.0 * rtmp1 * (r - r_l) / (rtmp2 * rtmp2 * rtmp2); 

	// composite force(multiplication rule)
	double force = ULennardJones * SmoothingGradient + ULennardJonesGradient * Smoothing;

	// apply force
	Vector<double, 3> forcevector = dir * force; // direction * force

	// add individual particle to particle force to sum
	p1.addForce(forcevector);
	
	// new function to avoid unnecessary object construction
	p2.substractForce(forcevector);
}

void Simulation::forceGravitationCalculator(void* data, Particle& p1, Particle& p2)
{
	SimulationDesc *desc = (SimulationDesc*)data;

	// using simple force calculation model
	double invdist = 1.0 / p1.x.distance(p2.x);
	
	// use for speed up
	double denominator = invdist * invdist * invdist; 

	double p1m = desc->materials[p1.type].mass;
	double p2m = desc->materials[p2.type].mass;

	double factor = p1m * p2m * denominator * desc->gravitational_constant;
	
	utils::Vector<double, 3> force = (p2.x - p1.x) * factor;
	
	// add individual particle to particle force to sum
	p1.addForce(force);
	p2.substractForce(force);
}	

void Simulation::forceLJCalculator(void* data, Particle& p1, Particle& p2)
{
#ifdef DEBUG
	// assert data is a valid pointer!
	if(!data)
	{
		LOG4CXX_ERROR(simulationLogger, "error: data is not a valid pointer!");
		return;
	}
#endif

	SimulationDesc *desc = (SimulationDesc*)data;

	// assert indices
	assert(p1.type >= 0);
	assert(p2.type >= 0);
	assert(p1.type < desc->materials.size());
	assert(p2.type < desc->materials.size());

	double eps24 = desc->getEpsilon24(p1.type, p2.type);
	double sigmaSq = desc->getSigmaSq(p1.type, p2.type);		

	// optimized version
	double invdistSq  = 1.0 / p1.x.distanceSq(p2.x);
				// cache this for better performance
	double temp1 = sigmaSq * invdistSq; 
	double pow6 = temp1 * temp1 * temp1;
	double pow12 = pow6 * pow6;

	double temp2 = pow6 - 2.0 * pow12;

	double prefactor = eps24 * invdistSq; // cache also here values...
	double factor = prefactor * temp2;

	utils::Vector<double, 3> force = (p2.x - p1.x) * factor;
	
	// add individual particle to particle force to sum
	p1.addForce(force);
	
	// new function to avoid unnecessary object construction
	p2.substractForce(force);
}

void Simulation::forceCalculatorMembranePull(void* data, Particle& p)
{
#ifdef DEBUG
	// assert data is a valid pointer!
	if(!data)
	{
		LOG4CXX_ERROR(simulationLogger, "error: data is not a valid pointer!");
		return;
	}
#endif

	SimulationDesc *desc = (SimulationDesc*)data;

	utils::Vector<double, 3> pull;
	pull[2] = 0.8;

	if (p.type == 1)
		p.addForce(pull);
}

void Simulation::gravityCalculator(void *data, Particle& p)
{
	SimulationDesc *desc = (SimulationDesc*)data;

	// add gravitational force, based on G = m * g
	utils::Vector<double, 3> grav_force;
	// act on the last dimension
	int dimensionToAffect = desc->dimensions - 1;

	// only the last dimension is affected...
	grav_force[dimensionToAffect] = desc->materials[p.type].mass * desc->gravitational_constant;

	p.addForce(grav_force);
}

void Simulation::calculateX() {
	
	assert(particles);

	// call particles.Iterate() on posCalculator
	particles->Iterate(posCalculator, (void*)&desc);
}

void Simulation::posCalculator(void* data, Particle& p) {
	#ifdef DEBUG
	// assert data is a valid pointer!
	if(!data)
	{
		LOG4CXX_ERROR(simulationLogger, "error: data is not a valid pointer!");
		return;
	}
#endif

	SimulationDesc *desc = (SimulationDesc*)data;
	
	// calc new pos
	// base calculation on Velocity-Störmer-Verlet-Algorithm
	// x_i ( t^{n+1} ) = x_i (t^n) + delta_t * p_i (t^n) + delta_t^2 * f_i (t^n) * 0.5 / m_i
	//p.x = p.x + desc->delta_t * p.v + desc->delta_t * desc->delta_t * p.getF() * 0.5 * (1.0 / p.m);
	// Cache here!
	p.x = p.x + desc->delta_t * p.v + desc->delta_t * desc->delta_t * p.getF() * desc->getHalfInvMass(p.type);
}

void Simulation::calculateV() {
	// call particles.Iterate() on velCalculator
	particles->Iterate(velCalculator, (void*)&desc);
}

void Simulation::velCalculator(void* data, Particle& p) {
	#ifdef DEBUG
	// assert data is a valid pointer!
	if(!data)
	{
		LOG4CXX_ERROR(simulationLogger, "error: data is not a valid pointer!");
		return;
	}
#endif

	SimulationDesc *desc = (SimulationDesc*)data;

	// calc new vel
	// base calculation on Velocity-Störmer-Verlet-Algorithm
	// v_i ( t^{n+1} ) = v_i(t^n) + dt * (F_i(t^n) + F_i(T^{n+1} ) ) / 2m_i
	//p.v = (p.v +  desc->delta_t * (p.getF() + p.getOldF() ) * 0.5 * (1.0 / desc->materials[p.type].mass));
	p.v = p.v +  desc->delta_t * (p.getF() + p.getOldF() ) * desc->getHalfInvMass(p.type);
}

void Simulation::plotParticles(int iteration) {

	assert(particles);
	
	// switch between VTK and XYZ output
	// depending on the value of desc.output_fmt
	switch(desc.output_fmt)
	{
	case SOF_VTK:
		{
			// VTK Output
			outputWriter::VTKWriter writer;
			writer.plotParticles(particles->getParticles(), desc.materials, desc.outname, iteration);
			break;
		}
	case SOF_XYZ:
		{
			// XYZ Output
			outputWriter::XYZWriter writer;
			writer.plotParticles(particles->getParticles(), desc.outname, iteration);
			break;
		}
	case SOF_TXT:
		{
			// set txt output
			TXTFile file;
			stringstream strstr;
			strstr<<desc.outname << "_" << (iteration < 10 ? "000" : (iteration < 100 ? "00" : ( iteration < 1000 ? "0" : "") )) << iteration << ".txt";
			file.writeFile(strstr.str().c_str(), particles->getParticles(), desc.materials);
		}
	default:
		{
		}
	}
}

void Simulation::kineticEnergyCalculator(void *data, Particle& p)
{
	SimulationDesc *desc = (SimulationDesc*)data;

	// calc according to
	// Task 1, Sheet 4
	desc->kineticEnergy += desc->materials[p.type].mass * p.v.L2NormSq() * 0.5;
}

void Simulation::totalMassCalculator(void *data, Particle& p)
{
	SimulationDesc *desc = (SimulationDesc*)data;

	desc->totalMass += desc->materials[p.type].mass;
}

void Simulation::setBrownianMotionCalculator(void *data, Particle& p)
{
	double brownianMotionFactor = ((SimulationDesc*)data)->brownianMotionFactor;
	unsigned int dim = ((SimulationDesc*)data)->dimensions;

	// apply brownian Motion
	MaxwellBoltzmannDistribution(p, brownianMotionFactor, dim);
}

void Simulation::applyTemperatureScalingFactor(void *data, Particle& p)
{
	double beta = *((double*)data);
	p.v = p.v * beta;
}

void Simulation::initializeThermostat()
{
	assert(particles);

	// calc kinetic Energy of the whole system
	// temperature is start temperature
	// currently we use a dimensionless method
	double EnergyProposed = 0.5 * desc.dimensions * particles->getParticleCount() *desc.temperature;

	// we need the total mass of all particles...
	
	desc.totalMass = 0.0;
	particles->Iterate(totalMassCalculator, (void*)&desc);
	double totalMass = desc.totalMass;
	
	// calculate inital velocity that is set for each particle
	desc.brownianMotionFactor = sqrt(2.0 * EnergyProposed / (desc.dimensions * totalMass));

	// set brownianMotion to every particle
	particles->Iterate(setBrownianMotionCalculator, (void*)&desc);

}

void Simulation::adjustThermostat()
{
	assert(particles);
	
	// inc temperature
	desc.temperature += desc.temperatureStepSize;

	// secure that temperature is not above target temperature
	if(desc.temperature > desc.targetTemperature)desc.temperature = desc.targetTemperature;
	// calc beta

	//get energy of the system
	desc.kineticEnergy = 0.0;	
	particles->Iterate(kineticEnergyCalculator, (void*)&desc);
	double EnergyAtTheMoment = desc.kineticEnergy;
	
	// proposed energy
	double EnergyProposed = 0.5 * desc.dimensions * particles->getParticleCount() *desc.temperature;

	// beta
	double beta = sqrt(EnergyProposed / EnergyAtTheMoment);

	particles->Iterate(applyTemperatureScalingFactor, (double*)&beta);

	
}

void Simulation::diffusionCalculator(void* data, Particle& p)
{
	// cast
	ThermostatisticalData *tdata = (ThermostatisticalData*)data;

	assert(p.id >= 0 && p.id < tdata->initialPositions.size());

	// add to sum
	tdata->diffusionValues.last()[1] += p.x.distanceSq(tdata->initialPositions[p.id]);
}

void Simulation::rdfCalculator(void* data, Particle& p1, Particle& p2)
{
	// cast
	ThermostatisticalData *tdata = (ThermostatisticalData*)data;

	// find right interval
	double distance = p1.x.distance(p2.x);

	// sort into right interval
	int index = distance / tdata->delta_t;

	if(index < 0)index = 0;
	if(index >= tdata->distanceCounter.size())index = tdata->distanceCounter.size() - 1;

	assert(index >= 0 && index < tdata->distanceCounter.size());

	tdata->distanceCounter[index]++;

}

void Simulation::calculateThermostatisticalData(const double t)
{
	// clear halo!
	if(particles->getType() == PCT_LINKEDCELL)((LinkedCellParticleContainer*)particles)->clearHaloParticles();
	
	// calc diffusion

	// add zero sum...
	this->thermostatistics->diffusionValues.push_back(utils::Vector<double, 2>(t, 0.0));
	// use iterate
	particles->Iterate(Simulation::diffusionCalculator, (void*)thermostatistics);

	// divide by N-1!!! this is a mistake in the formula, as 1/N is not an unbiased estimator!
	thermostatistics->diffusionValues.last()[1] /= (double)(thermostatistics->initialPositions.size() - 1);

	// calc rdf
	for(int i = 0; i < thermostatistics->distanceCounter.size(); i++)thermostatistics->distanceCounter[i] = 0;
	particles->IteratePairwise(Simulation::rdfCalculator, (void*)thermostatistics);

	// calc local density
	for(int i = 0; i < thermostatistics->distanceCounter.size(); i++)
	{
		double a = i * thermostatistics->delta_t;
		a = a * a * a;
		double b = (i+1) * thermostatistics->delta_t;
		b = b * b * b;
		double density = (double)thermostatistics->distanceCounter[i] / (4.0 * PI / 3.0 * (b - a));
		thermostatistics->rdfValues.push_back(utils::Vector<double, 2>(t, density));
	}
}

err_type Simulation::writeThermostatisticalDataToFile()
{
	using namespace utils;

	std::string filename = desc.outname +"_thermostatistics.csv";
	
	// open file
	FILE *pFile = fopen(filename.c_str(), "w");
	if(!pFile)
	{
		LOG4CXX_ERROR(generalOutputLogger, ">> error: failed to write thermostatistical data!");
		return E_FILEERROR;
	}

	// print header
	fprintf(pFile, "steps, point of time, diffusion");

	// add columns for rdf
	for(int j = 0; j < thermostatistics->rdfValues.size(); j++)fprintf(pFile, ", rdf%d", j);
	fprintf(pFile, "\n");

	// calc num diffusion pairs
	int count = thermostatistics->diffusionValues.size() / desc.iterationsPerTStatisticsCalculation;

	// go through and average
	for(int i = 0; i < count; i++)
	{
		Vector<double, 2> avg = Vector<double, 2>(0, 0);
		TFastArray<double> avgrdf;
		int numrdfsamples = thermostatistics->distanceCounter.size();
		avgrdf.setsize(numrdfsamples);		
		for(int k = 0; k < numrdfsamples; k++)avgrdf[k] = 0;

		for(int j = 0; j < desc.iterationsPerTStatisticsCalculation; j++)
		{
			avg = avg + thermostatistics->diffusionValues[i * desc.iterationsPerTStatisticsCalculation + j];

			for(int k = 0; k < numrdfsamples; k++)
			{
				avgrdf[k] += thermostatistics->rdfValues[(i * desc.iterationsPerTStatisticsCalculation + j) * numrdfsamples + k][1];
			}
		}
		avg[0] /= (double)desc.iterationsPerTStatisticsCalculation;
		avg[1] /= (double)desc.iterationsPerTStatisticsCalculation;

		for(int k = 0; k < numrdfsamples; k++)avgrdf[k] /= (double)desc.iterationsPerTStatisticsCalculation;

		// write to file
		fprintf(pFile, "%d, %lf, %lf", (int)((avg[0] - desc.start_time)/ desc.delta_t), avg[0], avg[1]);

		// add columns for rdf
		for(int j = 0; j < numrdfsamples; j++)fprintf(pFile, ", %lf", avgrdf[j]);
		fprintf(pFile, "\n");
	}

	fclose(pFile);

	return S_OK;
}

void Simulation::notifyViewer()
{
	using namespace std;
#ifndef ICE
	// notify viewer
		
	VParticle *particles = NULL;

	// wish to store 500 particles
	int count = this->particles->getParticleCount();

	count = Viewer::Instance().LockParticles(&particles, count);
	
	int pos = 0;

	// add particles from particle container
	vector<Particle> vp = this->particles->getParticles();
	if(!vp.empty())
		for(vector<Particle>::iterator it = vp.begin(); it != vp.end(); it++)
		{
			particles[pos].x = it->x[0];
			particles[pos].y = it->x[1];
			particles[pos].z = it->x[2];
			particles[pos].type = it->type;
			pos++;
		}

	Viewer::Instance().UnlockParticles(&particles);
#endif
}


void Simulation::idAssigner(void* data, Particle& p)
{
	int *id =(int*)data;

	// set id
	p.id = *id;

	// inc temp variable
	(*id)++;
}