//------------------------------------------------------------------------------------------------
// File LinkedCellParticleContainer.h
// contains class LinkedCelltParticleContainer
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------


#ifndef LINKEDCELL_PARTICLE_CONTAINER_H_
#define LINKEDCELL_PARTICLE_CONTAINER_H_

#include <vector>
#include "Logging.h"
#include "Particle.h"
#include "ParticleContainer.h"

/// a class that is used to store Particles and iterate over them
template<int dim> class LinkedCellParticleContainer : public ParticleContainer {

private:

	// TODO:
	// let's add the methods applyVelocity(...) and applyForce(...) as well (to the general interface ParticleContainer)
	// they do the same thing and take the same parameters as iterate and iteratePairwaise,
	// but in the case of this class, they also call additional functions to deal with the borders
	// and to reassign the particles to cells

	/// Cutoff distance
	double								cutoffDistance;

	/// vector containing X, Y, (Z) size of each cell
	utils::Vector<double, dim>			cellSize;

	/// vector containing count of cells in X, Y, (Z) direction
	utils::Vector<unsigned int, dim>	cellCount;

	/// get total cell count
	inline unsigned int					getCellCount()
	{
		unsigned int sum = 1;
		for(int i = 0; i < dim; i++)
			sum *= cellCount[i];
		return sum;
	}

	/// the Linked Cell algorithm establishes a grid, which is located in space
	/// via the frontLowerLeftCorner, (v beeing the corner)
	/// |-----|-----|-----|-----|-----|
	/// |     |     |     |     |     |
	/// |-----|-----|-----|-----|-----|
	/// |     |     |     |     |     |
	/// v-----|-----|-----|-----|-----|
	utils::Vector<double, dim>			frontLowerLeftCorner;

	/// a dynamic array, containing all Cells encoded in a 1D array
	std::vector<Particle>				*Cells;
	
	/// ???
	// Note that the halo acts only as a storage for the sake of completeness. The field could be dropped as well
	std::vector<Particle> halo;
	
	/// Array of Pairs of Cell Indices, e.g. (1, 2) is the pair adressing Cell 1 and Cell 2
	/// where 1, 2 are the index of the Cells array
	std::vector<utils::Vector<unsigned int, 2>>	cellPairs;

	/// updating cells every iteration is not very appropriate, better wait some
	/// iterations and perform then update
	/// this variable stores how many iterations the container should wait
	int iterationsPerParticleToCellReassignment;
	
	/// count of iterations so far, note this variable is incremented every time Iterate or (!) IteratePairwise is called!
	int	iterationCount;
	
	/// the distance below which reflective boundaries reflect
	double reflectiveBoundaryDistance;

	// each double in the Vector means something different
	// 0, int: the cell
	// 1, int: the axis
	// 2, double: the direction. -1.0 if the border is on the lower side, 1.0 on the positive side
	// 3, double: the border's coordinate on the given axis, in the given direction
	/// store what cells have reflective boundaries and of what type
	std::vector<utils::Vector<double, 4>>	reflectiveBoundaryCells;

	/// helper function to convert fast 2D indices to 1D based on cellCount
	/// note that indices should be asserted!
	inline unsigned int Index2DTo1D(unsigned int x, unsigned int y)
	{
		return x + cellCount[0] * y;
	}

	/// helper function to convert fast 3D indices to 1D based on cellCount
	/// note that indices should be asserted!
	inline unsigned int Index3DTo1D(unsigned int x, unsigned int y, unsigned int z)
	{
		return x + cellCount[0] * (y + z * cellCount[1]);
	}

	/// helper function to calculate the extent of the simulation area
	/// @return returns the extent of the Simulation area as a 3D vector, note that even if a 2D Grid is used a 3D Vector will be returned
	utils::Vector<double, 3> calcSimulationAreaExtent()
	{
		utils::Vector<double, 3> res;

		// calc based on number of cells in each dimension and their individual extent
		for(int i = 0; i < dim; i++)
		{
			res[i] = cellCount[i] * cellSize[i];
		}

		return res;
	}

	/// generate Index Pairs according to cellCount
	// TODO: write a testcase about this function
	void		generatePairs()	
	{
		// temp pair variable
		utils::Vector<unsigned int, 2> pair;
		
		// empty array?
		if(!cellPairs.empty())cellPairs.clear();

		// 2D neighborhood
		if (dim == 2)
		{
			// go through all x, y
			for (int x = 0; x < cellCount[0]; x++)
				for (int y = 0; y < cellCount[1]; y++)

					// TODO: comment
					for (int xp = 0; xp < 2; xp++)
						for (int yp = 0; yp < 2; yp++)
							if (x + xp < cellCount[0] && y + yp < cellCount[1])
							{
							
								// make pairs ((x,y), (x + xp, y + yp))
								pair[0] = Index2DTo1D(x, y);
								pair[1] = Index2DTo1D(x + xp, y + yp);
								
								cellPairs.push_back(pair);
							}
		}
		// 3D neighborhood
		else if (dim == 3)
		{
			//go through all x, y, z
			for (int x = 0; x < cellCount[0]; x++)
				for (int y = 0; y < cellCount[1]; y++)
					for (int z = 0; z < cellCount[2]; z++)

						// TODO: comment
						for (int xp = 0; xp < 2; xp++)
							for (int yp = 0; yp < 2; yp++)
								for (int zp = 0; zp < 2; zp++)
									if (x + xp < cellCount[0] && y + yp < cellCount[1] && z + zp < cellCount[2])
									{										
										// make pairs ((x,y,z), (x + xp, y + yp, z + zp))
										pair[0] = Index3DTo1D(x, y, z);
										pair[1] = Index3DTo1D(x + xp, y + yp, z + zp);
										cellPairs.push_back(pair);
									}
		}
	}
	
	/// define the reflective boundary cells
	void	SetReflectiveBoundaries(bool leftReflectiveBoundary, bool rightReflectiveBoundary,
									bool frontReflectiveBoundary, bool backReflectiveBoundary,
									// these two will be ignored in the two-dimensional case
									bool bottomReflectiveBoundary, bool topReflectiveBoundary)
	{
		// TODO: testing and debugging
		if (dim == 2) {
			if (leftReflectiveBoundary || rightReflectiveBoundary)
				for (int y = 0; y < cellCount[1]; y++) {
					if (leftReflectiveBoundary) {
						// temp variable
						utils::Vector<double, 4> vec;
						vec[0] = Index2DTo1D(0, y);
						vec[1] = 0;
						vec[2] = -1.0;
						vec[3] = frontLowerLeftCorner[0];
						reflectiveBoundaryCells.push_back(vec);
					}
					if (rightReflectiveBoundary) {
						// temp variable
						utils::Vector<double, 4> vec;
						vec[0] = Index2DTo1D(cellCount[0] - 1, y);;
						vec[1] = 0;
						vec[2] = 1.0;
						vec[3] = frontLowerLeftCorner[0] + calcSimulationAreaExtent()[0];
						reflectiveBoundaryCells.push_back(vec);
					}
				}
			if (frontReflectiveBoundary || backReflectiveBoundary)
				for (int x = 0; x < cellCount[0]; x++) {
					if (frontReflectiveBoundary) {
						// temp variable
						utils::Vector<double, 4> vec;
						vec[0] = Index2DTo1D(x, 0);
						vec[1] = 1;
						vec[2] = -1.0;
						vec[3] = frontLowerLeftCorner[1];
						reflectiveBoundaryCells.push_back(vec);
					}
					if (backReflectiveBoundary) {
						// temp variable
						utils::Vector<double, 4> vec;
						vec[0] = Index2DTo1D(x, cellCount[1] - 1);
						vec[1] = 1;
						vec[2] = 1.0;
						vec[3] = frontLowerLeftCorner[1] + calcSimulationAreaExtent()[1];
						reflectiveBoundaryCells.push_back(vec);
					}
				}
		}
	else if (dim == 3) {
			if (leftReflectiveBoundary || rightReflectiveBoundary)
				for (int y = 0; y < cellCount[1]; y++)
					for (int z = 0; z < cellCount[2]; z++) {
						if (leftReflectiveBoundary) {
							// temp variable
							utils::Vector<double, 4> vec;
							vec[0] = Index3DTo1D(0, y, z);
							vec[1] = 0;
							vec[2] = -1.0;
							vec[3] = frontLowerLeftCorner[0];
							reflectiveBoundaryCells.push_back(vec);
						}
						if (rightReflectiveBoundary) {
							// temp variable
							utils::Vector<double, 4> vec;
							vec[0] = Index3DTo1D(cellCount[0] - 1, y, z);
							vec[1] = 0;
							vec[2] = 1.0;
							vec[3] = frontLowerLeftCorner[0] + calcSimulationAreaExtent()[0];
							reflectiveBoundaryCells.push_back(vec);
						}
					}
		if (frontReflectiveBoundary || backReflectiveBoundary)
				for (int x = 0; x < cellCount[0]; x++)
					for (int z = 0; z < cellCount[2]; z++) {
						if (frontReflectiveBoundary) {
							// temp variable
							utils::Vector<double, 4> vec;
							vec[0] = Index3DTo1D(x, 0, z);
							vec[1] = 1;
							vec[2] = -1.0;
							vec[3] = frontLowerLeftCorner[1];
							reflectiveBoundaryCells.push_back(vec);
						}
						if (backReflectiveBoundary) {
							// temp variable
							utils::Vector<double, 4> vec;
							vec[0] = Index3DTo1D(x, cellCount[1] - 1, z);
							vec[1] = 1;
							vec[2] = 1.0;
							vec[3] = frontLowerLeftCorner[1] + calcSimulationAreaExtent()[1];
							reflectiveBoundaryCells.push_back(vec);
						}
					}
			if (bottomReflectiveBoundary || topReflectiveBoundary)
				for (int x = 0; x < cellCount[0]; x++)
					for (int y = 0; y < cellCount[1]; y++) {
						if (bottomReflectiveBoundary) {
							// temp variable
							utils::Vector<double, 4> vec;
							vec[0] = Index3DTo1D(x, y, 0);
							vec[1] = 2;
							vec[2] = -1.0;
							vec[3] = frontLowerLeftCorner[2];
							reflectiveBoundaryCells.push_back(vec);
						}
						if (topReflectiveBoundary) {
							// temp variable
							utils::Vector<double, 4> vec;
							vec[0] = Index3DTo1D(x, y, cellCount[2] - 1);
							vec[1] = 2;
							vec[2] = 1.0;
							vec[3] = frontLowerLeftCorner[2] + calcSimulationAreaExtent()[2];
							reflectiveBoundaryCells.push_back(vec);
						}
					}
		}
	}


	/// function to assign a vector of particles to formerly properly initialized cells
	/// note, that this function clears all arrays!
	/// @param particles a vector particles
	void	AssignParticles(const std::vector<Particle>& particles)
	{
		// assert
		assert(Cells);

		// clear if necessary, cells
		for(int i = 0; i < getCellCount(); i++)
		{
			if(!Cells[i].empty())Cells[i].clear();
		}

		// go through all particles in the vector...
		for (std::vector<Particle>::const_iterator it = particles.begin() ; it < particles.end(); it++)
		{
			Particle p = *it;
			AddParticle(p);		
		}

	}

	/// a method to check if all particles are properly assigned...
	void CheckAssignment()
	{
		//only 2D!
		if(dim != 2)return;

		// go through all Cells
		for(int i = 0; i < getCellCount(); i++)
		{
			//go through every individual cell
			if(Cells[i].empty())continue;

			bool falseAssignment = false;

			int xIndex = i % this->cellCount[0];
			int yIndex = i / cellCount[0];
			double xmin = this->frontLowerLeftCorner[0] + xIndex * cellSize[0];
			double xmax = xmin + cellSize[1];
			double ymin = this->frontLowerLeftCorner[1] + yIndex * cellSize[1];
			double ymax = ymin + cellSize[1];

			// go through particles and check if they are in boundaries
			for(std::vector<Particle>::iterator it = Cells[i].begin(); it != Cells[i].end(); )
			{
				Particle p = *it;

				if(p.x[0] < xmin)falseAssignment = true;
				if(p.x[0] > xmax)falseAssignment = true;
				if(p.x[1] < ymin)falseAssignment = true;
				if(p.x[1] > ymax)falseAssignment = true;
			}

			if(falseAssignment)std::cout<<"Error:!!!! False Assignment"<<std::endl;
		}
	}
	/// a method that reassigns the particles to the cells they belong to
	/// note that this method should only be called every 'iterationsPerParticleToCellReassignment'th call
	/// a good number for this must be determined experimentally for increased efficiency
	/// the default value of this variable for the LinkedCellAlgorithm is one
	void ReassignParticles()
	{
		int xIndex = 0, yIndex = 0, zIndex = 0;

		// go through all Cells...
		for(int i = 0; i < getCellCount(); i++)
		{
			// does cell contain any elements? - if not continue
			if(Cells[i].empty())continue;
			
			// go through every cell's particles...
			for(std::vector<Particle>::iterator it = Cells[i].begin(); it != Cells[i].end(); )
			{
				Particle p = *it;
			
				// 2D
				if (dim == 2) {
					// calc Index
					xIndex = (int)((p.x[0] - frontLowerLeftCorner[0]) / cellSize[0]);
					yIndex = (int)((p.x[1] - frontLowerLeftCorner[1]) / cellSize[1]);

					// is particle outside of its father cell?
					if(Index2DTo1D(xIndex, yIndex) != i)
					{
						// is particle contained in grid or halo?
						if (xIndex < 0 || yIndex < 0 ||
							xIndex >= cellCount[0] || yIndex >= cellCount[1])
							halo.push_back(p);
						else {
							int index = Index2DTo1D(xIndex, yIndex);
							Cells[index].push_back(p);
						}

						// remove particle from current cell (i-th cell)
						it = Cells[i].erase(it);

					}
				}
				// 3D
				else if (dim == 3) {
					// calc Index
					xIndex = (int)((p.x[0] - frontLowerLeftCorner[0]) / cellSize[0]);
					yIndex = (int)((p.x[1] - frontLowerLeftCorner[1]) / cellSize[1]);
					zIndex = (int)((p.x[2] - frontLowerLeftCorner[2]) / cellSize[2]);

					// is particle outside of its father cell?
					if(Index2DTo1D(xIndex, yIndex) != i)
					{
						// particle in halo?
						if (xIndex < 0 || yIndex < 0 || zIndex < 0 ||
							xIndex >= cellCount[0] || yIndex >= cellCount[1] || zIndex >= cellCount[2])
							halo.push_back(p);
						else {
							int index = Index3DTo1D(xIndex, yIndex, zIndex);
							Cells[index].push_back(p);
						}
					
						// remove particle from current cell (i-th cell)
						it = Cells[i].erase(it);

					}
				}

				// manually increment iterator, as erase may delete last element, and therefore vector::iterator inc will cause an error
				if(it != Cells[i].end())it++;
			}
		}

		//test
		this->CheckAssignment();
	}

	// Deprecated...
	
	//int iterationsToCellReassignmentLeft;


	

public:
	/// default constructor, set everthing to good values
	LinkedCellParticleContainer()
	{
		Cells = NULL;
		iterationCount = 0;
		cutoffDistance = 0.0;
		reflectiveBoundaryDistance = 0;
		iterationsPerParticleToCellReassignment = 2;
	}

	/// destructor
	~LinkedCellParticleContainer()
	{
		SAFE_DELETE_A(Cells);
	}

	/// a constructor that takes quite a lot of arguments
	/// @param particles
	/// TODO: comment params...
	LinkedCellParticleContainer(const std::vector<Particle>& particles,
								const double cutoffDistance,
								utils::Vector<double, dim> frontLowerLeftCorner,
								utils::Vector<double, dim> simulationAreaExtent,
								int iterationsPerParticleToCellReassignment,
								bool leftReflectiveBoundary, bool rightReflectiveBoundary,
								bool frontReflectiveBoundary, bool backReflectiveBoundary,
								// these two will be ignored in the two-dimensional case
								bool bottomReflectiveBoundary, bool topReflectiveBoundary,
								double sigma)
	{
		// set to zero, so Init doesn't crash...
		Cells = NULL;

		Init(particles, cutoffDistance, frontLowerLeftCorner, simulationAreaExtent, iterationsPerParticleToCellReassignment,
			leftReflectiveBoundary, rightReflectiveBoundary,
								frontReflectiveBoundary, backReflectiveBoundary,
								bottomReflectiveBoundary, topReflectiveBoundary, sigma);		
	}

	/// init function taking a lot of arguments
	void								Init(const std::vector<Particle>& particles,
											 const double cutoffDistance,
											 utils::Vector<double, dim> frontLowerLeftCorner,
											 utils::Vector<double, dim> simulationAreaExtent,
											 int iterationsPerParticleToCellReassignmentbool,
											 bool leftReflectiveBoundary, bool rightReflectiveBoundary,
											 bool frontReflectiveBoundary, bool backReflectiveBoundary,
											 // these two will be ignored in the two-dimensional case
											 bool bottomReflectiveBoundary, bool topReflectiveBoundary,
											 double sigma)
	{
		// member initialization
		this->iterationCount = 0;
		this->cutoffDistance = cutoffDistance;
		this->frontLowerLeftCorner = frontLowerLeftCorner;
		this->iterationsPerParticleToCellReassignment = iterationsPerParticleToCellReassignment;
		this->reflectiveBoundaryDistance = sigma * 1.1225;


		// calc number of cells in each dimension, note that we are rounding down
		// this is done because the number of cells that fit in the designated area may not be a natural number
		// therefore the number of cells is rounded down to give an exact fit
		for(int i = 0; i < dim; i++) {
			cellCount[i] = (int)(simulationAreaExtent[i] / cutoffDistance);
			// special case: not even one cell fits into the area completely
			if (cellCount[i] == 0)
				cellCount[i] = 1;
		}

		// set cell size: adapt the cellSize in each dimension
		// this can be slightly larger than the cutoffDistance if the cells don't fit
		// into the the area precisely and have to be enlarged, as just mentioned above
		for(int i = 0; i < dim; i++)cellSize[i] = (double)(simulationAreaExtent[i] / cellCount[i]);
		
		//delete if necessary
		SAFE_DELETE_A(Cells);

		// alloc mem
		Cells = new std::vector<Particle>[getCellCount()];

		// generate pairs
		generatePairs();
		
		// assign the particles to their initial cells
		AssignParticles(particles);

		// set the reflective boundary conditions
		SetReflectiveBoundaries(leftReflectiveBoundary, rightReflectiveBoundary,
								frontReflectiveBoundary, backReflectiveBoundary,
								bottomReflectiveBoundary, topReflectiveBoundary);
	}
	
	/// applies the reflective boundary condition to all cells that apply
	void ApplyReflectiveBoundaryConditions(void(*func)(void*, Particle&, Particle&), void *data) {
		for (std::vector<utils::Vector<double, 4>>::iterator it = reflectiveBoundaryCells.begin(); it != reflectiveBoundaryCells.end(); it++)
		{
			// TODO: this should obviously work by reference to the array, not copy it
			// does it do that right now?
			std::vector<Particle>& cell = cells[(int)(elem[0])];
			int axis = (int)(elem[1]);
			double direction = elem[2];
			double border = elem[3];
			for (std::vector<Particle>::iterator it = cell.begin() ; it < cell.end(); it++) {
				Particle& p = *it;
				double dist = direction * (border - p.x[axis]);
				// skip the particle if it is too far away from the border
				if (dist > reflectiveBoundaryDistance)
					continue;
				// create a temporary, virtual Particle
				Particle vp = new Particle(p);
				vp.x[axis] = border;
				(*func)(data, p, vp);
				// TODO: delete vp manually?
			}
		}
	}

	/// a method to add a Particle to the LinkedCellParticleContainer
	void AddParticle(const Particle& particle)
	{
		int xIndex = 0,yIndex = 0, zIndex = 0;

		assert(Cells);

		// 2D
		if(dim == 2) {
				// calc Index
				xIndex = (int)((particle.x[0] - frontLowerLeftCorner[0]) / cellSize[0]);
				yIndex = (int)((particle.x[1] - frontLowerLeftCorner[1]) / cellSize[1]);

				// is particle contained in grid or halo?
				if (xIndex < 0 || yIndex < 0 ||
					xIndex >= cellCount[0] || yIndex >= cellCount[1])
					halo.push_back(particle);
				else {
					int index = Index2DTo1D(xIndex, yIndex);
					Cells[index].push_back(particle);
				}
			}
		// 3D
		else if(dim == 3) {
			// calc Index
			xIndex = (int)((particle.x[0] - frontLowerLeftCorner[0]) / cellSize[0]);
			yIndex = (int)((particle.x[1] - frontLowerLeftCorner[1]) / cellSize[1]);
			zIndex = (int)((particle.x[2] - frontLowerLeftCorner[2]) / cellSize[2]);

			// particle in halo?
			if (xIndex < 0 || yIndex < 0 || zIndex < 0 ||
				xIndex >= cellCount[0] || yIndex >= cellCount[1] || zIndex >= cellCount[2])
				halo.push_back(particle);
			else {
				int index = Index3DTo1D(xIndex, yIndex, zIndex);
				Cells[index].push_back(particle);
			}
		}
	}

	/// a method that takes a void(*func)(void*, Particle) and uses it to iterate over all Particles
	/// @param data additional data given to func
	void Iterate(void(*func)(void*, Particle&), void *data) {
		
		assert(Cells);

		// go through all Cells...
		for(int i = 0; i < getCellCount(); i++)
		{
			for(std::vector<Particle>::iterator it = Cells[i].begin(); it != Cells[i].end(); it++)
			{
				Particle& p = *it;
				(*func)(data, p);				
			}
		}

		// inc counter for iterations, if needed reassign particles in cells...
		iterationCount++;
		if(iterationCount > this->iterationsPerParticleToCellReassignment)
		{
			ReassignParticles();
			iterationCount = 0;
		}
	}
	/// a method that takes a void(*func)(void*, Particle, Particle) and
	/// uses it to iterate over all pairs of Particles (each symmetrical pair is
	/// only taken once to reduce redundancy)
	/// @param data additional data given to func
	void IteratePairwise(void(*func)(void*, Particle&, Particle&), void *data) {
		
		// iterate over all elements of cellPairs. 
		for (std::vector<utils::Vector<unsigned int, 2>>::iterator it = cellPairs.begin(); it != cellPairs.end(); it++)
		{
			utils::Vector<unsigned int, 2> pair = *it;

			// are cells valid? - if not next iteration
			// TODO: doesn't this problem resolve itself if it's ignored? This is an unnecessary check that just costs time
			// Answer: no it doesn't! if we try to iterate over an empty vector, this will cause an error!
			if(Cells[pair[0]].empty() || Cells[pair[1]].empty())continue;

			// calc data for a pair (a, b) where a != b
			if(pair[0] != pair[1])
				for (std::vector<Particle>::iterator it1 = Cells[pair[0]].begin() ; it1 < Cells[pair[0]].end(); it1++)
					for (std::vector<Particle>::iterator it2 = Cells[pair[1]].begin() ; it2 < Cells[pair[1]].end(); it2++)
						{
							// call the function on the pair of Particles
							Particle& p1 = *it1;
							Particle& p2 = *it2;

							//is distance squared less than cutoff radius squared?
							if(p1.x.distanceSq(p2.x) < cutoffDistance * cutoffDistance)
								(*func)(data, p1, p2);
						}
			// calc data for a pair (a, a)
			else
				for (std::vector<Particle>::iterator it1 = Cells[pair[0]].begin() ; it1 < Cells[pair[0]].end(); it1++)
					for (std::vector<Particle>::iterator it2 = it1 + 1 ; it2 < Cells[pair[1]].end(); it2++)
					// make sure that a Particle is not paired with itself
					if (it1 != it2)			
					{
						// call the function on the pair of Particles
						Particle& p1 = *it1;
						Particle& p2 = *it2;
						//is distance squared less than cutoff radius squared?
						if(p1.x.distanceSq(p2.x) < cutoffDistance * cutoffDistance)
							(*func)(data, p1, p2);
					}
		}

		// inc counter for iterations, if needed reassign particles in cells...
		iterationCount++;
		if(iterationCount > this->iterationsPerParticleToCellReassignment)
		{
			ReassignParticles();
			iterationCount = 0;
		}
	}

	/// add particles from *.txt file
	void							AddParticlesFromFile(const char *filename)
	{
		//...
	}

	/// our new fileformat, replace later AddParticlesFromFile
	/// @return return true if file could be read
	bool							AddParticlesFromFileNew(const char *filename)
	{
		//	... 
		// nothing
		return true;
	}

	/// removes all particles
	void							Clear()
	{
		// delete halo's particles
		if(!halo.empty())halo.clear();

		assert(Cells);

		// go through all Cells
		for(int i = 0; i < getCellCount(); i++)
		{
			if(!Cells[i].empty())Cells[i].clear();
		}
	}

	/// are any particles contained?
	bool							IsEmpty()
	{
		// any particles in halo contained?
		if(!halo.empty())return false;

		// any cells contained?
		if(!Cells)return true;

		// go through all Cells
		for(int i = 0; i < getCellCount(); i++)
		{
			if(!Cells[i].empty())return false;
		}

		return true;
	}

	/// returns how many particles are contained in this container
	/// @return returns count of particles in this container
	unsigned int					getParticleCount()
	{
		unsigned int sum = 0;

		assert(Cells);

		//go through all cells...
		for(int i = 0; i < getCellCount(); i++)
		{
			if(!Cells[i].empty())sum += Cells[i].size();
		}

		return sum;
	}
	
	
	// Quick 'n' dirty hack
	std::vector<Particle> p;
	
	
	/// this method shall be later removed...
	/// returns ListParticleContainer's internal container
	const std::vector<Particle>&	getParticles()
	{
		// this method is only a quick hack, to test the algorithm...
		// for this reason we need a static variable to secure, the reference is valid

		

		// clear if necessary
		if(!p.empty())p.clear();

		
		/*
		// add halo...
		//p.insert(p.end(), halo.begin(), halo.end());

		assert(Cells);

		//go through cells
		for(int i = 0; i < getCellCount(); i++)
		{
			p.insert(p.end(), Cells[i].begin(), Cells[i].end());
		}

		*/


		// test wise
		int filledCells = 0;

		for(int i = 0; i < getCellCount(); i++)
		{
			if(!Cells[i].empty())
			{
					for(std::vector<Particle>::iterator it = Cells[i].begin(); it != Cells[i].end(); it++)
					{
						//p.insert(p.end(), Cells[i].begin(), Cells[i].end());
						Particle pt = *it;
						pt.type = i % 16; //4 colors!
						p.push_back(pt);
					}
					filledCells++;
			}
		}

		std::cout<<"filled cells: "<<filledCells<<" / "<<getCellCount()<<std::endl;
		return p;
	}
};

#endif 

