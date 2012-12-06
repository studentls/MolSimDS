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
#include <list>
#include "Logging.h"
#include "Particle.h"
#include "ParticleContainer.h"
#include "utils/utils.h"

// use flag System for Boundaries
// to combine flags use e.g. BC_RIGHT | BC_LEFT
#define BC_NONE			0
#define BC_OUTFLOW		0x1000
#define BC_ALL			0x3F		// right | left |...|back
#define BC_RIGHT		0x1
#define BC_LEFT			0x2
#define BC_TOP			0x4
#define BC_BOTTOM		0x8
#define BC_FRONT		0x10
#define BC_BACK			0x20
/// a class that is used to store Particles and iterate over them
class LinkedCellParticleContainer : public ParticleContainer {

private:
	
	/// dimension of the grid
	unsigned int						dim;

	/// Cutoff distance
	double								cutoffDistance;

	/// vector containing X, Y, (Z) size of each cell
	utils::Vector<double, 3>			cellSize;

	/// vector containing count of cells in X, Y, (Z) direction
	utils::Vector<unsigned int, 3>		cellCount;

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
	utils::Vector<double, 3>			frontLowerLeftCorner;

	/// a dynamic array, containing all Cells encoded in a 1D array
	std::vector<Particle>				*Cells;
	
	/// halo particle
	/// Note that the halo acts only as a storage for the sake of completeness.
	/// The field could be dropped as well, as no interaction with cells occurs
	std::vector<Particle>				halo;
	
	/// Array of Pairs of Cell Indices, e.g. (1, 2) is the pair adressing Cell 1 and Cell 2
	/// where 1, 2 are the index of the Cells array
	std::vector<utils::Vector<unsigned int, 2>>	cellPairs;
		
	/// the distance below which reflective boundaries reflect
	double reflectiveBoundaryDistance;

	/// boundary conditions encoded as flags
	unsigned int boundaryConditions;

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

public:
	/// default constructor, set everthing to good values
	LinkedCellParticleContainer()
	{
		Cells = NULL;
		cutoffDistance = 0.0;
		reflectiveBoundaryDistance = 0;
		dim = 0;
	}

	/// destructor
	~LinkedCellParticleContainer()
	{
		SAFE_DELETE_A(Cells);
	}

	/// a constructor that takes quite a lot of arguments
	/// @param dim dimension of the grid, shall be 2 or 3
	/// @param particles from which the grid will be constructed
	/// @param cutoffDistance in this simple implementation, essentially the grid meshwidth
	/// @param frontLowerLeftCorner offset of the grid
	/// @param simulationAreaExtent extent of the grid
	/// @param boundaryconditions use to specify flags like BC_LEFT, ...
	/// @param sigma used to calculate reflective Boundary distance
	/// TODO: comment params...
	LinkedCellParticleContainer(const unsigned int dim,
								const std::vector<Particle>& particles,
								const double cutoffDistance,
								utils::Vector<double, 3> frontLowerLeftCorner,
								utils::Vector<double, 3> simulationAreaExtent,
								const unsigned int boundaryConditions,
								double sigma)
	{
		// set to zero, so Init doesn't crash...
		Cells = NULL;
		this->dim = dim;
		this->cutoffDistance = 0.0;
		reflectiveBoundaryDistance = 0;


		// call init
		Init(particles, cutoffDistance, frontLowerLeftCorner, simulationAreaExtent, 
			boundaryConditions, sigma);		
	}

	/// init function taking a lot of arguments
	void								Init(const std::vector<Particle>& particles,
											 const double cutoffDistance,
											 utils::Vector<double, 3> frontLowerLeftCorner,
											 utils::Vector<double, 3> simulationAreaExtent,
											 const unsigned int boundaryConditions,
											 double sigma)
	{
		this->cutoffDistance = cutoffDistance;
		this->frontLowerLeftCorner = frontLowerLeftCorner;
		this->reflectiveBoundaryDistance = sigma * 1.1225;
		this->boundaryConditions = boundaryConditions;

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

		if(!Cells)LOG4CXX_ERROR(generalOutputLogger, "memory allocation for cells failed!");
		
		// assign the particles to their initial cells
		AssignParticles(particles);

		// set the reflective boundary conditions
		SetReflectiveBoundaries(boundaryConditions & BC_LEFT,
								boundaryConditions & BC_RIGHT,
								boundaryConditions & BC_FRONT,
								boundaryConditions & BC_BACK,
								boundaryConditions & BC_BOTTOM,
								boundaryConditions & BC_TOP);
	}
	
	/// applies the reflective boundary condition to all cells that apply
	void ApplyReflectiveBoundaryConditions(void(*func)(void*, Particle&, Particle&), void *data) {
		for (std::vector<utils::Vector<double, 4>>::iterator it = reflectiveBoundaryCells.begin(); it != reflectiveBoundaryCells.end(); it++)
		{
			// TODO: this should obviously work by reference to the array, not copy it
			// does it do that right now?
			utils::Vector<double, 4>& elem = *it;
			std::vector<Particle>& cell = Cells[(int)(elem[0])];
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
				// a copy of the current particle, but located at the border
				Particle vp(p);
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
	void Iterate(void(*func)(void*, Particle&), void *data);

	/// helper function, to get neighbours
	/// @return vector of neighbour cell indices...
	inline std::vector<unsigned int> getNeighbours(unsigned int index)
	{
		std::vector<unsigned int> neighbours;

		// only dim = 2 supported
		if(dim != 2)return neighbours;

		unsigned int x = index % cellCount[0];
		unsigned int y = index / cellCount[0];

		//beware of special cases!
		for(int xp = -1; xp <= 1; xp++)
			for(int yp = -1; yp <= 1; yp++)
			{
				int index = x +xp + (y + yp) * cellCount[0];

				//valid?
				if(index < 0 || index >= getCellCount())continue;

				neighbours.push_back(index);
			}

		return neighbours;
	}

	/// a method that takes a void(*func)(void*, Particle, Particle) and
	/// uses it to iterate over all pairs of Particles (each symmetrical pair is
	/// only taken once to reduce redundancy)
	/// @param data additional data given to func
	void IteratePairwise(void(*func)(void*, Particle&, Particle&), void *data);

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

		// add halo...
		//p.insert(p.end(), halo.begin(), halo.end());

		assert(Cells);

		//go through cells
		for(int i = 0; i < getCellCount(); i++)
		{
			//empty cell?
			if(Cells[i].empty())continue;
			
			for(std::vector<Particle>::const_iterator it = Cells[i].begin(); it != Cells[i].end(); it++)
			{
				p.push_back(*it);
			}
		}
		return p;
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
					if(Index3DTo1D(xIndex, yIndex, zIndex) != i)
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
	}

	/// method to identify container
	/// @return returns PCT_LIST
	ParticleContainerType			getType() {return PCT_LINKEDCELL;}

	/// method to return halo particles
	/// @return returns reference to halo particle vector
	std::vector<Particle>&			getHaloParticles()	{return this->halo;}

	/// method to clear halo
	void							clearHaloParticles()		{if(!halo.empty())halo.clear();}

	/// get all boundary particles
	std::vector<Particle>			getBoundaryParticles();
};

#endif 

