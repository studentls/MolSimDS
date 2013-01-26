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
// reflective boundaries
#define BC_NONE			0
#define BC_OUTFLOW		0x1000
#define BC_ALL			0x3F		// right | left |...|back
#define BC_RIGHT		0x1
#define BC_LEFT			0x2
#define BC_TOP			0x4
#define BC_BOTTOM		0x8
#define BC_FRONT		0x10
#define BC_BACK			0x20

// special periodic boundaries
#define BC_PERIODIC_XAXIS	0x100
#define BC_PERIODIC_YAXIS	0x200
#define BC_PERIODIC_ZAXIS	0x400

struct PeriodicBoundary
{
	unsigned int cell1;	/// index for first cell
	unsigned int cell2;	/// index for second cell
	double	xAxis;		/// distance on x - axis (0 or length of area on axis)
	double	yAxis;		/// distance on y - axis (0 or length of area on axis)
	double	zAxis;		/// distance on y - axis (0 or length of area on axis)
};


/// indices for Cells used for OpenMP parallelization
/// structure to hold a strip of Indices
class  IndexStrip : public utils::TFastArray<utils::Vector<unsigned int, 2> >
{
private:
	/// little helper function, to create fast a pair
	/// @param i1 first index
	/// @param i2 second index
	utils::Vector<unsigned int, 2> makePair(const unsigned int i1, const unsigned int i2)
	{
		utils::Vector<unsigned int, 2> res;
		res[0] = i1;
		res[1] = i2;

		return res;
	}

	/// little helper function, to create fast a 1D index
	/// @param x first index
	/// @param y second index
	/// @param cellCount 3D-vector with domain dimensions
	unsigned int Index2DTo1D(const unsigned int x, const unsigned int y, const utils::Vector<unsigned int, 3>		cellCount)
	{
		assert(x + cellCount[0] * y < cellCount[0] * cellCount[1]);

		return x + cellCount[0] * y;
	}
public:
	IndexStrip()
	{
		TFastArray<utils::Vector<unsigned int, 2> >();
	}

	/// method to construct indices for a vertical strip(for 2D)
	/// @param verticalCellCount amount for how many cells indices shall be constructed
	/// @param leftStart index of the left Cell. Indices will be constructed for verticalCellCount cells between leftStart and leftStart+1 Cells
	/// @param cellCount 3D Vector containing domain cell count in each dimension
	void	constructVerticalStripIndices(const int leftStart, const int verticalCellCount, const utils::Vector<unsigned int, 3>		cellCount)
	{
		for(int i = 0; i < verticalCellCount - 1; i++)
		{
			// - pair
			this->push_back(makePair(Index2DTo1D(leftStart, i, cellCount),
									 Index2DTo1D(leftStart + 1, i, cellCount)));
			// | left pair
			this->push_back(makePair(Index2DTo1D(leftStart, i, cellCount),
									 Index2DTo1D(leftStart, i + 1, cellCount)));
			//  right | pair
			this->push_back(makePair(Index2DTo1D(leftStart + 1, i, cellCount),
									 Index2DTo1D(leftStart + 1, i + 1, cellCount)));
			// \ pair
			this->push_back(makePair(Index2DTo1D(leftStart, i, cellCount),
									 Index2DTo1D(leftStart + 1, i + 1, cellCount)));
			// / pair
			this->push_back(makePair(Index2DTo1D(leftStart + 1, i, cellCount),
									 Index2DTo1D(leftStart, i + 1, cellCount)));
		}

		// final pair
		
		this->push_back(makePair(Index2DTo1D(leftStart, verticalCellCount - 1, cellCount),
									Index2DTo1D(leftStart + 1, verticalCellCount - 1, cellCount)));
	}
};

/// a class that is used to store Particles and iterate over them
/// utilizes Linked-Cell algorithm for improved performance (O(n) instead of O(n^2))
/// see this graph which compares LinkedCell with a brute-force algorithm:
/// \image html LinkedCellPerformance.png
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
	utils::Vector<double, 3>						frontLowerLeftCorner;

	/// a dynamic array, containing all Cells encoded in a 1D array
	/// note that we will encode in this array also a halo cell frame
	/// so the grid has size (x+2)(y+2)(z+2) after construction
	std::vector<Particle>							*Cells;
	
	
	/// Array of Pairs of Cell Indices, e.g. (1, 2) is the pair adressing Cell 1 and Cell 2
	/// where 1, 2 are the index of the Cells array
	/// note that here the space is very important, as g++ has problems otherwise parsing it
	std::vector<utils::Vector<unsigned int, 2> >	cellPairs;

	/// used for OpenMP
	utils::TFastArray<IndexStrip>							oddStrips;	/// contains index data for odd strips
	utils::TFastArray<IndexStrip>							evenStrips; /// contains index pairs for even strips

	/// array of indices of halo cells
	std::vector<unsigned int>						haloIndices;

	/// array of indices of boundary cells
	std::vector<unsigned int>						boundaryIndices;

	/// little helper function, to create fast a pair
	/// @param i1 first index
	/// @param i2 second index
	utils::Vector<unsigned int, 2> makePair(const unsigned int i1, const unsigned int i2)
	{
		utils::Vector<unsigned int, 2> res;
		res[0] = i1;
		res[1] = i2;

		return res;
	}

	/// little helper function, to create fast a triple
	/// @param i1 first index
	/// @param i2 second index
	/// @param i3 third index
	utils::Vector<unsigned int, 3> makeTriple(const unsigned int i1, const unsigned int i2, const unsigned int i3)
	{
		utils::Vector<unsigned int, 3> res;
		res[0] = i1;
		res[1] = i2;
		res[2] = i3;

		return res;
	}

	/// function, which generates pairs for linkedcellalgorithm
	void					generatePairs();

		
	/// the distance below which reflective boundaries reflect
	double reflectiveBoundaryDistance;

	/// boundary conditions encoded as flags
	unsigned int boundaryConditions;

	/// store groups of values that make it easier to apply periodicBoundaryConditions
	std::vector<PeriodicBoundary > periodicBoundaryGroups;

	/// stores at which axis a periodic boundary exists
	utils::Vector<bool, 3> periodicBoundaries;

	// store reflective boundaries as planes
	std::vector<utils::Plane> reflectiveBoundaries;

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

	/// helper function to convert 1D index to 2D indices
	/// @return pair of 2D indices (x, y)
	inline utils::Vector<unsigned int, 2> Index1DTo2D(const unsigned int index)
	{
		utils::Vector<unsigned int, 2> res;
		
		// make assertion
		assert(index < getCellCount());

		res[0] = index % cellCount[0];	// x
		res[1] = index / cellCount[0];  // y

		return res;
	}

	/// helper function to convert 1D index to 3D indices
	/// @return pair of 3D indices (x, y)
	inline utils::Vector<unsigned int, 3> Index1DTo3D(const unsigned int index)
	{
		utils::Vector<unsigned int, 3> res;
		
		// make assertion
		assert(index < getCellCount());

		
		res[0] = index % cellCount[0];	// x
		
		// y, z
		unsigned int temp = index / cellCount[0]; // == y + z * cellCount[1]

		res[1] = temp % cellCount[1]; // y
		res[2] = temp / cellCount[1]; // z
		
		return res;
	}

	/// get opposite HaloCell index
	/// @param index the index for which the opposite halo cell shall be found
	/// @return returns index of the opposite halo cell
	inline unsigned int getOppositeHaloCellIndex(unsigned int index)
	{
		unsigned int res = 0;

		assert(index < getCellCount());

		// extract x, y, z
		unsigned int x = dim == 2 ? Index1DTo2D(index)[0] : Index1DTo3D(index)[0];
		unsigned int y = dim == 2 ? Index1DTo2D(index)[1] : Index1DTo3D(index)[1];
		unsigned int z = dim == 2 ? 0 : Index1DTo3D(index)[2];

		unsigned int newx = x;
		unsigned int newy = y;
		unsigned int newz = z;

		unsigned int xmax = cellCount[0] - 1;
		unsigned int ymax = cellCount[1] - 1;
		unsigned int zmax = cellCount[2] - 1;

		// method is, if x, y, z are halo cells, then mirror indices ( that is for x e.g. newx = xmax - x   )

		if(dim == 2)
		{
			if(x == 0 || x == xmax)newx = xmax - x;
			if(y == 0 || y == ymax)newy = ymax - y;

			res = Index2DTo1D(newx, newy);
		}
		else if(dim == 3)
		{
			if(x == 0 || x == xmax)newx = xmax - x;
			if(y == 0 || y == ymax)newy = ymax - y;
			if(z == 0 || z == zmax)newz = zmax - z;

			res = Index3DTo1D(newx, newy, newz);
		}
		else LOG4CXX_ERROR(generalOutputLogger, "dimension error");

		return res;
	}

	/// helper function to calculate the extent of the simulation area
	/// @return returns the extent of the Simulation area as a 3D vector, note that even if a 2D Grid is used a 3D Vector will be returned
	utils::Vector<double, 3> calcSimulationAreaExtent()
	{
		utils::Vector<double, 3> res;

		// calc based on number of cells in each dimension and their individual extent
		for(int i = 0; i < dim; i++)
		{
			res[i] = this->frontLowerLeftCorner[i] + (cellCount[i] - 2) * cellSize[i];
		}

		return res;
	}
	
	/// set boundaries
	void SetReflectiveBoundaries();

	/// set periodic boundary cells
	/// for easier iteration later
	void SetPeriodicBoundaries(bool xAxis, bool yAxis, bool zAxis);

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

	/// function which will calculate all necessary index array
	void	calcIndices();

	/// function to calcualte Indices for Multithreading
	void	calcMTIndices();

	/// function to calculate indices of the r-th frame from the outside
	/// e.g. r = 0 will return indices of the halo frame
	///		 r = 1 the indices of the boundary cells
	/// @param r the r-th frame from the outside(starting by zero)
	/// @param out index conatiner where indices will be stored
	void	calcFrameIndices(std::vector<unsigned int> &out, const unsigned int r);

	/// applies the reflective boundary condition to all cells that apply
	/// @param func function pointer, to calculate interaction with boundary particle
	/// @data optional data given to func
	void ApplyReflectiveBoundaryConditions(void(*func)(void*, Particle&, Particle&), void *data);
	
	/// applies forces between particles on opposite sides of periodic boundaries
	/// @param func function pointer, to calculate interaction with boundary particle on the other side
	/// @data optional data given to func
	void ApplyPeriodicBoundaryConditions(void(*func)(void*, Particle&, Particle&), void *data);

	/// ensure that particles that leave a periodic boundary are correctly entered on the opposite side
	void ReassignHaloForPeriodicConditions();

	/// inline method to quickly calculate an index for reassignment
	/// @return returns index for grid including halo cells
	inline int PositionToIndex2D(const utils::Vector<double, 3>& pos)
	{
		int xIndex = 0, yIndex = 0;

		// calc Index
		double x = (pos[0] - frontLowerLeftCorner[0]) / cellSize[0];
		double y = (pos[1] - frontLowerLeftCorner[1]) / cellSize[1];
				
		// is particle contained in grid or halo?
		// note that the regular grid starts with (1, 1)
		
		if(x < 0)xIndex = 0;
		else xIndex =(int)x + 1;
		if(y < 0)yIndex = 0;
		else yIndex =(int)y + 1;
		
		if(xIndex > cellCount[0] - 2)xIndex = cellCount[0] - 1;
		if(yIndex > cellCount[1] - 2)yIndex = cellCount[1] - 1;

				
		// make index
		return Index2DTo1D(xIndex, yIndex);
	}

	/// inline method to quickly calculate an index for reassignment
	/// @return returns index for grid including halo cells
	inline int PositionToIndex3D(const utils::Vector<double, 3>& pos)
	{
		int xIndex = 0, yIndex = 0, zIndex = 0;

		// calc Index
		double x = (pos[0] - frontLowerLeftCorner[0]) / cellSize[0];
		double y = (pos[1] - frontLowerLeftCorner[1]) / cellSize[1];
		double z = (pos[2] - frontLowerLeftCorner[2]) / cellSize[2];
				
		// is particle contained in grid or halo?
		// note that the regular grid starts with (1, 1, 1)
		
		if(x < 0)xIndex = 0;
		else xIndex =(int)x + 1;
		if(y < 0)yIndex = 0;
		else yIndex =(int)y + 1;
		if(z < 0)zIndex = 0;
		else zIndex =(int)z + 1;
		
		if(xIndex > cellCount[0] - 2)xIndex = cellCount[0] - 1;
		if(yIndex > cellCount[1] - 2)yIndex = cellCount[1] - 1;
		if(zIndex > cellCount[2] - 2)zIndex = cellCount[2] - 1;

				
		// make index
		return Index3DTo1D(xIndex, yIndex, zIndex);
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
		

		// add +2 for every direction to create halo layer...
		for(int i = 0; i < dim; i++)cellCount[i] += 2;

		//delete if necessary
		SAFE_DELETE_A(Cells);

		// alloc mem
		Cells = new std::vector<Particle>[getCellCount()];

		if(!Cells)LOG4CXX_ERROR(generalOutputLogger, "memory allocation for cells failed!");
		
		// calc Indices (pairs, halo, boundary)
		calcIndices();

		// assign the particles to their initial cells
		AssignParticles(particles);

		// set the reflective boundary conditions
		SetReflectiveBoundaries();
		
		// set the periodic boundary conditions
		SetPeriodicBoundaries(this->boundaryConditions & BC_PERIODIC_XAXIS,
							  this->boundaryConditions & BC_PERIODIC_YAXIS,
							  this->boundaryConditions & BC_PERIODIC_ZAXIS);
	}

	/// a method to add a Particle to the LinkedCellParticleContainer
	/// @param particle the particle to add
	void AddParticle(const Particle& particle)
	{
		int index = 0;

		assert(Cells);

		// 2D
		if(dim == 2) {
				index = PositionToIndex2D(particle.x);
				Cells[index].push_back(particle);
			}
		// 3D
		else if(dim == 3) {
			index = PositionToIndex3D(particle.x);
			Cells[index].push_back(particle);
		}
	}

	/// applies boundary conditions if they are present(both periodic & reflective conditions)
	void	ApplyBoundaryConditions(void(*func)(void*, Particle&, Particle&), void *data)
	{
		if(boundaryConditions != 0)
			ApplyReflectiveBoundaryConditions(func, data);
		if(boundaryConditions & BC_PERIODIC_XAXIS ||
		   boundaryConditions & BC_PERIODIC_XAXIS ||
		   boundaryConditions & BC_PERIODIC_XAXIS)
		   ApplyPeriodicBoundaryConditions(func, data);
		
	}

	/// a method that takes a void(*func)(void*, Particle) and uses it to iterate over all Particles
	/// @param data additional data given to func
	/// @param func function pointer, to calculate interaction with boundary particle
	/// @data optional data given to func
	void Iterate(void(*func)(void*, Particle&), void *data);

	/// helper function, to get neighbours for a 2D grid
	/// @param index 1D index of a grid cell
	/// @return vector of neighbour cell indices...
	inline std::vector<unsigned int> getNeighbours2D(unsigned int index)
	{
		std::vector<unsigned int> neighbours;

		assert(dim == 2);

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

	/// helper function, to get neighbours for a 3D grid
	/// @param index 1D index of a grid cell
	/// @return vector of neighbour cell indices...
	inline std::vector<unsigned int> getNeighbours3D(unsigned int index)
	{
		std::vector<unsigned int> neighbours;

		assert(dim == 3);

		unsigned int x = index % cellCount[0];
		unsigned int z = index / (cellCount[0] * cellCount[1]);
		unsigned int y = (index - z * cellCount[0]*cellCount[1]) / cellCount[0];

		//beware of special cases!
		for(int xp = -1; xp <= 1; xp++)
			for(int yp = -1; yp <= 1; yp++)
				for(int zp = -1; zp <= 1; zp++)
				{
					int index = x +xp + (y + yp) * cellCount[0] + (z + zp) * cellCount[0] * cellCount[1];

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
	/// @param func function pointer, to calculate interaction with boundary particle
	/// @data optional data given to func
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
		assert(Cells);

		// go through all Cells and delete if necessary
		for(int i = 0; i < getCellCount(); i++)
		{
			if(!Cells[i].empty())Cells[i].clear();
		}
	}

	/// are any particles contained?
	bool							IsEmpty()
	{
		// any cells exist?
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
	/// this method also ensures that particles that leave a periodic boundary enter on the other side
	void ReassignParticles()
	{
		// ensure that particles that leave a periodic boundary enter on the other side
		ReassignHaloForPeriodicConditions();

		int index = 0;

		// go through all Cells...
		for(int i = 0; i < getCellCount(); i++)
		{
			// does cell contain any elements? - if not continue
			if(Cells[i].empty())continue;
			
			// go through every cell's particles...
			for(std::vector<Particle>::iterator it = Cells[i].begin(); it != Cells[i].end(); )
			{
				Particle& p = *it;
			
				// 2D
				if (dim == 2) {
					index = PositionToIndex2D(p.x);
					
					// is particle outside of its father cell?
					if(index != i)
					{
						// reassign
						Cells[index].push_back(p);
						
						// remove particle from current cell (i-th cell)
						it = Cells[i].erase(it);

					}
				}
				// 3D
				else if (dim == 3) {
					index = PositionToIndex3D(p.x);
					
					// is particle outside of its father cell?
					if(index != i)
					{
						// reassign			
						Cells[index].push_back(p);
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
	/// @return returns PCT_LINKEDCELL
	ParticleContainerType			getType() {return PCT_LINKEDCELL;}

	/// method to return halo particles
	/// @return returns a vector containing all halo particles
	std::vector<Particle>			getHaloParticles();

	/// method to clear halo
	void							clearHaloParticles();

	/// get all boundary particles
	/// @return new vector of boundary particles
	std::vector<Particle>			getBoundaryParticles();

	/// method to retrieve a Bounding Box, which surrounds all particles
	/// @return returns a BoundingBox, which defines extent and center of all particles in the container(bounding box)
	utils::BoundingBox				getBoundingBox()
	{
		using namespace utils;

		// easy task here, everything is stored!
		BoundingBox bb;

		bb.extent = this->calcSimulationAreaExtent();
		bb.center = bb.extent * 0.5 + frontLowerLeftCorner;
		
		return bb;
	}

	/// @return amount of cells in each direction
	utils::Vector<unsigned int, 3>	getCellExtent()	{return this->cellCount - utils::Vector<unsigned int, 3>(2);}

	// befriend with ParticleContainerTest in order to test private functions
	friend class					ParticleContainerTest;
};

#endif 

