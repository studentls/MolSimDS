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

	/// vector containg count of cells in X, Y, (Z) direction
	utils::Vector<unsigned int, dim>	cellCount;

	/// get total cell count
	inline unsigned int					getCellCount()
	{
		unsigned int sum = 0;
		for(int i = 0; i < dim; i++)
			sum += CellCount[i];
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

	/// a dynamic array, containg all Cells encoded in a 1D array
	std::vector<Particle>				*Cells;
	
	/// ???
	// Note that the halo acts only as a storage for the sake of completeness. The field could be dropped as well
	std::vector<Particle> halo;
	
	/// Array of Pairs of Cell Indices, e.g. (1, 2) is the pair adressing Cell 1 and Cell 2
	/// where 1, 2 are the index of the Cells array
	std::vector<utils::Vector<int, 2>>	cellPairs;

	/// updating cells every iteration is not very appropriate, better wait some
	/// iterations and perform then update
	/// this variable stores how many iterations the container should wait
	int iterationsPerParticleToCellReassignment;
	
	/// count of iterations so far, note this variable is incremented every time Iterate or (!) IteratePairwise is called!
	int	iterationCount;

	/// helper function to convert fast 2D indices to 1D based on cellCount
	inline unsigned int Index2DTo1D(unsigned int x, unsigned int y)
	{
		// correct index?
		assert(x + cellCount[0] < getCellCount());

		return x + cellCount[0] * y;
	}

	/// helper function to convert fast 3D indices to 1D based on cellCount
	inline unsigned int Index3DTo1D(unsigned int x, unsigned int y, unsigned int z)
	{
		// correct index?
		assert(x + y * (cellCount[0] + z * cellCount[1]) < getCellCount());

		return x + y * (cellCount[0] + z * cellCount[1]);
	}


	/// generate Index Pairs according to cellCount
	// TODO: write a testcase about this function
	void		generatePairs()	
	{
		// temp pair variable
		utils::Vector<int, 2> pair;
		
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
									if (x + xp < cellNumbers[0] && y + yp < cellNumbers[1] && z + zp < cellNumbers[2])
									{										
										// make pairs ((x,y,z), (x + xp, y + yp, z + zp))
										pair[0] = Index3DTo1D(x, y, z);
										pair[1] = Index3DTo1D(x + xp, y + yp, z + zp);
										cellPairs.push_back(pair);
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

	/// a method that reassigns the particles to the cells they belong to
	/// note that this method should only be called every 'iterationsPerParticleToCellReassignment'th call
	/// a good number for this must be determined experimentally for increased efficiency
	/// the default value of this variable for the LinkedCellAlgorithm is one
	void ReassignParticles()
	{
		int xIndex = 0,yIndex = 0, zIndex = 0;

		// go through all Cells...
		for(int i = 0; i < getCellCount(); i++)
		{
			// go through every cell's particles...
			for(std::vector<Particle>::iterator it = Cells[i].begin(); it != Cells[i].end(); it++)
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
						Cells[i].erase(it);

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
						Cells[i].erase(it);

					}
				}
			}
		}
	}

	// Deprecated...
	
	//int iterationsToCellReassignmentLeft;


	// Deprecated ...
	// each double in the Vector means something different
	// 0, int: the cell
	// 1, int: the axis
	// 2, double: the direction. -1.0 if the border is on the lower side, 1.0 on the positive side
	// 3, double: the border's coordinate on the given axis, in the given direction
	// TODO: question: is it possible to save an int in a slot that takes a double without manual conversion? Works in c#, unsure about c++
	///List<utils::Vector<double, 4>> reflectiveBoundaryCells;
	
	// TODO: initialize this
	//double reflectiveBoundaryDistance;
	// TODO: set reflectiveBoundaryCells; depends on how the value is taken as a constructor parameter

	

public:

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
								int iterationsPerParticleToCellReassignment
								)
	{
		// member initialization
		this->iterationCount = 0;
		this->cutoffDistance = cutoffDistance;
		this->frontLowerLeftCorner = frontLowerLeftCorner;
		this->iterationsPerParticleToCellReassignment = iterationsPerParticleToCellReassignment;

		// set cell size( is currently simply cutoffDimension)
		for(int i = 0; i < dim; i++)cellSize[i] = cutoffDistance;

		// calc number of cells in each dimension, note that we are rounding down
		for(int i = 0; i < dim; i++)cellCount[i] = simulationAreaExtent[i] / cellSize[i];
		
		// generate pairs
		generatePairs();
		
		// assign the particles to their initial cells
		AssignParticles(particles);
		
	}
	/*// applies the reflectove boundary condition to all cells that apply
	void ApplyReflectiveBoundaryConditions(void(*func)(void*, Particle&, Particle&), void *data) {
		// TODO: iterate over all elements of reflectiveBoundaryCells. This depends on the data type of reflectiveBoundaryCells
		for () {
			std::vector<Particle> cell = cells[(int)(elem[0])];
			int axis = (int)(elem[1]);
			double direction = elem[2];
			double border = elem[3];
			for (std::vector<Particle>::iterator it = cell.begin() ; it < cell.end(); it++) {
				Particle& p = *it;
				double dist = direction * (border - p.x[axis]);
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
	*/
	/// a method to add a Particle to the ListParticleContainer
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

			// calc data for every pair
			for (std::vector<Particle>::iterator it1 = Cells[pair[0]].begin() ; it1 < Cells[pair[0]].end(); it1++)
				for (std::vector<Particle>::iterator it2 = Cells[pair[1]].begin() ; it2 < Cells[pair[1]].end(); it2++)
					// make sure that a Particle is not paired with itself
					if (it1 != it2)
					{
						// call the function on the pair of Particles
						Particle& p1 = *it1;
						Particle& p2 = *it2;
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
	void							AddParticlesFromFile(const char *filename);

	/// our new fileformat, replace later AddParticlesFromFile
	/// @return return true if file could be read
	bool							AddParticlesFromFileNew(const char *filename);

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

		assert(Cells);

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

	/// this method shall be later removed...
	/// returns ListParticleContainer's internal container
	const std::vector<Particle>&	getParticles();
};

#endif 

