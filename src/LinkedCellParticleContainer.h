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

	// NOTE: my Visual Studio seems to be fucked for some reason
	// I don't get an intellisense list
	// therefore, some functions may be syntactic garbage and not work at all
	// this file just defines the logic, but is not compilable

	// TODO:
	// let's add the methods applyVelocity(...) and applyForce(...) as well (to the general interface ParticleContainer)
	// they do the same thing and take the same parameters as iterate and iteratePairwaise,
	// but in the case of this class, they also call additional functions to deal with the borders
	// and to reassign the particles to cells

	double cutoffDistance;
	double *cellSize;
	utils::Vector<double, dim> frontLowerLeftCorner;
	std::vector<Particle> particles;
	std::vector<Particle> *cells;
	int *cellNumbers;
	int cellCount;
	// Note that the halo acts only as a storage for the sake of completeness. The field could be dropped as well
	std::vector<Particle> halo;
	// TODO: what data type works as a simple list, but is relatively fast? Use it for cellPairs and reflectiveBoundaryCells
	List<utils::Vector<int, 2>> cellPairs;
	int iterationsPerParticleToCellReassignment;
	int iterationsToCellReassignmentLeft;
	// each double in the Vector means something different
	// 0, int: the cell
	// 1, int: the axis
	// 2, double: the direction. -1.0 if the border is on the lower side, 1.0 on the positive side
	// 3, double: the border's coordinate on the given axis, in the given direction
	// TODO: question: is it possible to save an int in a slot that takes a double without manual conversion? Works in c#, unsure about c++
	List<utils::Vector<double, 4>> reflectiveBoundaryCells;
	// TODO: initialize this
	double reflectiveBoundaryDistance;
	// TODO: set reflectiveBoundaryCells; depends on how the value is taken as a constructor parameter

public:

	/// a constructor that takes quite a lot of arguments
	LinkedCellParticleContainer(const std::vector<Particle>& particles,
								double cutoffDistance,
								utils::Vector<double, dim> frontLowerLeftCorner,
								utils::Vector<double, dim> simulationAreaExtent,
								int iterationsPerParticleToCellReassignment,
								bool leftReflectiveBoundary, bool rightReflectiveBoundary,
								bool frontReflectiveBoundary, bool backReflectiveBoundary,
								// these two will be ignored in the two-dimensional case
								bool bottomReflectiveBoundary, bool topReflectiveBoundary) {
		this.particles = particles;
		this.cutoffDistance = cutoffDistance;
		this.frontLowerLeftCorner = frontLowerLeftCorner;
		this.iterationsPerParticleToCellReassignment = iterationsPerParticleToCellReassignment;
		this.iterationsToCellReassignmentLeft = 1;
		cellSize = new int[dim]
		cellNumbers = new int[dim]
		cellCount = 1;
		for (int i = 0; i < dim; i++) {
			// round down if the number of cells in a dimension is not a round integer
			cellNumbers[i] = (int)(simulationAreaExtent[i] / cutoffDistance);
			cellSize[i] = cellNumbers[i] * cutoffDistance;
			// special case: a dimension is so small that only one cell fits in
			if (cellNumbers[i] == 0) {
				cellNumbers[i] = 1;
				cellSize[i] = simulationAreaExtent[i];
			}
			cellCount *= cellNumbers[i];
		}
		halo = new std::vector<Particle>();
		// TODO: initialize cellPairs correctly
		cellPairs = new List<utils::Vector<int, 2>>;
		if (dim == 2) {
			for (int x = 0; x < cellNumbers[0]; x++)
				for (int y = 0; y < cellNumbers[1]; y++)
					for (int xp = 0; xp < 2; xp++)
						for (int yp = 0; yp < 2; yp++)
							if (x + xp < cellNumbers[0] && y + yp < cellNumbers[1]) {
								utils::Vector<int, 2> pair = new utils::Vector<int, 2>;
								pair[0] = y * cellNumbers[0] + x;
								pair[1] = (y + yp) * cellNumbers[0] + (x + xp);
								cellPairs.add(pair);
							}
		}
		else if (dim == 3) {
			for (int x = 0; x < cellNumbers[0]; x++)
				for (int y = 0; y < cellNumbers[1]; y++)
					for (int z = 0; z < cellNumbers[2]; z++)
						for (int xp = 0; xp < 2; xp++)
							for (int yp = 0; yp < 2; yp++)
								for (int zp = 0; zp < 2; zp++)
									if (x + xp < cellNumbers[0] && y + yp < cellNumbers[1] && z + zp < cellNumbers[2]) {
										utils::Vector<int, 2> pair = new utils::Vector<int, 2>;
										pair[0] = (z * cellNumbers[1] + y) * cellNumbers[0] + x;
										pair[1] = ((z + zp) * cellNumbers[1] + (y + yp)) * cellNumbers[0] + (x + xp);
										cellPairs.add(pair);
									}
		}
		else
			// TODO: throw a sensible error message and log an error as well
			throw new Exception();
		// define the reflective boundary cells
		// TODO: this might need some testing and debugging, to be on the safe side
		reflectiveBoundaryCells = new List<utils::Vector<double, 4>>;
		if (dim == 2) {
			if (leftReflectiveBoundary || rightReflectiveBoundary)
				for (int y = 0; y < cellNumbers[1]; y++) {
					if (leftReflectiveBoundary) {
						utils::Vector<double, 4> vec = new utils::Vector<double, 4>();
						vec[0] = y * cellNumbers[0];
						vec[1] = 0;
						vec[2] = -1.0;
						vec[3] = frontLowerLeftCorner[0];
						reflectiveBoundaryCells.add(vec)
					}
					if (rightReflectiveBoundary) {
						utils::Vector<double, 4> vec = new utils::Vector<double, 4>();
						vec[0] = y * cellNumbers[0] + (cellNumbers[0] - 1);
						vec[1] = 0;
						vec[2] = 1.0;
						vec[3] = frontLowerLeftCorner[0] + simulationAreaExtent[0];
						reflectiveBoundaryCells.add(vec)
					}
				}
			if (frontReflectiveBoundary || backReflectiveBoundary)
				for (int x = 0; x < cellNumbers[0]; x++) {
					if (frontReflectiveBoundary) {
						utils::Vector<double, 4> vec = new utils::Vector<double, 4>();
						vec[0] = x;
						vec[1] = 1;
						vec[2] = -1.0;
						vec[3] = frontLowerLeftCorner[1];
						reflectiveBoundaryCells.add(vec)
					}
					if (backReflectiveBoundary) {
						utils::Vector<double, 4> vec = new utils::Vector<double, 4>();
						vec[0] = (cellNumbers[1] - 1) * cellNumbers[0] + x;
						vec[1] = 1;
						vec[2] = 1.0;
						vec[3] = frontLowerLeftCorner[1] + simulationAreaExtent[1];
						reflectiveBoundaryCells.add(vec)
					}
				}
		}
		else if (dim == 3) {
			if (leftReflectiveBoundary || rightReflectiveBoundary)
				for (int y = 0; y < cellNumbers[1]; y++)
					for (int z = 0; z < cellNumbers[2]; z++) {
						if (leftReflectiveBoundary) {
							utils::Vector<double, 4> vec = new utils::Vector<double, 4>();
							vec[0] = (z * cellNumbers[1] + y) * cellNumbers[0] + 0;
							vec[1] = 0;
							vec[2] = -1.0;
							vec[3] = frontLowerLeftCorner[0];
							reflectiveBoundaryCells.add(vec)
						}
						if (rightReflectiveBoundary) {
							utils::Vector<double, 4> vec = new utils::Vector<double, 4>();
							vec[0] = (z * cellNumbers[1] + y) * cellNumbers[0] + (cellNumbers[0] - 1);
							vec[1] = 0;
							vec[2] = 1.0;
							vec[3] = frontLowerLeftCorner[0] + simulationAreaExtent[0];
							reflectiveBoundaryCells.add(vec)
						}
					}
			if (frontReflectiveBoundary || backReflectiveBoundary)
				for (int x = 0; x < cellNumbers[0]; x++)
					for (int z = 0; z < cellNumbers[2]; z++) {
						if (frontReflectiveBoundary) {
							utils::Vector<double, 4> vec = new utils::Vector<double, 4>();
							vec[0] = (z * cellNumbers[1] + 0) * cellNumbers[0] + x;
							vec[1] = 1;
							vec[2] = -1.0;
							vec[3] = frontLowerLeftCorner[1];
							reflectiveBoundaryCells.add(vec)
						}
						if (backReflectiveBoundary) {
							utils::Vector<double, 4> vec = new utils::Vector<double, 4>();
							vec[0] = (z * cellNumbers[1] + (cellNumbers[1] - 1)) * cellNumbers[0] + x;
							vec[1] = 1;
							vec[2] = 1.0;
							vec[3] = frontLowerLeftCorner[1] + simulationAreaExtent[1];
							reflectiveBoundaryCells.add(vec)
						}
					}
			if (bottomReflectiveBoundary || topReflectiveBoundary)
				for (int x = 0; x < cellNumbers[0]; x++)
					for (int y = 0; y < cellNumbers[1]; y++) {
						if (bottomReflectiveBoundary) {
							utils::Vector<double, 4> vec = new utils::Vector<double, 4>();
							vec[0] = y * cellNumbers[0] + x;
							vec[1] = 2;
							vec[2] = -1.0;
							vec[3] = frontLowerLeftCorner[2];
							reflectiveBoundaryCells.add(vec)
						}
						if (topReflectiveBoundary) {
							utils::Vector<double, 4> vec = new utils::Vector<double, 4>();
							vec[0] = ((cellNumbers[2] -1) * cellNumbers[1] + y) * cellNumbers[0] + x;
							vec[1] = 2;
							vec[2] = 1.0;
							vec[3] = frontLowerLeftCorner[2] + simulationAreaExtent[2];
							reflectiveBoundaryCells.add(vec)
						}
					}
		}
		else
			// TODO: throw a sensible error message and log an error as well
			throw new Exception();
		// assign the particles to their initial cells
		ReassignParticles();
	}

	// applies the reflectove boundary condition to all cells that apply
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

	/// a method that reassigns the particles to the cells they belong to
	/// note that this method will only work every 'iterationsPerParticleToCellReassignment'th call
	/// a good number for this must be determined experimentally for increased efficiency
	/// the default value of this variable for the LinkedCellAlgorithm is one
	void ReassignParticles() {
		iterationsToCellReassignmentLeft--;
		if (iterationsToCellReassignmentLeft != 0)
			return;
		iterationsToCellReassignmentLeft = iterationsPerParticleToCellReassignment;
		// NOTE (TODO): explicitly delete the old cells? before creating new ones?
		Cells = new std::vector<Particle>[cellCount];
		for (int i = 0; i < cellCount; i++)
			Cells[i] = new std::vector<Particle>();
		for (std::vector<Particle>::iterator it = particles.begin() ; it < particles.end(); it++)
		{
			Particle& p = *it;
			if (dim == 2) {
				xIndex = (int)((p.x[0] - frontLowerLeftCorner[0]) / cellSize[0]);
				yIndex = (int)((p.y[1] - frontLowerLeftCorner[1]) / cellSize[1]);
				if (xIndex < 0 || yIndex < 0 ||
					xIndex >= cellNumbers[0] || yIndex >= cellNumbers[1])
					halo.push_back(p);
				else {
					int index = yIndex * cellNumbers[0] + xIndex;
					cells[index].push_back(p);
				}
			}
			else if (dim == 3) {
				xIndex = (int)((p.x[0] - frontLowerLeftCorner[0]) / cellSize[0]);
				yIndex = (int)((p.y[1] - frontLowerLeftCorner[1]) / cellSize[1]);
				zIndex = (int)((p.y[2] - frontLowerLeftCorner[2]) / cellSize[2]);
				if (xIndex < 0 || yIndex < 0 || zIndex < 0 ||
					xIndex >= cellNumbers[0] || yIndex >= cellNumbers[1] || zIndex >= cellNumbers[2])
					halo.push_back(p);
				else {
					int index = (zIndex * cellNumbers[1] + yIndex) * cellNumbers[0] + xIndex;
					cells[index].push_back(p);
				}
			}
			else
				// TODO: throw a sensible error message and log an error as well
				throw new Exception();
		}
	}

	/// a method to add a Particle to the ListParticleContainer
	void AddParticle(const Particle& particle);

	/// a method that takes a void(*func)(void*, Particle) and uses it to iterate over all Particles
	/// @param data additional data given to func
	void Iterate(void(*func)(void*, Particle&), void *data) {
		for (std::vector<Particle>::iterator it = particles.begin() ; it < particles.end(); it++)
		{
			Particle& p = *it;
			(*func)(data, p);
		}
	}
	
	/// a method that takes a void(*func)(void*, Particle, Particle) and uses it to iterate over all pairs of Particles (each symmetrical pair is only taken once to reduce redundancy)
	/// @param data additional data given to func
	void IteratePairwise(void(*func)(void*, Particle&, Particle&), void *data) {// iterate over all Particles
		// TODO: iterate over all elements of cellPairs. This depends on the data type of cellPairs
		for () {
			std::vector<Particle> cell1 = cells[elem[0]];
			std::vector<Particle> cell2 = cells[elem[1]];
			for (std::vector<Particle>::iterator it1 = cell1.begin() ; it1 < cell1.end(); it1++)
				for (std::vector<Particle>::iterator it2 = cell2.begin() ; it2 < cell2.end(); it2++)
					// make sure that a Particle is not paired with itself
					if (it1 != it2)
					{
						// call the function on the pair of Particles
						Particle& p1 = *it1;
						Particle& p2 = *it2;
						(*func)(data, p1, p2);
					}
		}
	}

	/// add particles from *.txt file
	void AddParticlesFromFile(const char *filename);

	/// our new fileformat, replace later AddParticlesFromFile
	/// @return return true if file could be read
	bool AddParticlesFromFileNew(const char *filename);

	/// removes all particles
	void Clear();

	/// are any particles contained?
	bool IsEmpty();

	/// returns how many particles are contained in this container
	/// @return returns count of particles in this container
	unsigned int					getParticleCount();

	// this method shall be later removed...
	///returns ListParticleContainer's internal container
	const std::vector<Particle>& getParticles();
	
	~LinkedCellParticleContainer()
	{
		SAFE_DELETE_A(Cells);
	}
};


// implement here...

template<typename T> void LinkedCellParticleContainer<T>::AddParticle(const Particle& particle)
{
	// functions have to look like this...
	// first template<typename T> ... LinkedCellParticleContainer<T>::...
}

#endif 

