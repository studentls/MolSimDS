//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File LinkedCellParticleContainer.cpp
// contains class LinkedCellParticleContainer
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------

#include <cstdio>
#include <cstring>

#include "LinkedCellParticleContainer.h"

void LinkedCellParticleContainer::Iterate(void(*func)(void*, Particle&), void *data)
{		
	assert(Cells);

	// go through all Cells...
	for(int i = 0; i < getCellCount(); i++)
	{
		if(Cells[i].empty())continue;

		for(std::vector<Particle>::iterator it = Cells[i].begin(); it != Cells[i].end(); it++)
		{
			Particle& p = *it;
			(*func)(data, p);				
		}
	}
}

void LinkedCellParticleContainer::IteratePairwise(void(*func)(void*, Particle&, Particle&), void *data)
{
	// Step 1: calc pairwise force between particles for each cell

	// go through cells
	for(unsigned int i = 0; i < getCellCount(); i++)
	{
		// for all particles in cell_i, calc forces in cell_i
		for(std::vector<Particle>::iterator it1 = Cells[i].begin(); it1 != Cells[i].end(); it1++)
		{
			for(std::vector<Particle>::iterator it2 = it1 + 1; it2 != Cells[i].end(); it2++)
			{
				Particle& p1 = *it1;
				Particle& p2 = *it2;

				func(data, p1, p2);
			}
		}
	}

	// Step 2: calc force with neighbouring cells
	// use stored pairs to iterate, these may be ordered to increase performance in a special way
	// ramaining an unsolved problem for simulations yet
	
	// pairs should exist!
	assert(!cellPairs.empty());

	for(std::vector<utils::Vector<unsigned int, 2> >::iterator it = cellPairs.begin(); it != cellPairs.end(); it++)
	{
		int i = (*it)[0];
		int j = (*it)[1];
		
		// only for one combination
		assert(i < j);
		
		//interact cell_i and cell_j

		// for all particles in cell_i
		for(std::vector<Particle>::iterator it1 = Cells[i].begin(); it1 != Cells[i].end(); it1++)
		{		
			// calc force, based on actio / reaction between cell_i and cell_j
			for(std::vector<Particle>::iterator it2 = Cells[j].begin(); it2 != Cells[j].end(); it2++)
			{
				Particle& p1 = *it1;
				Particle& p2 = *it2;

				func(data, p1, p2);
			}			
		}
	}
}

std::vector<Particle> LinkedCellParticleContainer::getBoundaryParticles()
{
	std::vector<Particle> particles;

	// case 2D
	if(dim == 2)
	{
		// |.........|
		// |         |
		// |.........|
		for(int x = 0; x < cellCount[0]; x++)
			for(int y = 0; y < cellCount[1]; y++)
			{
				// only boundaries are allowed
				if(x != 0 && y != 0 && x != cellCount[0] && y != cellCount[1])continue;

				// add to container
				int index = Index2DTo1D(x, y);
				particles.insert(particles.end(), Cells[index].begin(), Cells[index].end());
			}
	}

	// case 3D
	if(dim == 3)
	{
		// |.........|
		// |         |
		// |.........|
		for(int x = 0; x < cellCount[0]; x++)
			for(int y = 0; y < cellCount[1]; y++)
				for(int z = 0; z < cellCount[2]; z++)
				{
					// only boundaries are allowed
					if(x != 0 && y != 0 && z != 0 && x != cellCount[0] && y != cellCount[1] && z != cellCount[2])continue;

					// add to container
					int index = Index3DTo1D(x, y, z);
					particles.insert(particles.end(), Cells[index].begin(), Cells[index].end());
				}
	}

	return particles;
}

void LinkedCellParticleContainer::ApplyReflectiveBoundaryConditions(void(*func)(void*, Particle&, Particle&), void *data)
{
	// go through boundaries...
	for (std::vector<utils::Vector<double, 4> >::iterator it = reflectiveBoundaryCells.begin(); it != reflectiveBoundaryCells.end(); it++)
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

	// Annotation:
	// In Task 3 it is mentioned that each boundary shall specified if it reflects or outflows the particle
	// as there is currently no interaction with halo at all, the particles will be simply deleted if a global
	// outflow flag is set

	// if BC_OUTFLOW is set, delete halo particles
	if(boundaryConditions & BC_OUTFLOW)clearHaloParticles();
}

void	LinkedCellParticleContainer::generatePairs()
{
	// if necessary, delete pair container
	if(!cellPairs.empty())cellPairs.clear();

	// to create pairs, we go through all points in the grid and determine possible neighbours
	// let i1 be the 1D-index of the current point, and i2 a index of a possible neighbour
	// if i1 < i2 the pair will be added
	// note that boundaries shall be treated properly
	switch(dim)
	{
	case 2:
		{
			// go through grid points x, y
			for(int x = 0; x < cellCount[0]; x++)
				for(int y = 0; y < cellCount[1]; y++)
				{
					// calc index for point (x,y)
					unsigned int i1 = Index2DTo1D(x, y);

					// get neighbours
					std::vector<unsigned int> neighbours = getNeighbours2D(i1);

					// go through neighbours
					assert(!neighbours.empty());
					for(std::vector<unsigned int>::iterator it = neighbours.begin(); it < neighbours.end(); it++)
					{
						unsigned int i2 = *it;

						// add pair if i1 < i2
						if(i1 < i2)cellPairs.push_back(makePair(i1, i2));
					}
				}

			break;
		}
	case 3:
		{
			// go through grid points x, y, z
			for(int x = 0; x < cellCount[0]; x++)
				for(int y = 0; y < cellCount[1]; y++)
					for(int z = 0; z < cellCount[2]; z++)
					{
						// calc index for point (x,y, z)
						unsigned int i1 = Index3DTo1D(x, y, z);

						// get neighbours
						std::vector<unsigned int> neighbours = getNeighbours3D(i1);

						// go through neighbours
						assert(!neighbours.empty());
						for(std::vector<unsigned int>::iterator it = neighbours.begin(); it < neighbours.end(); it++)
						{
							unsigned int i2 = *it;

							// add pair if i1 < i2
							if(i1 < i2)cellPairs.push_back(makePair(i1, i2));
						}
					}

			break;
		}
	default:
		LOG4CXX_ERROR(generalOutputLogger, "failed to calculate pairs, as only dimensions 2, 3 are supported yet");
	}
}


void	LinkedCellParticleContainer::SetReflectiveBoundaries(bool leftReflectiveBoundary, bool rightReflectiveBoundary,
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
					this->reflectiveBoundaryCells.push_back(vec);
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
					this->reflectiveBoundaryCells.push_back(vec);
				}
				if (backReflectiveBoundary) {
					// temp variable
					utils::Vector<double, 4> vec;
					vec[0] = Index2DTo1D(x, cellCount[1] - 1);
					vec[1] = 1;
					vec[2] = 1.0;
					vec[3] = frontLowerLeftCorner[1] + calcSimulationAreaExtent()[1];
					this->reflectiveBoundaryCells.push_back(vec);
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
						this->reflectiveBoundaryCells.push_back(vec);
					}
					if (rightReflectiveBoundary) {
						// temp variable
						utils::Vector<double, 4> vec;
						vec[0] = Index3DTo1D(cellCount[0] - 1, y, z);
						vec[1] = 0;
						vec[2] = 1.0;
						vec[3] = frontLowerLeftCorner[0] + calcSimulationAreaExtent()[0];
						this->reflectiveBoundaryCells.push_back(vec);
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
						this->reflectiveBoundaryCells.push_back(vec);
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
