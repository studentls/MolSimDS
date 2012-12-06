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

	// go through cells
	for(unsigned int i = 0; i < getCellCount(); i++)
	{
		// for all particles in cell_i, calc forces with neighbours
		for(std::vector<Particle>::iterator it1 = Cells[i].begin(); it1 != Cells[i].end(); it1++)
		{

			// get neighbours of cell_i with little helper function
			std::vector<unsigned int> neighbours = getNeighbours(i);

			// go through neighbours
			for(std::vector<unsigned int>::iterator nt = neighbours.begin(); nt != neighbours.end(); nt++)
			{
				unsigned int j = *nt;

				// calc force, based on actio / reaction between cell_i and cell_j
				for(std::vector<Particle>::iterator it2 = Cells[j].begin(); it2 != Cells[j].end(); it2++)
				{
					// only for one combination
					if(i < j)
					{
						Particle& p1 = *it1;
						Particle& p2 = *it2;

						func(data, p1, p2);
					}
				}
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

	// 3D not yet supported...

	return particles;
}

void LinkedCellParticleContainer::ApplyReflectiveBoundaryConditions(void(*func)(void*, Particle&, Particle&), void *data)
{
	// go through boundaries...
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

	// if BC_OUTFLOW is set, delete halo particles
	if(boundaryConditions & BC_OUTFLOW)clearHaloParticles();
}