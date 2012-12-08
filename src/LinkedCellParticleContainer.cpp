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

		// if one cell is empty go to next pair
		if(Cells[i].empty() || Cells[j].empty())continue;

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

void LinkedCellParticleContainer::ApplyReflectiveBoundaryConditions(void(*func)(void*, Particle&, Particle&), void *data)
{
	using namespace std;
	using namespace utils;
	
	// go through boundaries (max. 6)
	for(vector<Boundary>::iterator bt = boundaries.begin(); bt != boundaries.end(); bt++)
	{
		// go through boundary cells
		for(vector<unsigned int>::iterator it = boundaryIndices.begin(); it != boundaryIndices.end(); it++)
		{
			if(Cells[*it].empty())continue;

			// go through boundary cell
			for(vector<Particle>::iterator pt = Cells[*it].begin(); pt != Cells[*it].end(); pt++)
			{
				double dist = bt->p.distance(pt->x);
				// skip the particle if it is too far away from the border
				if (dist > reflectiveBoundaryDistance)
					continue;
				// create a temporary, virtual Particle
				// a copy of the current particle, but located at the border
				/*Particle vp(*pt);
				vp.x[axis] = border;
				(*func)(data, p, vp);*/

			}
		}
	}
	
	//
	//// go through boundaries...
	//for (std::vector<utils::Vector<double, 4> >::iterator it = reflectiveBoundaryCells.begin(); it != reflectiveBoundaryCells.end(); it++)
	//{
	//	// TODO: this should obviously work by reference to the array, not copy it
	//	// does it do that right now?
	//	utils::Vector<double, 4>& elem = *it;
	//	std::vector<Particle>& cell = Cells[(int)(elem[0])];
	//	int axis = (int)(elem[1]);
	//	double direction = elem[2];
	//	double border = elem[3];
	//	for (std::vector<Particle>::iterator it = cell.begin() ; it < cell.end(); it++) {
	//		Particle& p = *it;
	//		double dist = direction * (border - p.x[axis]);
	//		// skip the particle if it is too far away from the border
	//		if (dist > reflectiveBoundaryDistance)
	//			continue;
	//		// create a temporary, virtual Particle
	//		// a copy of the current particle, but located at the border
	//		Particle vp(p);
	//		vp.x[axis] = border;
	//		(*func)(data, p, vp);
	//		// TODO: delete vp manually?
	//	}
	//}

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
	using namespace utils;

	Boundary boundary;
	boundary.type = BT_REFLECTIVE;

	// insert boundaries
	if(leftReflectiveBoundary)
	{		
		boundary.p = Plane(Vector(frontLowerLeftCorner[0], 0.0, 0.0), Vector(1.0, 0.0, 0.0));
		boundaries.push_back(boundary);
	}
	if(rightReflectiveBoundary)
	{		
		boundary.p = Plane(Vector(frontLowerLeftCorner[0] + cellCount[0] * cellSize[0], 0.0, 0.0), Vector(1.0, 0.0, 0.0));
		boundaries.push_back(boundary);
	}
	if(frontReflectiveBoundary)
	{		
		boundary.p = Plane(Vector(0.0, frontLowerLeftCorner[1], 0.0), Vector(0.0, 1.0, 0.0));
		boundaries.push_back(boundary);
	}
	if(backReflectiveBoundary)
	{		
		boundary.p = Plane(Vector(0.0, frontLowerLeftCorner[1] + cellCount[1] * cellSize[1], 0.0), Vector(0.0, 1.0, 0.0));
		boundaries.push_back(boundary);
	}
	if(bottomReflectiveBoundary)
	{		
		boundary.p = Plane(Vector(0, 0.0, frontLowerLeftCorner[2]), Vector(0.0, 0.0, 1.0));
		boundaries.push_back(boundary);
	}
	if(topReflectiveBoundary)
	{		
		boundary.p = Plane(Vector(0, 0.0, frontLowerLeftCorner[2] + cellCount[2] * cellSize[2]), Vector(0.0, 0.0, 1.0));
		boundaries.push_back(boundary);
	}

	// note: for performance reasons, maybe index lists can be now generated to make iteration faster
}

void LinkedCellParticleContainer::clearHaloParticles()
{
	// check if indices exist
	assert(!haloIndices.empty());

	// go through halo cells, and clear cells if not empty
	for(std::vector<unsigned int>::iterator it = haloIndices.begin(); it != haloIndices.end(); it++)
	{
		// assert index
		assert(*it < getCellCount());

		if(Cells[*it].empty())continue;

		Cells[*it].clear();
	}
}

std::vector<Particle> LinkedCellParticleContainer::getHaloParticles()
{
	std::vector<Particle> particles;

	// check if indices exist
	assert(!haloIndices.empty());

	// go through halo cells, and add to particles if not empty
	for(std::vector<unsigned int>::iterator it = haloIndices.begin(); it != haloIndices.end(); it++)
	{
		// assert index
		assert(*it < getCellCount());

		if(Cells[*it].empty())continue;

		particles.insert(particles.end(), Cells[*it].begin(), Cells[*it].end());
	}

	return particles;
}

std::vector<Particle> LinkedCellParticleContainer::getBoundaryParticles()
{
	std::vector<Particle> particles;

	// check if indices exist
	assert(!haloIndices.empty());

	// go through boundary cells, and add to particles if not empty
	for(std::vector<unsigned int>::iterator it = boundaryIndices.begin(); it != boundaryIndices.end(); it++)
	{
		// assert index
		assert(*it < getCellCount());

		if(Cells[*it].empty())continue;

		particles.insert(particles.end(), Cells[*it].begin(), Cells[*it].end());
	}

	return particles;
}

void LinkedCellParticleContainer::calcIndices()
{
	// first generate Pairs
	generatePairs();

	// second calculate indices for halo & boundary cells
	if(!haloIndices.empty())haloIndices.clear();
	if(!boundaryIndices.empty())boundaryIndices.clear();

	calcFrameIndices(haloIndices, 0);
	calcFrameIndices(boundaryIndices, 1);
}

void LinkedCellParticleContainer::calcFrameIndices(std::vector<unsigned int> &out, const unsigned int r)
{
	int m = cellCount[0];
	int n = cellCount[1];
	int o = cellCount[2];

	// getting the r-th layer of a grid

	// calc Boundary Cells' indices similiar to halo
	switch(dim)
	{
	case 2:
		{
			// grid |------------|
			//      |            |
			//      |------------|
			// note that construction ensures, that halo layer has at least 3 cells in each direction!

			// ------- upper
			for(int x = r + 1; x <= m - r - 2; x++)
				out.push_back(Index3DTo1D(x, r, 0));

			// handle special case
			if(r != n - r - 1)
				// ------- lower
				for(int x = r + 1; x <= m - r - 2; x++)
					out.push_back(Index3DTo1D(x, n - r - 1, 0));
			// |
			// |
			// | left
			for(int y = r; y <= n - r - 1; y++)
				out.push_back(Index3DTo1D(r, y, 0));
			
			// handle special case
			if(r != m - r - 1)
				//        |
				//        |
				// right  | 
				for(int y = r; y <= n - r - 1; y++)
					out.push_back(Index3DTo1D(m - r - 1, y, 0));

			break;
		}
	case 3:
		{
			// do it the same way as in 2D

			// for a cube there are 6 faces!
			// top face
			for(int x = r; x <= m - r - 1; x++)
				for(int y = r; y <= n - r - 1; y++)
					out.push_back(Index3DTo1D(x, y, r));	

			// handle special case
			if(r != o - r - 1)
				// bottom face
				for(int x = r; x <= m - r - 1; x++)
					for(int y = r; y <= n - r - 1; y++)
						out.push_back(Index3DTo1D(x, y, o - r - 1));

			
			//front face
			for(int x = r; x <= m - r - 1; x++)
				for(int z = r + 1; z <= o - r - 2; z++)
					out.push_back(Index3DTo1D(x, r, z));
			
			// handle special case
			if(r != n - r - 1)
				//back face
				for(int x = r; x <= m - r - 1; x++)
					for(int z = r + 1; z <= o - r - 2; z++)
						out.push_back(Index3DTo1D(x, n - r - 1, z));

			
			//left face
			for(int y = r + 1; y <= n - r - 2; y++)
				for(int z = r + 1; z <= o - r - 2; z++)
					out.push_back(Index3DTo1D(r, y, z));

			// handle special case
			if(r != m - r - 1)
				//right face
				for(int y = r + 1; y <= n - r - 2; y++)
					for(int z = r + 1; z <= o - r - 2; z++)
						out.push_back(Index3DTo1D(m - r - 1, y, z));

			break;
		}
	default:
		LOG4CXX_ERROR(generalOutputLogger, "failed to calculate pairs, as only dimensions 2, 3 are supported yet");
	}
}
