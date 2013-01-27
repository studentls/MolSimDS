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
#ifdef OPENMP
	//openMP version
	
#pragma omp parallel for num_threads(threadCount)
	for(int i = 0; i < getCellCount(); ++i)
	{
		if(Cells[i].empty())continue;

		for(int j = 0; j < Cells[i].size(); ++j)
		{
			Particle& p = Cells[i][j];
			(*func)(data, p);				
		}
	}
#else
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
#endif
}

void LinkedCellParticleContainer::IteratePairwise(void(*func)(void*, Particle&, Particle&), void *data)
{
#ifndef OPENMP
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
#else
	//Open MP version
	// Step 1: calc pairwise force between particles for each cell

	// go through cells
#pragma omp parallel for num_threads(threadCount)
	for(int i = 0; i < getCellCount(); i++)
	{
		// for all particles in cell_i, calc forces in cell_i
		for(int j = 0; j < Cells[i].size(); j++)
		{
			for(int k = j + 1; k < Cells[i].size(); k++)
			{
				Particle& p1 = Cells[i][j];
				Particle& p2 = Cells[i][k];

				func(data, p1, p2);
			}
		}
	}

	// Step 2: calc force with neighbouring cells
	// use stored pairs to iterate, these may be ordered to increase performance in a special way
	// ramaining an unsolved problem for simulations yet
	
	// pairs should exist!
	assert(oddStrips.size() != 0);
	assert(evenStrips.size() != 0);

#pragma omp parallel for num_threads(threadCount)
	// 2.1 go through oddStrips and calc forces
	for(int iStrip = 0; iStrip < oddStrips.size(); iStrip++)
	{
		// go through all pairs in the strip
		for(int iPair = 0; iPair < oddStrips[iStrip].size(); iPair++)
		{
			int i = oddStrips[iStrip][iPair][0];
			int j = oddStrips[iStrip][iPair][1];

			// if one cell is empty go to next pair
			if(Cells[i].empty() || Cells[j].empty())continue;

			// for all particles in cell_i
			for(int  k = 0; k < Cells[i].size(); k++)
			{		
				// calc force, based on actio / reaction between cell_i and cell_j
				for(int l = 0; l < Cells[j].size(); l++)
				{
					Particle& p1 = Cells[i][k];
					Particle& p2 = Cells[j][l];

					func(data, p1, p2);
				}			
			}
		}
	}

#pragma omp parallel for num_threads(threadCount)
	// 2.2 go through evenStrips and calc forces
	for(int iStrip = 0; iStrip < evenStrips.size(); iStrip++)
	{
		// go through all pairs in the strip
		for(int iPair = 0; iPair < evenStrips[iStrip].size(); iPair++)
		{
			int i = evenStrips[iStrip][iPair][0];
			int j = evenStrips[iStrip][iPair][1];

			// if one cell is empty go to next pair
			if(Cells[i].empty() || Cells[j].empty())continue;

			// for all particles in cell_i
			for(int  k = 0; k < Cells[i].size(); k++)
			{		
				// calc force, based on actio / reaction between cell_i and cell_j
				for(int l = 0; l < Cells[j].size(); l++)
				{
					Particle& p1 = Cells[i][k];
					Particle& p2 = Cells[j][l];

					func(data, p1, p2);
				}			
			}
		}
	}
	

#endif
}

void LinkedCellParticleContainer::ApplyReflectiveBoundaryConditions(void(*func)(void*, Particle&, Particle&), void *data)
{
	using namespace std;
	using namespace utils;
	
	
	
		// go through boundary cells
		for(vector<unsigned int>::iterator it = boundaryIndices.begin(); it != boundaryIndices.end(); it++)
		{
			if(Cells[*it].empty())continue;

			// go through boundary cell
			for(vector<Particle>::iterator pt = Cells[*it].begin(); pt != Cells[*it].end(); pt++)
			{
				
				// go through boundaries (max. 6)
				for(vector<Plane>::iterator bt = reflectiveBoundaries.begin(); bt != reflectiveBoundaries.end(); bt++)
				{
					double dist = bt->distance(pt->x);
				
					// skip the particle if it is too far away from the border
					if (abs(dist) > reflectiveBoundaryDistance)
						continue;

					    // create a temporary, virtual Particle
						// a copy of the current particle, but located on the boundary
						Particle vp(*pt);
						vp.x = vp.x - dist * bt->n;
						(*func)(data, *pt, vp);								
			}
		}
	}

	// Annotation:
	// In Task 3 it is mentioned that each boundary shall specified if it reflects or outflows the particle
	// as there is currently no interaction with halo at all, the particles will be simply deleted if a global
	// outflow flag is set
	
	// TODO: test if periodic conditions are there, cause clearHaloParticles fails in this case!

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

					if(neighbours.empty())
					{
						LOG4CXX_ERROR(generalOutputLogger, "failure!");
					}
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



void LinkedCellParticleContainer::SetPeriodicBoundaries(bool xAxis, bool yAxis, bool zAxis)
{
	// clear container if not empty
	if(!periodicBoundaryGroups.empty())periodicBoundaryGroups.clear();

	// set axises
	periodicBoundaries[0] = xAxis;
	periodicBoundaries[1] = yAxis;
	periodicBoundaries[2] = zAxis;

	// get extent of simulation area
	utils::Vector<double, 3> SimulationAreaExtent = calcSimulationAreaExtent();

	if (dim == 2)
	{
		for (int x1 = 0; x1 < cellCount[0]; x1++)
			for (int y1 = 0; y1 < cellCount[1]; y1++)
				for (int x2 = 0; x2 < cellCount[0]; x2++)
					for (int y2 = 0; y2 < cellCount[1]; y2++)
					{
						int index1 = Index2DTo1D(x1, y1);
						int index2 = Index2DTo1D(x2, y2);

						// ensure that every pair is only taken once
						if (index1 <= index2)
							continue;

						// for every possibility of each dimension
						// being either a direct neighbour
						// or a neighbour through the periodic boundary
						for (int xOpp = 0; xOpp < 2; xOpp++)
							for (int yOpp = 0; yOpp < 2; yOpp++)
							{
								// presume that relation, then check later if it is actually true
								// this is necessary because the relation can't be determined directly
								// if there are too few cells, so that there are several valid relations
								// for each pair of cells
								bool xOpposite = (xOpp == 0);
								bool yOpposite = (yOpp == 0);
								// skip those where the cells are direct neighbours in every dimension
								if (!xOpposite && !yOpposite)// && !zOpposite)
									continue;
								// for every axis, skip the pair if
								// the above presumption regarding the positions
								// turns out to be wrong for the pair of particles
								if (xOpposite && !((x1 == 1 && x2 == cellCount[0] - 2) || (x2 == 1 && x1 == cellCount[0] - 2)) ||
									!xOpposite && !(x1 == x2 + 1 || x1 == x2 - 1 || x1 == x2))
									continue;
								if (yOpposite && !((y1 == 1 && y2 == cellCount[1] - 2) || (y2 == 1 && y1 == cellCount[1] - 2)) ||
									!yOpposite && !(y1 == y2 + 1 || y1 == y2 - 1 || y1 == y2))
									continue;
								
								// if the algorithm gets till here, it is a pair of indirect neighbours. store it.
								PeriodicBoundary b;
								b.cell1 = index1;
								b.cell2 = index2;
								b.xAxis = xOpposite ? (SimulationAreaExtent[0] * ((x1 == 1) ? 1.0 : -1.0)) : 0;
								b.yAxis = yOpposite ? (SimulationAreaExtent[1] * ((y1 == 1) ? 1.0 : -1.0)) : 0;
								b.zAxis = 0.0;

								periodicBoundaryGroups.push_back(b);
							}
					}
	}
	else if (dim == 3)
	{
		for (int x1 = 0; x1 < cellCount[0]; x1++)
			for (int y1 = 0; y1 < cellCount[1]; y1++)
				for (int z1 = 0; z1 < cellCount[2]; z1++)
					for (int x2 = 0; x2 < cellCount[0]; x2++)
						for (int y2 = 0; y2 < cellCount[1]; y2++)
							for (int z2 = 0; z2 < cellCount[2]; z2++)
							{
								int index1 = Index2DTo1D(x1, y1);
								int index2 = Index2DTo1D(x2, y2);

								// ensure that every pair is only taken once
								if (index1 <= index2)
									continue;

								// for every possibility of each dimension
								// being either a direct neighbour
								// or a neighbour through the periodic boundary
								for (int xOpp = 0; xOpp < 2; xOpp++)
									for (int yOpp = 0; yOpp < 2; yOpp++)
										for (int zOpp = 0; zOpp < 2; zOpp++)
										{
											// presume that relation, then check later if it is actually true
											// this is necessary because the relation can't be determined directly
											// if there are too few cells, so that there are several valid relations
											// for each pair of cells
											bool xOpposite = (xOpp == 0);
											bool yOpposite = (yOpp == 0);
											bool zOpposite = (zOpp == 0);
											// skip those where the cells are direct neighbours in every dimension
											if (!xOpposite && !yOpposite && !zOpposite)
												continue;
											// for every axis, skip the pair if
											// the above presumption regarding the positions
											// turns out to be wrong for the pair of particles
											if (xOpposite && !((x1 == 1 && x2 == cellCount[0] - 2) || (x2 == 1 && x1 == cellCount[0] - 2)) ||
												!xOpposite && !(x1 == x2 + 1 || x1 == x2 - 1 || x1 == x2))
												continue;
											if (yOpposite && !((y1 == 1 && y2 == cellCount[1] - 2) || (y2 == 1 && y1 == cellCount[1] - 2)) ||
												!yOpposite && !(y1 == y2 + 1 || y1 == y2 - 1 || y1 == y2))
												continue;
											if (zOpposite && !((z1 == 1 && z2 == cellCount[2] - 2) || (z2 == 1 && z1 == cellCount[2] - 2)) ||
												!zOpposite && !(z1 == z2 + 1 || z1 == z2 - 1 || z1 == z2))
												continue;
								
											// if the algorithm gets till here, it is a pair of indirect neighbours. store it.
											PeriodicBoundary b;
											b.cell1 = index1;
											b.cell2 = index2;
											b.xAxis = xOpposite ? (SimulationAreaExtent[0] * ((x1 == 1) ? 1.0 : -1.0)) : 0;
											b.yAxis = yOpposite ? (SimulationAreaExtent[1] * ((y1 == 1) ? 1.0 : -1.0)) : 0;
											b.zAxis = zOpposite ? (SimulationAreaExtent[2] * ((z1 == 1) ? 1.0 : -1.0)) : 0;
											periodicBoundaryGroups.push_back(b);
										}
							}
	}
}

void LinkedCellParticleContainer::ApplyPeriodicBoundaryConditions(void(*func)(void*, Particle&, Particle&), void *data)
{
	for (std::vector<PeriodicBoundary>::iterator it = periodicBoundaryGroups.begin(); it != periodicBoundaryGroups.end(); it++)
	{
		PeriodicBoundary& b = *it;

		std::vector<Particle>& cell1 = Cells[b.cell1];
		std::vector<Particle>& cell2 = Cells[b.cell2];
		double axis1 = b.xAxis;
		double axis2 = b.yAxis;
		double axis3 = b.zAxis;

		// interact cells
		// to achieve this, temporarily shift coordinates
		for (std::vector<Particle>::iterator it = cell2.begin() ; it < cell2.end(); it++)
		{
			Particle& p2 = *it;

			// start shift
			p2.x[0] -= axis1;
			p2.x[1] -= axis2;
			p2.x[2] -= axis3;

			for (std::vector<Particle>::iterator it = cell1.begin() ; it < cell1.end(); it++)
			{
				Particle& p1 = *it;
				func(data, p1, p2);
			}

			// end shift
			p2.x[0] += axis1;
			p2.x[1] += axis2;
			p2.x[2] += axis3;
		}
	}
}

void LinkedCellParticleContainer::ReassignHaloForPeriodicConditions()
{
	// skip this function if no periodic boundaries are used
	if (!periodicBoundaries[0] && !periodicBoundaries[1] && !periodicBoundaries[2])
		return;

	utils::Vector<double, 3> SimulationAreaExtent = calcSimulationAreaExtent();

	// go through all Cells...
	for(int i = 0; i < getCellCount(); i++)
	{
		if(Cells[i].empty())continue;

		for(std::vector<Particle>::iterator it = Cells[i].begin(); it != Cells[i].end(); it++)
		{
			Particle& p = *it;
			for (int i = 0; i < dim; i++)
				if (periodicBoundaries[i])
				{
					if (p.x[i] < this->frontLowerLeftCorner[i])
						p.x[i] += SimulationAreaExtent[i];
					else if (p.x[i] > this->frontLowerLeftCorner[i] + SimulationAreaExtent[i])
						p.x[i] -= SimulationAreaExtent[i];
				}
		}
	}
}

void	LinkedCellParticleContainer::SetReflectiveBoundaries()
{
	using namespace utils;


	Plane boundary;

	// insert boundaries
	// additionally create normals of planes such as all normals point into the cuboidal simulation area

	if(boundaryConditions & BC_LEFT)
	{		
		boundary.constructFromPoint(Vector<double, 3>(frontLowerLeftCorner[0], 0.0, 0.0), Vector<double, 3>(1.0, 0.0, 0.0));		
		reflectiveBoundaries.push_back(boundary);
	}
	if(boundaryConditions & BC_RIGHT)
	{		
		boundary.constructFromPoint(Vector<double, 3>(frontLowerLeftCorner[0] + (cellCount[0] - 2) * cellSize[0], 0.0, 0.0), Vector<double, 3>(-1.0, 0.0, 0.0));
		reflectiveBoundaries.push_back(boundary);
	}
	if(boundaryConditions & BC_FRONT)
	{		
		boundary.constructFromPoint(Vector<double, 3>(0.0, frontLowerLeftCorner[1], 0.0), Vector<double, 3>(0.0, 1.0, 0.0));
		reflectiveBoundaries.push_back(boundary);
	}
	if(boundaryConditions & BC_BACK)
	{		
		boundary.constructFromPoint(Vector<double, 3>(0.0, frontLowerLeftCorner[1] + (cellCount[1] - 2) * cellSize[1], 0.0), Vector<double, 3>(0.0, -1.0, 0.0));
		reflectiveBoundaries.push_back(boundary);
	}

	if(dim > 2)
	{	
		if(boundaryConditions & BC_BOTTOM)
		{		
			boundary.constructFromPoint(Vector<double, 3>(0, 0.0, frontLowerLeftCorner[2]), Vector<double, 3>(0.0, 0.0, 1.0));
			reflectiveBoundaries.push_back(boundary);
		}
		if(boundaryConditions & BC_TOP)
		{		
			boundary.constructFromPoint(Vector<double, 3>(0, 0.0, frontLowerLeftCorner[2] + (cellCount[2] - 2) * cellSize[2]), Vector<double, 3>(0.0, 0.0, -1.0));
			reflectiveBoundaries.push_back(boundary);
		}
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

	// if OpenMP calc indices
#ifdef OPENMP
	calcMTIndices();
#endif
}

void LinkedCellParticleContainer::calcMTIndices()
{
	// only for 2D at the moment!
	if(dim == 2)
	{

		// (1) calc Indices for even stripes
		for(int i = 0; i < cellCount[0] - 1; i += 2)
		{
			IndexStrip strip;
			strip.constructVerticalStripIndices(i, cellCount[1], cellCount);

			evenStrips.push_back(strip);
		}

		// (2) calc Indices for odd stripes
		for(int i = 1; i < cellCount[0] - 1; i += 2)
		{
			IndexStrip strip;
			strip.constructVerticalStripIndices(i, cellCount[1], cellCount);

			oddStrips.push_back(strip);
		}
	}
	else
	{
		//...
		LOG4CXX_INFO(generalOutputLogger, "3D MT not yet implemented!");
	}
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
