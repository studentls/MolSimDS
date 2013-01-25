/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <vector>
#include "utils/Vector.h"


/// a class that is used to represent a particle and all its attributes and physical circumstances
class Particle {

private:
	/// the force effective on this particle
	utils::Vector<double, 3> f;

	/// the force which was effective on this particle
	utils::Vector<double, 3> old_f;
	
	/// the indices of the direct neighbors of the Particle, if in a membrane
	/// note that for performance reasons only the neighbors with a lower index are stored
	/// to reduce redundancy
	std::vector<int> directNeighbors;

	/// the indices of the diagonal neighbors of the Particle, if in a membrane
	/// note that for performance reasons only the neighbors with a lower index are stored
	/// to reduce redundancy
	std::vector<int> diagonalNeighbors;

public:
	/// a constructor that takes and sets its type
	Particle(int type = 0);

	/// a constructor that clones an existing Particle
	Particle(const Particle& other);

	/// a full constructor that takes the position and velocity (both as 3-tuple vectors) as well as the mass and type
	Particle(
			/// for visualization, we need always 3 coordinates
			/// -> in case of 2d, we use only the first and the second
			const utils::Vector<double, 3>& x_arg,
	        const utils::Vector<double, 3>& v_arg,
	        double m_arg,
	        int type = 0
	);

	/// destructor
	~Particle()  		{}
	
	/// the position of the particle
	utils::Vector<double, 3> x;

	/// the velocity of the particle
	utils::Vector<double, 3> v;

	/// the mass of this particle
	double m;

	/// type of the particle. The use of this value depends on particularities of the implementation. Currently not in use
	int type;

	/// get the current force acting on the Particle
	utils::Vector<double, 3>& getF();

	/// get the force that acted on the Particle last iteration
	utils::Vector<double, 3>& getOldF();

	/// get the indices of the direct neighbors of the Particle
	std::vector<int>& getDirectNeighbors();

	/// set the indices of the direct neighbors of the Particle
	void setDirectNeighbors(std::vector<int>& n);

	/// get the indices of the diagonal neighbors of the Particle
	std::vector<int>& getDiagonalNeighbors();

	/// set the indices of the diagonal neighbors of the Particle
	void setDiagonalNeighbors(std::vector<int>& n);
	
	/// reset the force acting on the Particle
	void resetForce();

	/// adds force acting on the Particle
	inline void addForce(const utils::Vector<double, 3>& force)
	{
		//f = f + force;
		// direct calculation for better performance
		f[0] += force[0];
		f[1] += force[1];
		f[2] += force[2];
	}

	/// substracts force acting on the Particle
	inline void substractForce(const utils::Vector<double, 3>& force)
	{
		//f = f - force;
		// direct calculation for better performance
		f[0] -= force[0];
		f[1] -= force[1];
		f[2] -= force[2];
	}

	/// compares two particles by their type, position, velocity, force and mass
	bool operator == (Particle& other);

	/// assignment operator
	Particle& operator = (const Particle& p)
	{
		f		= p.f;
		old_f	= p.old_f;
		x		= p.x;
		v		= p.v;
		type	= p.type;
		m		= p.m;
		directNeighbors = p.directNeighbors;
		diagonalNeighbors = p.diagonalNeighbors;

		return *this;
	}

	/// print information about this Particle
	std::string toString();
};

std::ostream& operator<<(std::ostream& stream, Particle& p);

#endif 
