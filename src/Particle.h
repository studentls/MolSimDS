/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "utils/Vector.h"


/// a class that is used to represent a particle and all its attributes and physical circumstances
class Particle {

private:
	/// the force effective on this particle
	utils::Vector<double, 3> f;

	/// the force wich was effective on this particle
	utils::Vector<double, 3> old_f;

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
	        int type = 0
	);

	/// destructor
	~Particle()  		{}
	
	/// the position of the particle
	utils::Vector<double, 3> x;

	/// the velocity of the particle
	utils::Vector<double, 3> v;

	/// type of the particle, used to index a particle
	short type;

	/// get the current force acting on the Particle
	utils::Vector<double, 3> getF() const {return f;}

	/// get the force that acted on the Particle last iteration
	utils::Vector<double, 3> getOldF() const	{return old_f;}
	
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
	inline Particle& operator = (const Particle& p)
	{
		f		= p.f;
		old_f	= p.old_f;
		x		= p.x;
		v		= p.v;
		type	= p.type;

		return *this;
	}

	/// print information about this Particle
	std::string toString();
};

std::ostream& operator<<(std::ostream& stream, Particle& p);

#endif 
