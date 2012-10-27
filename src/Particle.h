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

	/// a full constructor that takes the position and velocity (both as 3-tuple vectors)
	/// as well as the mass and type
	Particle(
			/// for visualization, we need always 3 coordinates
			/// -> in case of 2d, we use only the first and the second
			utils::Vector<double, 3> x_arg,
	        utils::Vector<double, 3> v_arg,
	        double m_arg,
	        int type = 0
	);

	/// the position of the particle
	utils::Vector<double, 3> x;

	/// the velocity of the particle
	utils::Vector<double, 3> v;

	/// the mass of this particle
	double m;

	/// type of the particle
	/// The use of this value depends on particularities of the implementation
	/// currently not in use
	int type;

	/// get the current force acting on the Particle
	utils::Vector<double, 3>& getF();

	/// get the force that acted on the Particle last iteration
	utils::Vector<double, 3>& getOldF();

	/// change the force acting on the Particle
	void changeForce(utils::Vector<double, 3> force);

	/// overwrites the equality operator of Particles
	bool operator==(Particle& other);

	/// print information about this Particle
	std::string toString();
};

std::ostream& operator<<(std::ostream& stream, Particle& p);

#endif /* PARTICLE_H_ */
