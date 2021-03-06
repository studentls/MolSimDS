/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "Particle.h"

#include <sstream>
#include <iostream>

Particle::Particle(const int type_arg) {
	type = type_arg;
	f = 0.0;
	old_f = 0.0;
	v = 0.0;
	id = -1;
	directNeighbors = std::vector<int>();
	diagonalNeighbors = std::vector<int>();
}

Particle::Particle(const Particle& other) {
	// copy every attribute of the Particle to create a new one
	x = other.x;
	v = other.v;
	f = other.f;
	old_f = other.old_f;
	type = other.type;
	id = other.id;
	directNeighbors = other.directNeighbors;
	diagonalNeighbors = other.diagonalNeighbors;
}

// Todo: maybe use initializater list instead of copy?
Particle::Particle(const	utils::Vector<double, 3>& x_arg,
	        const utils::Vector<double, 3>& v_arg,
	        const int type_arg
) {
    x = x_arg;
    v = v_arg;
    type = type_arg;
    f = 0.0;
    old_f = 0.0;
	id = -1;
	directNeighbors = std::vector<int>();
	diagonalNeighbors = std::vector<int>();
}

std::vector<int>& Particle::getDirectNeighbors() {
	return directNeighbors;
}

void Particle::setDirectNeighbors(std::vector<int>& n) {
	directNeighbors = n;
}

std::vector<int>& Particle::getDiagonalNeighbors() {
	return diagonalNeighbors;
}

void Particle::setDiagonalNeighbors(std::vector<int>& n) {
	diagonalNeighbors = n;
}

void Particle::resetForce() {
	// set the old_force to the current one
	old_f = f;
	// and replace the current force with 0
	f = 0.0;
}

std::string Particle::toString() {
	std::stringstream stream;
	stream << "Particle: X:" << x <<  " v: " << v << " f: " << f << " old_f: " << old_f << " type: " << type;
	return stream.str();
}

bool Particle::operator == (Particle& other) {
	// two Particles are considered identical iff
	// all of their values are equal
	// excepting 'neighbors', since that value is irrelevant for most purposes
	if ( (x == other.x) && (v == other.v) && (f == other.f) &&
			(type == other.type) && (old_f == other.old_f)) {
		return true;
	}

	return false;
}

std::ostream& operator<<(std::ostream& stream, Particle& p) {
	stream << p.toString();
	return stream;
}
