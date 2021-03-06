/*
 * Vector
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#ifndef VECTOR_
#define VECTOR_

#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <cassert> 


// defines for axises
#define AXIS_X		0
#define AXIS_Y		1
#define AXIS_Z		2
#define AXIS_XY		10
#define AXIS_XZ		20
#define AXIS_YZ		30
#define AXIS_XYZ	100

//fixed issue, namespace was used inapproriately
namespace utils {

template <typename type, int length>
class Vector;



/** Global operators first */

template <typename type, int length>
std::ostream& operator << (std::ostream& , const Vector<type, length>& );

template <typename type, int length>
Vector<type, length> operator*(double scalar, const Vector<type, length>& v) {
	return v * scalar;
}

/** Vector class definition */

template <typename type, int length>
class Vector {
private:
	type content[length];

	friend Vector operator* <type,length>(double scalar, const Vector& v);
	
	friend std::ostream& operator << <>(std::ostream& stream, const Vector<type, length>& v);

public:
	
	Vector() {
		for (int i = 0; i < length; i++) {
			content[i] = 0;
		}
	}

	Vector(const type& arg) {
		for (int i = 0; i < length; i++) {
			content[i] = arg;
		}
	}

	// useful constructors
	Vector(const type& x_1, const type& x_2)
	{
		assert(length >= 2);

		content[0] = x_1;
		content[1] = x_2;

		for(int i = 2; i < length; i++)content[i] = 0;
	}

	Vector(const type& x_1, const type& x_2, const type& x_3)
	{
		assert(length >= 3);

		content[0] = x_1;
		content[1] = x_2;
		content[2] = x_3;

		for(int i = 3; i < length; i++)content[i] = 0;
	}

	// improve performance by reducing senseless copy calls...
	Vector(type *args) {

		// check if valid and correct length!
#ifdef DEBUG
		assert(args);
		//assert(args[length - 1]); // valid array?
#endif
		for (int i = 0; i < length; i++) {
			content[i] = args[i];
		}
	}

	Vector(const Vector& other) {
		for (int i = 0; i < length; i++) {
			content[i] = other[i];
		}
	}

	Vector operator + (const Vector& rhs) const {
		 type result[length];

		 for (int i = 0; i < length; i++) {
			result[i] = this->content[i] + rhs.content[i];
		 }
		 return Vector(result);
	}

	Vector operator - (const Vector& rhs) const {
		type result[length];

		for (int i = 0; i < length; i++) {
			result[i] = this->content[i] - rhs.content[i];
		}
		return Vector(result);
	}

	Vector operator*(const double scalar) const{
		type result[length];

		for (int i = 0; i < length; i++) {
			result[i] = this->content[i] * scalar;
		}
		return Vector(result);
	}

	double operator * (const Vector& rhs) const {
		double sum = 0;
		for (int i = 0; i < length; i++) {
			sum += (this->content[i] * rhs.content[i]);
		}
		return sum;
	}

	double L2Norm() const {
		double square_sum = 0;
		for (int i = 0; i < length; i++) {
			square_sum += (this->content[i] * this->content[i]);
		}
		return sqrt(square_sum);
	}

	/// @return Norm squared, so dot product of v with v
	double L2NormSq() const {
		double square_sum = 0;
		for (int i = 0; i < length; i++) {
			square_sum += (this->content[i] * this->content[i]);
		}
		return square_sum;
	}

	///
	/// @param to distance between caller vector and vector to
	/// @return returns distance of vector to given vector
	///
	inline double distance(const Vector& to)
	{
		Vector r = to - *this;
		return r.L2Norm();
	}

	///
	/// @param to distance squared between caller vector and vector to,
	/// can be useful as sqrt is a very expensive operation
	/// @return returns distance of vector to given vector
	///
	inline double distanceSq(const Vector& to)
	{
		Vector r = to - *this;

		// dot product
		return r * r;
	}

	/// normalize vector, changes components so euclidean norm of the vector is now 1
	/// @return returns normalized vector
	Vector normalize()
	{
		double norm = L2Norm();
		for(int i = 0; i < length; i++)
			content[i] /= norm;
		return *this;
	}


	bool equals(const Vector& rhs) const {
		for (int i = 0; i < length; i++) {
			if (rhs.content[i] != this->content[i]) {
				return false;
			}
		}
		return true;
	}

	Vector& operator=(const Vector& rhs) {
		if(this != &rhs) {
			for (int i = 0; i < length; i++) {
				content[i] = rhs.content[i];
			}
		}
		return *this;
	}

	Vector& operator=(double rhs) {
		for (int i = 0; i < length; i++) {
			content[i] = rhs;
		}
		return *this;
	}

	type& operator[](int i) {

		//this function was totally unsafe...
		assert(i >= 0 && i < length);

		return content[i];
	}

	const type& operator[](int i) const {
			return content[i];
	}

	bool operator==(const Vector& rhs) const {
		for (int i = 0; i < length; i++) {
			if(content[i] != rhs.content[i]) {
				return false;
			}
		}
		return true;
	}

	std::string toString() const {
		std::stringstream str;
		str << *this;
		return str.str();
	}
};

template <typename type, int length>
std::ostream& operator << (std::ostream& stream, const Vector<type, length>& v) {

	stream << "[";
	for (int i = 0; i < length; i++) {
		stream << v.content[i] << ";";
	}
	stream << "]";
	return stream;
}

}


//typedef for easier understanding
typedef utils::Vector<double, 3> Vec3;
typedef utils::Vector<double, 2> Vec2;

namespace utils
{

/// simple struct to encapsulate size/area information
struct BoundingBox
{
	Vector<double, 3> extent; /// extent of the BB
	Vector<double, 3> center; /// center of the BB in absolute coordinates

	BoundingBox()	{}
};

}

#endif /* VECTOR_ */
