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
	}

	Vector(type arg) {
		for (int i = 0; i < length; i++) {
			content[i] = arg;
		}
	}

	Vector(type args[length]) {
		for (int i = 0; i < length; i++) {
			content[i] = args[i];
		}
	}

	Vector(const Vector& other) {
		for (int i = 0; i < length; i++) {
			content[i] = other[i];
		}
	}

	Vector operator+(const Vector& rhs) const {
		type result[length];

		for (int i = 0; i < length; i++) {
			result[i] = this->content[i] + rhs.content[i];
		}
		return Vector(result);
	}

	Vector operator-(const Vector& rhs) const {
		type result[length];

		for (int i = 0; i < length; i++) {
			result[i] = this->content[i] - rhs.content[i];
		}
		return Vector(result);
	}

	Vector operator*(double scalar) const{
		type result[length];

		for (int i = 0; i < length; i++) {
			result[i] = this->content[i] * scalar;
		}
		return Vector(result);
	}

	double L2Norm() const {
		double square_sum = 0;
		for (int i = 0; i < length; i++) {
			square_sum += (this->content[i] * this->content[i]);
		}
		return sqrt(square_sum);
	}

	///
	/// @param to distance between caller vector and vector to
	/// @return returns distance of vector to given vector
	///
	double distance(const Vector& to)
	{
		Vector r = to - *this;
		return r.L2Norm();
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

#endif /* VECTOR_ */
