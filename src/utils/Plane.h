//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File Plane.h
// includes a class to handle a plane in 3D-space
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------



#ifndef PLANE_HEADER_
#define PLANE_HEADER_

#include "utils.h"

namespace utils
{
	class Plane
	{
	public:
		// encoded in Hessian Normalform, all x � R� are in the plane, iff
		// n*x - d = 0
		// and for simplicity reasons ||n||_2 = 1 always
		Vector<double, 3>	n; // normal vector, always normalized!
		double				d; // distance

		// constructors
		Plane():n(0, 0, 0), d(0)												{}
		Plane(const Vector<double, 3>& _n, const double _d):n(_n), d(_d)		{n.normalize();}
		Plane(const Plane& p):n(p.n), d(p.d)									{}

		/// construct Plane from Point and normal vector
		/// sets internal values to constructed plane
		Plane constructFromPoint(const Vector<double, 3> point, const Vector<double, 3>& normal)
		{
			n = normal;
			n.normalize();
			d = -(n*point);

			return *this;
		}

		Plane& operator =(const Plane& p)										{n = p.n; d = p.d; return *this;}

		Vector<double, 3> normalize()	{return n.normalize();}

		double				absdistance(const Vector<double, 3>& point)	{return abs(n * point + d);/*vector dot product, note that n is normalized*/}
		double				distance(const Vector<double, 3>& point)	{return n * point + d;/*vector dot product, note that n is normalized*/}
		
	
	};
}

#endif