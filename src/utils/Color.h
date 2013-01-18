//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File Color.h
// contains class to handle colors
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------


#ifndef COLOR_HEADER_
#define COLOR_HEADER_

#include <cmath>
#include <iostream>
#include <sstream>
#include <cassert>

namespace utils
{

static const double g_ColConv = 0.003921568627450980392156862745098;

/// Basic Color class
class Color
{
public:
	float a;
	float r;
	float g;
	float b;

	Color():r(0.0),g(0.0),b(0.0),a(1.0) {}
	Color(const float _r, const float _g, const float _b):r(_r), g(_g), b(_b), a(1.0)	{} 
	Color(const float _r, const float _g, const float _b, const float _a):r(_r), g(_g), b(_b), a(_a)	{} 
	Color(const int _r, const int _g, const int _b)
	{
		r = (float)_r / 255.0;
		g = (float)_g / 255.0;
		b = (float)_b / 255.0;
		a = 1.0;
	}
	Color(const Color& c):r(c.r), g(c.g), b(c.b), a(c.a)	{}

	/// value constructor out of unsigned long value
	Color(const unsigned long& col): r(g_ColConv * (float)(unsigned char)((col & 0x00ff0000) >> 16)),
								  g(g_ColConv * (float)(unsigned char)((col & 0x0000ff00) >> 8)),
								  b(g_ColConv * (float)(unsigned char)((col & 0x000000ff))),
								  a(g_ColConv * (float)(unsigned char)((col & 0xff000000) >> 24))	{}

	/// convert to standard 32 Bit representation in ARGB format
	operator unsigned long ()
	{
		
		return((a >= 1.0 ? 255 : a <= 0.0 ? 0 : (unsigned long)(a * 255.0)) << 24 |
			   (r >= 1.0 ? 255 : r <= 0.0 ? 0 : (unsigned long)(r * 255.0)) << 16 |
			   (g >= 1.0 ? 255 : g <= 0.0 ? 0 : (unsigned long)(g * 255.0)) << 8 |
			   (b >= 1.0 ? 255 : b <= 0.0 ? 0 : (unsigned long)(b * 255.0)));
				
	}

	int getRed()const 	{return (r >= 1.0 ? 255 : r <= 0.0 ? 0 : (unsigned long)(r * 255.0)); } 
	int getGreen()const	{return (g >= 1.0 ? 255 : g <= 0.0 ? 0 : (unsigned long)(g * 255.0)); }
	int getBlue()const	{return (b >= 1.0 ? 255 : b <= 0.0 ? 0 : (unsigned long)(b * 255.0)); }
	
	
	void operator =		(const Color& c) {r =  c.r; g =  c.g; b =  c.b; a = c.a;}

	inline std::string toString()
	{
		std::stringstream s;
		s<<"("<<r<<","<<g<<","<<b<<","<<a<<")";
		return s.str();
	}

	//to def!
	static const Color white;
	static const Color black;
	static const Color yellow; 
	static const Color green;
	static const Color blue;
};


//operators

/// multiply (alpha value will not be multiplied)
inline Color operator * (const Color& c, const float d)			{return Color(c.r * d, c.g * d, c.b * d); } 
inline Color operator * (const float d, const Color& c)			{return Color(c.r * d, c.g * d, c.b * d); } 
inline Color operator * (const Color& c, const Color& k)		{return Color(c.r * k.r, c.g * k.g, c.b * k.b);}
inline Color operator / (const Color& c, const float d)			{float invd = 1.0f / d; return Color(c.r * invd, c.g * invd, c.b * invd); } 
inline Color operator + (const Color& a, const Color& b)		{return Color(a.r + b.r, a.g + b.g, a.b + b.b);}
inline Color operator - (const Color& a, const Color& b)		{return Color(a.r - b.r, a.g - b.g, a.b - b.b);}
inline Color operator - (const Color& c)						{return Color(-c.r, -c.g, -c.b);}

// Little Helper for Reordering
inline unsigned long ARGBToABGR(unsigned long col)
{
	return ((col & 0x000000ff) << 16) | (col & 0xff00ff00) | ((col & 0x00ff0000) >> 16);
}

}
#endif
