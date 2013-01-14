//------------------------------------------------------------------------------------------------
// (C) 2013 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File Viewer.h
// contains class to visualize output
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------
#ifndef VIEWER_HEADER_
#define VIEWER_HEADER_

#include "utils/Base.h"
#include "utils/Vector.h"
#include "Simulation.h"

// include GLFW library
//windows
#ifdef WINDOWS
#include <gl/glfw.h>
#else
//Mac/Linux
#include <GL/glfw.h>
#endif

// boost
#include <boost/thread.hpp>
#include <boost/date_time.hpp>

/// based on OpenGL
/// use singleton pattern to easily communicate
class Viewer
{
private:
	/// BackBuffer Height/Width
	int	width;
	int	height;

	/// Mouse cursor position
	utils::Vector<int, 2> mouseOldPos;
	utils::Vector<int, 2> mousePos;

	/// viewport translation
	utils::Vector<float, 2> translation;

	/// variable indicates if viewer is running
	bool			isRunning;

	/// Draw function
	void			Draw();

	/// draws the grid for the particles...
	/// @param xcount cells in x direction
	/// @param ycount cells in y direction
	void			DrawInftyGrid(float x, float y, float w, float h, int xcount, int ycount, float factor = 0.06125f);

	/// clears everything up
	void			Shutdown();

	/// callback functions
	static void GLFWCALL reshape(int w, int h)
	{
		glViewport(0, 0, w, h);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		if (w <= h)
			gluOrtho2D (0.0, 1.0, 0.0, 1.0*(GLfloat)h/(GLfloat)w);
		else
			gluOrtho2D (0.0, 1.0*(GLfloat)w/(GLfloat)h, 0.0, 1.0);
	}

	/// update function
	inline void	updateMousePos(const int x, const int y)	{mouseOldPos = mousePos; mousePos[0] = x; mousePos[1] = y;}

	/// event function for left mouse button pushed
	void	OnMouseLDown();

	/// wrapper
	static void GLFWCALL mousePosFun(int x, int y)
	{
		Viewer::Instance().updateMousePos(x, y);
    
		if(glfwGetMouseButton(GLFW_MOUSE_BUTTON_LEFT))
		{
			Viewer::Instance().OnMouseLDown();
		}  
	}

	/// for singleton, constructor is private
	Viewer():width(0), height(0), isRunning(false)	{}

public:
	/// instance method
	static Viewer& Instance()
	{
		static Viewer theoneandonly;
		return theoneandonly;
	}
	
	/// starts viewer
	void	InitAndDisplay();

	/// run message loop
	/// @return returns true if loop is left without any errors
	bool	MessageLoop();

	/// @return returns if viewer is running
	bool	IsRunning()	{return isRunning;}
};

#endif