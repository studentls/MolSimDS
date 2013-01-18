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
#include "utils/Color.h"
#include "Simulation.h"

// include GLFW library
//windows
#ifdef WINDOWS
// for point sprites the GL extensions are needed
#include <gl/glew.h>
#include <gl/glfw.h>
#else
//Mac/Linux
#include <GL/glfw.h>
#endif

// boost
#include <boost/thread.hpp>
#include <boost/date_time.hpp>

#define COLOR_COUNT 10

/// simple struct for opengl particles
struct VParticle
{
	float x;
	float y;
	float z;
	short type;
};

/// based on OpenGL
/// use singleton pattern to easily communicate
class Viewer
{
private:
	/// mutex for multithread access to array
	boost::mutex			mutex;

	/// BackBuffer Height/Width
	int						width;
	int						height;

	/// OpenGL id of the used point sprite
	GLuint					pointSpriteID;

	/// Mouse cursor position
	utils::Vector<int, 2>	mouseOldPos;
	utils::Vector<int, 2>	mousePos;
	int						mouseOldWheel;
	int						mouseWheel;

	/// viewport translation
	utils::Vector<float, 2> translation;

	/// zoomfactor	
	float					zoomFactor;

	/// array for opengl particles
	VParticle				*particles;
	unsigned int			particle_size;	/// array size
	unsigned int			particle_count;	/// count of valid particles in array

	/// color array			
	utils::Color			colArray[COLOR_COUNT];

	/// struct to store grid information
	struct Grid 
	{
		float x;
		float y;
		float w;
		float h;
		int xcount;		/// how many cells in x direction?
		int ycount;		/// how many cells in y direction?
	};

	/// store Grid information
	Grid					grid;

	/// variable indicates if viewer is running
	bool					isRunning;

	/// Draw function
	void					Draw();

	/// Process function, handle events
	void					Process();

	/// draws the grid for the particles...
	/// @param xcount cells in x direction
	/// @param ycount cells in y direction
	void					DrawInftyGrid(float x, float y,
										  float w, float h,
										  int xcount, int ycount,
										  bool crosses = true,
										  float factor = 0.06125f);

	/// clears everything up
	void					Shutdown();

	/// callback functions
	static void GLFWCALL	reshape(int w, int h)
	{
		// set data
		Viewer::Instance().width  = w;
		Viewer::Instance().height = h;

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

	inline void updateMouseWheel(const int pos)				{mouseOldWheel = mouseWheel; mouseWheel = pos;}
	
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

	/// wrapper
	static void GLFWCALL mouseWheelFun(int pos)
	{
		Viewer::Instance().updateMouseWheel(pos);
	}

	/// helper function to create a nice looking point sprite
	void fillcircleRGBA(unsigned char *data, int w, int h, int centerx, int centery, int radius, float feather, const utils::Color& col);

	/// helper to generate some nice colors
	void generateColors();

	/// for singleton, constructor is private
	Viewer():width(0), height(0), particle_size(0), particle_count(0), particles(NULL), isRunning(false)
	{
		grid.x = grid.y = 0.25f;
		grid.w = grid.h = 0.6f;
		grid.xcount = grid.ycount = 5;

		zoomFactor = 1.0f;
	}

public:
	/// public destructor
	~Viewer()
	{
		// delete GL Textures
		if(particles)glDeleteTextures(1, &pointSpriteID);

		SAFE_DELETE(particles);
	}
	
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

	/// @return size of interior particle array
	int		getParticleArraySize()	{return particle_size;}

	/// @return count of displayed particles
	int		getParticleCount()	{return particle_count;}

	/// set grid values
	void	setGrid(float x, float y, float w, float h, int xcount, int ycount)
	{
		grid.x		= x;
		grid.y		= y;
		grid.w		= w;
		grid.h		= h;
		grid.xcount = xcount;
		grid.ycount = ycount;
	}


	/// set size of array
	/// @param max_count maximum count of particles that can be saved
	void	setParticleArraySize(const int max_count)
	{
		// delete array if necessary
		if(particles)SAFE_DELETE_A(particles);

		particle_size = max_count;
		particle_count = 0;

		particles = new VParticle[max_count];
	}

	/// retrieve pointer to particle array
	/// be careful! this method shall be only called in combination with UnlockParticles!
	/// @param count desired count of particles to store
	/// @param out pointer to the particle array, where particles can be stored
	/// @return max possible count of particles that can be stored
	inline int		LockParticles(VParticle **out, const int count)
	{
		particle_count = min(count, particle_size);
		
		// first lock mutex
		mutex.lock();

		// out
		*out = particles;

		return particle_count;
	}

	/// unlock particle array, for security reasons a pointer is passed, which shall be NULLED afterwards
	/// @param out pointer to the particle array retrieved from LockParticles, where particles have been stored
	inline void	UnlockParticles(VParticle **out)
	{
		// mutex unlock
		mutex.unlock();

		// null pointer
		*out = NULL;
	}
};

#endif