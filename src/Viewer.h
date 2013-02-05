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

// for HPC errors will occur
#ifndef ICE

// include GLFW library
//windows
#ifdef WINDOWS
#define NOMINMAX
// for point sprites the GL extensions are needed
#include <gl/glew.h>
#include <gl/glfw.h>

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

#else
//Mac/Linux
#include <GL/glfw.h>
#endif


/// some states for the viewer
#define STATE_VELOCITY  0x0
#define STATE_POSITIONS 0x1
#define STATE_COUNT 2

// define appropriate lib to use for thread synchronization
#ifdef USE_BOOST

// boost
#include <boost/thread.hpp>
#include <boost/date_time.hpp>
#endif

#define COLOR_COUNT 10

/// simple struct for opengl particles
struct VParticle
{
	float x;
	float y;
	float z;

	float vx;
	float vy;
	float vz;

	short type;
};

/// holds a color and position
struct GradientEntry
{
	utils::Color c;
	float		pos;
};

/// simple gradient class
class Gradient
{
private:
	utils::TFastArray<GradientEntry> entries;
public:

	utils::Color getColor(const float position)
	{
		assert(entries.size() >= 2);
		
		int i = 1;

		if(position <= 0.0f)return entries.first().c;
		if(position >= 1.0f)return entries.last().c;

		// go through entries
		while(i < entries.size())
		{
			if(position < entries[i].pos)break;
			i++;
		}

		// fix float problems
		i = utils::min((unsigned int)i, entries.size() -1);

		assert(i < entries.size());

		// found i index of greater entry
		float dist = entries[i].pos - entries[i-1].pos;

		// interpolate colors
		float interpos = (position - entries[i-1].pos) / dist;

		return (1.0f - interpos) * entries[i - 1].c + interpos * entries[i].c;
	}

	void	createDefaultGradient()
	{
		// put some nice default colors
		GradientEntry entry;
		entry.c = utils::Color(0xFF53007C);
		entry.pos = 0.0f;
		entries.push_back(entry);
		entry.c = utils::Color(0xFF174477);
		entry.pos = 0.25f;
		entries.push_back(entry);
		entry.c = utils::Color(0xFF4A97D6);
		entry.pos = 0.75f;
		entries.push_back(entry);
		entry.c = utils::Color(0xFFAECDE5);
		entry.pos = 1.0f;
		entries.push_back(entry);	
		
	}
};

/// based on OpenGL
/// use singleton pattern to easily communicate
class Viewer
{
private:

	/// mutex for multithread access to array
#ifdef USE_BOOST
	boost::mutex			mutex;
#else 
	GLFWmutex				mutex;
#endif

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
	utils::Vector<float, 3> minVel;
	utils::Vector<float, 3> maxVel;
	Gradient				gradient; /// some nice colors
	int						state; /// state of the viewer, what will be rendered?

	/// array to form a key up event
	bool					oldKeys[512];
	bool					curKeys[512];

	inline bool				keyUp(const int key)	{return oldKeys[key] && !curKeys[key];}

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

	/// calculate min/max velocity
	void					calcMinMaxValues();

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

		// init gradient
		gradient.createDefaultGradient();

		// per default velocity
		state = STATE_VELOCITY;

		// false keys
		for(int i = 0; i < 512; i++)
		{
			oldKeys[i] = curKeys[i] = false;
		}

		// Init GLFW Lib
		if(glfwInit() != GL_TRUE)
			Shutdown();

		// initialize mutex if needed
#ifndef USE_BOOST
		mutex = glfwCreateMutex();
#endif
	}

public:
	/// public destructor
	~Viewer()
	{
		// delete GL Textures
		if(particles)glDeleteTextures(1, &pointSpriteID);

		SAFE_DELETE(particles);

#ifndef USE_BOOST
		glfwDestroyMutex(mutex);
#endif
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
	inline int		LockParticles(VParticle **out, const int unsigned count)
	{
#ifdef WINDOWS
		particle_count = utils::min(count, particle_size);
#else
		particle_count = std::min(count, particle_size);
#endif
#ifdef USE_BOOST
		// first lock mutex
		mutex.lock();
#else
		glfwLockMutex(mutex);
#endif

		// out
		*out = particles;

		return particle_count;
	}

	/// unlock particle array, for security reasons a pointer is passed, which shall be NULLED afterwards
	/// @param out pointer to the particle array retrieved from LockParticles, where particles have been stored
	inline void	UnlockParticles(VParticle **out)
	{
#ifdef USE_BOOST
		// mutex unlock
		mutex.unlock();
#else
		glfwUnlockMutex(mutex);
#endif

		// calc min max values
		calcMinMaxValues();

		// null pointer
		*out = NULL;
	}
};

#endif

#endif
