//------------------------------------------------------------------------------------------------
// (C) 2013 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File Viewer.cpp
// contains class to visualize output
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------

#include "Viewer.h"

// not for HPC
#ifndef ICE

// Little Helper for Reordering
unsigned long ARGBToABGR(unsigned long col)
{
	return ((col & 0x000000ff) << 16) | (col & 0xff00ff00) | ((col & 0x00ff0000) >> 16);
}


void	Viewer::InitAndDisplay()
{
	// now viewer runs
	isRunning = true;

	// open window, and set basic attributes
    glfwOpenWindow(0, 0, 8, 8, 8, 8, 8, 8, GLFW_WINDOW);
    glfwSetWindowTitle("molsim");
    glfwEnable(GLFW_STICKY_KEYS);
    glfwSwapInterval(1);
    glfwSetTime(0.0);
    
    glfwEnable(GLFW_MOUSE_CURSOR);
    
    glfwGetWindowSize(&width, &height);
    glfwSetWindowSizeCallback(reshape);
    glfwSetMousePosCallback(mousePosFun);
	glfwSetMouseWheelCallback(mouseWheelFun);
    
	// set mouse position
	int x = 0, y = 0;
	glfwGetMousePos(&y, &y);
	updateMousePos(x, y);
	updateMousePos(x, y);
    
	// set mouse wheel
	int pos = 0;
	pos = glfwGetMouseWheel();
	updateMouseWheel(pos);
	updateMouseWheel(pos);

    // init point sprite
	glGenTextures(1, &pointSpriteID);
    glBindTexture(GL_TEXTURE_2D, pointSpriteID);
    
    // generate some data
    unsigned char *data = NULL;    
    int w = 128, h = 128;    
    data = new unsigned char[w * h * 4]; // use rgba format
    
    memset(data, 0, w * h * 4);
  
	// fill with a nice feathered circle
	fillcircleRGBA(data, w, h, 32*2, 32*2, 15*2, 2*10.0f, utils::Color(1.0f, 1.0f, 1.0f));
           
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, data); 
	//glGenerateMipmap(GL_TEXTURE_2D);  // generate mipmaps, glu is obsolete!
	
	delete [] data;
        
    glBindTexture(GL_TEXTURE_2D, pointSpriteID);


	// generate some nice colors for visualization
	generateColors();
}

bool	Viewer::MessageLoop()
{
	//loop
    while(isRunning)
    {
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

		// process any data
		Process();

#ifdef USE_BOOST
		// lock
		mutex.lock();
#else
		glfwLockMutex(mutex);
#endif

		Draw();

#ifdef USE_BOOST
		mutex.unlock();
#else
		glfwUnlockMutex(mutex);
#endif
        
        glfwSwapBuffers();
        
		// this should be made thread safe
        isRunning = !glfwGetKey(GLFW_KEY_ESC) && glfwGetWindowParam(GLFW_OPENED);
    }

	return true;
}

void Viewer::Process()
{
	// handle mouse wheel delta
	float factor = 0.15f;
	int deltaw = mouseWheel - mouseOldWheel;
	mouseOldWheel = mouseWheel;

	zoomFactor *= powf(1.0f - factor, (float)deltaw);


	// sync keys
	for(int i = 0; i < 512; i++)oldKeys[i] = curKeys[i];

	// new key states
	for(int i = 0; i < GLFW_KEY_SPECIAL+69; i++)curKeys[i] = (glfwGetKey(i) == GLFW_PRESS);


	// change state if key up is registered
	if(keyUp(GLFW_KEY_SPACE))state = (state + 1) % STATE_COUNT;
}

void Viewer::Draw()
{
	float w = (float)width;
	float h = (float)height;

	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	// scale factor
	float s = zoomFactor;

	if (w <= h)
		gluOrtho2D (0.0 * s, s * 1.0, s * 0.0, s * 1.0*(GLfloat)h/(GLfloat)w);
	else
		gluOrtho2D (0.0 * s, s * 1.0*(GLfloat)w/(GLfloat)h, s *0.0, s*1.0);	
	

	glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();		

    // translate by translation vector
    glTranslatef(translation[0] * s, translation[1] * s, 0.0f);

	
	// fit to width
	// scale a bit
	float aspectratio = (float)width / (float)height;
	float scalefactor = aspectratio / (grid.w);
	glScalef(scalefactor, scalefactor, 1.0f);

	// draw grid	
	bool drawCrosses = zoomFactor <= 0.1 ? true : false; // draw crosses only for small zoomfactor
	int xcount = zoomFactor <= 1.0 ? grid.xcount : 1;	 // draw grid only when visibility is good
	int ycount = zoomFactor <= 1.0 ? grid.ycount : 1;
	DrawInftyGrid(grid.x, grid.y, grid.w, grid.h, xcount, ycount, drawCrosses);
    
 //  	// point sprites for a nice visualization
	//glEnable(GL_TEXTURE_2D);
 //   glBindTexture(GL_TEXTURE_2D, pointSpriteID);

	//glEnable (GL_BLEND);
 //   glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//glEnable(GL_COLOR_MATERIAL);
 //   glEnable(GL_POINT_SPRITE);
 //   glBindTexture(GL_POINT_SPRITE, pointSpriteID);
 //   glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
	//glTexEnvi(GL_POINT_SPRITE, GL_TEXTURE_ENV_MODE, GL_MODULATE );
	//glTexEnvf(GL_POINT_SPRITE, GL_TEXTURE_ENV_MODE, GL_BLEND);

	//// set correct size
	//glPointSize(10.0f);

 //   glBegin(GL_POINTS);
 //   if(particles)
 //   for(int i = 0; i < particle_count; i++)
 //   {
 //       glColor3f(1.0f, 0.0f, 0.0f);
 //       glVertex3f(particles[i].x, particles[i].y, 0.0f);
 //   }    
 //   glEnd();    
 //   glDisable(GL_POINT_SPRITE);

	// use quads to visualize sprites
	glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, pointSpriteID);

	glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	float fSize = 3.0f;
	float fSize2 = fSize * 0.5f; // half size

    glBegin(GL_QUADS);
    if(particles)
    for(int i = 0; i < particle_count; i++)
    {
      /*  glColor3f(1.0f, 0.0f, 0.0f);
        glVertex3f(particles[i].x, particles[i].y, 0.0f);*/

		float x = particles[i].x;
		float y = particles[i].y;
		float vx = particles[i].vx;
		float vy = particles[i].vy;

		utils::Color col;

		// render color due to state
		if(state == STATE_POSITIONS)
		{
			// color according to array
			col = colArray[particles[i].type % COLOR_COUNT];
		}
		else if(state == STATE_VELOCITY)
		{
			// interpolate color
			
			/*float minLengthSq = minVel[0] * minVel[0] + minVel[1] * minVel[1];
			float maxLengthSq = maxVel[0] * maxVel[0] + maxVel[1] * maxVel[1];*/

			float l = sqrtf(vx * vx + vy*vy);

			float minLengthSq = 0.0f;
			float maxLengthSq = 30.0f;

			float t = l / (maxLengthSq - minLengthSq);

			if(t > 1.0f)t = 1.0f;
			if(t < 0.0f)t = 0.0f;

			// interpolate color from gradient
			col = gradient.getColor(1.0f - t);
		}


		glColor3f(col.r, col.g, col.b);


		glTexCoord2d(0.0,0.0); glVertex3f(x - fSize2, y - fSize2, 0.0);
		glTexCoord2d(1.0,0.0); glVertex3f(x + fSize2, y - fSize2, 0.0f);
		glTexCoord2d(1.0,1.0); glVertex3f(x + fSize2, y + fSize2, 0.0f);
		glTexCoord2d(0.0,1.0); glVertex3f(x - fSize2, y + fSize2, 0.0f);

    }    
    glEnd();    

    glDisable(GL_TEXTURE_2D);

    glColor3f(1.0f, 1.0f, 1.0f);
    glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);
    glFlush();
}

void Viewer::Shutdown()
{
	glfwTerminate();
}

void Viewer::OnMouseLDown()
{
    glfwGetWindowSize(&width, &height);
        
    float factor = 1.0 / (float)height;
        
    translation[0] += factor * (mousePos[0] - mouseOldPos[0]);
	translation[1] -= factor * (mousePos[1] - mouseOldPos[1]);
}

void Viewer::DrawInftyGrid(float x, float y, float w, float h, int xcount, int ycount, bool crosses, float factor)
{
   	float dx = w / (float)(xcount);
	float dy = h / (float)(ycount);

    glColor3f(0.8f, 0.8f, 0.8f);
    glBegin(GL_LINES);
    
    for(int i = 0; i < xcount + 1; i++)
    {
        glVertex3f(x + (float)i * dx, y, 0.0f);
        glVertex3f(x + (float)i * dx, y + h, 0.0f);
    }
    
    for(int i = 0; i < ycount + 1; i++)
    {
        glVertex3f(x, y + i * dy, 0.0f);
        glVertex3f(x + w, y + i * dy, 0.0f);
    }
    
    
    //draw infty lines...
    
    float posinfty = 99999.9f;
    float neginfty = -99999.9f; // set later to viewport
    
    glColor3f(0.95f, 0.95f, 0.95f);
    
    // x - direction
    for(int i = 0; i < xcount + 1; i++)
    {
        glVertex3f(x + dx * i, y, 0.0f);
        glVertex3f(x + dx * i, y + neginfty, 0.0f);
        glVertex3f(x + dx * i, y + h, 0.0f);
        glVertex3f(x + dx * i, y + h + posinfty, 0.0f);

    }
    
    // y - direction
    for(int i = 0; i < ycount + 1; i++)
    {
        glVertex3f(x, y + i * dy, 0.0f);
        glVertex3f(x + neginfty, y + i * dy, 0.0f);
        glVertex3f(x + w, y + i * dy, 0.0f);
        glVertex3f(x + w + posinfty, y + i * dy, 0.0f);
    }
    
    
    
    glColor3f(0.0f, 0.0f, 0.0f);
    

    // draw crosses?
	if(crosses)
	{
		for(int i = 0; i < xcount + 1; i++)
			for(int j = 0; j < ycount + 1; j++)
			{
				if(i > 0)glVertex3f(x + i * dx - dx * factor, y + j * dy, 0.0f);
				else glVertex3f(x + i * dx, y + j * dy, 0.0f);
				if(i < xcount)glVertex3f(x + i * dx + dx * factor, y + j * dy, 0.0f);
				else glVertex3f(x + i * dx, y + j * dy, 0.0f);
				if(j > 0)glVertex3f(x + i * dx, y + j * dy - dy *factor, 0.0f);
				else glVertex3f(x + i * dx, y + j * dy, 0.0f);
				if(j < ycount)glVertex3f(x + i * dx, y + j * dy + dy *factor, 0.0f);
				else glVertex3f(x + i * dx, y + j * dy, 0.0f);
			}
	}
    glEnd();
}

void Viewer::fillcircleRGBA(unsigned char *data, int w, int h, int centerx, int centery, int radius, float feather, const utils::Color& col)
{
    
    float innerradiusSq = (radius - feather) * (radius - feather);
    float outerradiusSq = (radius + feather) * (radius + feather);
    
	unsigned char r = (unsigned char)col.getRed();
	unsigned char g = (unsigned char)col.getGreen();
	unsigned char b = (unsigned char)col.getBlue();

    //loop
    for(int x = 0; x < w; x++)
        for(int y = 0; y < h; y++)
        {
            float distSq = (x - centerx) * (x - centerx) + (y - centery) * (y - centery);
            int index = 4 * (x + y * w);
            // inside inner radius?
            if(distSq < innerradiusSq)
            {
                data[index + 0] = r;
                data[index + 1] = g;
                data[index + 2] = b;
                data[index + 3] = 255;
            }
            // inside outer radius(between inner and outer => feather)
            else if(distSq < outerradiusSq)
            {
                float factor = 0.0f;
                
                float distfrominnerradius = sqrtf(distSq) - sqrtf(innerradiusSq);
                
                factor = distfrominnerradius / ( 2.0f * feather);
                
                if(factor > 1.0f)factor = 1.0f;
                else if(factor < 0.0f)factor = 0.0f;
                
                factor = 1.0f - factor;
                
                data[index + 0] = r;
                data[index + 1] = g;
                data[index + 2] = b;
                data[index + 3] = 255 * factor;

            }
            else
            {
                data[index + 0] = 0;
                data[index + 1] = 0;
                data[index + 2] = 0;
                data[index + 3] = 0;
            }
	}
    
}

void Viewer::generateColors()
{
	using namespace utils;

	//manually set some nice colors
	colArray[0] = Color(0.0f, 0.0f, 0.0f);
	colArray[1] = Color(0xff2c87c7);//Color(1.0f, 0.0f, 0.0f);
	colArray[2] = Color(0xff26adde);//Color(0.0f, 1.0f, 0.0f);	
	colArray[3] = Color(0xff2459d4);//Color(0.0f, 0.0f, 1.0f);	
	colArray[4] = Color(0xff21bed4);	
	colArray[5] = Color(0.0f, 0.0f, 0.0f);	
	colArray[6] = Color(0.0f, 0.0f, 0.0f);		
	colArray[7] = Color(0.0f, 0.0f, 0.0f);	
	colArray[8] = Color(0.0f, 0.0f, 0.0f);	
	colArray[9] = Color(0.0f, 0.0f, 0.0f);
}

void Viewer::calcMinMaxValues()
{
	// set to infty
	for(int i = 0; i < 3; i++)minVel[i] = +99999.9f;
	for(int i = 0; i < 3; i++)maxVel[i] = -99999.9f;

	// go through particles
	for(int i = 0; i < this->particle_count; i++)
	{
		minVel[0] = utils::min(minVel[0], particles[i].vx);
		minVel[1] = utils::min(minVel[1], particles[i].vy);
		minVel[2] = utils::min(minVel[2], particles[i].vz);

		maxVel[0] = utils::max(maxVel[0], particles[i].vx);
		maxVel[1] = utils::max(maxVel[1], particles[i].vy);
		maxVel[2] = utils::max(maxVel[2], particles[i].vz);
	}
}

#endif
