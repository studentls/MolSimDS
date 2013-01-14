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


void	Viewer::InitAndDisplay()
{
	// now viewer runs
	isRunning = true;
	
	// Init GLFW Lib
	 if(glfwInit() != GL_TRUE)
       Shutdown();

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
    
	// set mouse position
	int x = 0, y = 0;
	glfwGetMousePos(&y, &y);
	updateMousePos(x, y);
	updateMousePos(x, y);
    
    
    //setupgl();
    //
    //pmutex.lock();
    //createparticles();
    //pmutex.unlock();*/
    
  
    
  //  releasegl();
}

bool	Viewer::MessageLoop()
{
	//loop
    while(isRunning)
    {
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
        
		Draw();
        
        glfwSwapBuffers();
        
		// this should be made thread safe
        isRunning = !glfwGetKey(GLFW_KEY_ESC) && glfwGetWindowParam(GLFW_OPENED);
    }

	return true;
}

void Viewer::Draw()
{
	glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

	// translate by translation vector
    glTranslatef(translation[0], translation[1], 0.0f);
    
    DrawInftyGrid(0.25f, 0.25f, 0.6f, 0.5f, 5, 5);
    
    //pmutex.lock();
    //plotparticles();
    //pmutex.unlock();
    //
    //
    
    glColor3f(1.0f, 1.0f, 1.0f);
    glDisable(GL_TEXTURE_2D);
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

void Viewer::DrawInftyGrid(float x, float y, float w, float h, int xcount, int ycount, float factor)
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
    
    float posinfty = 999.9f;
    float neginfty = -999.9f; // set later to viewport
    
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
    
    // draw crosses0
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
    
    glEnd();
}