//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File main.cpp
// contains main function
//------------------------------------------------------------------------------------------------


#include "Logging.h"
#include "MolSim.h"


#include <gl\GL.h>
#include <gl\glut.h>
#include <gl\GLU.h>

//make it global for test reasons..
MolSim *molsim;


float Colors[] = {	1.0f, 1.0f, 0.0f, 1.0f,
					0.0f, 0.0f, 1.0f, 1.0f,
					1.0f, 0.0f, 0.0f, 0.0f};

void draw()
{
	std::vector<Particle> particles = molsim->getParticles();

	glutPostRedisplay();

	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPointSize(3.0f);
	glColor3f(0.3, 0.3, 0.3);
	//show grid
	for(int x = 0; x < 180/ 3 - 1; x++)
		for(int y = 0; y < 90 / 3 - 1; y++)
		{
			float fx = 3*x / (180.0f / 2.0f) - 1.0f;
			float fy = 3*y / (90.0f / 2.0f) -1.0f;
			float fx2 = 3*(x+1) / (180.0f / 2.0f) - 1.0f;
			float fy2 = 3*(y+1) / (90.0f / 2.0f) -1.0f;
			glBegin(GL_LINES);
			glVertex3f(fx, fy, 0.0f);
			glVertex3f(fx2, fy, 0.0f);
			glVertex3f(fx, fy, 0.0f);
			glVertex3f(fx, fy2, 0.0f);
			glEnd();
		}

	glBegin(GL_POINTS);
	glColor3f(1.0, 1.0, 1.0);

	for(std::vector<Particle>::iterator it = particles.begin(); it != particles.end(); it++)
	{
		//change color according to particle's container...
		float d = 1.0f / ((int)(it->type / 2));
		switch(it->type % 8)
		{
			case 0:
				glColor3f(0.0f * d , 0.0f * d, 0.0f * d);
				break;
			case 1:
				glColor3f(1.0f * d, 0.0f * d, 1.0f * d);
				break;
			case 2:
				glColor3f(1.0f * d, 1.0f * d, 0.0f * d);
				break;
			case 3:
				glColor3f(0.0f * d, 1.0f * d, 1.0f * d);
				break;

			case 4:
				glColor3f(1.0f * d , 0.0f * d, 0.0f * d);
				break;
			case 5:
				glColor3f(0.0f * d, 0.0f * d, 1.0f * d);
				break;
			case 6:
				glColor3f(.0f * d, 1.0f * d, 1.0f * d);
				break;
			case 7:
				glColor3f(0.0f * d, 1.0f * d, 0.4f * d);
				break;
		}
		float x = (*it).x[0];
		float y = (*it).x[1];
		 x = x / (180.0f / 2.0f) - 1.0f;
		 y = y / (90.0f / 2.0f) -1.0f;
		glVertex3f(x, y, 0.0f);
	}

	glEnd();

	glutSwapBuffers();
}

void idle()
{
	molsim->Step();
}

// entry point
int main(int argc, char* argsv[]) {


	configureLoggers();

	/*MolSim **/molsim = new MolSim();

	if(FAILED(molsim->Init(argc, argsv)))
	{
		LOG4CXX_ERROR(generalOutputLogger, ">> Initialization of Molecular Simulator failed");
		LOG4CXX_ERROR(generalOutputLogger, ">> quitting program...");
		return 0;
	}


	//for test reasons
	glutInit(&argc, argsv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize(200*4,90*4);
    glutCreateWindow("MolSimDS");

	glutDisplayFunc(draw);
	glutIdleFunc(idle);
	glutMainLoop();



	molsim->Run();

	molsim->Release();

	SAFE_DELETE(molsim);

	return 0;
}
