#include "camera.h"
#include <GL\glew.h>
#include <GL\glut.h>


Camera::Camera(int width, int height)
	:mEye(0.0f, -2.0f, 0.0f)
	, mTarget(0.0f, 0.0f, 0.0f)
	, mUp(0.0f, 0.0f, 1.0f)
{
	mWindowWidth = width;
	mWindowHeight = height;
}

void Camera::SetCamera()
{
	glViewport(0, 0, mWindowWidth, mWindowHeight);
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(mEye.x, mEye.y, mEye.z,
		mTarget.x, mTarget.y, mTarget.z,
		mUp.x, mUp.y, mUp.z);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); /*Set the stack top matrix to Identiy*/
	//glOrtho(-2.0, 2.0,-2.0, 2.0,-2.0, 2.0);
	gluPerspective(75.0f,
		float(mWindowWidth)/float(mWindowHeight),
		0.1f, 1000.0f);

}

Camera::~Camera()
{

}