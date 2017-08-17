#include "inteSimVis.h"

#include <GL/glew.h>
#include <GL/glut.h>

#include "Core\Math\vector2D.h"
#include "Core\Math\vector3D.h"
#include "Core\Math\matrix4x4.h"

#include "Core\Utility.h"
#include "Core\perf.h"

#include <iostream>

#include "Space\uniformGrid.h"
#include "Space\nestedGrid.h"
#include "Space\uniformGridMath.h"
#include "Sim\Vorton\vorton.h"
#include "Sim\Vorton\vorticityDistribution.h"

static InteSimVis* sInstance = 0;
static const float timeStep = 1.0f / 30.0f;

static const int WIDTH = 1280;
static const int HEIGHT = 720;
static int sMousePrevX = -999;
static int sMousePrevY = -999;
static int sMouseButtons[3] = { 0, 0, 0 };

static const size_t vortonStride = sizeof(Vorton);
static const size_t tracerStride = sizeof(Particle);

static void DrawCoordinates(void)
{
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);

	glMatrixMode(GL_MODELVIEW);
	glBegin(GL_LINES);
	glColor4f(1.0f, 0.0f, 0.0f, 0.5f);
	glVertex3i(0, 0, 0);
	glVertex3i(1, 0, 0);

	glColor4f(0.0f, 1.0f, 0.0f, 0.5f);
	glVertex3i(0, 0, 0);
	glVertex3i(0, 1, 0);

	glColor4f(0.0f, 0.0f, 1.0f, 0.5f);
	glVertex3i(0, 0, 0);
	glVertex3i(0, 0, 1);
	glEnd();

	glEnable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);
}

static void GlutDisplayCallback(void)
{
	//Set Camera
	//Set Material
	//Draw
	sInstance->mCamera.SetCamera();

	QUERY_PERFORMANCE_ENTER;
	const VortonSim & rVortonSim = sInstance->mFluidBodySim.GetVortonSim();
	sInstance->mTracerRenderer.SetParticleData((char*) & rVortonSim.GetTracers()[0]);
	sInstance->mTracerRenderer.Render(sInstance->mTimeNow, timeStep, rVortonSim.GetTracers().size());

	DrawCoordinates();

	glutSwapBuffers();
	QUERY_PERFORMANCE_EXIT(InteSiVis_Render);
}

static void GlutReshapeCallback(int width, int height)
{

	/*Tells which part of the window the final window appears in
		AKA: FITS THE RENDERED IMAGE INTO THE WINDOW, normalized to view port size
		glViewPort(ldx, ldy, with, height)
		ld: left down corner of the image
		ldx 0, ldy 0: left down corner of the image is the left down 
			corner of the window 
		width, height: we use the width, height of the resized window, so the 
			whole window will show the rendered image
	*/
	sInstance->mCamera.SetViewPort(width, height);
	glViewport(0, 0, width, height);
}

static void GlutIdleCallback()
{
	/*Animation part, update all Particles
	*/
	QUERY_PERFORMANCE_ENTER;
	sInstance->mFluidBodySim.Update(timeStep, sInstance->mFrame);
	QUERY_PERFORMANCE_EXIT(InteSiVis_FluidBodySim_Update);


	++sInstance->mFrame;
	sInstance->mTimeNow += timeStep;

	GlutDisplayCallback();
}

static void GlutSpecialKeyCallback(int key, int x, int y)
{
	Camera& cam = sInstance->mCamera;
	Vector3D vTarget = cam.GetTarget();

	switch (key)
	{
	case GLUT_KEY_F1: sInstance->InitialConditions(1); break;
	case GLUT_KEY_F2: sInstance->InitialConditions(2); break;
	case GLUT_KEY_F3: sInstance->InitialConditions(3); break;
	case GLUT_KEY_F4: sInstance->InitialConditions(4); break;
	case GLUT_KEY_F5: sInstance->InitialConditions(5); break;
	case GLUT_KEY_F6: sInstance->InitialConditions(6); break;
	case GLUT_KEY_F7: sInstance->InitialConditions(0); break;

	case GLUT_KEY_UP	: vTarget.z += 0.1f; cam.SetTarget(vTarget);break;
	case GLUT_KEY_DOWN	: vTarget.z -= 0.1f; cam.SetTarget(vTarget); break;
	case GLUT_KEY_LEFT	: vTarget.x += 0.1f; cam.SetTarget(vTarget); break;
	case GLUT_KEY_RIGHT	: vTarget.x -= 0.1f; cam.SetTarget(vTarget); break;
	default:
		break;
	}

	glutPostRedisplay();
}

void GlutKeyboardHandler(unsigned char key, int x, int y)
{
	Camera & cam = sInstance->mCamera;
	float azimuth, elevation, radius;
	cam.GetOrbit(azimuth, elevation, radius);
	switch (key)
	{
	case ',': radius += 0.1f; break;
	case '.': radius -= 0.1;  break;
	case 27: exit(0); break;
	default:
		return;
	}
	cam.SetOrbit(azimuth, elevation, radius);
}

void GlutMouseHandler(int button, int state, int x, int y)
{
	if (GLUT_DOWN == state)
	{
		sMouseButtons[button] = 1;
	}
	else
	{
		sMouseButtons[button] = 0;
	}

	sMousePrevX = x;
	sMousePrevY = y;
}

void GlutMouseMotionHandler(int x, int y)
{
	if (sMousePrevX == -999)
	{
		sMousePrevX = x;
		sMousePrevY = y;
		return;
	}

	static float fMouseMotionSensitivity = 0.005f;
	const float dx = CLAMP(fMouseMotionSensitivity * float(x - sMousePrevX), -100.0f, 100.0f);
	const float dy = CLAMP(fMouseMotionSensitivity * float(y - sMousePrevY), -100.0f, 100.0f);

	if (sMouseButtons[0])
	{
		/*left click dragging*/
		Camera & cam = sInstance->mCamera;
		float azimuth, elevation, radius;
		// Obtain previous camera parameters.
		cam.GetOrbit(azimuth, elevation, radius);
		// Avoid gimbal lock by limiting elevation angle to avoid the poles.
		elevation = CLAMP(elevation - dy, 4.f * FLT_EPSILON, PI * (1.0f - 4.f * FLT_EPSILON));
		// Set new camera parameters based on how much mouse moved.
		cam.SetOrbit(azimuth - dx, elevation, radius);
	}
	else if (sMouseButtons[2])
	{
		/*right click dragging mouse*/
		Camera & cam = sInstance->mCamera;
		float azimuth, elevation, radius;
		// Obtain previous camera parameters.
		cam.GetOrbit(azimuth, elevation, radius);
		// Set new camera parameters based on how much mouse moved.
		cam.SetOrbit(azimuth, elevation, radius - (dx + dy));

	}
	else if (sMouseButtons[1])
	{
		Camera & cam = sInstance->mCamera;
		const Vector3D & rEye = cam.GetEye();
		const Vector3D & rTarget = cam.GetTarget();
		// Extract world space direction vectors associated with view (used to compute camera-facing coordinates).
		// Note that these vectors are the unit vectors of the inverse of the view matrix.
		// They are the world-space unit vectors of the view transformation.
		GLfloat viewMatrixData[4][4];
		glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat *)& viewMatrixData);

		Vector3D viewRight(viewMatrixData[0][0], viewMatrixData[1][0], viewMatrixData[2][0]);
		Vector3D viewUp(viewMatrixData[0][1], viewMatrixData[1][1], viewMatrixData[2][1]);
		Vector3D viewForward(viewMatrixData[0][2], viewMatrixData[1][2], viewMatrixData[2][2]);
		Vector3D delta = 2.0f * (-viewRight * dx + viewUp * dy);

		// This ought to use the view to change position.
		cam.SetEye(Vector3D(rEye.x + delta.x, rEye.y + delta.y, rEye.z + delta.z));
		cam.SetTarget(Vector3D(rTarget.x + delta.x, rTarget.y + delta.y, rTarget.z + delta.z));
	}

	sMousePrevX = x;
	sMousePrevY = y;
}

void GlutEntryHandler(int state)
{
	if (GLUT_LEFT == state)
	{   // Focus left this window so act like mouse buttons got released.
		sMouseButtons[0] = sMouseButtons[1] = sMouseButtons[2];
	}
}

InteSimVis::InteSimVis(float viscosity, float density)
	: mFluidBodySim(viscosity, density)
	, mCamera( WIDTH, HEIGHT )
	, mVortonRenderer( 0, vortonStride, 0, 0 )
	, mTracerRenderer( 0, tracerStride, 0, 0 )
	, mRenderWindow(0)
	, mStatusWindow(0)
	, mFrame(0)
	, mTimeNow(0.0)
	, mInitialized(false)
{
	sInstance = this;
	InitialConditions(2);
}

InteSimVis::~InteSimVis()
{
	sInstance = nullptr;
}

void	InteSimVis::InitialConditions(unsigned ic)
{
	mScenario = ic;
	sInstance->mFrame	= 0;
	sInstance->mTimeNow = 0;

	mFluidBodySim.Clear();
	static const float      fRadius = 1.0f;
	static const float      fThickness = 1.0f;
	static const float      fMagnitude = 20.0f;
	static const unsigned   numCellsPerDim = 10;
	static const unsigned   numVortonsMax = numCellsPerDim * numCellsPerDim * numCellsPerDim;
	unsigned                numTracersPer = 6;

	std::vector<Vorton> & rVortons = mFluidBodySim.GetVortonSim().GetVortons();

	switch (ic)
	{
	case 0: // vortex ring -- vorticity in [0,1]
		AssignVorticity(rVortons, 2.0f * fMagnitude, numVortonsMax, VortexRing(fRadius, fThickness, Vector3D(1.0f, 0.0f, 0.0f)));
		//mCamera.SetTarget(Vector3D(10.f, 0.0f, 0.0f));
		mCamera.SetEye(Vector3D(0.0f, -5.0f, 0.0f));
	case 1: // Jet ring -- vorticity in [0, 1]
		AssignVorticity(rVortons, fMagnitude, 2048, JetRing(fRadius, fThickness, Vector3D(1.0f, 0.0f, 0.0f)));
		mCamera.SetTarget(Vector3D(10.f, 0.f, 0.f));
		mCamera.SetEye(Vector3D(10.f, -10.0f, 0.0f));
		numTracersPer = 6;
		break;
	case 2: // Projectile
		AssignVorticity(rVortons, 0.125f * FLT_EPSILON, 2048, VortexNoise(Vector3D(4.0f * fThickness, 0.5 * fThickness, 0.5 * fThickness)));
		// Add a sphere
		mFluidBodySim.GetSpheres().push_back(RbSphere(Vector3D(-2.0f, 0.0f, -0.1f), Vector3D(15.0f, 0.0f, 0.0f), 0.2f, 0.2f));
		mCamera.SetTarget(Vector3D(0.0f, 0.0f, 0.0f));
		mCamera.SetEye(Vector3D(0.0f, -2.0f, 0.0f));
		numTracersPer = 3;
		break;
	case 3: // Projectile - side view
		AssignVorticity(rVortons, 0.125f * FLT_EPSILON, 2048, VortexNoise(Vector3D(4.0f * fThickness, 0.5f * fThickness, 0.5f * fThickness)));
		// Add a sphere
		mFluidBodySim.GetSpheres().push_back(RbSphere(Vector3D(-2.0f, 0.0f, -0.1f), Vector3D(15.0f, 0.0f, 0.0f), 0.2f, 0.2f));
		mCamera.SetTarget(Vector3D(1.0f, 0.0f, 0.0f));
		mCamera.SetEye(Vector3D(2.0f, 0.0f, 0.0f));
		numTracersPer = 3;
		break;
	case 4: // Projectile - slow velocity
		AssignVorticity(rVortons, 0.125f * FLT_EPSILON, 2048, VortexNoise(Vector3D(4.0f * fThickness, 0.5 * fThickness, 0.5 * fThickness)));
		// Add a sphere
		mFluidBodySim.GetSpheres().push_back(RbSphere(Vector3D(-2.0f, 0.0f, -0.1f), Vector3D(5.0f, 0.0f, 0.0f), 0.2f, 0.2f));
		mCamera.SetTarget(Vector3D(0.0f, 0.0f, 0.0f));
		mCamera.SetEye(Vector3D(0.0f, -2.0f, 0.0f));
		numTracersPer = 3;
		break;
	case 5:
		AssignVorticity(rVortons, 0.125f * FLT_EPSILON, 2048, VortexNoise(Vector3D(fThickness, fThickness, fThickness)));
		// Add a sphere
		mFluidBodySim.GetSpheres().push_back(RbSphere(Vector3D(0.f, 0.0f, 0.0f), Vector3D(0.f, 0.0f, 0.0f), 100.0f, 0.2f));
		mFluidBodySim.GetSpheres()[0].ApplyImpulsiveTorque(Vector3D(0.0f, 0.0f, 25.0f)); // Make sphere spin
		mCamera.SetTarget(Vector3D(0.0f, 0.0f, 0.0f));
		mCamera.SetEye(Vector3D(-1.5f, 0.0f, 0.0f));
		numTracersPer = 4;
		break;
	case 6:
		AssignVorticity(rVortons, 0.125f * FLT_EPSILON, 2048, VortexNoise(Vector3D(fThickness, fThickness, fThickness)));
		// Add a sphere
		mFluidBodySim.GetSpheres().push_back(RbSphere(Vector3D(0.f, 0.0f, 0.0f), Vector3D(0.f, 0.0f, 0.0f), 100.0f, 0.2f));
		mFluidBodySim.GetSpheres()[0].ApplyImpulsiveTorque(Vector3D(0.0f, 0.0f, 25.0f)); // Make sphere spin
		mCamera.SetTarget(Vector3D(0.0f, 0.0f, 0.0f));
		mCamera.SetEye(Vector3D(0.0f, 0.0f, 1.5f));
		numTracersPer = 4;
		break;
	default:
		break;
	}

	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHT2);
	glDisable(GL_LIGHT3);
	glDisable(GL_LIGHT4);
	glDisable(GL_LIGHT5);
	glDisable(GL_LIGHT6);
	glDisable(GL_LIGHT7);

	mFluidBodySim.Initialize(numTracersPer);
}

void	InteSimVis::InitDevice(int argc, char ** argv)
{
	glutInit(&argc, argv);
	
	unsigned int mode = GLUT_RGBA
		| GLUT_DEPTH //depth buffer
		| GLUT_DOUBLE //double buffered (GL_FRONT/GL_BACK)
		| GLUT_ACCUM //accumulation buffer (add the drawing buffer(GL_FRONT/BACK) to accumulatin buffer)
		| GLUT_STENCIL //stencil buffer
		| GLUT_MULTISAMPLE; //multisampling support
	
	glutInitDisplayMode(mode);

	glutInitWindowPosition(0, 0);
	glutInitWindowSize(WIDTH, HEIGHT);
	
	mRenderWindow = glutCreateWindow("Fluid Simulation");

	glutDisplayFunc(GlutDisplayCallback);
	glutReshapeFunc(GlutReshapeCallback);
	glutIdleFunc(GlutIdleCallback);
	glutSpecialFunc(GlutSpecialKeyCallback);
	glutKeyboardFunc(GlutKeyboardHandler);
	glutMotionFunc(GlutMouseMotionHandler);
	glutMouseFunc(GlutMouseHandler);
	glutEntryFunc(GlutEntryHandler);

	glutSetWindow(mRenderWindow);
	glutMainLoop();
	
}

int main(int argc, char ** argv)
{

	//InteSimVis inte;
	//inte.InitDevice(argc, argv);
	//UniformGrid<unsigned>::UnitTest();
	//GridMathUnitTest();
	//Vorton::UnitTest();
	//VortDistributionTest(1.0f, 2048, VortexRing(1.0f, 0.5f, Vector3D(0.0f, 0.0f, 1.0f)));
	//VortDistributionTest(1.0f, 2048, JetRing(1.5f, 0.3f, Vector3D(0.0f, 0.0f, 1.0f)));
	//VortDistributionTest(1.0f, 2048, VortexNoise( Vector3D(0.5, 0.5, 0.5) ));
	//VortDistributionTest(1.0f, 2048, VortexTube(1.0f, 0.0f, 10.0f, 2, 0));
	//VortDistributionTest(1.0f, 2048, VortexTube(1.0f, 0.0f, 10.0f, 2, -1));

	InteSimVis inte( 0.05f, 1.0f );
	inte.InitDevice( argc, argv );
	return 0;
}

