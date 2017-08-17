#ifndef INTE_SIM_VIS_H
#define INTE_SIM_VIS_H

#include "Render\camera.h"
#include "Render\particleRenderer.h"
#include "Sim\fluidBodySim.h"

class InteSimVis
{
public:
	InteSimVis(float viscosity, float density);
	~InteSimVis();
	
	void	InitDevice(int argc, char ** argv);
	void	InitialConditions(unsigned ic);

	FluidBodySim		mFluidBodySim;		///< Simulation of fluid and rigid bodies
	Camera				mCamera;			///< Camera for rendering
	ParticleRenderer	mTracerRenderer;	///< Renderer for tracers
	ParticleRenderer	mVortonRenderer;	///< Renderer for vortons
	int					mRenderWindow;		///< Identifier for render window
	int                 mStatusWindow;		///< Identifier for status window
	unsigned            mFrame;				///< Frame counter
	double              mTimeNow;			///< Current virtual time
	bool                mInitialized;		///< Whether this application has been initialized
	int                 mScenario;			///< Which scenario is being simulated now
	
	
};

#endif