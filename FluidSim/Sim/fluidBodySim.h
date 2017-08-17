#ifndef FLUID_BODY_SIM_H
#define FLUID_BODY_SIM_H

#include "RigidBody\rbSphere.h"
#include "Vorton\vortonSim.h"

class FluidBodySim
{
public:

	FluidBodySim( float viscosity, float density)
		:mVortonSim( viscosity, density)
	{}

	void						Initialize( unsigned numTracersPerCellCubeRoot );
	void						Update(float timeStep, unsigned uFrame);
	VortonSim &					GetVortonSim(void)    { return mVortonSim; }
	std::vector< RbSphere > &	GetSpheres(void)      { return mSpheres; }
	void						Clear(void)
	{
		mVortonSim.Clear();
		mSpheres.clear();
	}

private:

	void RemoveEmbeddedParticles();
	void SolveBoundaryConditions();

	VortonSim					mVortonSim;
	std::vector< RbSphere >		mSpheres;

};

#endif