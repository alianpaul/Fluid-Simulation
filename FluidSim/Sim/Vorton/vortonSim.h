#ifndef VORTON_SIM_H
#define VORTON_SIM_H

#include <vector>

#include "particle.h"
#include "vorton.h"
#include "Space\nestedGrid.h"


#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/tick_count.h"

/*
	Dynamic simulation of a fluid, using tiny vortex elements
	This implements a PORTION of a fluid simulation, and 
	effectively NEGLECTS boundary conditions, boundary conditions 
	is in another module
*/

class VortonSim
{
public:
	VortonSim( float viscosity = 0.0f , float density = 1.0f )
		: mMinCorner(FLT_MAX, FLT_MAX, FLT_MAX)
		, mMaxCorner(-mMinCorner)
		, mViscosity(viscosity)
		, mCirculationInitial(0.0f, 0.0f, 0.0f)
		, mLinearImpulseInitial(0.0f, 0.0f, 0.0f)
		, mAverageVorticity(0.0f, 0.0f, 0.0f)
		, mFluidDensity(density)
		, mMassPerParticle(0.0f)
	{
	}

	void								Initialize( unsigned numTracersPerCellCubeRoot );
	std::vector< Vorton >  &			GetVortons( )                  { return mVortons; }
	const std::vector< Vorton >  &		GetVortons( ) const            { return mVortons; }
	std::vector< Particle > &			GetTracers( )                  { return mTracers; }
	const std::vector< Particle > &		GetTracers( ) const            { return mTracers; }
	/*
		If the order doesn't matter, copy the tail to the objective delete position.
		and then pop back.
	*/
	void								KillTracer(size_t iTracer)
	{
		mTracers[iTracer] = mTracers[mTracers.size() - 1];
		mTracers.pop_back();
	}

	const Vector3D						GetTracerCenterOfMass( ) const;
	const UniformGrid< Vector3D > &		GetVelocityGrid( ) const { return mVelGrid; };
	const float &						GetMassPerParticle() const { return mMassPerParticle; }
	void								Update( float timeStep, unsigned uFrame );
	void								Clear()
	{
		mVortons.clear();
		mInfluenceTree.Clear();
		mVelGrid.Clear();
		mTracers.clear();
	}


private:
	void		AssignVortonsFromVorticity(	const UniformGrid< Vector3D > & vortGrid );
	void		ConservedQuantities( Vector3D & vCirculation, Vector3D & vLinearImpulse );
	void		FindBoundingBox();
	void		MakeBaseVortonGrid();
	void		AggregateClusters( unsigned uParentLayer );
	void		CreateInfluenceTree();
	Vector3D	ComputeVelocity( const Vector3D & vPosition, const unsigned idxParent[3], unsigned iLayer );
	Vector3D	ComputeVelocityBruteForce( const Vector3D & vPosition );
	void		ComputeVelocityGridSlice( size_t izStart, size_t izEnd);
	void		ComputeVelocityGrid();
	void		StretchAndTiltVortons( const float & timeStep, const unsigned & uFrame );
	void		ComputeAverageVorticity(void);
	void		DiffuseVorticityGlobally( const float & timeStep, const unsigned & uFrame );
	void		DiffuseVorticityPSE( const float & timeStep, const unsigned & uFrame );
	void		AdvectVortons( const float & timeStep );

	void		InitializePassiveTracers( unsigned multiplier );
	void		AdvectTracersSlice(const float & timeStep, const unsigned & uFrame, size_t izStart, size_t izEnd);
	void		AdvectTracers(const float & timeStep, const unsigned & uFrame);


	std::vector< Vorton >		mVortons				;	//Dynamic array of tiny vortex elements
	NestedGrid< Vorton >		mInfluenceTree			;	//Influence tree
	UniformGrid< Vector3D >		mVelGrid				;	//Uniform grid of velocity values
	Vector3D					mMinCorner				;	//Minimal corner of axis-aligned bounding box
	Vector3D					mMaxCorner				;	//Maximal corner of axis-aligned bounding box
	float						mViscosity				;	//Viscosity. Used to compute viscous diffusion
	Vector3D					mCirculationInitial		;	//Initial circulation, which should be conserved when viscosity is zero
	Vector3D					mLinearImpulseInitial	;	//Initial linear impulse
	Vector3D					mAverageVorticity		;	//Hack, average vorticity used to compute a kind of viscous vortex diffusion
	float						mFluidDensity			;	//Uniform density of fluid
	float						mMassPerParticle		;	//Mass of each fluid particle
	std::vector< Particle >		mTracers				;	//Passive tracer particles

	tbb::task_scheduler_init	tbb_init;

	friend class VortonSim_AdvectTracers_TBB;

};

#endif