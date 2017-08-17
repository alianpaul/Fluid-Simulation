#include "fluidBodySim.h"
#include "Core\perf.h"

void FluidBodySim::Initialize( unsigned numTracersPerCellCubeRoot )
{
	RemoveEmbeddedParticles(); //Remove The Invalid vortons
	mVortonSim.Initialize(numTracersPerCellCubeRoot);
	RemoveEmbeddedParticles(); //Remove The Invalid vortons
}

/*
	Only be called in Initialize to remove the vortons and particles inside
	RigidBody
*/

void FluidBodySim::RemoveEmbeddedParticles()
{
	const unsigned numBodies = mSpheres.size();
	for (unsigned uBody = 0; uBody < numBodies; ++uBody)
	{
		RbSphere & rSphere = mSpheres[ uBody ];
		
		//VortonSim not provide function to erase vorton
		//We need to get the vector vorton container, and use container erase function to
		//delete in valid vorton
		std::vector< Vorton > &		rVortons = mVortonSim.GetVortons();
		for (std::vector< Vorton >::iterator iVorton = rVortons.begin(); iVorton != rVortons.end(); )
		{
			Vorton &		rVorton			= *iVorton;
			const Vector3D	vSphereToVorton = rVorton.mPosition - rSphere.mPosition;
			const float		fSphereToVorton = vSphereToVorton.norm();
			
			if (fSphereToVorton < (rVorton.mRadius + rSphere.mRadius))
			{
				iVorton = rVortons.erase( iVorton );
			}
			else
			{
				++ iVorton;
			}

		}

		//VortonSim has a mem func KillTracer to delete particle at some offset
		for (unsigned uTracer = 0; uTracer < mVortonSim.GetTracers().size();)
		{
			Particle &		rTracer			= mVortonSim.GetTracers()[uTracer];
			const Vector3D	vSphereToTracer = rTracer.mPosition - rSphere.mPosition;
			const float		fSphereToTracer = vSphereToTracer.norm();
			if (fSphereToTracer < (rTracer.mSize + rSphere.mRadius))
			{
				mVortonSim.KillTracer( uTracer );
			}
			else
			{
				++ uTracer;
			}
		}
	}
}

void FluidBodySim::SolveBoundaryConditions()
{
	const unsigned numBodies	= mSpheres.size();
	const unsigned numVortons	= mVortonSim.GetVortons().size();
	const unsigned numTracers	= mVortonSim.GetTracers().size();

	const float & rMassPerParticle = mVortonSim.GetMassPerParticle();

	//For each Body
	for (unsigned uBody = 0; uBody < numBodies; ++uBody)
	{
		
		RbSphere & rSphere	= mSpheres[uBody];

		//For Each Vortons that near the body
		for (unsigned uVorton = 0; uVorton < numVortons; ++uVorton)
		{
			Vorton & rVorton	= mVortonSim.GetVortons()[uVorton];
			const Vector3D	vSphereToVorton		= rVorton.mPosition - rSphere.mPosition;
			const float		fSphereToVorton		= vSphereToVorton.norm();
			// This boundary thickness compensates for low discretization resolution,
			// by spreading the influence of the body surface to just outside the body,
			// deeper into the fluid.  This also has an effect somewhat like
			// instantaneous viscous diffusion, in the immediate vicinity of
			// the boundary.  It should be kept as small as possible,
			// but must be at least 1.  A value of 1 means only vortons
			// colliding with the body receive influence.  A value of 2 seems
			// most appropriate since that is the size of a grid cell, so
			// 2 essentially means vortons within a grid cell receive influence.
			// So a value in [1,2] seems appropriate. But values over 1.2 trap
			// vortons inside the body, because the "bend" can draw vortons back
			// toward the body.
			// Note, the larger fBndThkFactor is, the more vortons get influenced,
			// which drives the simulation to instability and also costs more CPU
			// time due to the increased number of vortons involved.
			const float fBndThkFactor		= 1.2f; // Thickness of boundary, in vorton radii.
			const float fBoundaryThickness	= fBndThkFactor * rVorton.mRadius; // Thickness of boundary, i.e. region within which body sheds vorticity into fluid.

			if ( fSphereToVorton < (rSphere.mRadius + fBoundaryThickness) )
			{
				//Vortons that interact with the rigid body.

				//STEP 1 Find the desired velocity that should be 
				//induced by this vorton
				//fv: Ambient Flow velocity, rv: Rigid body velocity, vv: velocity induced by the old voritcity 
				// fv - vv - rv: velocity respect to rigid body without the influence of this vorton
				// The new velocity induced by this vorton should be minus
				const Vector3D vSphereToVortonDir	= vSphereToVorton / fSphereToVorton;
				const Vector3D vSphereToContactPt	= vSphereToVortonDir * rSphere.mRadius;
				const Vector3D vContactPt			= rSphere.mPosition + vSphereToContactPt;
				//fv
				Vector3D vVelAmbientAtContactPt;
				mVortonSim.GetVelocityGrid().Interpolate( vVelAmbientAtContactPt, vContactPt );
				//rv
				const Vector3D vVelDueToRotAtConPt	= cross( rSphere.mAngVelocity, vSphereToContactPt );
				const Vector3D vVelBodyAtContactPt	= rSphere.mVelocity + vVelDueToRotAtConPt;
				//vv
				Vector3D vVelDueToVort;
				rVorton.AccumulateVelocity(vVelDueToVort, vContactPt);

				const Vector3D vVelFlowRelBodyAtConPt(vVelAmbientAtContactPt - vVelDueToVort - vVelBodyAtContactPt);

				const Vector3D vVorticityOld = rVorton.mVorticity;

				//STEP 2 Find where to place the vorton
				//Dir
				//Normal is vSphereToVortonDir
				const Vector3D vVelDir	= vVelAmbientAtContactPt.unit();
				const Vector3D vVortDir = cross(vVelDir, vSphereToVortonDir);
					  Vector3D vBendDir	= cross(vVortDir, vVelDir);
				vBendDir.normalize();
				//Length
				//If the vorton is collide with the rigid body, move it out(tangent to the body)
				//If Not, The distance to the surface of the body don't change.
				const float fConPtToVort	= fSphereToVorton - rSphere.mRadius;
				const float fBendDist		= (fConPtToVort < rVorton.mRadius) ? rVorton.mRadius : fConPtToVort;
				const Vector3D vBend		= vBendDir * fBendDist;
				rVorton.mPosition			= vContactPt + vBend;

				//STEP 3 Update the vorticity of this vorton
				rVorton.AssignByVelocity( vContactPt, -vVelFlowRelBodyAtConPt );
#define DELAY_SHEDDING 1
#if DELAY_SHEDDING 
				const float fGain = 0.1f;
				const float fOneMinusGain = 1.0f - fGain;
				rVorton.mVorticity = fGain * rVorton.mVorticity + fOneMinusGain * vVorticityOld;
#endif
				//Just like particle
				//Update the Vorton's velocity and Body's velocity according to the momentum conservation
				const Vector3D vAngVelDiff		= rVorton.mVorticity - vVorticityOld; //???
				const float	   fMassOfVorton    = 0.3f * rMassPerParticle;
				rSphere.ApplyImpulsiveTorque(vAngVelDiff * fMassOfVorton);

				const Vector3D vVelChange = rVorton.mVelocity - vVelBodyAtContactPt;
				rSphere.ApplyImpulse(vVelChange * rMassPerParticle);
				rVorton.mVelocity = vVelBodyAtContactPt;
			}

		}


		//For Each Tracers that collide with rigie body
		for (unsigned uTracer = 0; uTracer < numTracers; ++uTracer)
		{

			Particle &		rTracer			= mVortonSim.GetTracers()[ uTracer ];
			const Vector3D	vSphereToTracer = rTracer.mPosition - rSphere.mPosition;
			const float		fSphereToTracer = vSphereToTracer.norm();
			if (fSphereToTracer < (rTracer.mSize + rSphere.mRadius))
			{
				//The Tracer is collide with the rigid body
				//Move the tracer to the fSphereToTracer = rTracer.mSize + rSphere.mRadius
				const float		distRescale			= (rTracer.mSize + rSphere.mRadius) / fSphereToTracer * (1.0f + FLT_EPSILON); 
				const Vector3D	vSphereToTracerNew	= vSphereToTracer * distRescale;
				rTracer.mPosition					= rSphere.mPosition + vSphereToTracerNew;

				//Update the tracer's velocity, Now the tracer is at the boundary of the 
				//rigid body, The tracer's velocity is 0 respect to rigid body. So we 
				//set tracer's velocity = rigid body velocity at that point
				const Vector3D vVelDueToRotation	= cross(rSphere.mAngVelocity, vSphereToTracerNew);
				const Vector3D vVelNew				= rSphere.mVelocity + vVelDueToRotation;
				
				//The Rigid body change the rTracer's velocity
				//Due to the momentum conservation, a reaction momentum
				//will act on the rigid body.
				const Vector3D vVelChange = rTracer.mVelocity - vVelNew;
				rSphere.ApplyImpulse(vVelChange * rMassPerParticle);

				rTracer.mVelocity = vVelNew;

			}
		}

	}

}

/*
	TODO: Interact with the bodies.
*/
void FluidBodySim::Update( float timeStep, unsigned uFrame )
{
	QUERY_PERFORMANCE_ENTER;
	mVortonSim.Update(timeStep, uFrame);
	QUERY_PERFORMANCE_EXIT(FluidBodySim_VortonSim_Update);

	QUERY_PERFORMANCE_ENTER;
	SolveBoundaryConditions();
	QUERY_PERFORMANCE_EXIT(FluidBodySim_SolveBoundaryConditions);

	QUERY_PERFORMANCE_ENTER;
	RbSphere::UpdateSystem(mSpheres, timeStep, uFrame);
	QUERY_PERFORMANCE_EXIT(FluidBodySim_RbSphere_UpdateSystem);
}
