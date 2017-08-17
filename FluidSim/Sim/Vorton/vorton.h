#ifndef VORTON_H
#define VORTON_H

#include "Core\Math\vector3D.h"
#include "Core\Utility.h"

#include <math.h>
#include <float.h>

static const float sfFourPi        = 4.0f * PI;
static const float sfOneOverFourPi = 1.0f / sfFourPi;
static const float sfAvoidSingularity = powf(FLT_MIN, 1.0f / 3.0f);
//static const float sfAvoidSingularity_FLT = FLT_MIN;

#define VORTON_ACCUMULATE_VELOCITY_private( vVelocity , vPosQuery , mPosition , mVorticity , mRadius )      \
{                                                                                                           \
    const Vector3D      vNeighborToSelf     = vPosQuery - mPosition ;                                       \
    const float         radius2             = mRadius * mRadius ;                                           \
    const float         dist2               = vNeighborToSelf.norm() + sfAvoidSingularity ;                  \
    const float         oneOverDist         = finvsqrtf( dist2 ) ;                                          \
    const Vector3D      vNeighborToSelfDir  = vNeighborToSelf * oneOverDist ;                               \
    /* If the reciprocal law is used everywhere then when 2 vortices get close, they tend to jettison. */   \
    /* Mitigate this by using a linear law when 2 vortices get close to each other. */                      \
    const float         distLaw             = ( dist2 < radius2 )                                           \
                                                ?   /* Inside vortex core */                                \
                                                ( oneOverDist / radius2 )                                   \
                                                :   /* Outside vortex core */                               \
                                                ( oneOverDist / dist2 ) ;                                   \
    vVelocity +=  sfOneOverFourPi * ( 8.0f * radius2 * mRadius ) * cross(mVorticity, vNeighborToSelf) * distLaw ;   \
}

#define VORTON_ACCUMULATE_VELOCITY( vVelocity , vPosQuery , vorton ) VORTON_ACCUMULATE_VELOCITY_private( vVelocity , vPosQuery , vorton.mPosition , vorton.mVorticity , vorton.mRadius )


class Vorton
{
public:
	Vorton()
		: mPosition(0.0f, 0.0f, 0.0f)
		, mVorticity(0.0f, 0.0f, 0.0f)
		, mRadius(0.0f)
		, mVelocity(0.0f, 0.0f, 0.0f)
	{

	}

	Vorton( const Vector3D & vPos, const Vector3D & vVort, float fRadius = 0.0f)
		: mPosition(vPos)
		, mVorticity(vVort)
		, mRadius(fRadius)
	{

	}

	/*
		Attention:
			copy ctor ÖÐ£¬Ã»ÓÐcopy velocity, WHY??
	*/
	Vorton( const Vorton & rhs)
		: mPosition(rhs.mPosition)
		, mVorticity(rhs.mVorticity)
		, mRadius(rhs.mRadius)
	{

	}

	Vorton( Vorton && rhs)
		: mPosition(rhs.mPosition)
		, mVorticity(rhs.mVorticity)
		, mRadius(rhs.mRadius)
	{

	}

	~Vorton()
	{
	};

	/*
		Calculate the velocity at point vPosQuery induced by this Vorton
		Attention:
			By the Original Formular, The velocity approaches INF when dist between mPosition
			and vPosQuery equals 0. 

			volumeElement's size is 8 * radius * radius * radius

			if dist < radius
			 1																	1
			--- * volumeElement's size * cross(vorticity, relativePos) * ---------------
			4pi															 dist * raduis * raduis

			if dist > radius
			 1																	1
			--- * volumeElement's size * cross(vorticity, relativePos) * ---------------
			4pi															 dist * dist * dist

	*/
	void AccumulateVelocity( Vector3D & vVelocity, const Vector3D & vPosQuery ) const
	{
		const Vector3D		vSelfToPosQuery = vPosQuery - mPosition;
		const float			fRadius2		= mRadius * mRadius;
		const float			fDist2			= vSelfToPosQuery.norm2() + sfAvoidSingularity;
		const float			fOneOverDist	= 1.0f / (vSelfToPosQuery.norm() + sfAvoidSingularity);
		const float			fDistLaw		= (fDist2 < fRadius2) 
											  ? (fOneOverDist / fRadius2) 
											  : (fOneOverDist / fDist2);

		vVelocity += sfOneOverFourPi * (8.0f * mRadius * fRadius2) * cross(mVorticity, vSelfToPosQuery) * fDistLaw;

	}

	/*
		Calculate the vorticity of this vorton to abtain the velocity at query position.
		Attention:
			cross( Vorticity, Dist ) = Velocity
			cross( Dist, Velocity ) = new_Vorticity
			The new_Vorticity is not the same with Vorticity (except that Vorticity is orthogonal with Dist)
			new_Vorticity is in the Vorticity Dist plane.
			
	*/

	void AssignByVelocity( const Vector3D & vQueryPosition, const Vector3D & vVelocity)
	{
		const Vector3D  vPosRelative = vQueryPosition - mPosition;
		const float dist = vPosRelative.norm();
		mVorticity = sfFourPi * dist * cross(vPosRelative, vVelocity) / (8.0f * mRadius * mRadius * mRadius);
	}

	static void UnitTest();

	Vector3D    mPosition;   ///< Position (in world units) of center of vortex particle
	Vector3D    mVorticity;   ///< Vorticity of vortex particle
	float		mRadius;	///< Radius of vortex particle
	Vector3D    mVelocity;   ///< Velocity of this vorton -- used to cache value obtained during advection, to optimize collision response.
};


#endif