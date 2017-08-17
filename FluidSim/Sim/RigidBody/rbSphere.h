#ifndef RB_SPHERE_H
#define RB_SPHERE_H

#include "rigidBody.h"

class RbSphere : public RigidBody
{
public:
	typedef RigidBody Parent;

	RbSphere()
		: RigidBody()
		, mRadius( 0.0f )
	{
	}

	RbSphere(const Vector3D & vPos, const Vector3D & vVelocity, const float & fMass, const float & fRadius)
		: RigidBody( vPos, vVelocity, fMass)
		, mRadius(fRadius)
	{
		mInertiaInv = Matrix3x3::identity() * ( 5.0f * mInverseMass / (2.0f * fRadius * fRadius) );
	}

	static void UpdateSystem(std::vector< RbSphere > & rbSpheres, float timeStep, unsigned uFrame);

	float mRadius;

};

#endif