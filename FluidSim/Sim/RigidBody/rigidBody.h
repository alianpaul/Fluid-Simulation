#ifndef RIGID_BODY_H
#define RIGID_BODY_H

#include "Core\Math\vector3D.h"
#include "Core\Math\matrix3x3.h"
#include <vector>

class RigidBody
{
public:

	RigidBody()
		: mPosition(0.0f, 0.0f, 0.0f)
		, mVelocity(0.0f, 0.0f, 0.0f)
		, mOrientation(0.0f, 0.0f, 0.0f)
		, mAngVelocity(0.0f, 0.0f, 0.0f)
		, mInverseMass(0.0f)
		, mInertiaInv(Matrix3x3::identity())
		, mForce(0.0f, 0.0f, 0.0f)
		, mTorque(0.0f, 0.0f, 0.0f)
		, mMomentum(0.0f, 0.0f, 0.0f)
		, mAngMomentum(0.0f, 0.0f, 0.0f)
	{
	}

	RigidBody(const Vector3D &	vPos, const Vector3D &	vVelocity, const float &  fMass)
		: mPosition(vPos)
		, mVelocity(vVelocity)
		, mOrientation(0.0f, 0.0f, 0.0f)
		, mAngVelocity(0.0f, 0.0f, 0.0f)
		, mInverseMass(1.0f / fMass)
		, mInertiaInv(Matrix3x3::identity() * mInverseMass)
		, mForce(0.0f, 0.0f, 0.0f)
		, mTorque(0.0f, 0.0f, 0.0f)
		, mMomentum(vVelocity * fMass)
		, mAngMomentum(0.0f, 0.0f, 0.0f)
	{
	}

	/* Apply a force to a rigid body at a given location
	*/
	void	ApplyForce(const Vector3D & vForce, const Vector3D & vPosition)
	{
		mForce += vForce;
		const Vector3D vPosRelBody = vPosition - mPosition;
		mTorque += cross(vPosRelBody, vForce); // By the Defination of torque
	}


	/* Apply an impulse to a rigid body through its center-of-mass (i.e. without applying a torque
	*/
	void	ApplyImpulse(const Vector3D & vImpulse)
	{
		mMomentum += vImpulse;
		mVelocity = mInverseMass * mMomentum; //Update the linear velocity accordingly
	}

	/* Apply an impulse to a rigid body at a given location
	*/
	void	ApplyImpulse(const Vector3D & vImpulse, const Vector3D & vPosition)
	{
		mMomentum += vImpulse;                         // Apply impulse
		mVelocity = mInverseMass * mMomentum;        // Update linear velocity accordingly
		const Vector3D vPosRelBody = vPosition - mPosition;
		ApplyImpulsiveTorque(cross(vPosRelBody, vImpulse));
	}

	/* Apply an impulsive torque to a rigid body
	*/
	void	ApplyImpulsiveTorque(const Vector3D & vImpulsiveTorque)
	{
		mAngMomentum += vImpulsiveTorque;          // Apply impulsive torque
		mAngVelocity = mInertiaInv * mAngMomentum; // Update angular velocity accordingly
	}

	void	Update(const float & timeStep)
	{
		mMomentum += mForce * timeStep;
		mVelocity = mInverseMass * mMomentum;
		mPosition += mVelocity * timeStep;

		mAngMomentum += mTorque * timeStep;
		mAngVelocity = mInverseMass *  mAngMomentum;
		mOrientation += mAngVelocity * timeStep;

		mForce = mTorque = Vector3D(0.0f, 0.0f, 0.0f);
	}

	static void UpdateSystem(std::vector< RigidBody > & rigidBodies, float timeStep, unsigned uFrame);

	Vector3D    mPosition;   ///< Position (in world units) of center of vortex particle
	Vector3D    mVelocity;   ///< Linear velocity of sphere
	Vector3D    mOrientation;   ///< Orientation of sphere in axis-angle form
	Vector3D    mAngVelocity;   ///< Angular velocity of sphere

protected:
	float		mInverseMass;   ///< Reciprocal of the mass of this body
	Matrix3x3   mInertiaInv;   ///< Inverse of inertial tensor

private:
	Vector3D    mForce;   ///< Total force applied to this body for a single frame.
	Vector3D    mTorque;   ///< Total torque applied to this body for a single frame.
	Vector3D    mMomentum;   ///< Linear momentum of sphere
	Vector3D    mAngMomentum;   ///< Angular momentum of sphere
};

#endif