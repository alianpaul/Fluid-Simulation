#ifndef PARTICLE_H
#define PARTICLE_H

#include "Core\Math\vector3D.h"

class Particle
{
public:

	Particle()
		: mPosition(0.0f, 0.0f, 0.0f)
		, mVelocity(0.0f, 0.0f, 0.0f)
		, mOrientation(0.0f, 0.0f, 0.0f)
		, mAngularVelocity(0.0f, 0.0f, 0.0f)
		, mMass(0.0f)
		, mSize(0.0f)
		, mBirthTime(0)
	{	
	}

	Vector3D	mPosition			;	//Position of center of particle
	Vector3D	mVelocity			;   //Velocity of the particle
	Vector3D	mOrientation		;   //Orientation of the particle, in axis angle???
	Vector3D	mAngularVelocity	;   //Angular velocity of the particle ???
	float		mMass				;   //Mass
	float		mSize				;	//Size
	int			mBirthTime			;	//Birth time of the particle, in "ticks"


};

#endif