#ifndef CAMERA_H
#define CAMERA_H

#include "Core\Math\vector3D.h"
#include <math.h>

class Camera
{
public:
	Camera(int width, int height);
	~Camera();
	
	Camera(const Camera&) = delete;
	Camera(Camera&&) = delete;
	Camera& operator=(const Camera&) = delete;
	Camera& operator=(Camera&&) = delete;

	void SetCamera();
	
	const int &		GetWidth() const { return mWindowWidth; }
	void			SetWidth(int width) { mWindowWidth = width; }

	const int &		GetHeight() const { return mWindowHeight; }
	void			SetHeight(int height){ mWindowHeight = height; }

	void			SetViewPort(int w, int h){ mWindowWidth = w; mWindowHeight = h; }

	const Vector3D &	GetTarget() const { return mTarget; }
	void				SetTarget(const Vector3D& vTarget) { mTarget = vTarget; }

	const Vector3D &	GetEye() const { return mEye; }
	void				SetEye(const Vector3D& vEye){ mEye = vEye; }

	void	GetOrbit(float & azimuth, float & elevation, float & radius)
	{
		Vector3D vRelativeDistance = mEye - mTarget;
		radius = vRelativeDistance.norm();
		azimuth = atan2f(vRelativeDistance.y, vRelativeDistance.x);
		elevation = acosf(vRelativeDistance.z / radius);
	}

	void SetOrbit(float azimuth, float elevation, float radius)
	{
		Vector3D vRelativeDistance = 
			radius * Vector3D(sin(elevation) * cos(azimuth),
			sin(elevation) * sin(azimuth), 
			cos(elevation));

		mEye = mTarget + vRelativeDistance;
	}

private:
	Vector3D mEye;
	Vector3D mTarget;
	Vector3D mUp;
	int mWindowWidth, mWindowHeight; //Used by glViewport()
};
#endif