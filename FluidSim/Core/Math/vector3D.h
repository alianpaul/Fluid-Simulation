#ifndef VEC3D_H
#define VEC3D_H

#include <iosfwd>
#include <cmath>
#include "Core\Utility.h"

class Vector3D {
public:

	// components
	float x, y, z;

	/**
	* Constructor.
	* Initializes tp vector (0,0,0).
	*/
	Vector3D() : x(0.0), y(0.0), z(0.0) { }

	/**
	* Constructor.
	* Initializes to vector (x,y,z).
	*/
	Vector3D(float x, float y, float z) : x(x), y(y), z(z) { }

	/**
	* Constructor.
	* Initializes to vector (c,c,c)
	*/
	Vector3D(float c) : x(c), y(c), z(c) { }

	/**
	* Constructor.
	* Initializes from existing vector
	*/
	Vector3D(const Vector3D& v) : x(v.x), y(v.y), z(v.z) { }

	// returns reference to the specified component (0-based indexing: x, y, z)
	inline float& operator[] (const int& index) {
		return (&x)[index];
	}

	// returns const reference to the specified component (0-based indexing: x, y, z)
	inline const float& operator[] (const int& index) const {
		return (&x)[index];
	}

	inline bool operator==(const Vector3D& v) const {
		return v.x == x && v.y == y && v.z == z;
	}

	// negation
	inline Vector3D operator-(void) const {
		return Vector3D(-x, -y, -z);
	}

	// addition
	inline Vector3D operator+(const Vector3D& v) const {
		return Vector3D(x + v.x, y + v.y, z + v.z);
	}

	// subtraction
	inline Vector3D operator-(const Vector3D& v) const {
		return Vector3D(x - v.x, y - v.y, z - v.z);
	}

	// right scalar multiplication
	inline Vector3D operator*(const float& c) const {
		return Vector3D(x * c, y * c, z * c);
	}

	// scalar division
	inline Vector3D operator/(const float& c) const {
		const float rc = 1.0 / c;
		return Vector3D(rc * x, rc * y, rc * z);
	}

	// addition / assignment
	inline void operator+=(const Vector3D& v) {
		x += v.x; y += v.y; z += v.z;
	}

	// subtraction / assignment
	inline void operator-=(const Vector3D& v) {
		x -= v.x; y -= v.y; z -= v.z;
	}

	// scalar multiplication / assignment
	inline void operator*=(const float& c) {
		x *= c; y *= c; z *= c;
	}

	// scalar division / assignment
	inline void operator/=(const float& c) {
		(*this) *= (1. / c);
	}

	/**
	* Returns Euclidean length.
	*/
	inline float norm(void) const {
		return sqrt(x*x + y*y + z*z);
	}

	/**
	* Returns Euclidean length squared.
	*/
	inline float norm2(void) const {
		return x*x + y*y + z*z;
	}

	/**
	* Returns unit vector.
	*/
	inline Vector3D unit(void) const {
		float rNorm = 1. / sqrt(x*x + y*y + z*z);
		return Vector3D(rNorm*x, rNorm*y, rNorm*z);
	}

	/**
	* Divides by Euclidean length.
	*/
	inline void normalize(void) {
		(*this) /= norm();
	}

	bool Resembles(const Vector3D & rhs) const
	{
		const Vector3D vDiff = (*this) - rhs;
		static const float svDiffTol = 1.0e-6f;
		return vDiff.norm() < svDiffTol;
	}

}; // class Vector3D

// left scalar multiplication
inline Vector3D operator* (const float& c, const Vector3D& v) {
	return Vector3D(c * v.x, c * v.y, c * v.z);
}

// dot product (a.k.a. inner or scalar product)
inline float dot(const Vector3D& u, const Vector3D& v) {
	return u.x*v.x + u.y*v.y + u.z*v.z;
}

// cross product
inline Vector3D cross(const Vector3D& u, const Vector3D& v) {
	return Vector3D(u.y*v.z - u.z*v.y,
		u.z*v.x - u.x*v.z,
		u.x*v.y - u.y*v.x);
}

inline Vector3D RandomSpread(const Vector3D& vSpread){
	return Vector3D(RandomSpread(vSpread.x),
		RandomSpread(vSpread.y),
		RandomSpread(vSpread.z));
}

// prints components
std::ostream& operator<<(std::ostream& os, const Vector3D& v);

#endif