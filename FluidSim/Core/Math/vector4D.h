#ifndef VEC4D_H
#define VEC4D_H

#include "vector3D.h"
#include <iosfwd>

class Vector4D {
public:

	// components
	float x, y, z, w;

	/**
	* Constructor.
	* Initializes tp vector (0,0,0, 0).
	*/
	Vector4D() : x(0.0), y(0.0), z(0.0), w(0.0) { }

	/**
	* Constructor.
	* Initializes to vector (x,y,z,w).
	*/
	Vector4D(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) { }

	/**
	* Constructor.
	* Initializes to vector (x,y,z,0).
	*/
	Vector4D(float x, float y, float z) : x(x), y(y), z(z), w(0.0) { }


	/**
	* Constructor.
	* Initializes to vector (c,c,c,c)
	*/
	Vector4D(float c) : x(c), y(c), z(c), w(c) { }

	/**
	* Constructor.
	* Initializes from existing vector4D.
	*/
	Vector4D(const Vector4D& v) : x(v.x), y(v.y), z(v.z), w(v.w) { }

	/**
	* Constructor.
	* Initializes from existing vector3D.
	*/
	Vector4D(const Vector3D& v) : x(v.x), y(v.y), z(v.z), w(0.0) { }

	/**
	* Constructor.
	* Initializes from existing vector3D and w value.
	*/
	Vector4D(const Vector3D& v, float w) : x(v.x), y(v.y), z(v.z), w(w) { }

	// returns reference to the specified component (0-based indexing: x, y, z)
	inline float& operator[] (const int& index) {
		return (&x)[index];
	}

	// returns const reference to the specified component (0-based indexing: x, y, z)
	inline const float& operator[] (const int& index) const {
		return (&x)[index];
	}

	// negation
	inline Vector4D operator-(void) const {
		return Vector4D(-x, -y, -z, -w);
	}

	// addition
	inline Vector4D operator+(const Vector4D& v) const {
		return Vector4D(x + v.x, y + v.y, z + v.z, w + v.w);
	}

	// subtraction
	inline Vector4D operator-(const Vector4D& v) const {
		return Vector4D(x - v.x, y - v.y, z - v.z, w - v.w);
	}

	// right scalar multiplication
	inline Vector4D operator*(const float& c) const {
		return Vector4D(x * c, y * c, z * c, w * c);
	}

	// scalar division
	inline Vector4D operator/(const float& c) const {
		const float rc = 1.0 / c;
		return Vector4D(rc * x, rc * y, rc * z, rc * w);
	}

	// addition / assignment
	inline void operator+=(const Vector4D& v) {
		x += v.x; y += v.y; z += v.z; z += v.w;
	}

	// subtraction / assignment
	inline void operator-=(const Vector4D& v) {
		x -= v.x; y -= v.y; z -= v.z; w -= v.w;
	}

	// scalar multiplication / assignment
	inline void operator*=(const float& c) {
		x *= c; y *= c; z *= c; w *= c;
	}

	// scalar division / assignment
	inline void operator/=(const float& c) {
		(*this) *= (1. / c);
	}

	/**
	* Returns Euclidean distance metric extended to 4 dimensions.
	*/
	inline float norm(void) const {
		return sqrt(x*x + y*y + z*z + w*w);
	}

	/**
	* Returns Euclidean length squared.
	*/
	inline float norm2(void) const {
		return x*x + y*y + z*z + w*w;
	}

	/**
	* Returns unit vector. (returns the normalized copy of this vector.)
	*/
	inline Vector4D unit(void) const {
		float rNorm = 1. / sqrt(x*x + y*y + z*z + w*w);
		return Vector4D(rNorm*x, rNorm*y, rNorm*z);
	}

	/**
	* Divides by Euclidean length.
	* This vector will be of unit length i.e. "normalized" afterwards.
	*/
	inline void normalize(void) {
		(*this) /= norm();
	}

	/**
	* Converts this vector to a 3D vector ignoring the w component.
	*/
	Vector3D to3D();

	/**
	* Converts this vector to a 3D vector by dividing x, y, and z by w.
	*/
	Vector3D projectTo3D();

}; // class Vector4D

// left scalar multiplication
inline Vector4D operator* (const float& c, const Vector4D& v) {
	return Vector4D(c * v.x, c * v.y, c * v.z, c*v.w);
}

// dot product (a.k.a. inner or scalar product)
inline float dot(const Vector4D& u, const Vector4D& v) {
	return u.x*v.x + u.y*v.y + u.z*v.z + u.w*v.w;;
}

// prints components
std::ostream& operator<<(std::ostream& os, const Vector4D& v);

#endif