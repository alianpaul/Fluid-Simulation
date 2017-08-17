#ifndef VEC2D_H
#define VEC2D_H

#include <iosfwd>
#include <cmath>

class Vector2D
{
public:

	// components
	float x, y;

	/**
	* Constructor.
	* Initializes to vector (0,0).
	*/
	Vector2D() : x(0.0), y(0.0) { }

	/**
	* Constructor.
	* Initializes to vector (a,b).
	*/
	Vector2D(float x, float y) : x(x), y(y) { }

	/**
	* Constructor.
	* Copy constructor. Creates a copy of the given vector.
	*/
	Vector2D(const Vector2D& v) : x(v.x), y(v.y) { }

	// additive inverse
	inline Vector2D operator-(void) const {
		return Vector2D(-x, -y);
	}

	// addition
	inline Vector2D operator+(const Vector2D& v) const {
		Vector2D u = *this;
		u += v;
		return u;
	}

	// subtraction
	inline Vector2D operator-(const Vector2D& v) const {
		Vector2D u = *this;
		u -= v;
		return u;
	}

	// right scalar multiplication
	inline Vector2D operator*(float r) const {
		Vector2D vr = *this;
		vr *= r;
		return vr;
	}

	// scalar division
	inline Vector2D operator/(float r) const {
		Vector2D vr = *this;
		vr /= r;
		return vr;
	}

	// add v
	inline void operator+=(const Vector2D& v) {
		x += v.x;
		y += v.y;
	}

	// subtract v
	inline void operator-=(const Vector2D& v) {
		x -= v.x;
		y -= v.y;
	}

	// scalar multiply by r
	inline void operator*=(float r) {
		x *= r;
		y *= r;
	}

	// scalar divide by r
	inline void operator/=(float r) {
		x /= r;
		y /= r;
	}

	/**
	* Returns norm.
	*/
	inline float norm(void) const {
		return sqrt(x*x + y*y);
	}

	/**
	* Returns norm squared.
	*/
	inline float norm2(void) const {
		return x*x + y*y;
	}

	/**
	* Returns unit vector parallel to this one.
	*/
	inline Vector2D unit(void) const {
		return *this / this->norm();
	}


}; // clasd Vector2D

// left scalar multiplication
inline Vector2D operator*(float r, const Vector2D& v) {
	return v*r;
}

// inner product
inline float dot(const Vector2D& v1, const Vector2D& v2) {
	return v1.x*v2.x + v1.y*v2.y;
}

// cross product
inline float cross(const Vector2D& v1, const Vector2D& v2) {
	return v1.x*v2.y - v1.y*v2.x;
}

// prints components
std::ostream& operator<<(std::ostream& os, const Vector2D& v);

#endif