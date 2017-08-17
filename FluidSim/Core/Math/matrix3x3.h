#ifndef MAT3X3_H
#define MAT3X3_H

#include <iosfwd>
#include "vector3D.h"

class Matrix3x3 {

public:

	// The default constructor.
	Matrix3x3(void) { }

	// Constructor for row major form data.
	// Transposes to the internal column major form.
	// REQUIRES: data should be of size 9 for a 3 by 3 matrix..
	Matrix3x3(float * data)
	{
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				// Transpostion happens within the () query.
				(*this)(i, j) = data[i * 3 + j];
			}
		}
	}

	/**
	* Sets all elements to val.
	*/
	void zero(float val = 0.0);

	/**
	* Returns the determinant of A.
	*/
	float det(void) const;

	/**
	* Returns the Frobenius norm of A.
	*/
	float norm(void) const;

	/**
	* Returns the 3x3 identity matrix.
	*/
	static Matrix3x3 identity(void);

	/**
	* Returns a matrix representing the (left) cross product with u.
	*/
	static Matrix3x3 crossProduct(const Vector3D& u);

	/**
	* Returns the ith column.
	*/
	Vector3D& column(int i);
	const Vector3D& column(int i) const;

	/**
	* Returns the transpose of A.
	*/
	Matrix3x3 T(void) const;

	/**
	* Returns the inverse of A.
	*/
	Matrix3x3 inv(void) const;

	// accesses element (i,j) of A using 0-based indexing
	float& operator()(int i, int j);
	const float& operator()(int i, int j) const;

	// accesses the ith column of A
	Vector3D& operator[](int i);
	const Vector3D& operator[](int i) const;

	// increments by B
	void operator+=(const Matrix3x3& B);

	// returns -A
	Matrix3x3 operator-(void) const;

	// returns A-B
	Matrix3x3 operator-(const Matrix3x3& B) const;

	Matrix3x3 operator+(const Matrix3x3& B) const;

	// returns c*A
	Matrix3x3 operator*(float c) const;

	// returns A*B
	Matrix3x3 operator*(const Matrix3x3& B) const;

	// returns A*x
	Vector3D operator*(const Vector3D& x) const;

	// divides each element by x
	void operator/=(float x);

protected:

	// column vectors
	Vector3D entries[3];

}; // class Matrix3x3

// returns the outer product of u and v
Matrix3x3 outer(const Vector3D& u, const Vector3D& v);

// returns c*A
Matrix3x3 operator*(float c, const Matrix3x3& A);

// prints entries
std::ostream& operator<<(std::ostream& os, const Matrix3x3& A);

#endif