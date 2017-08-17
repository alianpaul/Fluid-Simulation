#ifndef UNIFORM_GRID_H
#define UNIFORM_GRID_H

#include <GL/glew.h>
#include <GL/glut.h>

#include "Core\Math\vector3D.h"
#include "Core\Utility.h"
#include <math.h>
#include <float.h>
#include <vector>

class UniformGridGeometry
{
public:
	UniformGridGeometry()
		: mMinCorner(0.0f, 0.0f, 0.0f)
		, mGridExtent(0.0f, 0.0f, 0.0f)
		, mCellExtent(0.0f, 0.0f, 0.0f)
		, mCellPerExtent(0.0f, 0.0f, 0.0f)
	{
		mNumPoints[0] = mNumPoints[1] = mNumPoints[2] = 0;
	}

	UniformGridGeometry(unsigned uNumElements , const Vector3D & vMin , const Vector3D & vMax, bool bPowerOf2)
	{
		DefineShape(uNumElements, vMin, vMax, bPowerOf2);
	}
	
	void	DefineShape(unsigned uNumElements, const Vector3D & vMin, const Vector3D & vMax, bool bPowerOf2)
	{
		mMinCorner = vMin;
		static const float Nudge = 1.0f + FLT_EPSILON;
		mGridExtent = (vMax - vMin) * Nudge;


		int numDims = 3;
		Vector3D vSizeEffective(GetExtent());
		if (mGridExtent.x == 0.0f)
		{
			vSizeEffective.x = 1.0f;
			--numDims;
		}
		else if (mGridExtent.y  == 0.0f)
		{
			vSizeEffective.y = 1.0f;
			--numDims;
		}
		else if (mGridExtent.z == 0.0f)
		{
			vSizeEffective.z = 1.0f;
			--numDims;
		}

		float fSize = vSizeEffective.x * vSizeEffective.y * vSizeEffective.z;
		float fCellSideSize = powf(fSize / float(uNumElements) , 1.0f / float(numDims));

		int numCells[3] = {	MAX2(1, unsigned(mGridExtent.x / fCellSideSize + 0.5)),
							MAX2(1, unsigned(mGridExtent.y / fCellSideSize + 0.5)),
							MAX2(1, unsigned(mGridExtent.z / fCellSideSize + 0.5))};

		if (bPowerOf2)
		{
			numCells[0] = NearestPowerOfTwo(numCells[0]);
			numCells[1] = NearestPowerOfTwo(numCells[1]);
			numCells[2] = NearestPowerOfTwo(numCells[2]);
		}

		/*???*/
		while (numCells[0] * numCells[1] * numCells[2] >= uNumElements * 8)
		{
			numCells[0] = MAX2(1, numCells[0] / 2);
			numCells[1] = MAX2(1, numCells[1] / 2);
			numCells[2] = MAX2(1, numCells[2] / 2);
		}

		mNumPoints[0] = numCells[0] + 1;
		mNumPoints[1] = numCells[1] + 1;
		mNumPoints[2] = numCells[2] + 1;

		PrecomputeSpacing();
	}

	void Decimate(const UniformGridGeometry & src, int iDecimation)
	{
		mMinCorner = src.mMinCorner;
		mGridExtent = src.mGridExtent;

		//cells = src.cells/iDecimation
		mNumPoints[0] = src.GetNumCells(0) / iDecimation + 1;
		mNumPoints[1] = src.GetNumCells(1) / iDecimation + 1;
		mNumPoints[2] = src.GetNumCells(2) / iDecimation + 1;

		mNumPoints[0] = MAX2(2, GetNumPoints(0));
		mNumPoints[1] = MAX2(2, GetNumPoints(1));
		mNumPoints[2] = MAX2(2, GetNumPoints(2));

		PrecomputeSpacing();
	}

	void CopyShape(const UniformGridGeometry & src)
	{
		Decimate(src, 1);
	}

	const Vector3D &	GetExtent() const	{ return mGridExtent; }

	unsigned			GetNumCells(unsigned idx) const	{ return mNumPoints[idx] - 1; }

	unsigned			GetNumPoints(unsigned idx) const { return mNumPoints[idx]; }

	unsigned			GetGridCapacity() const { return GetNumPoints(0) * GetNumPoints(1) * GetNumPoints(2); }

	const Vector3D &	GetCellSpacing() const { return mCellExtent; }

	const Vector3D &	GetCellPerExtent() const { return mCellPerExtent; }

	const Vector3D &	GetMinCorner() const { return mMinCorner; }

	/*
		the indices of the Grid Points
	*/
	void				IndicesOfPosition(unsigned indices[3], const Vector3D & vPosition) const
	{
		Vector3D vRelDis(vPosition - GetMinCorner());
		indices[0] = unsigned(vRelDis.x * GetCellPerExtent().x); //Down 0.9 -> 0
		indices[1] = unsigned(vRelDis.y * GetCellPerExtent().y);
		indices[2] = unsigned(vRelDis.z * GetCellPerExtent().z);
	}

	/*
		the actual contents array is 1-D,
		this function map the 3-D point to the corresponding grid point indices,
		then map the grid point indices to the 1-D array offset
		the offset of point = grid point before it
	*/
	unsigned			OffsetOfPosition(const Vector3D & vPosition) const
	{
		unsigned indices[3];
		IndicesOfPosition(indices, vPosition);
		unsigned offset = indices[0] + GetNumPoints(0) * (indices[1] + GetNumPoints(1) * indices[2]);
		//indices[0] + GetNumPoints(0)*indices[1] + GetNumPoints(0)*GetNumPoints(1)*indices[2] 
		return offset;
	}

	void				PositionFromIndices(Vector3D & vPosition, const unsigned indices[3]) const
	{
		vPosition.x = GetMinCorner().x + float(indices[0]) * GetCellSpacing().x;
		vPosition.y = GetMinCorner().y + float(indices[1]) * GetCellSpacing().y;
		vPosition.z = GetMinCorner().z + float(indices[2]) * GetCellSpacing().z;
	}

	void				IndicesFromOffset(unsigned indices[3], const unsigned & offset) const
	{
		indices[2] = offset / (GetNumPoints(0) * GetNumPoints(1));
		indices[1] = (offset - indices[2] * GetNumPoints(0) * GetNumPoints(1)) / GetNumPoints(0);
		indices[0] = offset - GetNumPoints(0) * (indices[1] + GetNumPoints(1) * indices[2]);
	}

	void				PositionFromOffset(Vector3D & vPosition, const unsigned & offset)
	{
		unsigned indices[3];
		IndicesFromOffset(indices, offset);
		vPosition.x = GetMinCorner().x + float(indices[0]) * GetCellSpacing().x;
		vPosition.y = GetMinCorner().y + float(indices[1]) * GetCellSpacing().y;
		vPosition.z = GetMinCorner().z + float(indices[2]) * GetCellSpacing().z;
	}


	/*For Debug, draw the grid out
	*/
	void				Render()
	{
		glDisable(GL_DEPTH_TEST);
		glDisable(GL_LIGHTING);

		glMatrixMode(GL_MODELVIEW);
		glBegin(GL_LINES);
		glColor4f(0.5f, 0.5f, 0.5f, 0.5f);

		for (int idx = 0; idx < mNumPoints[0]; idx++)
		{
			for (int idy = 0; idy < mNumPoints[1]; idy++)
			{
				glVertex3f(
					mMinCorner.x + idx * mCellExtent.x,
					mMinCorner.y + idy * mCellExtent.y,
					mMinCorner.z);
				glVertex3f(
					mMinCorner.x + idx * mCellExtent.x,
					mMinCorner.y + idy * mCellExtent.y,
					mMinCorner.z + mGridExtent.z);
			}
		}

		for (int idz = 0; idz < mNumPoints[2]; idz++)
		{
			for (int idy = 0; idy < mNumPoints[1]; idy++)
			{
				glVertex3f(
					mMinCorner.x,
					mMinCorner.y + idy * mCellExtent.y,
					mMinCorner.z + idz * mCellExtent.z);
				glVertex3f(
					mMinCorner.x + mGridExtent.x,
					mMinCorner.y + idy * mCellExtent.y,
					mMinCorner.z + idz * mCellExtent.z);
			}
		}

		for (int idz = 0; idz < mNumPoints[2]; idz++)
		{
			for (int idx = 0; idx < mNumPoints[0]; idx++)
			{
				glVertex3f(
					mMinCorner.x + idx * mCellExtent.x,
					mMinCorner.y,
					mMinCorner.z + idz * mCellExtent.z);
				glVertex3f(
					mMinCorner.x + idx * mCellExtent.x,
					mMinCorner.y + mGridExtent.y,
					mMinCorner.z + idz * mCellExtent.z);
			}
		}

		glEnd();

		glEnable(GL_LIGHTING);
		glEnable(GL_DEPTH_TEST);
	}
protected:
	void PrecomputeSpacing()
	{
		mCellExtent.x = GetExtent().x / float( GetNumCells(0) );
		if (GetExtent().x == 0.0f) mCellPerExtent.x = 1.0f / FLT_MIN;
		else mCellPerExtent.x = float(GetNumCells(0)) / GetExtent().x;

		mCellExtent.y = GetExtent().y / float( GetNumCells(1) );
		if (GetExtent().y == 0.0f) mCellPerExtent.y = 1.0f / FLT_MIN;
		else mCellPerExtent.y = float(GetNumCells(1)) / GetExtent().y;

		mCellExtent.z = GetExtent().z / float( GetNumCells(2) );
		if (GetExtent().z == 0.0f) mCellPerExtent.z = 1.0f / FLT_MIN;
		else mCellPerExtent.z = float(GetNumCells(2)) / GetExtent().z;
	}

	unsigned GetOffsetOfPenultimateCell(void)
	{
		return (GetNumPoints(0) - 2) + GetNumPoints(0) * ((GetNumPoints(1) - 2) + GetNumPoints(1) * (GetNumPoints(2) - 2));
	}

	unsigned OffsetFromIndices(const unsigned indices[3]) const
	{
		return indices[0] + GetNumPoints(0) * (indices[1] + GetNumPoints(1) * indices[2]);
	}

	void Clear(void)
	{
		mMinCorner =
		mGridExtent =
		mCellExtent =
		mCellPerExtent = Vector3D(0.0f, 0.0f, 0.0f);
		mNumPoints[0] = mNumPoints[1] = mNumPoints[2] = 0;
	}

	Vector3D	mMinCorner; //the axis-aligned Bounding Box's min corner(in world)
	Vector3D	mGridExtent; //the Size of the Bounding Box's
	Vector3D	mCellExtent; //the Size of the Cell in the Bounding Box
	Vector3D	mCellPerExtent; //the Reciprocal of mCellExtent
	unsigned	mNumPoints[3]; //the Number of the Grid Points on the axis 0-x 1-y 2-z
};


template <class ItemT> class UniformGrid : public UniformGridGeometry
{
public:
	
	typedef UniformGridGeometry Parent;

	UniformGrid() : UniformGridGeometry() {}

	UniformGrid(unsigned uNumElements, const Vector3D & vMin, const Vector3D & vMax, bool bPowerOf2)
		: UniformGridGeometry(uNumElements, vMin, vMax, bPowerOf2)
	{

	}

	explicit UniformGrid(const UniformGridGeometry & rhs)
		: UniformGridGeometry(rhs)
	{

	}

	/*
		UniformGrid copy ctor does not copy the mContents!!!
		massive memory moves is inefficient
		UniformGrid will also be pushed onto a vector in NestedGrid
		this copy ctor will be call
		Base class part must be copied(with Base class part copy ctor, we can not
		access it directly, because the member is private, so we use the copy ctor)
		Deliberately miss the Derived class member for efficiency
	*/
	UniformGrid(const UniformGrid & rhs)
		: UniformGridGeometry( rhs )
	{

	}

	/*
		one Grid Point stores one Content
	*/

		  ItemT &	operator[](const unsigned offset) 			{ return mContents[offset]; }
	const ItemT &	operator[](const unsigned offset) const		{ return mContents[offset]; }
		  ItemT &	operator[](const Vector3D & vPosition)		{ return mContents[OffsetOfPosition(vPosition)]; }

	void Init(void)
	{
		mContents.resize( GetGridCapacity() );
	}

	void DefineShape(unsigned uNumElements, const Vector3D & vMin, const Vector3D & vMax, bool bPowerOf2)
	{
		mContents.clear();
		UniformGridGeometry::DefineShape(uNumElements, vMin, vMax, bPowerOf2);
	}

	size_t Size() const { return mContents.size(); }

	void Decimate(const UniformGridGeometry & src, int iDecimation)
	{
		UniformGridGeometry::Decimate(src, iDecimation);
	}

	void ComputeStatistics(ItemT & min, ItemT & max) const
	{
		max = min = (*this)[0];
		unsigned numCells = GetGridCapacity();
		for (unsigned offset = 1; offset < numCells; ++offset)
		{
			const ItemT & rVal = (*this)[offset];
			min = MIN2( min, rVal);
			max = MAX2( max, rVal);
		}
	}

	void Interpolate(ItemT & vResult, const Vector3D & vPosition) const
	{
		

		unsigned indices[3];
		IndicesOfPosition(indices, vPosition);
		Vector3D vMinCorner;
		PositionFromIndices(vMinCorner, indices);

		unsigned numXY = GetNumPoints(0) * GetNumPoints(1);
		unsigned offsetX0Y0Z0 = OffsetOfPosition( vMinCorner );
		unsigned offsetX1Y0Z0 = offsetX0Y0Z0 + 1;
		unsigned offsetX0Y1Z0 = offsetX0Y0Z0 + GetNumPoints(0);
		unsigned offsetX1Y1Z0 = offsetX0Y0Z0 + GetNumPoints(0) + 1;
		unsigned offsetX0Y0Z1 = offsetX0Y0Z0 + numXY;
		unsigned offsetX1Y0Z1 = offsetX0Y0Z0 + numXY + 1;
		unsigned offsetX0Y1Z1 = offsetX0Y0Z0 + numXY + GetNumPoints(0);
		unsigned offsetX1Y1Z1 = offsetX0Y0Z0 + numXY + GetNumPoints(0) + 1;
		
		//Debug

		Vector3D vDiff		= vPosition - vMinCorner;
		Vector3D vT			= Vector3D(vDiff.x * GetCellPerExtent().x, vDiff.y * GetCellPerExtent().y, vDiff.z * GetCellPerExtent().z);
		Vector3D vOneMinusT = Vector3D(1.0f, 1.0f, 1.0f) - vT;

		vResult = vOneMinusT.x * vOneMinusT.y * vOneMinusT.z * (*this)[offsetX0Y0Z0]
				+		  vT.x * vOneMinusT.y * vOneMinusT.z * (*this)[offsetX1Y0Z0]
				+ vOneMinusT.x *         vT.y * vOneMinusT.z * (*this)[offsetX0Y1Z0]
				+		  vT.x *	     vT.y * vOneMinusT.z * (*this)[offsetX1Y1Z0]
				+ vOneMinusT.x * vOneMinusT.y *			vT.z * (*this)[offsetX0Y0Z1]
				+         vT.x * vOneMinusT.y * 		vT.z * (*this)[offsetX1Y0Z1]
				+ vOneMinusT.x *         vT.y * 		vT.z * (*this)[offsetX0Y1Z1]
				+         vT.x *	     vT.y * 		vT.z * (*this)[offsetX1Y1Z1];

	}

	/*Given a position and item value in this position, 
	calculate it's cell 8 grid point itemValue
	and add to their origin value*/
	void Insert(const Vector3D & vPosition, const ItemT & item)
	{
		unsigned numXY = GetNumPoints(0) * GetNumPoints(1);
		unsigned offsetX0Y0Z0 = OffsetOfPosition(vPosition);
		unsigned offsetX1Y0Z0 = offsetX0Y0Z0 + 1;
		unsigned offsetX0Y1Z0 = offsetX0Y0Z0 + GetNumPoints(0);
		unsigned offsetX1Y1Z0 = offsetX0Y0Z0 + GetNumPoints(0) + 1;
		unsigned offsetX0Y0Z1 = offsetX0Y0Z0 + numXY;
		unsigned offsetX1Y0Z1 = offsetX0Y0Z0 + numXY + 1;
		unsigned offsetX0Y1Z1 = offsetX0Y0Z0 + numXY + GetNumPoints(0);
		unsigned offsetX1Y1Z1 = offsetX0Y0Z0 + numXY + GetNumPoints(0) + 1;

		unsigned indices[3];
		IndicesFromOffset(indices, offsetX0Y0Z0);
		Vector3D vMinCorner;
		PositionFromIndices(vMinCorner, indices);
		Vector3D vDiff = vPosition - vMinCorner;
		Vector3D vT = Vector3D(vDiff.x * GetCellPerExtent().x, vDiff.y * GetCellPerExtent().y, vDiff.z * GetCellPerExtent().z);
		Vector3D vOneMinusT = Vector3D(1.0f, 1.0f, 1.0f) - vT;

		(*this)[offsetX0Y0Z0] += vOneMinusT.x * vOneMinusT.y * vOneMinusT.z * item;
		(*this)[offsetX1Y0Z0] +=		 vT.x * vOneMinusT.y * vOneMinusT.z * item;
		(*this)[offsetX0Y1Z0] += vOneMinusT.x *         vT.y * vOneMinusT.z * item;
		(*this)[offsetX1Y1Z0] +=         vT.x *         vT.y * vOneMinusT.z * item;
		(*this)[offsetX0Y0Z1] += vOneMinusT.x * vOneMinusT.y *         vT.z * item;
		(*this)[offsetX1Y0Z1] +=         vT.x * vOneMinusT.y *         vT.z * item;
		(*this)[offsetX0Y1Z1] += vOneMinusT.x *         vT.y *         vT.z * item;
		(*this)[offsetX1Y1Z1] +=         vT.x *         vT.y *         vT.z * item;
	}


	void			Clear(void)
	{
		mContents.clear();
		Parent::Clear();
	}

	static void		UnitTest(void);

	void			GenerateBrickOfBytes(const char * strFilenameBase, unsigned uFrame) const;

private:
	std::vector<ItemT>	 mContents;
};

#endif