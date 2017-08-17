#include "uniformGridMath.h"


void ComputeCurlFromJacobian(UniformGrid<Vector3D> & curl, const UniformGrid<Matrix3x3> & jacobian)
{

	const unsigned dims[3] = { jacobian.GetNumPoints(0), jacobian.GetNumPoints(1), jacobian.GetNumPoints(2) };
	const unsigned numXY = dims[0] * dims[1];
	unsigned index[3];

	for ( index[2] = 0; index[2] < dims[2]; ++index[2] )
	{
		const unsigned offsetZ = numXY * index[2];
		for ( index[1] = 0; index[1] < dims[1]; ++index[1] )
		{
			const unsigned offsetYZ = offsetZ + index[1] * dims[0];
			for (index[0] = 0; index[0] < dims[0]; ++index[0])
			{
				const unsigned offsetXYZ = offsetYZ + index[0];
				const Matrix3x3 & rj     = jacobian[offsetXYZ];
				Vector3D & vc            = curl[offsetXYZ];
				vc = Vector3D(rj(2, 1) - rj(1, 2), rj(0, 2) - rj(2, 0), rj(1, 0) - rj(0, 1));
			}
		}
	}

}


void ComputeJacobian(UniformGrid<Matrix3x3> & jacobian, const UniformGrid<Vector3D> & vec)
{

	const Vector3D reciprocalSpacing = vec.GetCellPerExtent();
	const Vector3D halfreciprocalSpacing = reciprocalSpacing * 0.5;

	const unsigned dims[3] = { vec.GetNumPoints(0), vec.GetNumPoints(1), vec.GetNumPoints(2) };
	const unsigned dimsMinus1[3] = { dims[0] - 1, dims[1] - 1, dims[2] - 1 };
	const unsigned numXY = dims[0] * dims[1];

	unsigned index[3];
	/*
		Compute the Jacobian of grid points inside the 6 faces
		Attention: 
				Although UniformGrid has a member function called OffsetFromIndices,
			But it's not efficient to be use when we iterate through all the grid points inside
			the uniform grid. Because everytime we get the offset, we need to do 4 operations,
			no matter which grid point we are in.
				For example, in the inner most loop when we increate the x-index, the offset only
			needs to be added by 1, not all 4 operations everytime. So it saves a lot we manually maintain
			the offset when we iterate through the grid points.

	*/

#define ASSIGN_Z_OFFSETS								\
	const unsigned offsetZ0 = index[2] * numXY;			\
	const unsigned offsetZM = (index[2] - 1) * numXY;	\
	const unsigned offsetZP = (index[2] + 1) * numXY;

#define ASSIGN_YZ_OFFSETS												\
	const unsigned offsetY0Z0 = offsetZ0 +  index[1]      * dims[0];	\
	const unsigned offsetYMZ0 = offsetZ0 + (index[1] - 1) * dims[0];	\
	const unsigned offsetYPZ0 = offsetZ0 + (index[1] + 1) * dims[0];	\
	const unsigned offsetY0ZM = offsetZM +  index[1]      * dims[0];	\
	const unsigned offsetY0ZP = offsetZP +  index[1]      * dims[0];

#define ASSIGN_XYZ_OFFSETS										\
	const unsigned offsetX0Y0Z0 = offsetY0Z0 + index[0];		\
	const unsigned offsetXMY0Z0 = offsetY0Z0 + index[0] - 1;	\
	const unsigned offsetXPY0Z0 = offsetY0Z0 + index[0] + 1;	\
	const unsigned offsetX0YMZ0 = offsetYMZ0 + index[0];		\
	const unsigned offsetX0YPZ0 = offsetYPZ0 + index[0];		\
	const unsigned offsetX0Y0ZM = offsetY0ZM + index[0];		\
	const unsigned offsetX0Y0ZP = offsetY0ZP + index[0];		

	
	for (index[2] = 1; index[2] < dimsMinus1[2]; ++index[2])
	{
		
		ASSIGN_Z_OFFSETS;
		for (index[1] = 1; index[1] < dimsMinus1[1]; ++index[1])
		{
			
			ASSIGN_YZ_OFFSETS;
			for (index[0] = 1; index[0] < dimsMinus1[0]; ++index[0])
			{

				ASSIGN_XYZ_OFFSETS;
				Matrix3x3 & rMatrix = jacobian[offsetX0Y0Z0];
				// d vec / d x
				rMatrix[0] = (vec[offsetXPY0Z0] - vec[offsetXMY0Z0]) * halfreciprocalSpacing.x;
				// d vec / d y
				rMatrix[1] = (vec[offsetX0YPZ0] - vec[offsetX0YMZ0]) * halfreciprocalSpacing.y;
				// d vec / d z
				rMatrix[2] = (vec[offsetX0Y0ZP] - vec[offsetX0Y0ZM]) * halfreciprocalSpacing.z;
				
			}
		}
	}

	/*
		The outter six faces;
		Attention:
				Using COMPUTE_FINITE_DIFF, We can do one iterate over all grid points and calculate the Jacobian Matrix.
			But its not efficient, because for most of the points inside extent (1, dimMinus1 -1, 1, dimsMinus1 - 1, 1, dimsMinus1 - 1)
			can calculate their Jabobian directly. There is no need to swith on the 3 cases(+, -, inside+ -).
				So to optimize the efficiency.We first iterate through the inside extent like above. Then we calculate each face's points Jacobian.
			like below. We will calculate the points on the edge twice, but it will not affect the efficiency every much.
		*/

#define COMPUTE_FINITE_DIFF \
	Matrix3x3 & rMatrix = jacobian[offsetX0Y0Z0]; \
	if( index[0] == 0)				   { rMatrix[0] = (vec[offsetXPY0Z0] - vec[offsetX0Y0Z0]) * reciprocalSpacing.x; } \
	else if(index[0] == dimsMinus1[0]) { rMatrix[0] = (vec[offsetX0Y0Z0] - vec[offsetXMY0Z0]) * reciprocalSpacing.x; } \
	else						       { rMatrix[0] = (vec[offsetXPY0Z0] - vec[offsetXMY0Z0]) * halfreciprocalSpacing.x;} \
	if( index[1] == 0)				   { rMatrix[1] = (vec[offsetX0YPZ0] - vec[offsetX0Y0Z0]) * reciprocalSpacing.y; } \
	else if(index[1] == dimsMinus1[1]) { rMatrix[1] = (vec[offsetX0Y0Z0] - vec[offsetX0YMZ0]) * reciprocalSpacing.y; } \
	else						       { rMatrix[1] = (vec[offsetX0YPZ0] - vec[offsetX0YMZ0]) * halfreciprocalSpacing.y;} \
	if( index[2] == 0)				   { rMatrix[2] = (vec[offsetX0Y0ZP] - vec[offsetX0Y0Z0]) * reciprocalSpacing.z; } \
	else if(index[2] == dimsMinus1[2]) { rMatrix[2] = (vec[offsetX0Y0Z0] - vec[offsetX0Y0ZM]) * reciprocalSpacing.z; } \
	else						       { rMatrix[2] = (vec[offsetX0Y0ZP] - vec[offsetX0Y0ZM]) * halfreciprocalSpacing.x;} 

	//-X
	index[0] = 0;
	for (index[2] = 0; index[2] < dims[2]; ++index[2])
	{
		ASSIGN_Z_OFFSETS;
		for (index[1] = 0; index[1] < dims[1]; ++index[1])
		{
			ASSIGN_YZ_OFFSETS;
			{
				ASSIGN_XYZ_OFFSETS;
				COMPUTE_FINITE_DIFF;
			}
		}

	}

	//+X
	index[0] = dimsMinus1[0];
	for (index[2] = 0; index[2] < dims[2]; ++index[2])
	{
		ASSIGN_Z_OFFSETS;
		for (index[1] = 0; index[1] < dims[1]; ++index[1])
		{
			ASSIGN_YZ_OFFSETS;
			{
				ASSIGN_XYZ_OFFSETS;
				COMPUTE_FINITE_DIFF;
			}
		}

	}

	//-Y
	index[1] = 0;
	for (index[2] = 0; index[2] < dims[2]; ++index[2])
	{
		ASSIGN_Z_OFFSETS;
		{
			
			ASSIGN_YZ_OFFSETS;
			for (index[0] = 0; index[0] < dims[0]; ++index[0])
			{
				ASSIGN_XYZ_OFFSETS;
				COMPUTE_FINITE_DIFF;
			}
		}
	}

	//+Y
	index[1] = dimsMinus1[1];
	for (index[2] = 0; index[2] < dims[2]; ++index[2])
	{
		ASSIGN_Z_OFFSETS;
		{

			ASSIGN_YZ_OFFSETS;
			for (index[0] = 0; index[0] < dims[0]; ++index[0])
			{
				ASSIGN_XYZ_OFFSETS;
				COMPUTE_FINITE_DIFF;
			}
		}
	}

	//-Z
	index[2] = 0;
	{
		
		ASSIGN_Z_OFFSETS;
		for (index[1] = 0; index[1] < dims[1]; ++index[1])
		{
			ASSIGN_YZ_OFFSETS;
			for (index[0] = 0; index[0] < dims[0]; ++index[0])
			{
				ASSIGN_XYZ_OFFSETS;
				COMPUTE_FINITE_DIFF;
			}
		}
	
	}

	//+Z
	index[2] = dimsMinus1[2];
	{
		ASSIGN_Z_OFFSETS;
		for (index[1] = 0; index[1] < dims[1]; ++index[1])
		{
			ASSIGN_YZ_OFFSETS;
			for (index[0] = 0; index[0] < dims[0]; ++index[0])
			{
				ASSIGN_XYZ_OFFSETS;
				COMPUTE_FINITE_DIFF;
			}
		}
	}

#undef COMPUTE_FINITE_DIFF
#undef ASSIGN_XYZ_OFFSETS
#undef ASSIGN_YZ_OFFSETS
#undef ASSIGN_X_OFFSETS

}

