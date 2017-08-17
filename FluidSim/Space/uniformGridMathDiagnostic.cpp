#include "uniformGridMath.h"

void GridMathUnitTest()
{
	/*
		The velocity field is:
		V(x, y, z) = [y, -x, 0];
		the curl should be [0, 0, -2] at all points
	*/

	const unsigned			num = 1024;
	Vector3D				vRange    (2.0f, 3.0f, 5.0f);
	Vector3D				vMin      (-0.5f * vRange);
	UniformGrid<Vector3D>	ugVelocity(num, -0.5f * vRange, 0.5f * vRange, true);
	ugVelocity.Init();
	const unsigned			numXY = ugVelocity.GetNumPoints(0) * ugVelocity.GetNumPoints(1);
	const unsigned			dims[3] = { ugVelocity.GetNumPoints(0), ugVelocity.GetNumPoints(1), ugVelocity.GetNumPoints(2) };
	unsigned				index[3];
	Vector3D				vCellExtent = ugVelocity.GetCellSpacing();

	for (index[2] = 0; index[2] < ugVelocity.GetNumPoints(2); ++index[2])
	{
		const unsigned offsetZ = index[2] * numXY;
		      Vector3D vPos    = vMin + index[2] * vCellExtent.z;

		for (index[1] = 0; index[1] < ugVelocity.GetNumPoints(1); ++index[1])
		{

			const unsigned offsetYZ = index[1] * dims[0] + offsetZ;
			               vPos.y   = vMin.y + index[1] * vCellExtent.y;

			for (index[0] = 0; index[0] < ugVelocity.GetNumPoints(0); ++index[0])
			{
				const unsigned offsetXYZ = index[0] + offsetYZ;
							   vPos.x    = vMin.x + index[0] * vCellExtent.x;

				ugVelocity[offsetXYZ].x = vPos.y;
				ugVelocity[offsetXYZ].y = -vPos.x;
				ugVelocity[offsetXYZ].z = 0;
			}
		}
	}

	UniformGrid<Vector3D>	ugCurl(ugVelocity);
	ugCurl.Init();
	UniformGrid<Matrix3x3>  ugJacobian(ugVelocity);
	ugJacobian.Init();

	ComputeJacobian(ugJacobian, ugVelocity);
	ComputeCurlFromJacobian(ugCurl, ugJacobian);

	/*
		Check if all points is [0, 0, -2]
	*/
	const Vector3D ans(0.0f, 0.0f, -2.0f);
	for (index[2] = 0; index[2] < ugVelocity.GetNumPoints(2); ++index[2])
	{
		const unsigned offsetZ = index[2] * numXY;
	
		for (index[1] = 0; index[1] < ugVelocity.GetNumPoints(1); ++index[1])
		{

			const unsigned offsetYZ = index[1] * dims[0] + offsetZ;

			for (index[0] = 0; index[0] < ugVelocity.GetNumPoints(0); ++index[0])
			{
				const unsigned offsetXYZ = index[0] + offsetYZ;
				_ASSERT(ugCurl[offsetXYZ].Resembles(ans));
			}
		}
	}
}