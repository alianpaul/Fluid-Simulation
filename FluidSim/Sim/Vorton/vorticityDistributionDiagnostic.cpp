#include "vorticityDistribution.h"
#include "Space\uniformGrid.h"

void VortDistributionTest(float fMagnitude, unsigned numVortonsMax, const IVorticityDistribution & vorticityDistribution)
{
	const Vector3D			vCenter(0.0f, 0.0f, 0.0f);
	const Vector3D			vDomainExtent = vorticityDistribution.GetDomainSize();
	const Vector3D			vMin = vCenter - 0.5f * vDomainExtent;
	const Vector3D			vMax = vCenter + 0.5f * vDomainExtent;
	UniformGrid<Vorton>		ugVortons(numVortonsMax, vMin, vMax, true);
	unsigned				numCells[3] = { ugVortons.GetNumCells(0), ugVortons.GetNumCells(1), ugVortons.GetNumCells(2) };
	const unsigned			numX = ugVortons.GetNumPoints(0);
	const unsigned			numXY = ugVortons.GetNumPoints(0) * ugVortons.GetNumPoints(1);

	const float				oneOverN[3] = { 1.0f / float(numCells[0]), 1.0f / float(numCells[1]), 1.0f / float(numCells[2]) };
	const Vector3D			vCellExtent(vDomainExtent.x * oneOverN[0], vDomainExtent.y * oneOverN[1], vDomainExtent.z * oneOverN[2]);
	const Vector3D			vNoise(0.0f * vCellExtent); //?? Original version is 0.0f * vCellExtent
	float					vortonRadius = powf(vCellExtent.x * vCellExtent.y * vCellExtent.z, 1.0f / 3.0f) * 0.5f;
	

	ugVortons.Init();
	Vector3D	position;
	unsigned	idx[3];
	for (idx[2] = 0; idx[2] < numCells[2]; ++idx[2])
	{
		const unsigned offsetZ = idx[2] * numXY;
		position.z				= vMin.z + (float(idx[2]) + 0.25) * vCellExtent.z;
		
		for (idx[1] = 0; idx[1] < numCells[1]; ++idx[1])
		{
			const unsigned offsetYZ = offsetZ + idx[1] * numX;
			position.y = vMin.y + (float(idx[1]) + 0.25) * vCellExtent.y;

			for (idx[0] = 0; idx[0] < numCells[0]; ++idx[0])
			{
				const unsigned offsetXYZ = offsetYZ + idx[0];
				position.x = vMin.x + (float(idx[0]) + 0.25) * vCellExtent.x;
				position += RandomSpread(vNoise);

				Vector3D vorticity;
				vorticityDistribution.AssignVorticity(vorticity, position, vCenter);

				Vorton vorton(position, vorticity * fMagnitude, vortonRadius);
				ugVortons[offsetXYZ] = vorton;
			}
		}
	}

	ugVortons.GenerateBrickOfBytes("test", 0);
}
