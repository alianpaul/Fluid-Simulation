#include "vorticityDistribution.h"
#include "Space\uniformGrid.h"

/*! \brief A very small number, between FLT_EPSILON and FLT_MIN.
*/
static const float sTiny = expf(0.5f * (logf(FLT_EPSILON) + logf(FLT_MIN)));

/*
	Return vortons in the domain,
	The domain and the vorticity of the vortons is 
	determined by IVorticityDistribution.
*/

void AssignVorticity(std::vector<Vorton> & vortons, float fMagnitude, unsigned numVortonsMax, const IVorticityDistribution & vorticityDistribution)
{
	const Vector3D			vCenter			(0.0f, 0.0f, 0.0f);
	const Vector3D			vDomainExtent	= vorticityDistribution.GetDomainSize();
	const Vector3D			vMin			= vCenter - 0.5f * vDomainExtent;
	const Vector3D			vMax			= vCenter + 0.5f * vDomainExtent;
	UniformGridGeometry		skeleton		(numVortonsMax, vMin, vMax, true);
	unsigned				numCells[3]		= { skeleton.GetNumCells(0), skeleton.GetNumCells(1), skeleton.GetNumCells(2) };
	//Attention, the numCells initialized by UniformGridGeometry maybe bigger than numVortonsMax ( < uNumElements * 8)
	//So we need to ensure that it's small than numVortonsMax

	while ( numCells[0] * numCells[1] * numCells[2] > numVortonsMax )
	{
		numCells[0] = MAX2(1, numCells[0] / 2);
		numCells[1] = MAX2(1, numCells[1] / 2);
		numCells[2] = MAX2(1, numCells[2] / 2);
	}

	const float				oneOverN[3]		= { 1.0f / float(numCells[0]), 1.0f / float(numCells[1]), 1.0f / float(numCells[2]) };
	const Vector3D			vCellExtent		(vDomainExtent.x * oneOverN[0], vDomainExtent.y * oneOverN[1], vDomainExtent.z * oneOverN[2]);
	const Vector3D			vNoise			(0.0f * vCellExtent); //?? Original version is 0.0f * vCellExtent
	float					vortonRadius	= powf(vCellExtent.x * vCellExtent.y * vCellExtent.z, 1.0f / 3.0f) * 0.5f;
	if (0.0f == vDomainExtent.z)
	{
		vortonRadius = powf(vCellExtent.x * vCellExtent.y, 0.5f) * 0.5f;
	}

	Vector3D	position;
	unsigned	idx[3];
	for (idx[2] = 0; idx[2] < numCells[2]; ++ idx[2] )
	{
		position.z = vMin.z + ( float(idx[2]) + 0.25 ) * vCellExtent.z;
		for (idx[1] = 0; idx[1] < numCells[1]; ++idx[1])
		{
			position.y = vMin.y + ( float(idx[1]) + 0.25 ) * vCellExtent.y;
			for (idx[0] = 0; idx[0] < numCells[0]; ++idx[0])
			{
				position.x	= vMin.x + ( float(idx[0]) + 0.25 ) * vCellExtent.x;
				position	+= RandomSpread(vNoise);

				Vector3D vorticity;
				vorticityDistribution.AssignVorticity(vorticity, position, vCenter);

				Vorton vorton(position, vorticity * fMagnitude, vortonRadius);
				if (vorticity.norm() > sTiny)
				{
					vortons.push_back(vorton);
				}
			}
		}
	}
}