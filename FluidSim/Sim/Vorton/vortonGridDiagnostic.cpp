#include "Space\uniformGrid.h"
#include "vorton.h"

/*
	Overwrite the original function in the tempale
	Min store the minPosition and minVorticity (in 3 dimensions), max vice versa
*/
void UniformGrid<Vorton>::ComputeStatistics( Vorton & min, Vorton & max ) const
{
	const Vector3D	vMax ( FLT_MAX, FLT_MAX, FLT_MAX );
	min = Vorton(vMax, vMax);
	max = Vorton(-vMax, -vMax);
	const unsigned numGridPoints = Size();
	for (unsigned offset = 0; offset < numGridPoints; ++offset)
	{
		const Vorton &		rVorton = (*this)[offset];
		const Vector3D &	rPosition = rVorton.mPosition;
		min.mPosition.x		= MIN2(min.mPosition.x, rPosition.x);
		min.mPosition.y		= MIN2(min.mPosition.y, rPosition.y);
		min.mPosition.z		= MIN2(min.mPosition.z, rPosition.z);
		max.mPosition.x		= MAX2(max.mPosition.x, rPosition.x);
		max.mPosition.y		= MAX2(max.mPosition.y, rPosition.y);
		max.mPosition.z		= MAX2(max.mPosition.z, rPosition.z);
		const Vector3D &	rVorticity = rVorton.mVorticity;
		min.mVorticity.x	= MIN2(min.mVorticity.x, rVorticity.x);
		min.mVorticity.y	= MIN2(min.mVorticity.y, rVorticity.y);
		min.mVorticity.z	= MIN2(min.mVorticity.z, rVorticity.z);
		max.mVorticity.x	= MAX2(max.mVorticity.x, rVorticity.x);
		max.mVorticity.y	= MAX2(max.mVorticity.y, rVorticity.y);
		max.mVorticity.z	= MAX2(max.mVorticity.z, rVorticity.z);
		
	}
}

void UniformGrid<Vorton>::GenerateBrickOfBytes( const char * strFilenameBase, unsigned uFrame ) const
{
	//Prepare the file name
	//strFilenameBase_position strFilenameBase_vorticity
	char strDataFilenames[2][256];
	sprintf_s(strDataFilenames[0], "%s_%d_position.txt", strFilenameBase, uFrame);
	sprintf_s(strDataFilenames[1], "%s_%d_vorticity.txt", strFilenameBase, uFrame);

	FILE * pDataFiles[2];
	fopen_s(&pDataFiles[0], strDataFilenames[0], "w");
	fopen_s(&pDataFiles[1], strDataFilenames[1], "w");

	if (!pDataFiles[0]) { _ASSERT(0); }
	if (!pDataFiles[1]) { _ASSERT(0); }

	//Find the min and max vorticity
	//We can use this range to normalize the vorticity, easy to visualize
	Vorton min, max;
	ComputeStatistics( min, max );

	Vector3D	vPosMin		= min.mPosition;
	Vector3D	vPosMax		= max.mPosition;

	//OutPut the MinPos and MaxPos
	//fprintf( pDataFiles[0], "MIN: %f %f %f \n", vPosMin.x, vPosMin.y, vPosMin.z );
	//fprintf( pDataFiles[0], "MAX: %f %f %f \n", vPosMax.x, vPosMax.y, vPosMax.z);

	const unsigned numGridPoints = Size();
	for (unsigned offset = 0; offset < numGridPoints; ++offset)
	{
		const Vorton &		rVorton			= (*this)[offset];
		const Vector3D &	rVorticity		= rVorton.mVorticity;
		const Vector3D &	rPosition		= rVorton.mPosition;

		fprintf( pDataFiles[0], "%f %f %f \n", rPosition.x, rPosition.y, rPosition.z );
		fprintf( pDataFiles[1], "%f %f %f \n", rVorticity.x, rVorticity.y, rVorticity.z);
	}

	fclose(pDataFiles[0]);
	fclose(pDataFiles[1]);

	return;
}