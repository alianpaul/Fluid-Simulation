#include "vorton.h"
/*
	Check the Vorton' member functin:
		AccumulateVelocity and AssignByVelocity works fine.
*/
void Vorton::UnitTest()
{
	fprintf(stderr, "Vorton::UnitTest------------\n");

	static const float		radiusVort			= 0.125f;
	static const unsigned	numPointsPerSide	= 4;
	static const float		gridSize			= float(numPointsPerSide) * 2.0f * radiusVort;

	//for each vorton v at position (xv + r, yv + r, zv + r)
	for (unsigned iz = 0; iz < numPointsPerSide; ++ iz)
	{
		const float zv = gridSize * (float(iz) / float(numPointsPerSide - 1) - 0.5f);

		for (unsigned iy = 0; iy < numPointsPerSide; ++ iy)
		{
			const float yv = gridSize * (float(iy) / float(numPointsPerSide - 1) - 0.5f);

			for (unsigned ix = 0; ix < numPointsPerSide; ++ ix)
			{
				const float xv = gridSize * (float(ix) / float(numPointsPerSide - 1) - 0.5f);

				//for each vorton's  query point at position (xq, yq, zq)
				for (unsigned jz = 0; jz < numPointsPerSide; ++ jz)
				{
					const float zq = gridSize * (float(jz) / float(numPointsPerSide - 1) - 0.5f);

					for (unsigned jy = 0; jy < numPointsPerSide; ++ jy)
					{
						const float yq = gridSize * (float(jy) / float(numPointsPerSide - 1) - 0.5f);

						for (unsigned jx = 0; jx < numPointsPerSide; ++ jx)
						{
							const float xq = gridSize * (float(jx) / float(numPointsPerSide - 1) - 0.5f);

							//for different kind of vorticity (xw, yw, zw)
							for (unsigned kz = 0; kz < numPointsPerSide; ++ kz)
							{
								const float zw = gridSize * (float(kz) / float(numPointsPerSide - 1) - 0.5f);
								
								for (unsigned ky = 0; ky < numPointsPerSide; ++ ky)
								{
									const float yw = gridSize * (float(ky) / float(numPointsPerSide - 1) - 0.5f);

									for (unsigned kx = 0; kx < numPointsPerSide; ++kx)
									{
										const float xw = gridSize * (float(kx) / float(numPointsPerSide - 1) - 0.5f);

										//Vorton at position (xv + r, yv + r, zv + r)
										const Vector3D	positionOfVorton(Vector3D(xv, yv, zv) + Vector3D(1.0f, 1.0f, 1.0f) * radiusVort);
										const Vorton	vort0(positionOfVorton, Vector3D(xw, yw, zw), radiusVort);
										
										//Compute the Velocity at query point(xq, yq, zq)
											  Vector3D	velDueToVort0( 0.0f, 0.0f, 0.0f );
										const Vector3D	posQuery( xq, yq, zq );
										vort0.AccumulateVelocity(velDueToVort0, posQuery);

										//vorticity and velocity should be orthogonal
										_ASSERT( dot( vort0.mVorticity, velDueToVort0 ) - 0.0f <= FLT_EPSILON);
										if (velDueToVort0.norm2() > FLT_MIN)
										{

											_ASSERT(dot((posQuery - vort0.mPosition),velDueToVort0) < FLT_EPSILON);

											//Do the reverse, calculate the vorticity to reach the velocity at 
											//queryPos, from this vorticity, we calculate again the velocity at queryPos
											//And check the velocity consistency
											Vorton vortTest(positionOfVorton, Vector3D(0.0f, 0.0f, 0.0f), radiusVort);
											vortTest.AssignByVelocity( posQuery, velDueToVort0);
											//This vorticity is not the same as teh vort0's vorticity, BUT what we care is that if it can induce
											//the same velocity at queryPos
											Vector3D velDueToVortTest(0.0f, 0.0f, 0.0f);
											vortTest.AccumulateVelocity(velDueToVortTest, posQuery);

											_ASSERT(velDueToVort0.Resembles(velDueToVortTest));

										}

									}
								}
							}
						}
					}
				}
			}
		}
	}


	fprintf(stderr, "Vorton::UnitTest END------------\n");
}