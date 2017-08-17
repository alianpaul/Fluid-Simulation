
#include "vortonSim.h"
#include "vortonClusterAux.h"
#include "Space\uniformGridMath.h"
#include "Core\perf.h"

#include <future>
#include <thread>

static unsigned gNumberOfProcessors = 8;

using namespace std;

void UpdateBoundingBox(Vector3D & vMinCorner, Vector3D & vMaxCorner, const Vector3D & vPoint)
{
	vMinCorner.x = MIN2(vMinCorner.x, vPoint.x);
	vMinCorner.y = MIN2(vMinCorner.y, vPoint.y);
	vMinCorner.z = MIN2(vMinCorner.z, vPoint.z);
	vMaxCorner.x = MAX2(vMaxCorner.x, vPoint.x);
	vMaxCorner.y = MAX2(vMaxCorner.y, vPoint.y);
	vMaxCorner.z = MAX2(vMaxCorner.z, vPoint.z);
}

//public:
/*
	Initialize:
	mCirculationInitial, mLinearImpulseInitial
	mAverageVorticity
	mInfluenceTree (For InitializePassiveTracers)
	mTracers
	mMassPerParticle
*/
void VortonSim::Initialize( unsigned numTracersPerCellCubeRoot )
{
	ConservedQuantities( mCirculationInitial, mLinearImpulseInitial );
	ComputeAverageVorticity();
	CreateInfluenceTree();
	InitializePassiveTracers(numTracersPerCellCubeRoot);

	//Calculate the mMassPerParticleUse by the mFluidDensity and  the volume size.
	{
		float domainVolume = mInfluenceTree[0].GetExtent().x * mInfluenceTree[0].GetExtent().y * mInfluenceTree[0].GetExtent().z;
		if (0.0f == mInfluenceTree[0].GetExtent().z)
			domainVolume = mInfluenceTree[0].GetExtent().x * mInfluenceTree[0].GetExtent().y;
		const float		totalMass			= domainVolume * mFluidDensity;
		const unsigned	numTracersPerCell	= numTracersPerCellCubeRoot * numTracersPerCellCubeRoot * numTracersPerCellCubeRoot;
		mMassPerParticle					= totalMass / float(numTracersPerCell * mInfluenceTree[0].GetGridCapacity());
	}
}

const Vector3D VortonSim::GetTracerCenterOfMass() const
{
	Vector3D vCoM(0.0f, 0.0f, 0.0f);
	for (unsigned iTracer = 0; iTracer < mTracers.size(); ++iTracer)
	{
		const Particle & rPcl = mTracers[ iTracer ];
		vCoM += rPcl.mPosition;
	}
	vCoM /= mTracers.size();
	return vCoM;
}

void VortonSim::Update(float timeStep, unsigned uFrame)
{
	QUERY_PERFORMANCE_ENTER;
	CreateInfluenceTree();
	QUERY_PERFORMANCE_EXIT(VortonSim_CreateInfluenceTree);


	QUERY_PERFORMANCE_ENTER;
	ComputeVelocityGrid();
	QUERY_PERFORMANCE_EXIT(VortonSim_ComputeVelocityGrid);

	QUERY_PERFORMANCE_ENTER;
	StretchAndTiltVortons(timeStep, uFrame);
	QUERY_PERFORMANCE_EXIT(VortonSim_StretchAndTiltVortons);

	QUERY_PERFORMANCE_ENTER;
	//DiffuseVorticityPSE(timeStep, uFrame);
	DiffuseVorticityGlobally(timeStep, uFrame);
	QUERY_PERFORMANCE_EXIT(VortonSim_DiffuseVorticityPSE);

	QUERY_PERFORMANCE_ENTER;
	AdvectVortons(timeStep);
	QUERY_PERFORMANCE_EXIT(VortonSim_AdvectVortons);

	QUERY_PERFORMANCE_ENTER;
	AdvectTracers(timeStep, uFrame);
	QUERY_PERFORMANCE_EXIT(VortonSim_AdvectTracers);
}

//private:
/*
	Assign each vorton's vorticity from a vorticity grid
	Attention:
		If the vorticity at the point is big enough, we assign a vorton at
		that point.
		This function will loop over the grid point by z, y, x indices.
		We rarely use function OffsetFromIndices and PoistionFromIndices because it's inefficient.
*/
void VortonSim::AssignVortonsFromVorticity(const UniformGrid< Vector3D > & vortGrid)
{
	mVortons.clear();

	const UniformGridGeometry &	ug			  = vortGrid;
	const float					fVortonRaduis = powf( ug.GetCellSpacing().x * ug.GetCellSpacing().y * ug.GetCellSpacing().z, 1.0f / 3.0f) * 0.5f;
	/*UNDO We don't use Nudge*/
	const Vector3D				vMin		  = ug.GetMinCorner();
	const Vector3D				vSpacing      = ug.GetCellSpacing();
	const unsigned				numPoints[3]  = { ug.GetNumPoints(0), ug.GetNumPoints(1), ug.GetNumPoints(2) };
	const unsigned				numXY         = ug.GetNumPoints(0) * ug.GetNumPoints(1);

	unsigned idx[3];
	for (idx[2] = 0; idx[2] < numPoints[2]; ++ idx[2])
	{
		Vector3D vPositionOfGridPoint ;
		vPositionOfGridPoint.z = vMin.z + float( idx[2] ) * vSpacing.z;
		const unsigned offsetZ = idx[2] * numXY;

		for (idx[1] = 0; idx[1] < numPoints[1]; ++ idx[1])
		{
			vPositionOfGridPoint.y  = vMin.y + float( idx[1] ) * vSpacing.y;
			const unsigned offsetYZ = offsetZ + idx[1] * numPoints[0];

			for (idx[0] = 0; idx[0] < numPoints[0]; ++idx[0])
			{
				vPositionOfGridPoint.x   = vMin.x + float(idx[0]) * vSpacing.x;
				const unsigned offsetXYZ = offsetYZ + idx[0];
				const Vector3D & vVorticity = vortGrid[offsetXYZ];
				//Only if the vorticity here is big enough,We create a new vorton
				if (vVorticity.norm2() > FLT_EPSILON)
				{
					Vorton vorton(vPositionOfGridPoint, vVorticity, fVortonRaduis);
					mVortons.push_back(vorton);
				}
			}
		}
	}
}

/*
	Compute the total circulation and linear impulse of all vortons in this simulation
	Attention:
		circulation of a vorton:
			vorticity of the vorton weighted by the volumeElement;
		vCirculation:
			integral of all vorton's circulation;
		vLinearImpulse:
			integral of all vorton's circulation weighed by it's position

*/
void VortonSim::ConservedQuantities(Vector3D & vCirculation, Vector3D & vLinearImpulse)
{
	vCirculation = vLinearImpulse = Vector3D(0.0f, 0.0f, 0.0f);
	for ( unsigned iVorton = 0; iVorton < mVortons.size(); ++ iVorton )
	{
		const Vorton & vorton = mVortons[iVorton];
		const float	volumeElement = POW3(vorton.mRadius) * 8.0f;

		vCirculation   += vorton.mVorticity * volumeElement;

		vLinearImpulse += cross(vorton.mPosition, vorton.mVorticity) * volumeElement;
	}
}

/*
	Find the axis-aligned bounding box for the vortons.
	Initilze mMinCorner,mMaxCorner member variables.
	Attention:
		Use the auxilary function UpdateBoundingBox,
		The result contain all the vortons and particles.

		The main reason we use Influence Tree is because is because we want the 
		velocity gird.
*/


void VortonSim::FindBoundingBox()
{
	mMinCorner = Vector3D(FLT_MAX, FLT_MAX, FLT_MAX);
	mMaxCorner = -mMinCorner;

	//Vortons
	for (unsigned iVorton = 0; iVorton < mVortons.size(); ++iVorton)
	{
		const Vorton & vorton = mVortons[iVorton];
		UpdateBoundingBox( mMinCorner, mMaxCorner, vorton.mPosition );
	}

	//Particles
	for (unsigned iTracer = 0; iTracer < mTracers.size(); ++iTracer)
	{
		const Particle & particle = mTracers[iTracer];
		UpdateBoundingBox(mMinCorner, mMaxCorner, particle.mPosition);
	}

	//prevent round-off errors
	const Vector3D extent( mMaxCorner - mMinCorner );
	const Vector3D nudge(extent * FLT_EPSILON);
	mMinCorner -= nudge;
	mMinCorner += nudge;

}



/*
	Create base layer of vorton influence tree
	Attention:
		Each CELL corresponds on average to a single supervorton.Each cell may contain
		multiple vortons or zero.
		This base grid is not for Eulerian operations (derivatives of vorticity). Because this routine
		associates each vortex with a single corner of the grid cell that contains it.
		Make sure theoretically conserved quantities behave like expected.
		Method assumes the influence tree skeleton has been created.

		We aggregate the vortons in the grid cell into a supervorton.
		The aggregation of position is weigted by the subvorton's vorticity(magnitude).
			Use VortonClusterAux to sum all the vorticity. Use this value to normlize.
			WHY we use VortonClusterAux??
				What VortonClusterAux store is just a simple float. So why don't we just use
				the UniformGrid<float>? The reason is that use this value to nomalize the 
				aggregated position(/VortonClusterAux). So it can not be 0. By the ctor of
				VortonClusterAux, we can achieve that easily.
		The aggregation of vorticity is simple sum.
		The raduis is the same with sub vortons
*/
void VortonSim::MakeBaseVortonGrid()
{
	UniformGrid< VortonClusterAux > ugAux( mInfluenceTree[0] );
	ugAux.Init();

	const unsigned numVortons = mVortons.size();
	for (unsigned iVorton = 0; iVorton < numVortons; ++ iVorton)
	{
		const Vorton &		vorton		=	mVortons[iVorton];
		const unsigned		offset		=	mInfluenceTree[0].OffsetOfPosition( vorton.mPosition );
		const float			vortMag		=	vorton.mVorticity.norm();
		Vorton &			superVorton	=	mInfluenceTree[0][offset];
		VortonClusterAux &	vortonAux	=	ugAux[offset];
		
		superVorton.mPosition	+=	vorton.mPosition * vortMag;
		superVorton.mVorticity	+=	vorton.mVorticity;
		superVorton.mRadius		 =	vorton.mRadius;
		vortonAux.mVortNormSum	+=	vortMag;
	}

	//Normalize to get the weighted aggregation of position.
	unsigned		numPoints[3] = { ugAux.GetNumPoints(0), ugAux.GetNumPoints(1), ugAux.GetNumPoints(2) };
	const unsigned	numXY		 = ugAux.GetNumPoints(0) * ugAux.GetNumPoints(1);
	
	unsigned idx[3];
	for (idx[2] = 0; idx[2] < numPoints[2]; ++idx[2])
	{
		const unsigned offsetZ = idx[2] * numXY;
		for (idx[1] = 0; idx[1] < numPoints[1]; ++idx[1])
		{
			const unsigned offsetYZ = offsetZ + idx[1] * numPoints[0];
			for (idx[0] = 0; idx[0] < numPoints[0]; ++idx[0])
			{
				const unsigned offsetXYZ = offsetYZ + idx[0];
				VortonClusterAux & VortonAux = ugAux[offsetXYZ];
				if (VortonAux.mVortNormSum != FLT_MIN)
				{
					Vorton & superVorton	= mInfluenceTree[0][offsetXYZ];
					superVorton.mPosition	/= ugAux[offsetXYZ].mVortNormSum;
				}
			}
		}
	}

}

/*
	Aggregate vorton clusters from a child layer into parent layer.
	Attention:
		The argument is the destination parent layer.
		The aggregation approach is the same with MakeBaseVortonGrid
		Loop through the Cells not the Points!!!!
*/
void VortonSim::AggregateClusters(unsigned uParentLayer)
{
	UniformGrid< Vorton > & rParentLayer = mInfluenceTree[ uParentLayer ];
	const unsigned *		pClusterDims = mInfluenceTree.GetDecimations( uParentLayer );

	const unsigned numCells[3]	= {rParentLayer.GetNumCells(0), rParentLayer.GetNumCells(1), rParentLayer.GetNumCells(2)};
	const unsigned numXY		= rParentLayer.GetNumPoints(0) * rParentLayer.GetNumPoints(1);
	unsigned idxParent[3];
	for (idxParent[2] = 0; idxParent[2] < numCells[2]; ++idxParent[2])
	{
		const unsigned offsetZ = idxParent[2] * numXY;
		for (idxParent[1] = 0; idxParent[1] < numCells[1]; ++idxParent[1])
		{
			const unsigned offsetYZ = offsetZ + idxParent[1] * rParentLayer.GetNumPoints(0);
			for (idxParent[0] = 0; idxParent[0] < numCells[0]; ++idxParent[0])
			{
				const unsigned	offsetXYZ		= offsetYZ + idxParent[0];
				Vorton &		rVortonParent	= rParentLayer[offsetXYZ];
				
				//Get the Vortons in the child layer
				unsigned	clusterMinIndices[3];
				mInfluenceTree.GetChildClusterMinCornerIndex( clusterMinIndices, pClusterDims, idxParent );

				UniformGrid< Vorton > &		rChildLayer = mInfluenceTree[ uParentLayer - 1 ];
				const unsigned				numXChild	= rChildLayer.GetNumPoints( 0 );
				const unsigned				numXYChild	= rChildLayer.GetNumPoints( 0 ) * rChildLayer.GetNumPoints( 1 );
				VortonClusterAux vortAux;

				//Loop through each child
				unsigned increment[3];
				for ( increment[2] = 0; increment[2] < pClusterDims[2]; ++ increment[2] )
				{
					const unsigned offsetZChild = ( clusterMinIndices[2] + increment[2] ) * numXYChild;
					for ( increment[1] = 0; increment[1] < pClusterDims[1]; ++increment[1] )
					{
						const unsigned offsetYZChild = offsetZChild + ( clusterMinIndices[1] + increment[1] ) * numXChild;
						for (increment[0] = 0; increment[0] < pClusterDims[0]; ++increment[0])
						{
							const unsigned	offsetXYZChild	= offsetYZChild + ( clusterMinIndices[0] + increment[0] );
							Vorton &		rVortonChild	= rChildLayer[offsetXYZChild];
							const float		vortMag			= rVortonChild.mVorticity.norm();

							rVortonParent.mVorticity	+= rVortonChild.mVorticity;
							rVortonParent.mPosition		+= rVortonChild.mPosition * vortMag;
							if (rVortonChild.mRadius != 0.0f)
							{
								rVortonParent.mRadius = rVortonChild.mRadius;
							}
							vortAux.mVortNormSum += vortMag;

						}
					}
				}

				//Normalize the Vorton parent position
				rVortonParent.mPosition /= vortAux.mVortNormSum;
			}
		}
	}
}

/*
	Create nested grid vorticity influence tree.
*/
void	VortonSim::CreateInfluenceTree()
{
	//Initialize the InfluenceTree skeleton
	FindBoundingBox();
	const unsigned numVortons = mVortons.size();
	UniformGrid< Vorton > ugSkeleton(numVortons, mMinCorner, mMaxCorner, true);

	mInfluenceTree.Initialize( ugSkeleton );

	//Initialize the base layer vorton grid
	MakeBaseVortonGrid();

	//Initialize all the parent layer vorton grid
	const unsigned numLayers = mInfluenceTree.GetDepth();
	for (unsigned uParentLayer = 1; uParentLayer < numLayers; ++uParentLayer)
	{
		AggregateClusters(uParentLayer);
	}
}

/*
	Compute the velocity at a given point due to the Vorton 
	at position idxParent[3] in layer iLayer.
	Attention:
		Once we have the GridCell idxParent[3] at layer iLayer,
		We loop through all it's child cells in the child layer,
		And check if the query position in the child cell. 
		If in, recursive call(also means iLayer >= 1)
		Not in, calculate this child cells vorton influence.

	How to improve accuracy?
		Only if the query position is in the Vorton's Cell, We use the cell' child to do Calculate.
		The consequency of this method maybe that we involved too few vortons.When the query position is close
		to the cell, we also use the cell aggregated vorton,Which induce a lot of error.
		In ensuring complexity of calculating the velocity with nested grid, How to improve accuracy?
		INVOLVE MORE VORTONS IN THE CALCULATION.
		When do the checking of query position, we can expand the cell extent a little bit.
		In this way, When the query position is near the cell(not in), we will also use this cell's children
		not the cell itself.

*/
Vector3D	VortonSim::ComputeVelocity(const Vector3D & vPosition, const unsigned idxParent[3], unsigned iLayer)
{
	UniformGrid< Vorton > &		rChildLayer		= mInfluenceTree[iLayer - 1];
	const unsigned *			pClusterDims	= mInfluenceTree.GetDecimations( iLayer );
	const unsigned				numXChild		= rChildLayer.GetNumPoints(0);
	const unsigned				numXYChild		= numXChild * rChildLayer.GetNumPoints(1);
	const Vector3D &			vGridMinCorner	= rChildLayer.GetMinCorner();
	const Vector3D &			vSpacing		= rChildLayer.GetCellSpacing();
	unsigned					clusterMinIndices[3];
	mInfluenceTree.GetChildClusterMinCornerIndex(clusterMinIndices, pClusterDims, idxParent);
		  Vector3D				velocityAccumulator(0.0f, 0.0f, 0.0f);

	static const float	marginFactor	= 0.0001f; //reasonable value is [0.0001f, 4.0],The bigger, the more accurate,but more complex.
	const Vector3D		margin			= marginFactor * vSpacing;


	unsigned increment[3];
	for ( increment[2] = 0; increment[2] < pClusterDims[2]; ++ increment[2] )
	{
		unsigned idxChild[3];
		Vector3D vCellMinCorner, vCellMaxCorner;

		idxChild[2] = clusterMinIndices[2] + increment[2];
		unsigned offsetZ = idxChild[2] * numXYChild;
		vCellMinCorner.z = vGridMinCorner.z + float( idxChild[2] ) * vSpacing.z;
		vCellMaxCorner.z = vGridMinCorner.z + float( idxChild[2] + 1 ) * vSpacing.z;

		for ( increment[1] = 0; increment[1] < pClusterDims[1]; ++ increment[1] )
		{

			idxChild[1] = clusterMinIndices[1] + increment[1];
			unsigned offsetYZ = offsetZ + idxChild[1] * numXChild;
			vCellMinCorner.y = vGridMinCorner.y + float( idxChild[1] ) * vSpacing.y;
			vCellMaxCorner.y = vGridMinCorner.y + float( idxChild[1] + 1) * vSpacing.y;

			for ( increment[0] = 0; increment[0] < pClusterDims[0]; ++ increment[0] )
			{

				idxChild[0] = clusterMinIndices[0] + increment[0];
				unsigned offsetXYZ = offsetYZ + idxChild[0];
				vCellMinCorner.x = vGridMinCorner.x + float( idxChild[0] ) * vSpacing.x;
				vCellMaxCorner.x = vGridMinCorner.x + float( idxChild[0] + 1 ) * vSpacing.x;

				//Not the BaseLayer Cell and query position is in the cell
				if (   ( iLayer > 1 )
					&& ( vPosition.x >= vCellMinCorner.x - margin.x)
					&& ( vPosition.y >= vCellMinCorner.y - margin.y)
					&& ( vPosition.z >= vCellMinCorner.z - margin.z)
					&& ( vPosition.x <  vCellMaxCorner.x + margin.x)
					&& ( vPosition.y <  vCellMaxCorner.y + margin.y)
					&& ( vPosition.z <  vCellMaxCorner.z + margin.z)
					)
				{
					velocityAccumulator += ComputeVelocity( vPosition, idxChild, iLayer - 1 );
				}
				else
				{
					const Vorton & rVortonChild = rChildLayer[ offsetXYZ ];
					//rVortonChild.AccumulateVelocity( velocityAccumulator, vPosition );
					VORTON_ACCUMULATE_VELOCITY(velocityAccumulator, vPosition, rVortonChild);
				}

			}
		}
	}

	return velocityAccumulator;
}

/*
	Compute velocity at vPosition
	Use it to compare the previous methos.
*/
Vector3D	VortonSim::ComputeVelocityBruteForce(const Vector3D & vPosition)
{
	Vector3D	velocityAccumulator(0.0f, 0.0f, 0.0f);

	for (unsigned iVorton = 0; iVorton < mVortons.size(); ++iVorton)
	{
		const Vorton & rVorton = mVortons[ iVorton ];
		rVorton.AccumulateVelocity( velocityAccumulator, vPosition );
	}

	return velocityAccumulator;
}

/*
	Compute velocity due to the vortons, for a subset of points in a uniform grid.
	A horizontal slice marked by izStart and izEnd.
	Attention:
		This method assumes that the mInfluenceTree is initialized and mVelGrid is Created.
		There two ways of computing velocity:
			From the mInfluenceTree.(ComputeVelocity)
			From the orignal vortons in mVortons.(ComputeVelocityBruteForce) For compare the accuracy.
		We use Macros to switch between these two modes.
*/
void    VortonSim::ComputeVelocityGridSlice(size_t izStart, size_t izEnd)
{
	const unsigned		dims[3]		= { mVelGrid.GetNumPoints(0), mVelGrid.GetNumPoints(1), mVelGrid.GetNumPoints(2) };
	const unsigned		numXY		= dims[0] * dims[1];
	const Vector3D		vMinCorner	= mVelGrid.GetMinCorner();
	const float			nudge		= 1.0f - 2.0f * FLT_EPSILON;
	const Vector3D		vSpacing	= mVelGrid.GetCellSpacing() * nudge;

	unsigned idx[3];
	for ( idx[2] = izStart; idx[2] < izEnd; ++ idx[2] )
	{
		Vector3D vPosition;
		vPosition.z = vMinCorner.z + float(idx[2]) * vSpacing.z;

		const unsigned offsetZ = idx[2] * numXY;
		for ( idx[1] = 0; idx[1] < dims[1]; ++ idx[1] )
		{
			vPosition.y = vMinCorner.y + float(idx[1]) * vSpacing.y;

			const unsigned offsetYZ = offsetZ + idx[1] * dims[0];
			for ( idx[0] = 0; idx[0] < dims[0]; ++ idx[0] )
			{
				vPosition.x = vMinCorner.x + float(idx[0]) * vSpacing.x;
				const unsigned offsetXYZ = offsetYZ + idx[0];

				//From the mInfluenceTree.
				unsigned iLayer			= mInfluenceTree.GetDepth();
				unsigned idxParent[3]	= { 0, 0, 0 };
				mVelGrid[offsetXYZ] = ComputeVelocity(vPosition, idxParent, iLayer - 1);

				//mVelGrid[offsetXYZ] = ComputeVelocityBruteForce(vPosition);
			}
		}
	}

}

/*
	Compute velocity due to vortons.
	Attention:
		For the initial version, we don't use the multi-thread.
*/
void	VortonSim::ComputeVelocityGrid()
{
	mVelGrid.Clear();
	mVelGrid.CopyShape(mInfluenceTree[0]);
	mVelGrid.Init();

	const unsigned numZ = mVelGrid.GetNumPoints( 2 );

	//THREAD
	thread* threadsWorker[THREAD_NUM];
	const unsigned uStride = numZ / THREAD_NUM;
	for (size_t iThread = 0; iThread < THREAD_NUM; iThread++)
	{
		unsigned iStart = iThread * uStride;
		unsigned iEnd	= (iThread == THREAD_NUM - 1) ? numZ : iStart + uStride;
		threadsWorker[iThread] = new thread(&VortonSim::ComputeVelocityGridSlice, this, iStart, iEnd);
	}

	for (unsigned iThread = 0; iThread < THREAD_NUM; ++iThread)
	{
		threadsWorker[iThread]->join();
		delete threadsWorker[iThread];
	}
}

/*
	Stretch and tilt vortons using velocity field.
	-timeStep: amount of time to advance simulation
	-uFrame: frame counter
	Attetion:
		stretchTilt: vorton's vorticity transformed by velocity's gradient(Jacobian).
		vorton's velocity gradient is interpolated with grid.
*/

void    VortonSim::StretchAndTiltVortons(const float & timeStep, const unsigned & uFrame)
{
	//if Domain is 2D, stretching & tilting does not occur.
	if (	(mVelGrid.GetExtent().x == 0.0f) 
		||	(mVelGrid.GetExtent().y == 0.0f)
		||	(mVelGrid.GetExtent().z == 0.0f))
	{
		return;
	}

	UniformGrid< Matrix3x3 > velocityJacobianGrid( mVelGrid );
	velocityJacobianGrid.Init();
	ComputeJacobian( velocityJacobianGrid, mVelGrid );

	for ( unsigned iVorton = 0; iVorton < mVortons.size(); ++iVorton )
	{
		Vorton &	rVorton = mVortons[iVorton];
		Matrix3x3	velJac;
		velocityJacobianGrid.Interpolate( velJac, rVorton.mPosition );

		/*
			stretchTilt = rVorton.mVorticity * velJac ; //original
			??? WHY should calculate in this way. ????
		*/
		const Vector3D	vort = rVorton.mVorticity;
		const Vector3D	stretchTilt( dot( velJac[0], vort ), dot( velJac[1], vort ), dot( velJac[2], vort) );
		//const Vector3D	stretchTilt(velJac * vort);

		//Update the vorticity by stretchTilt
		//0.5f fudge factor for stability
		rVorton.mVorticity += 0.5f * stretchTilt * timeStep;
	}
	
}

/*
	Compute the average vorticity of all vortons in this simulation.
	Attention:
		This is used to compute a hacky, non-physical approximation to
        viscous vortex diffusion.
*/
void    VortonSim::ComputeAverageVorticity(void)
{
	// Zero accumulators.
	mAverageVorticity = Vector3D(0.0f, 0.0f, 0.0f);
	const size_t numVortons = mVortons.size();
	for (unsigned iVorton = 0; iVorton < numVortons; ++iVorton)
	{   // For each vorton in this simulation...
		const Vorton &  rVorton = mVortons[iVorton];
		mAverageVorticity += rVorton.mVorticity;
	}
	mAverageVorticity /= float(numVortons);
}

/*
	Diffuse vorticity globally.
	Approximation of viscous diffusion.
	Attention:
		Bring this vorton's vorticity closer to the average.
		By this way, we can effectively exchange vorticity between vortons.
*/
void    VortonSim::DiffuseVorticityGlobally(const float & timeStep, const unsigned & uFrame)
{
	const Vector3D	vAvgVorticity = mAverageVorticity ;
	mAverageVorticity = Vector3D(0.0f, 0.0f, 0.0f); //Store the new calculated average vorticity.
	const unsigned numVortons = mVortons.size();

	for (unsigned iVorton = 0; iVorton < numVortons; ++ iVorton)
	{
		Vorton &	rVorton = mVortons[ iVorton ];
		Vector3D &	rVorticitySelf = rVorton.mVorticity;

		mAverageVorticity += rVorticitySelf;

		const Vector3D vortDiff		= rVorticitySelf - vAvgVorticity;
		const Vector3D exchange		= mViscosity * timeStep * vortDiff;
		rVorticitySelf -= exchange;
	}

	mAverageVorticity /= float( numVortons );
}

/*
	Diffuse vorticity using a particle strength exchange method.
	Attention:
		We use the mInfluenceTree[0] resolution to partition vortons.
		The cell of this structure stores the Vortons offset in this cell. 
		Exchange Steps:
		for vortons in this cell
		1. Exchange with vortons in the same cell (X0Y0Z0)
		2. Exchange with vortons in the +z cell (X0Y0ZP)
		3. Exchange with vortons in the +y cell (X0YPZ0)
		4. Exchange with vortons in the +x cell (XPY0Z0)

*/
void    VortonSim::DiffuseVorticityPSE(const float & timeStep, const unsigned & uFrame)
{
	UniformGrid< std::vector< unsigned > > ugVortRef( mInfluenceTree[0] );
	ugVortRef.Init() ;

	for ( unsigned iVorton = 0; iVorton < mVortons.size(); ++ iVorton )
	{
		Vorton & rVorton	= mVortons[ iVorton ];
		ugVortRef[ rVorton.mPosition ].push_back( iVorton );
	}

	const unsigned numCells[3]	= { ugVortRef.GetNumCells(0), ugVortRef.GetNumCells(1), ugVortRef.GetNumCells(2) };
	const unsigned numXY		= ugVortRef.GetNumPoints(0) * ugVortRef.GetNumPoints(1);
	const unsigned numX			= ugVortRef.GetNumPoints(0);

	unsigned idx[3];
	for (idx[2] = 0; idx[2] < numCells[2]; ++idx[2])
	{
		const unsigned offsetZ0 =	idx[2] * numXY;
		const unsigned offsetZP = ( idx[2] + 1 )* numXY;

		for (idx[1] = 0; idx[1] < numCells[1]; ++idx[1])
		{
			const unsigned offsetY0Z0 = offsetZ0 +  idx[1] * numX;
			const unsigned offsetY0ZP = offsetZP +  idx[1] * numX;
			const unsigned offsetYPZ0 = offsetZ0 + (idx[1] + 1) * numX;

			for (idx[0] = 0; idx[0] < numCells[0]; ++idx[0])
			{
				const unsigned offsetX0Y0Z0 = offsetY0Z0 + idx[0];
				//For each vorton in this cell.
				for ( unsigned ivHere = 0; ivHere < ugVortRef[offsetX0Y0Z0].size(); ++ ivHere )
				{
					
					const unsigned &	rVortOffsetHere	= ugVortRef[offsetX0Y0Z0][ivHere];
					Vorton &			rVortonHere		= mVortons[rVortOffsetHere];
					Vector3D &			rVorticityHere	= rVortonHere.mVorticity;

					//Exchange vorticity with vortons in the same cell.(not exchanged before).
					for ( unsigned ivThere = ivHere + 1; ivThere < ugVortRef[offsetX0Y0Z0].size(); ++ivThere )
					{
						const unsigned &	rVortOffsetThere	= ugVortRef[offsetX0Y0Z0][ivThere];
						Vorton &			rVortonThere		= mVortons[rVortOffsetThere];
						Vector3D &			rVorticityThere		= rVortonThere.mVorticity;

						const Vector3D vortDiff = rVorticityHere - rVorticityThere;
						const Vector3D exchange = 2.0f * mViscosity * timeStep * vortDiff;
						rVorticityHere	-= exchange;
						rVorticityThere += exchange;
					}

					//Exchange vorticity with vortons in XPY0Z0.(all vortons in XPY0Z0)
					{
						const unsigned offsetXPY0Z0 = offsetY0Z0 + idx[0] + 1;
						for (unsigned ivThere = 0; ivThere < ugVortRef[offsetXPY0Z0].size(); ++ ivThere)
						{
							const unsigned &	rVortOffsetThere	= ugVortRef[offsetXPY0Z0][ivThere];
							Vorton &			rVortonThere		= mVortons[rVortOffsetThere];
							Vector3D &			rVorticityThere		= rVortonThere.mVorticity;

							const Vector3D vortDiff = rVorticityHere - rVorticityThere;
							const Vector3D exchange = mViscosity * timeStep * vortDiff;
							rVorticityHere -= exchange;
							rVorticityThere += exchange;
						}
					}

					//Exchange vorticity with vortons in X0YPZ0.(all vortons in X0YPZ0)
					{
						const unsigned offsetX0YPZ0 = offsetYPZ0 + idx[0];
						for (unsigned ivThere = 0; ivThere < ugVortRef[offsetX0YPZ0].size(); ++ivThere)
						{
							const unsigned &	rVortOffsetThere = ugVortRef[offsetX0YPZ0][ivThere];
							Vorton &			rVortonThere = mVortons[rVortOffsetThere];
							Vector3D &			rVorticityThere = rVortonThere.mVorticity;

							const Vector3D vortDiff = rVorticityHere - rVorticityThere;
							const Vector3D exchange = mViscosity * timeStep * vortDiff;
							rVorticityHere -= exchange;
							rVorticityThere += exchange;
						}
					}

					//Exchange vorticity with vortons in X0Y0ZP
					{
						const unsigned offsetX0Y0ZP = offsetY0ZP + idx[0];
						for (unsigned ivThere = 0; ivThere < ugVortRef[offsetX0Y0ZP].size(); ++ivThere)
						{
							const unsigned &	rVortOffsetThere = ugVortRef[offsetX0Y0ZP][ivThere];
							Vorton &			rVortonThere = mVortons[rVortOffsetThere];
							Vector3D &			rVorticityThere = rVortonThere.mVorticity;

							const Vector3D vortDiff = rVorticityHere - rVorticityThere;
							const Vector3D exchange = mViscosity * timeStep * vortDiff;
							rVorticityHere -= exchange;
							rVorticityThere += exchange;
						}
					}

					//Dissipate vorticity, reduces the vorticity of each vorton.
					//Vorticity dissipates analogously to energy.
					rVorticityHere -= mViscosity * timeStep * rVorticityHere;
				}

			}
		}
	}

}

/*
	Advect vortons using velocity field mVelGrid (calclated by ComputeVelocityGrid)
*/

void    VortonSim::AdvectVortons(const float & timeStep)
{
	for (unsigned iVorton = 0; iVorton < mVortons.size(); ++ iVorton)
	{
		Vorton &	rVorton		= mVortons[ iVorton ];
		Vector3D	velocity;
		mVelGrid.Interpolate(velocity, rVorton.mPosition);
		rVorton.mPosition	+= velocity * timeStep;
		rVorton.mVelocity	= velocity;
	}
}

/*
	Initialize passive tracers
	Attetion:
		The tracer(mTracer) is initialized according to the mInfluenceTree Base Grid.
		For interior each cells:
			we generate multiplier*multiplier*multiplier particles and store
			into mTracer. 
*/

void    VortonSim::InitializePassiveTracers(unsigned multiplier)
{
	//For interior cells in the mInfluenceTree base grid.
	const unsigned begin[3]		= { 1 * mInfluenceTree[0].GetNumCells(0) / 8,
									1 * mInfluenceTree[0].GetNumCells(1) / 8,
									1 * mInfluenceTree[0].GetNumCells(2) / 8 };
	const unsigned end[3]		= { 7 * mInfluenceTree[0].GetNumCells(0) / 8,
									7 * mInfluenceTree[0].GetNumCells(1) / 8,
									7 * mInfluenceTree[0].GetNumCells(2) / 8 };
	const Vector3D &	rMinCorner	= mInfluenceTree[0].GetMinCorner();
	const Vector3D &	rSpacing	= mInfluenceTree[0].GetCellSpacing();
	const Vector3D		vNoise		= rSpacing / float(multiplier);
	const float			pclSize = 2.0f * powf(rSpacing.x * rSpacing.y * rSpacing.z, 2.0f / 3.0f) / float(multiplier);

	const unsigned nt[3] = { multiplier, multiplier, multiplier };
	
	unsigned idx[3];
	for (idx[2] = begin[2]; idx[2] < end[2]; ++idx[2])
	{
		

		for (idx[1] = begin[1]; idx[1] < end[1]; ++idx[1])
		{
			

			for (idx[0] = begin[0]; idx[0] < end[0]; ++idx[0])
			{
				Vector3D vPosMinCorner;
				mInfluenceTree[0].PositionFromIndices(vPosMinCorner, idx);

				Particle pcl;
				pcl.mVelocity			= Vector3D( 0.0f, 0.0f, 0.0f );
				pcl.mOrientation		= Vector3D( 0.0f, 0.0f, 0.0f );
				pcl.mAngularVelocity	= Vector3D( 0.0f, 0.0f, 0.0f );
				pcl.mMass = 1.0f;
				pcl.mSize = pclSize;
				pcl.mBirthTime = 0;

				unsigned it[3];
				for (it[2] = 0; it[2] < nt[2]; ++it[2])
				for (it[1] = 0; it[1] < nt[1]; ++it[1])
				for (it[0] = 0; it[0] < nt[0]; ++it[0])
				{
					Vector3D vShift(float(it[0]) * rSpacing.x / float(nt[0]), 
									float(it[1]) * rSpacing.y / float(nt[1]),
									float(it[2]) * rSpacing.z / float(nt[2]));
					pcl.mPosition = vPosMinCorner + vShift + RandomSpread( vNoise );
					mTracers.push_back( pcl );
				}

			}
		}
	}


}

void	VortonSim::AdvectTracersSlice(const float & timeStep, const unsigned & uFrame, unsigned itStart, unsigned itEnd)
{
	for (unsigned iParticle = itStart; iParticle < itEnd; ++ iParticle)
	{
		Particle &	rTracer = mTracers[ iParticle ];
		Vector3D	velocity;
		mVelGrid.Interpolate( velocity, rTracer.mPosition );
		rTracer.mPosition	+= velocity * timeStep;
		rTracer.mVelocity	=	velocity;
	}
}
/*
	Advect passive tracers
	Attention:
		Need to be multi-thread.
*/

/*! \brief Function object to advect passive tracer particles using Threading Building Blocks
*/
class VortonSim_AdvectTracers_TBB
{
	VortonSim * mVortonSim;    ///< Address of VortonSim object
	const float & mTimeStep;
	const unsigned & mFrame;
public:
	void operator() (const tbb::blocked_range<size_t> & r) const
	{   // Advect subset of tracers.
		mVortonSim->AdvectTracersSlice(mTimeStep, mFrame, r.begin(), r.end());
	}
	VortonSim_AdvectTracers_TBB(VortonSim * pVortonSim, const float & timeStep, const unsigned & uFrame)
		: mVortonSim(pVortonSim)
		, mTimeStep(timeStep)
		, mFrame(uFrame)
	{}
};

void    VortonSim::AdvectTracers(const float & timeStep, const unsigned & uFrame)
{
	//THREAD
	/*
	thread* threadsWorker[THREAD_NUM];
	const unsigned uStride = mTracers.size() / THREAD_NUM;
	for (unsigned iThread = 0; iThread < THREAD_NUM; ++iThread)
	{
		unsigned	iStart	= iThread * uStride;
		unsigned	iEnd	= (iThread == THREAD_NUM - 1) ? mTracers.size() : iStart + uStride;
		threadsWorker[iThread] = new thread(&VortonSim::AdvectTracersSlice, this, timeStep, uFrame, iStart, iEnd);
	}
	for (unsigned iThread = 0; iThread < THREAD_NUM; ++iThread)
	{
		threadsWorker[iThread]->join();
		delete threadsWorker[iThread];
	}
	*/
	const unsigned uStride = mTracers.size() / gNumberOfProcessors;
	tbb::parallel_for(tbb::blocked_range<size_t>(0, mTracers.size(), uStride), VortonSim_AdvectTracers_TBB(this, timeStep, uFrame));
}