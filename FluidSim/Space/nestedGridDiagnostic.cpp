#include "Space\nestedGrid.h"

/*static*/ void NestedGrid<unsigned>::UnitTest()
{
	static unsigned const num = 1024;
	static const Vector3D vRange(2.0f, 3.0f, 5.0f); // range of random positions
	Vector3D vMin(-0.5f * vRange);   // Minimum coordinate of UniformGrid.
	Vector3D vMax(0.5f * vRange);   // Maximum coordinate of UniformGrid.
	UniformGrid<unsigned> ug(num, vMin, vMax, true);
	NestedGrid<unsigned> ng0(ug);

	/*Check the root layer, has the same geometry with ug*/
	UniformGrid<unsigned> & ngRoot = ng0[0];
	_ASSERT(& ug != & ngRoot);
	_ASSERT(ug.GetMinCorner() == ngRoot.GetMinCorner());
	_ASSERT(ug.GetNumPoints(0) == ngRoot.GetNumPoints(0));
	_ASSERT(ug.GetNumPoints(1) == ngRoot.GetNumPoints(1));
	_ASSERT(ug.GetNumPoints(1) == ngRoot.GetNumPoints(2));
	_ASSERT(ug.GetExtent() == ngRoot.GetExtent());

	for (unsigned iLayer = ng0.GetDepth() - 1 ; iLayer > 0 ; -- iLayer)
	{
		/*Check the decimation of each Layer*/
		unsigned decimations[3];
		ng0.ComputeDecimations(decimations, iLayer);
		_ASSERT(decimations[0] >= 1 && decimations[0] <= 2);
		_ASSERT(decimations[1] >= 1 && decimations[1] <= 2);
		_ASSERT(decimations[2] >= 1 && decimations[2] <= 2);

		/*Check the realationship between decimation and NumOfCells is corrent*/
		UniformGrid<unsigned> & layerParent = ng0[ iLayer    ];
		UniformGrid<unsigned> & layerChild	= ng0[ iLayer - 1];

		_ASSERT(layerChild.GetNumCells(0) == layerParent.GetNumCells(0) * decimations[0]);
		_ASSERT(layerChild.GetNumCells(1) == layerParent.GetNumCells(1) * decimations[1]);
		_ASSERT(layerChild.GetNumCells(2) == layerParent.GetNumCells(2) * decimations[2]);

		/*Check the realtionship between decimation and cell extent*/
		_ASSERT(layerParent.GetCellSpacing().x == layerChild.GetCellSpacing().x * decimations[0]);
		_ASSERT(layerParent.GetCellSpacing().y == layerChild.GetCellSpacing().y * decimations[1]);
		_ASSERT(layerParent.GetCellSpacing().z == layerChild.GetCellSpacing().z * decimations[2]);

		unsigned indexOfParentCell[3];
		
		for (indexOfParentCell[2] = 0; indexOfParentCell[2] < layerParent.GetNumCells(2); ++indexOfParentCell[2])
		{
			//Vector3D vCellMinCornerParent;
			//Vector3D vCellMaxCornerParent;

			for (indexOfParentCell[1] = 0; indexOfParentCell[1] < layerParent.GetNumCells(1); ++indexOfParentCell[1])
			{
				for (indexOfParentCell[0] = 0; indexOfParentCell[0] < layerParent.GetNumCells(0); ++indexOfParentCell[0])
				{
					unsigned indexIntoChild[3];
					ng0.GetChildClusterMinCornerIndex(indexIntoChild, decimations, indexOfParentCell);

					/*Check indexIntoChild is valid*/
					_ASSERT((indexIntoChild[0] >= 0) && (indexIntoChild[0] >= indexOfParentCell[0]) && (indexIntoChild[0] < layerChild.GetNumCells(0)));
					_ASSERT((indexIntoChild[1] >= 0) && (indexIntoChild[1] >= indexOfParentCell[1]) && (indexIntoChild[1] < layerChild.GetNumCells(1)));
					_ASSERT((indexIntoChild[2] >= 0) && (indexIntoChild[2] >= indexOfParentCell[2]) && (indexIntoChild[2] < layerChild.GetNumCells(2)));

					/*Check minCornerCell maxCornerCell is the same with minCorner maxCorner of Parent*/
					Vector3D vPositonMinCornerChild;
					layerChild.PositionFromIndices(vPositonMinCornerChild, indexIntoChild);
					Vector3D vPositionMinCornerParent;
					layerParent.PositionFromIndices(vPositionMinCornerParent, indexOfParentCell);
					_ASSERT(vPositionMinCornerParent.Resembles(vPositonMinCornerChild));

					Vector3D vPositionMaxCornerChild;
					vPositionMaxCornerChild.x = vPositonMinCornerChild.x + float(decimations[0]) * layerChild.GetCellSpacing().x;
					vPositionMaxCornerChild.y = vPositonMinCornerChild.y + float(decimations[1]) * layerChild.GetCellSpacing().y;
					vPositionMaxCornerChild.z = vPositonMinCornerChild.z + float(decimations[2]) * layerChild.GetCellSpacing().z;
					Vector3D vPositionMaxCornerParent = vPositionMinCornerParent + layerParent.GetCellSpacing();
					_ASSERT(vPositionMaxCornerChild.Resembles(vPositionMaxCornerParent));
				}
			}
		}

	}
}