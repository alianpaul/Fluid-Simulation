#ifndef NESTED_GRID
#define NESTED_GRID

#include <GL/glew.h>
#include <GL/glut.h>

#include "Core\Math\vector3D.h"
#include "Core\Utility.h"
#include "uniformGrid.h"
#include <math.h>
#include <float.h>
#include <vector>

template<class ItemT> class NestedGrid
{
public:
	typedef UniformGrid< ItemT > Layer;

	NestedGrid()
		:mDecimations(0)
	{

	}

	NestedGrid(const Layer & src)
		:mDecimations(0)
	{
		Initialize(src);
	}

	~NestedGrid()
	{
		delete[] mDecimations;
	}

	void Initialize(const Layer & src)
	{
		mLayers.clear();
		unsigned numLayers = PrecomputeNumLayers(src);
		mLayers.reserve(numLayers);
		AddLayer(src, 1);
		unsigned index = 0;
		while (mLayers[index].GetGridCapacity() > 8)
		{
			AddLayer(mLayers[index], 2);
			++index;
		}

		PrecomputeDecimations();
	}

	/*
		在mLayers 中新加入一个UniformGrid，该Layer为 layerTemplate 以 iDecimate 为参数做Decimate
		并且初始化好UniformGrid 中的container(Init())
	*/

	void AddLayer(const UniformGridGeometry & layerTemplate, unsigned iDecimation)
	{
		mLayers.push_back(Layer());
		mLayers.back().Decimate(layerTemplate, iDecimation);
		mLayers.back().Init();
	}

	unsigned GetDepth() const { return mLayers.size(); }

		  Layer & operator[](unsigned index)	   { return mLayers[index]; }
	const Layer & operator[](unsigned index) const { return mLayers[index]; }
	      
	const unsigned * GetDecimations(unsigned index){ return mDecimations[index]; }

	void GetChildClusterMinCornerIndex(unsigned clusterMinIndices[3], const unsigned decimations[3], const unsigned indicesOfParentCell[3])
	{
		clusterMinIndices[0] = decimations[0] * indicesOfParentCell[0];
		clusterMinIndices[1] = decimations[1] * indicesOfParentCell[1];
		clusterMinIndices[2] = decimations[2] * indicesOfParentCell[2];
	}

	void Clear()
	{
		for (unsigned index = 0; index < GetDepth(); ++index)
		{
			mLayers[index].Clear();
		}

		mLayers.clear();
	}

	static void UnitTest();

private:
	unsigned PrecomputeNumLayers(const Layer & src)
	{
		unsigned numLayers = 1;
		unsigned numPoints[3] = { src.GetNumPoints(0), src.GetNumPoints(1), src.GetNumPoints(2) };
		unsigned size = numPoints[0] * numPoints[1] * numPoints[2];
		while (size > 8)
		{
			numPoints[0] = MAX2(2, (numPoints[0] - 1) / 2 + 1);
			numPoints[1] = MAX2(2, (numPoints[1] - 1) / 2 + 1);
			numPoints[2] = MAX2(2, (numPoints[2] - 1) / 2 + 1);
			size = numPoints[0] * numPoints[1] * numPoints[2];
			++numLayers;
		}

		return numLayers;
	}

	/*
	计算当前ParentLayer 比其前的 ChildLayer，在各位维度上的真实的Decimation
	*/
	void ComputeDecimations( unsigned decimation[3], unsigned iParentLayer)
	{
		const Layer & parent = mLayers[iParentLayer    ];
		const Layer & child  = mLayers[iParentLayer - 1];
		decimation[0] = child.GetNumCells(0) / parent.GetNumCells(0);
		decimation[1] = child.GetNumCells(1) / parent.GetNumCells(1);
		decimation[2] = child.GetNumCells(2) / parent.GetNumCells(2);
	}

	void PrecomputeDecimations(void)
	{
		unsigned numLayers = GetDepth();
		mDecimations = new unsigned[numLayers][3];
		for (unsigned index = 1; index < numLayers; ++index)
		{
			ComputeDecimations( mDecimations[index], index );
		}

		//src is the bottomest child, which has no Decimation
		mDecimations[0][0] = mDecimations[0][1] = mDecimations[0][2];
	}


	std::vector<Layer> mLayers;
	unsigned (*mDecimations)[3];
};

#endif