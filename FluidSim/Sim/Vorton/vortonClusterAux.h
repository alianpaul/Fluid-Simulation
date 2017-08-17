#ifndef VORTON_CLUSTER_AUX_H
#define VORTON_CLUSTER_AUX_H

#include <float.h>

class VortonClusterAux
{
public:
	VortonClusterAux() : mVortNormSum(FLT_MIN) { }

	float	mVortNormSum;
};

#endif