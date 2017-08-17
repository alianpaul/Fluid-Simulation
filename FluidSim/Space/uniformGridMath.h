#ifndef UNIFROM_GRID_MATH
#define UNIFORM_GRID_MATH

#include "Core\Math\matrix3x3.h"
#include "uniformGrid.h"


extern void ComputeJacobian( UniformGrid<Matrix3x3> & jacobian,  const UniformGrid<Vector3D> & vec);

extern void ComputeCurlFromJacobian(UniformGrid<Vector3D> & curl, const UniformGrid<Matrix3x3> & jacobian);

extern void GridMathUnitTest();

#endif