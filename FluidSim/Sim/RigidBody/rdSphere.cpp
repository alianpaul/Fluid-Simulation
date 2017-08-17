#include "rbSphere.h"

void RbSphere::UpdateSystem(std::vector< RbSphere > & rbSpheres, float timeStep, unsigned uFrame)
{
	const size_t numBodies = rbSpheres.size();
	for (unsigned uBody = 0; uBody < numBodies; ++uBody)
	{   // For each body in the simulation...
		RbSphere & rBody = rbSpheres[uBody];

		// Update body physical state
		rBody.Update(timeStep);
	}
}