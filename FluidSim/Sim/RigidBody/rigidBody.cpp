
#include "rigidBody.h"

void RigidBody::UpdateSystem( std::vector< RigidBody > & rigidBodies, float timeStep, unsigned uFrame )
{
	for (unsigned uBody = 0; uBody < rigidBodies.size(); ++uBody)
	{
		RigidBody & rbody = rigidBodies[uBody];
		rbody.Update(timeStep);
	}
}