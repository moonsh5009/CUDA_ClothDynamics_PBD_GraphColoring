#include "CollisionSolver.h"
#include "SelfCollisionSolver.h"

bool CollisionSolver::ResolveSelfCollision(
	const ObjParam& clothParams,
	const DPrefixArray<uint>& clothEdges,
	BVH* clothBvh, RTriangle* RTri,
	const REAL delta, const REAL dt) 
{
	return SelfCollisionSolver::ResolveSelfCollision(clothParams, clothEdges, clothBvh, RTri, delta, dt);
}