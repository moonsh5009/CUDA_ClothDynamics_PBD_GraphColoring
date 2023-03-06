#ifndef __SELF_COLLISION_SOLVER_H__
#define __SELF_COLLISION_SOLVER_H__

#pragma once
#include "CollisionSolver.h"

//#define BVTT_SHARED_SIZE		1028
#define BVTT_SHARED_SIZE		2048
//#define BVTT_SHARED_SIZE		2560
//#define BVTT_SHARED_SIZE		4096

#define CCD_ITERATION			10000
#define RIGID_ITERATION			100

class SelfCollisionSolver {
public:
	SelfCollisionSolver() {}
	~SelfCollisionSolver() {}
public:
	static void ResolveSelfCollisionProximity(
		ContactElems& ceParam,
		const ObjParam& clothParams,
		BVHParam& clothBvh, RTriParam& RTri,
		Dvector<uint2>& NodeCEs,
		Dvector<uint>& iNodeCEs,
		Dvector<uint>& icurrs,
		Dvector<uint>& isEnds,
		const REAL delta, const REAL dt);
	static void ResolveSelfCollisionCCD(
		ContactElems& ceParam,
		const ObjParam& clothParams,
		BVHParam& clothBvh, RTriParam& RTri,
		Dvector<uint2>& NodeCEs,
		Dvector<uint>& iNodeCEs,
		Dvector<uint>& icurrs,
		Dvector<uint>& isEnds,
		const REAL delta, const REAL dt);
	static uint ResolveRigidImpactZone(
		ContactElemParam& ceParams,
		const ObjParam& clothParams,
		const DPrefixArray<uint>& d_edges,
		BVH* clothBvh,
		const REAL delta, const REAL dt);
public:
	static bool ResolveSelfCollision(
		const ObjParam& clothParams,
		const DPrefixArray<uint>& clothEdges,
		BVH* clothBvh, RTriangle* RTri,
		const REAL delta, const REAL dt);
};

#endif