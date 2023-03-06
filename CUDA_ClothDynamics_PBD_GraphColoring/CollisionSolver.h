#ifndef __COLLISION_SOLVER_H__
#define __COLLISION_SOLVER_H__

#pragma once
#include "BVH.h"

#define COLLISION_TESTTIMER

struct ContactElem {
	bool _isfv;
	uint _i[4];
	REAL _info;
};
struct ContactElemParam {
	ContactElem	*_elems;
	uint		*d_size;
	uint		h_size;
};
class ContactElems {
public:
	Dvector<ContactElem> _elems;
	Dvector<uint>		d_size;
	uint				h_size;
public:
	ContactElems() { d_size.resize(1); }
	~ContactElems() {}
public:
	inline void resize(void) {
		_elems.resize(h_size);
	}
	inline void resize(uint size) {
		_elems.resize(size);
		h_size = size;
	}
	inline ContactElemParam param(void) {
		ContactElemParam p;
		p._elems = _elems._list;
		p.d_size = d_size._list;
		p.h_size = h_size;
		return p;
	}
};
struct ContactElem_CMP
{
	__host__ __device__
		bool operator()(const ContactElem& a, const ContactElem& b) {
		if (a._info != b._info)
			return a._info < b._info;
		if (a._i[0] != b._i[0])
			return a._i[0] < b._i[0];
		if (a._i[1] != b._i[1])
			return a._i[1] < b._i[1];
		if (a._i[2] != b._i[2])
			return a._i[2] < b._i[2];
		return a._i[3] < b._i[3];
	}
};

class CollisionSolver {
public:
	CollisionSolver() {}
	~CollisionSolver() {}
//public:
//	static bool ResolveObstacleCollision(
//		const ObjectParams& clothParams,
//		const ObjectParams& obsParams,
//		BVH* clothBvh,
//		BVH* obsBvh,
//		const double delta, const double friction, const double dt);
public:
	static bool ResolveSelfCollision(
		const ObjParam& clothParams,
		const DPrefixArray<uint>& clothEdges,
		BVH* clothBvh, RTriangle* RTri,
		const REAL delta, const REAL dt);
//	static bool ResolveBoundaryCollision(
//		const ObjectParams& clothParams,
//		const double3& minBoundary, const double3& maxBoundary,
//		const double delta, const double friction, const double dt);
};

#endif