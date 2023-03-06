#ifndef __BVH_TREE__
#define __BVH_TREE__

#pragma once
#include "Constraint.h"
#include "RTriangle.h"

#define REFIT_BLOCKSIZE		512u

using namespace std;

struct TriInfo {
	uint _face;
	uint _id;
	REAL _pos;
};

struct BVHNode
{
	uint _level;
	uint _face;
	uint _path;
	AABB _aabb;
};
struct BVHParam
{
public:
	uint* _levels;
	uint* _faces;
	REAL* _mins[3];
	REAL* _maxs[3];
public:
	uint _maxLevel;
	uint _size;
	uint _pivot;
	uint _numFaces;
};

class BVH
{
public:
	Dvector<uint> _levels;
	Dvector<uint> _faces;
	Dvector<REAL> _mins[3];
	Dvector<REAL> _maxs[3];
	uint _test;
public:
	uint _maxLevel;
	uint _size;
	uint _pivot;
	uint _numFaces;
public:
	BVH() {};
	virtual ~BVH() {}
public:
	void initBVHTreeDevice(uint numFaces) {
		_maxLevel = Log2(numFaces - 1u << 1u);
		_pivot = (1u << _maxLevel) - numFaces;
		_size = (1u << _maxLevel + 1u) - 1u - (_pivot << 1u);
		_pivot = ((1u << (_maxLevel - 1u)) - _pivot) << 1u;
		_numFaces = numFaces;

		_faces.resize(_numFaces);
		_levels.resize(_size);
		_mins[0].resize(_size);
		_mins[1].resize(_size);
		_mins[2].resize(_size);
		_maxs[0].resize(_size);
		_maxs[1].resize(_size);
		_maxs[2].resize(_size);
	}
	inline void free(void) {
		_faces.freeChk();
		_levels.freeChk();
		_mins[0].freeChk();
		_mins[1].freeChk();
		_mins[2].freeChk();
		_maxs[0].freeChk();
		_maxs[1].freeChk();
		_maxs[2].freeChk();
	}
	inline BVHParam param(void) {
		BVHParam p;
		p._levels = _levels._list;
		p._faces = _faces._list;
		p._mins[0] = _mins[0]._list;
		p._mins[1] = _mins[1]._list;
		p._mins[2] = _mins[2]._list;
		p._maxs[0] = _maxs[0]._list;
		p._maxs[1] = _maxs[1]._list;
		p._maxs[2] = _maxs[2]._list;
		p._maxLevel = _maxLevel;
		p._size = _size;
		p._pivot = _pivot;
		p._numFaces = _numFaces;
		return p;
	}
public:
	void build(Dvector<uint>& fs, Dvector<REAL>& ns, const DPrefixArray<uint>& nbFs);
	void refit(ObjParam obj, REAL delta, const REAL dt, bool isCCD);
	void draw(const AABB& aabb);
	void draw(void);
};

inline __device__ void getBVHAABB(AABB& aabb, BVHParam& bvh, uint ind) {
	aabb._min.x = bvh._mins[0][ind];
	aabb._min.y = bvh._mins[1][ind];
	aabb._min.z = bvh._mins[2][ind];
	aabb._max.x = bvh._maxs[0][ind];
	aabb._max.y = bvh._maxs[1][ind];
	aabb._max.z = bvh._maxs[2][ind];
}
inline __device__ void updateBVHAABB(BVHParam& bvh, const AABB& aabb, uint ind) {
	bvh._mins[0][ind] = aabb._min.x;
	bvh._mins[1][ind] = aabb._min.y;
	bvh._mins[2][ind] = aabb._min.z;
	bvh._maxs[0][ind] = aabb._max.x;
	bvh._maxs[1][ind] = aabb._max.y;
	bvh._maxs[2][ind] = aabb._max.z;
}
inline __device__ void RefitBVHLeaf(
	uint fid, AABB& aabb, const ObjParam& params,
	const REAL delta, const REAL dt, const bool isCCD)
{
	fid *= 3u;
	uint ino0 = params._fs[fid + 0u];
	uint ino1 = params._fs[fid + 1u];
	uint ino2 = params._fs[fid + 2u];
	ino0 *= 3u;
	ino1 *= 3u;
	ino2 *= 3u;
	REAL3 p0, p1, p2;
	p0.x = params._ns[ino0 + 0u]; p0.y = params._ns[ino0 + 1u]; p0.z = params._ns[ino0 + 2u];
	p1.x = params._ns[ino1 + 0u]; p1.y = params._ns[ino1 + 1u]; p1.z = params._ns[ino1 + 2u];
	p2.x = params._ns[ino2 + 0u]; p2.y = params._ns[ino2 + 1u]; p2.z = params._ns[ino2 + 2u];
	setAABB(aabb, p0, delta);
	addAABB(aabb, p1, delta);
	addAABB(aabb, p2, delta);
	if (isCCD) {
		REAL3 v0, v1, v2;
		v0.x = params._vs[ino0 + 0u]; v0.y = params._vs[ino0 + 1u]; v0.z = params._vs[ino0 + 2u];
		v1.x = params._vs[ino1 + 0u]; v1.y = params._vs[ino1 + 1u]; v1.z = params._vs[ino1 + 2u];
		v2.x = params._vs[ino2 + 0u]; v2.y = params._vs[ino2 + 1u]; v2.z = params._vs[ino2 + 2u];
		p0 += v0 * dt; p1 += v1 * dt; p2 += v2 * dt;
		addAABB(aabb, p0, delta);
		addAABB(aabb, p1, delta);
		addAABB(aabb, p2, delta);
	}
}
inline __host__ __device__ __forceinline__ int getBVHIndex(uint path, uint level) {
	return (1u << level) - 1u + path;
}
inline __host__ __device__ __forceinline__ int getBVHChild(uint path, uint level) {
	return (1u << level + 1u) - 1u + (path << 1u);
}inline __device__ void setBVHNode(const BVHParam& bvh, const BVHNode& node, uint ind) {
	uint ileaf = bvh._size - bvh._numFaces;
	if (ind >= ileaf) 
		bvh._faces[ind - ileaf] = node._face;
	bvh._levels[ind] = node._level;
	bvh._mins[0][ind] = node._aabb._min.x;
	bvh._mins[1][ind] = node._aabb._min.y;
	bvh._mins[2][ind] = node._aabb._min.z;
	bvh._maxs[0][ind] = node._aabb._max.x;
	bvh._maxs[1][ind] = node._aabb._max.y;
	bvh._maxs[2][ind] = node._aabb._max.z;
}
inline __device__ void getBVHNode(BVHNode& node, const BVHParam& bvh, uint ind) {
	uint ileaf = bvh._size - bvh._numFaces;
	if (ind >= ileaf) 
		node._face = bvh._faces[ind - ileaf];
	node._level = bvh._levels[ind];
	node._aabb._min.x = bvh._mins[0][ind];
	node._aabb._min.y = bvh._mins[1][ind];
	node._aabb._min.z = bvh._mins[2][ind];
	node._aabb._max.x = bvh._maxs[0][ind];
	node._aabb._max.y = bvh._maxs[1][ind];
	node._aabb._max.z = bvh._maxs[2][ind];
}

#endif