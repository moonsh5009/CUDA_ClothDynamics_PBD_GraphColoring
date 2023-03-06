#include "BVH.h"
#include "DeviceManager.cuh"

//-------------------------------------------------------------------
__global__ void initBVHInfo_kernel(
	uint* fs, REAL* ns, TriInfo* infos, BVHParam bvh)
{
	uint id = blockDim.x * blockIdx.x + threadIdx.x;
	if (id >= bvh._numFaces)
		return;

	if (id == 0u) {
		bvh._mins[0][0] = DBL_MAX;
		bvh._mins[1][0] = DBL_MAX;
		bvh._mins[2][0] = DBL_MAX;
		bvh._maxs[0][0] = -DBL_MAX;
		bvh._maxs[1][0] = -DBL_MAX;
		bvh._maxs[2][0] = -DBL_MAX;
		bvh._levels[0] = 0u;
	}

	uint ino = id * 3u;
	uint ino0 = fs[ino + 0u];
	uint ino1 = fs[ino + 1u];
	uint ino2 = fs[ino + 2u];
	REAL3 p0, p1, p2;
	ino0 *= 3u; ino1 *= 3u; ino2 *= 3u;
	p0.x = ns[ino0 + 0u]; p0.y = ns[ino0 + 1u]; p0.z = ns[ino0 + 2u];
	p1.x = ns[ino1 + 0u]; p1.y = ns[ino1 + 1u]; p1.z = ns[ino1 + 2u];
	p2.x = ns[ino2 + 0u]; p2.y = ns[ino2 + 1u]; p2.z = ns[ino2 + 2u];
	p0 += p1 + p2;

	ino = bvh._size - bvh._numFaces + id;
	bvh._mins[0][ino] = p0.x;
	bvh._mins[1][ino] = p0.y;
	bvh._mins[2][ino] = p0.z;

	TriInfo info;
	info._face = id;
	info._id = 0u;
	infos[id] = info;
}
__global__ void InitMinMaxKernel(BVHParam bvh)
{
	__shared__ REAL s_mins[3][MAX_BLOCKSIZE];
	__shared__ REAL s_maxs[3][MAX_BLOCKSIZE];
	uint id = blockDim.x * blockIdx.x + threadIdx.x;
	if (id >= bvh._numFaces) {
		s_mins[0][threadIdx.x] = DBL_MAX;
		s_mins[1][threadIdx.x] = DBL_MAX;
		s_mins[2][threadIdx.x] = DBL_MAX;
		s_maxs[0][threadIdx.x] = -DBL_MAX;
		s_maxs[1][threadIdx.x] = -DBL_MAX;
		s_maxs[2][threadIdx.x] = -DBL_MAX;
		return;
	}

	uint ibvh = bvh._size - bvh._numFaces + id;
	s_mins[0][threadIdx.x] = bvh._mins[0][ibvh];
	s_mins[1][threadIdx.x] = bvh._mins[1][ibvh];
	s_mins[2][threadIdx.x] = bvh._mins[2][ibvh];
	s_maxs[0][threadIdx.x] = s_mins[0][threadIdx.x];
	s_maxs[1][threadIdx.x] = s_mins[1][threadIdx.x];
	s_maxs[2][threadIdx.x] = s_mins[2][threadIdx.x];
	for (uint s = BLOCKSIZE >> 1u; s > 32u; s >>= 1u) {
		__syncthreads();
		if (threadIdx.x < s) {
			for (uint i = 0u; i < 3u; i++) {
				if (s_mins[i][threadIdx.x] > s_mins[i][threadIdx.x + s])
					s_mins[i][threadIdx.x] = s_mins[i][threadIdx.x + s];
				if (s_maxs[i][threadIdx.x] < s_maxs[i][threadIdx.x + s])
					s_maxs[i][threadIdx.x] = s_maxs[i][threadIdx.x + s];
			}
		}
	}
	__syncthreads();

	if (threadIdx.x < 32u) {
		warpMin(s_mins[0], threadIdx.x);
		warpMin(s_mins[1], threadIdx.x);
		warpMin(s_mins[2], threadIdx.x);
		warpMax(s_maxs[0], threadIdx.x);
		warpMax(s_maxs[1], threadIdx.x);
		warpMax(s_maxs[2], threadIdx.x);
		if (threadIdx.x == 0) {
			atomicMin_REAL(bvh._mins[0], s_mins[0][threadIdx.x]);
			atomicMin_REAL(bvh._mins[1], s_mins[1][threadIdx.x]);
			atomicMin_REAL(bvh._mins[2], s_mins[2][threadIdx.x]);
			atomicMax_REAL(bvh._maxs[0], s_maxs[0][threadIdx.x]);
			atomicMax_REAL(bvh._maxs[1], s_maxs[1][threadIdx.x]);
			atomicMax_REAL(bvh._maxs[2], s_maxs[2][threadIdx.x]);
		}
	}
}
__global__ void updateBVHInfo_kernel(
	TriInfo* infos, BVHParam bvh)
{
	uint id = blockDim.x * blockIdx.x + threadIdx.x;
	if (id >= bvh._numFaces)
		return;

	TriInfo info = infos[id];
	uint ino = info._id;
	//if (ino >= bvh._size - bvh._numFaces)
	//	return;

	REAL3 min, max;
	min.x = bvh._mins[0][ino];
	min.y = bvh._mins[1][ino];
	min.z = bvh._mins[2][ino];
	max.x = bvh._maxs[0][ino];
	max.y = bvh._maxs[1][ino];
	max.z = bvh._maxs[2][ino];
	max -= min;

	uint elem = 0u;
	if (max.x < max.y) {
		max.x = max.y;
		elem = 1u;
	}
	if (max.x < max.z)
		elem = 2u;

	ino = bvh._size - bvh._numFaces + info._face;
	info._pos = bvh._mins[elem][ino];
	infos[id] = info;
}
__global__ void updateMinMax_kernel(
	BVHParam bvh, uint level, uint size)
{
	uint id = blockDim.x * blockIdx.x + threadIdx.x;
	if (id >= size)
		return;

	uint ino;
	ino = (1u << level) - 1u + id;
	//if (ino >= bvh._size - bvh._numFaces)
	//	return;

	uint istart = bvh._levels[ino];
	uint iend;
	uint half;
	if (id != size - 1u) 
		iend = bvh._levels[ino + 1u];
	else 
		iend = bvh._numFaces;
	half = iend - istart;

	uint minHalf = 1u << bvh._maxLevel - level - 2u;
	uint maxHalf = 1u << bvh._maxLevel - level - 1u;
	if (level > bvh._maxLevel - 2u) {
		minHalf = 0;
		if (level > bvh._maxLevel - 1u)
			maxHalf = 0;
	}
	if (half < minHalf + maxHalf)
		half -= minHalf;
	else
		half = maxHalf;

	REAL3 min, max;
	min.x = bvh._mins[0][ino];
	min.y = bvh._mins[1][ino];
	min.z = bvh._mins[2][ino];
	max.x = bvh._maxs[0][ino];
	max.y = bvh._maxs[1][ino];
	max.z = bvh._maxs[2][ino];

	ino = (ino << 1u) + 1u;
	bvh._levels[ino] = istart;
	bvh._mins[0][ino] = min.x;
	bvh._mins[1][ino] = min.y;
	bvh._mins[2][ino] = min.z;
	bvh._maxs[0][ino] = max.x;
	bvh._maxs[1][ino] = max.y;
	bvh._maxs[2][ino] = max.z;
	ino++;
	bvh._levels[ino] = istart + half;
	bvh._mins[0][ino] = min.x;
	bvh._mins[1][ino] = min.y;
	bvh._mins[2][ino] = min.z;
	bvh._maxs[0][ino] = max.x;
	bvh._maxs[1][ino] = max.y;
	bvh._maxs[2][ino] = max.z;
}
__global__ void subdivBVH_kernel(
	TriInfo* infos, BVHParam bvh)
{
	uint id = blockDim.x * blockIdx.x + threadIdx.x;
	if (id >= bvh._numFaces)
		return;

	TriInfo info = infos[id];
	uint ino = info._id;
	//if (ino >= bvh._size - bvh._numFaces)
	//	return;

	REAL3 min, max;
	min.x = bvh._mins[0][ino];
	min.y = bvh._mins[1][ino];
	min.z = bvh._mins[2][ino];
	max.x = bvh._maxs[0][ino];
	max.y = bvh._maxs[1][ino];
	max.z = bvh._maxs[2][ino];
	max -= min;

	uint elem = 0u;
	if (max.x < max.y) {
		max.x = max.y;
		elem = 1u;
	}
	if (max.x < max.z)
		elem = 2u;

	ino = (ino << 1u) + 2u;
	uint pivot = bvh._levels[ino];
	if (id < pivot) {
		info._id = ino - 1u;
		if (id == pivot - 1u)
			bvh._maxs[elem][info._id] = info._pos;
	}
	else {
		info._id = ino;
		if (id == pivot)
			bvh._mins[elem][info._id] = info._pos;
	}
	infos[id] = info;
}
__global__ void buildBVH_kernel(
	TriInfo* infos, BVHParam bvh)
{
	uint id = blockDim.x * blockIdx.x + threadIdx.x;
	if (id >= bvh._size)
		return;

	bvh._levels[id] = Log2(id + 1);
	
	/*if (id >= bvh._size - bvh._numFaces) {
		uint fid = id - (bvh._size - bvh._numFaces);
		TriInfo info = infos[fid];
		bvh._faces[fid] = info._face;
	}*/
	uint ileaf = bvh._size - bvh._numFaces;
	if (id >= ileaf) {
		uint fid;
		id -= ileaf;
		fid = id + bvh._pivot;
		if (fid >= bvh._numFaces)
			fid -= bvh._numFaces;

		TriInfo info = infos[fid];
		bvh._faces[id] = info._face;
	}
}

__global__ void RefitLeafBVHKernel(
	ObjParam params, BVHParam bvh, uint num,
	const REAL delta, const REAL dt, const bool isCCD)
{
	uint id = blockDim.x * blockIdx.x + threadIdx.x;
	if (id >= num)
		return;

	uint ind = num - 1u + id;
	uint ileaf = bvh._size - bvh._numFaces;
	AABB node;
	if (ind < ileaf) {
		uint ichild = (ind << 1u) + 1u;
		AABB lchild, rchild;
		uint lfid = bvh._faces[ichild - ileaf];
		uint rfid = bvh._faces[ichild + 1 - ileaf];
		RefitBVHLeaf(lfid, lchild, params, delta, dt, isCCD);
		RefitBVHLeaf(rfid, rchild, params, delta, dt, isCCD);
		updateBVHAABB(bvh, lchild, ichild);
		updateBVHAABB(bvh, rchild, ichild + 1);

		setAABB(node, lchild);
		addAABB(node, rchild);
	}
	else {
		uint fid = bvh._faces[ind - ileaf];
		RefitBVHLeaf(fid, node, params, delta, dt, isCCD);
	}
	updateBVHAABB(bvh, node, ind);
}
__global__ void RefitBVHKernel(
	ObjParam params, BVHParam bvh, uint num,
	const REAL delta, const REAL dt, const bool isCCD)
{
	uint id = blockDim.x * blockIdx.x + threadIdx.x;
	if (id >= num)
		return;

	uint ind = num - 1 + id;
	uint ichild = (ind << 1) + 1;
	AABB parent, lchild, rchild;
	getBVHAABB(lchild, bvh, ichild);
	getBVHAABB(rchild, bvh, ichild + 1);

	setAABB(parent, lchild);
	addAABB(parent, rchild);
	updateBVHAABB(bvh, parent, ind);
}
__global__ void RefitNodeBVHKernel(
	ObjParam params, BVHParam bvh, uint level,
	const REAL delta, const REAL dt, const bool isCCD)
{
	uint id = blockDim.x * blockIdx.x + threadIdx.x;
	int currLev = level;
	uint ind0, ind, ichild;
	AABB parent, lchild, rchild;
	while (currLev > 5) {
		ind0 = (1 << currLev--);
		if (id >= ind0--)
			return;

		ind = ind0 + id;
		ichild = (ind << 1) + 1;
		getBVHAABB(lchild, bvh, ichild);
		getBVHAABB(rchild, bvh, ichild + 1);

		setAABB(parent, lchild);
		addAABB(parent, rchild);
		updateBVHAABB(bvh, parent, ind);
		__syncthreads();
	}
	while (currLev >= 0) {
		ind0 = (1 << currLev--);
		if (id < ind0--) {
			ind = ind0 + id;
			ichild = (ind << 1) + 1;
			getBVHAABB(lchild, bvh, ichild);
			getBVHAABB(rchild, bvh, ichild + 1);

			setAABB(parent, lchild);
			addAABB(parent, rchild);
			updateBVHAABB(bvh, parent, ind);
		}
	}
}