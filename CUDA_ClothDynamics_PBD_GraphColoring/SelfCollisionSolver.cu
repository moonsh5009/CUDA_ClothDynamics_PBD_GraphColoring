#include "CollisionManager.cuh"
#include "SelfCollisionSolver.h"

//-------------------------------------------------------------------------
__device__ void getNumSelfContactElementsProximity(
	const uint lRTri, const uint rRTri,
	const uint* ino, const REAL3* pi,
	const uint* jno, const REAL3* pj,
	REAL delta, uint& num)
{
	uint i, j, i1, j1;
	REAL pene;
	for (i = 0u; i < 3u; i++) {
		if (RTriVertex(rRTri, i)) {
			if (isSelfContactTV_Proximity(
				ino[0], ino[1], ino[2], jno[i], pi[0], pi[1], pi[2], pj[i], delta, &pene)) {
				num++;
			}
		}
		if (RTriVertex(lRTri, i)) {
			if (isSelfContactTV_Proximity(
				jno[0], jno[1], jno[2], ino[i], pj[0], pj[1], pj[2], pi[i], delta, &pene)) {
				num++;
			}
		}

		i1 = (i + 1u) % 3u;
		if (RTriEdge(lRTri, i)) {
			for (j = 0u; j < 3u; j++) {
				j1 = (j + 1u) % 3u;
				if (RTriEdge(rRTri, j)) {
					if (isSelfContactEE_Proximity(
						ino[i], ino[i1], jno[j], jno[j1], pi[i], pi[i1], pj[j], pj[j1], delta, &pene)) {
						num++;
					}
				}
			}
		}
	}
}
__device__ void getNumSelfContactElementsCCD(
	const uint lRTri, const uint rRTri,
	const uint* ino, const REAL3* pi, const REAL3* qi,
	const uint* jno, const REAL3* pj, const REAL3* qj,
	uint& num)
{
	REAL t;
	uint i, j, i1, j1;
	for (i = 0u; i < 3u; i++) {
		if (RTriVertex(rRTri, i)) {
			if (isSelfContactTV_CCD(
				ino[0], ino[1], ino[2], jno[i], pi[0], pi[1], pi[2], pj[i], qi[0], qi[1], qi[2], qj[i], &t)) {
				num++;
			}
		}
		if (RTriVertex(lRTri, i)) {
			if (isSelfContactTV_CCD(
				jno[0], jno[1], jno[2], ino[i], pj[0], pj[1], pj[2], pi[i], qj[0], qj[1], qj[2], qi[i], &t)) {
				num++;
			}
		}

		i1 = (i + 1u) % 3u;
		if (RTriEdge(lRTri, i)) {
			for (j = 0u; j < 3u; j++) {
				j1 = (j + 1u) % 3u;
				if (RTriEdge(rRTri, j)) {
					if (isSelfContactEE_CCD(
						ino[i], ino[i1], jno[j], jno[j1], pi[i], pi[i1], pj[j], pj[j1], qi[i], qi[i1], qj[j], qj[j1], &t)) {
						num++;
					}
				}
			}
		}
	}
}
__device__ void getSelfContactElementsProximity(
	const uint lRTri, const uint rRTri,
	const uint* ino, const REAL3* pi,
	const uint* jno, const REAL3* pj,
	REAL delta,
	ContactElem* ces, uint& num)
{
	uint i, j, i1, j1;
	REAL pene;
	num = 0u;
	for (i = 0u; i < 3u; i++) {
		if (RTriVertex(rRTri, i))
			if (isSelfContactTV_Proximity(
				ino[0], ino[1], ino[2], jno[i], pi[0], pi[1], pi[2], pj[i], delta, &pene))
				makeSelfCE(ces[num++], true, pene, ino[0], ino[1], ino[2], jno[i]);
		if (RTriVertex(lRTri, i))
			if (isSelfContactTV_Proximity(
				jno[0], jno[1], jno[2], ino[i], pj[0], pj[1], pj[2], pi[i], delta, &pene))
				makeSelfCE(ces[num++], true, pene, jno[0], jno[1], jno[2], ino[i]);

		i1 = (i + 1u) % 3u;
		if (RTriEdge(lRTri, i)) {
			for (j = 0u; j < 3u; j++) {
				j1 = (j + 1u) % 3u;
				if (RTriEdge(rRTri, j))
					if (isSelfContactEE_Proximity(
						ino[i], ino[i1], jno[j], jno[j1], pi[i], pi[i1], pj[j], pj[j1], delta, &pene))
						makeSelfCE(ces[num++], false, pene, ino[i], ino[i1], jno[j], jno[j1]);
			}
		}
	}
}
__device__ void getSelfContactElementsCCD(
	const uint lRTri, const uint rRTri,
	const uint* ino, const REAL3* pi, const REAL3* qi,
	const uint* jno, const REAL3* pj, const REAL3* qj,
	ContactElem* ces, uint& num)
{
	REAL t;
	uint i, j, i1, j1;
	num = 0;
	for (i = 0u; i < 3u; i++) {
		if (RTriVertex(rRTri, i))
			if (isSelfContactTV_CCD(
				ino[0], ino[1], ino[2], jno[i], pi[0], pi[1], pi[2], pj[i], qi[0], qi[1], qi[2], qj[i], &t))
				makeSelfCE(ces[num++], true, t, ino[0], ino[1], ino[2], jno[i]);
		if (RTriVertex(lRTri, i))
			if (isSelfContactTV_CCD(
				jno[0], jno[1], jno[2], ino[i], pj[0], pj[1], pj[2], pi[i], qj[0], qj[1], qj[2], qi[i], &t))
				makeSelfCE(ces[num++], true, t, jno[0], jno[1], jno[2], ino[i]);

		i1 = (i + 1u) % 3u;
		if (RTriEdge(lRTri, i)) {
			for (j = 0u; j < 3u; j++) {
				j1 = (j + 1u) % 3u;
				if (RTriEdge(rRTri, j))
					if (isSelfContactEE_CCD(
						ino[i], ino[i1], jno[j], jno[j1], pi[i], pi[i1], pj[j], pj[j1], qi[i], qi[i1], qj[j], qj[j1], &t))
						makeSelfCE(ces[num++], false, t, ino[i], ino[i1], jno[j], jno[j1]);
			}
		}
	}
}
//-------------------------------------------------------------------------
__device__ __forceinline__ void getNumSelfContactElementsProximity_device(
	const ObjParam& param,
	uint f0, uint f1,
	uint RTriA, uint RTriB,
	REAL delta, uint& num)
{
	REAL3 pi[3], pj[3];
	uint ino[3], jno[3];
	getIfromParam(param._fs, f0, ino);
	getIfromParam(param._fs, f1, jno);
	getXfromParam(param._ns, ino, pi);
	getXfromParam(param._ns, jno, pj);

	getNumSelfContactElementsProximity(
		RTriA, RTriB, ino, pi, jno, pj, delta, num);
}
__device__ __forceinline__ void getNumSelfContactElementsCCD_device(
	const ObjParam& param,
	uint f0, uint f1,
	uint RTriA, uint RTriB,
	REAL dt, uint& num)
{
	REAL3 pi[3], qi[3], pj[3], qj[3];
	uint ino[3], jno[3];

	getIfromParam(param._fs, f0, ino);
	getIfromParam(param._fs, f1, jno);
	getXfromParam(param._ns, ino, pi);
	getXfromParam(param._ns, jno, pj);
	getXfromParam(param._vs, ino, qi);
	getXfromParam(param._vs, jno, qj);
	qi[0] = pi[0] + qi[0] * dt;
	qi[1] = pi[1] + qi[1] * dt;
	qi[2] = pi[2] + qi[2] * dt;
	qj[0] = pj[0] + qj[0] * dt;
	qj[1] = pj[1] + qj[1] * dt;
	qj[2] = pj[2] + qj[2] * dt;

	getNumSelfContactElementsCCD(
		RTriA, RTriB, ino, pi, qi, jno, pj, qj, num);
}
__device__ __forceinline__ void getSelfContactElementsProximity_device(
	const ObjParam& param,
	const ContactElemParam& ceParam,
	uint f0, uint f1,
	const uint RTriA, const uint RTriB,
	REAL delta)
{
	ContactElem ces[15];
	uint ceSize;

	REAL3 pi[3], pj[3];
	uint ino[3], jno[3];
	getIfromParam(param._fs, f0, ino);
	getIfromParam(param._fs, f1, jno);
	getXfromParam(param._ns, ino, pi);
	getXfromParam(param._ns, jno, pj);

	getSelfContactElementsProximity(
		RTriA, RTriB, ino, pi, jno, pj, delta, ces, ceSize);
	if (ceSize > 0)
		addCE(ceParam, ces, ceSize);
}
__device__ __forceinline__ void getSelfContactElementsCCD_device(
	const ObjParam& param,
	const ContactElemParam& ceParam,
	uint f0, uint f1,
	const uint RTriA, const uint RTriB,
	REAL dt)
{
	ContactElem ces[15];
	uint ceSize;

	REAL3 pi[3], qi[3], pj[3], qj[3];
	uint ino[3], jno[3];

	getIfromParam(param._fs, f0, ino);
	getIfromParam(param._fs, f1, jno);
	getXfromParam(param._ns, ino, pi);
	getXfromParam(param._ns, jno, pj);
	getXfromParam(param._vs, ino, qi);
	getXfromParam(param._vs, jno, qj);
	qi[0] = pi[0] + qi[0] * dt;
	qi[1] = pi[1] + qi[1] * dt;
	qi[2] = pi[2] + qi[2] * dt;
	qj[0] = pj[0] + qj[0] * dt;
	qj[1] = pj[1] + qj[1] * dt;
	qj[2] = pj[2] + qj[2] * dt;

	getSelfContactElementsCCD(
		RTriA, RTriB, ino, pi, qi, jno, pj, qj, ces, ceSize);
	if (ceSize > 0)
		addCE(ceParam, ces, ceSize);
}
//-------------------------------------------------------------------------
__global__ void getNumSelfCEProximity_kernel(
	ObjParam param,
	BVHParam bvh, RTriParam RTri,
	uint* CEsize,
	const REAL delta, const REAL dt)
{
	/*__shared__ uint s_bvtts[2][BVTT_SHARED_SIZE];
	__shared__ uint s_bvttNum;
	__shared__ uint s_endNum;
	uint id = blockDim.x * blockIdx.x + threadIdx.x;

	BVHNode node, nodeA, nodeB;
	uint RTriA, RTriB, ind, icomp;
	uint ibvtt, i;
	if (threadIdx.x == 0u) {
		i = blockDim.x * (blockIdx.x + 1u);
		if (i > param._numFaces)
			s_endNum = i - param._numFaces;
		else s_endNum = 0u;
		s_bvttNum = 0u;
	}
	__syncthreads();

	uint ileaf = bvh._size - param._numFaces;
	uint path, level, lastPath;

	bool isLoop = true;
	bool goback = false;

	if (id >= param._numFaces)
		isLoop = false;
	else {
		path = id;
		level = bvh._maxLevel;
		if (id >= bvh._pivot) {
			path = id - bvh._pivot + (bvh._pivot >> 1u);
			level--;
		}
		i = getBVHIndex(path, level);
		getBVHNode(node, bvh, i);
		{
			uint half = bvh._numFaces >> 1u;
			lastPath = (id + half);
			if (!(bvh._numFaces & 1u) && id >= half)
				lastPath--;
			if (lastPath >= bvh._numFaces) {
				goback = true;
				lastPath -= bvh._numFaces;
			}

			if (lastPath >= bvh._pivot)
				lastPath += lastPath - bvh._pivot;
		}

		if ((path & 1u) == 0u)
			path |= 1u;
		else {
			if (path + 1u >> level) {
				if (goback) {
					path = 0u;
					level = 1u;
					goback = false;
				}
			}
			else {
				do {
					level--;
					path >>= 1u;
				} while (path & 1u);
				path |= 1u;
			}
		}
	}
	uint num = 0u;
	do {
		__syncthreads();
		if (isLoop) {
			icomp = getBVHIndex(path, level);
			getBVHAABB(nodeA._aabb, bvh, icomp);
			bool isIntersect = intersect(node._aabb, nodeA._aabb);
			if (isIntersect && icomp < ileaf) {
				level++;
				path <<= 1u;
			}
			else {
				if (isIntersect) {
					ibvtt = atomicAdd(&s_bvttNum, 1u);
					s_bvtts[0][ibvtt] = node._face;
					s_bvtts[1][ibvtt] = bvh._faces[icomp - ileaf];
				}

				if ((path & 1u) == 0u)
					path |= 1u;
				else {
					if (path + 1u >> level) {
						if (goback) {
							path = 0u;
							level = 1u;
							goback = false;
						}
						else isLoop = false;
					}
					else {
						do {
							level--;
							path >>= 1u;
						} while (path & 1u);
						path |= 1u;
					}
					if (!goback) {
						if (path > (lastPath >> bvh._maxLevel - level))
							isLoop = false;
					}
				}
			}
			if (!isLoop) atomicAdd(&s_endNum, 1u);
		}
		__syncthreads();
		if (s_endNum >= blockDim.x || s_bvttNum >= BVTT_SHARED_SIZE - blockDim.x) {
			for (i = threadIdx.x; i < s_bvttNum; i += blockDim.x) {
				RTriA = RTri._info[s_bvtts[0][i]];
				RTriB = RTri._info[s_bvtts[1][i]];
				getNumSelfContactElementsProximity_device(
					param, s_bvtts[0][i], s_bvtts[1][i], RTriA, RTriB, delta, num);
			}
			__syncthreads();
			if (threadIdx.x == 0u)
				s_bvttNum = 0u;
		}
	} while (s_endNum < blockDim.x);

	s_bvtts[0][threadIdx.x] = num;
	for (i = blockDim.x >> 1u; i > 32u; i >>= 1u) {
		__syncthreads();
		if (threadIdx.x < i)
			s_bvtts[0][threadIdx.x] += s_bvtts[0][threadIdx.x + i];
	}
	__syncthreads();
	if (threadIdx.x < 32u) {
		warpSum(s_bvtts[0], threadIdx.x);
		if (threadIdx.x == 0u)
			atomicAdd(CEsize, s_bvtts[0][0]);
	}*/
	/*extern __shared__ uint s_sums[];
	uint id = blockDim.x * blockIdx.x + threadIdx.x;

	uint compfaces[3000];
	uint icompfaces = 0u;

	BVHNode node, comp;
	uint icomp, RTriA, RTriB, i;

	uint ileaf = bvh._size - param._numFaces;
	uint path, level, lastPath;

	bool isLoop = true;
	bool goback = false;

	if (id >= param._numFaces) {
		s_sums[threadIdx.x] = 0u;
		return;
	}

	path = id;
	level = bvh._maxLevel;
	if (id >= bvh._pivot) {
		path = id - bvh._pivot + (bvh._pivot >> 1u);
		level--;
	}
	i = getBVHIndex(path, level);
	getBVHNode(node, bvh, i);

	{
		uint half = bvh._numFaces >> 1u;
		lastPath = (id + half);
		if (!(bvh._numFaces & 1u) && id >= half)
			lastPath--;
		if (lastPath >= bvh._numFaces) {
			goback = true;
			lastPath -= bvh._numFaces;
		}

		if (lastPath >= bvh._pivot)
			lastPath += lastPath - bvh._pivot;
	}

	if ((path & 1u) == 0u)
		path |= 1u;
	else {
		if (path + 1u >> level) {
			if (goback) {
				path = 0u;
				level = 1u;
				goback = false;
			}
		}
		else {
			do {
				level--;
				path >>= 1u;
			} while (path & 1u);
			path |= 1u;
		}
	}
	if (!goback) {
		if (path >= (lastPath >> bvh._maxLevel - level))
			isLoop = false;
	}

	uint num = 0u;

	uint ino[3], jno[3];
	REAL3 pi[3], pj[3];
	RTriA = RTri._info[node._face];
	getIfromParam(param._fs, node._face, ino);
	getXfromParam(param._ns, ino, pi);
	while (isLoop) {
		icomp = getBVHIndex(path, level);
		getBVHAABB(comp._aabb, bvh, icomp);
		bool isIntersect = intersect(node._aabb, comp._aabb);
		if (isIntersect && icomp < ileaf) {
			level++;
			path <<= 1u;
		}
		else {
			if (isIntersect)
				compfaces[icompfaces++] = bvh._faces[icomp - ileaf];

			if ((path & 1u) == 0u)
				path |= 1u;
			else {
				if (path + 1u >> level) {
					if (goback) {
						path = 0u;
						level = 1u;
						goback = false;
					}
				}
				else {
					do {
						level--;
						path >>= 1u;
					} while (path & 1u);
					path |= 1u;
				}
			}
			if (!goback) {
				if (path >= (lastPath >> bvh._maxLevel - level))
					isLoop = false;
			}
		}
		if (icompfaces >= 3000u) {
			for (i = 0u; i < icompfaces; i++) {
				RTriB = RTri._info[compfaces[i]];
				getIfromParam(param._fs, compfaces[i], jno);
				getXfromParam(param._ns, jno, pj);
				getNumSelfContactElementsProximity(
					RTriA, RTriB, ino, pi, jno, pj, delta, num);
			}
			icompfaces = 0u;
		}
	}

	for (i = 0u; i < icompfaces; i++) {
		RTriB = RTri._info[compfaces[i]];
		getIfromParam(param._fs, compfaces[i], jno);
		getXfromParam(param._ns, jno, pj);
		getNumSelfContactElementsProximity(
			RTriA, RTriB, ino, pi, jno, pj, delta, num);
	}

	s_sums[threadIdx.x] = num;
	for (i = blockDim.x >> 1u; i > 32u; i >>= 1u) {
		__syncthreads();
		if (threadIdx.x < i)
			s_sums[threadIdx.x] += s_sums[threadIdx.x + i];
	}
	__syncthreads();
	if (threadIdx.x < 32u) {
		warpSum(s_sums, threadIdx.x);
		if (threadIdx.x == 0u)
			atomicAdd(CEsize, s_sums[0]);
	}*/
	__shared__ uint s_bvtts[2][BVTT_SHARED_SIZE];
	__shared__ uint s_bvttNum;
	__shared__ uint s_endNum;
	uint id = blockDim.x * blockIdx.x + threadIdx.x;

	BVHNode node, nodeA, nodeB;
	uint RTriA, RTriB, ind, icomp;
	uint ibvtt, i;
	if (threadIdx.x == 0u) {
		i = blockDim.x * (blockIdx.x + 1u);
		if (i > param._numFaces)
			s_endNum = i - param._numFaces;
		else s_endNum = 0u;
		s_bvttNum = 0u;
	}
	__syncthreads();

	uint ileaf = bvh._size - param._numFaces;
	uint path, level;

	bool isLoop = true;

	if (id >= param._numFaces)
		isLoop = false;
	else {
		path = id;
		level = bvh._maxLevel;
		if (id >= bvh._pivot) {
			path = id - bvh._pivot + (bvh._pivot >> 1u);
			level--;
		}
		i = getBVHIndex(path, level);
		getBVHNode(node, bvh, i);
		path = 0u;
		level = 1u;
	}
	uint num = 0u;
	do {
		__syncthreads();
		if (isLoop) {
			icomp = getBVHIndex(path, level);
			getBVHAABB(nodeA._aabb, bvh, icomp);
			bool isIntersect = intersect(node._aabb, nodeA._aabb);
			if (isIntersect && icomp < ileaf) {
				level++;
				path <<= 1u;
			}
			else {
				if (isIntersect) {
					uint f = bvh._faces[icomp - ileaf];
					if (node._face < f) {
						ibvtt = atomicAdd(&s_bvttNum, 1u);
						s_bvtts[0][ibvtt] = node._face;
						s_bvtts[1][ibvtt] = f;
					}
				}

				if ((path & 1u) == 0u)
					path |= 1u;
				else {
					if (path + 1u >> level)
						isLoop = false;
					else {
						do {
							level--;
							path >>= 1u;
						} while (path & 1u);
						path |= 1u;
					}
				}
			}
			if (!isLoop) atomicAdd(&s_endNum, 1u);
		}
		__syncthreads();
		if (s_endNum >= blockDim.x || s_bvttNum >= BVTT_SHARED_SIZE - blockDim.x) {
			for (i = threadIdx.x; i < s_bvttNum; i += blockDim.x) {
				RTriA = RTri._info[s_bvtts[0][i]];
				RTriB = RTri._info[s_bvtts[1][i]];
				getNumSelfContactElementsProximity_device(
					param, s_bvtts[0][i], s_bvtts[1][i], RTriA, RTriB, delta, num);
			}
			__syncthreads();
			if (threadIdx.x == 0u)
				s_bvttNum = 0u;
		}
	} while (s_endNum < blockDim.x);

	s_bvtts[0][threadIdx.x] = num;
	for (i = blockDim.x >> 1u; i > 32u; i >>= 1u) {
		__syncthreads();
		if (threadIdx.x < i)
			s_bvtts[0][threadIdx.x] += s_bvtts[0][threadIdx.x + i];
	}
	__syncthreads();
	if (threadIdx.x < 32u) {
		warpSum(s_bvtts[0], threadIdx.x);
		if (threadIdx.x == 0u)
			atomicAdd(CEsize, s_bvtts[0][0]);
	}
}
__global__ void getNumSelfCECCDtest_kernel(
	ObjParam param,
	BVHParam bvh, RTriParam RTri,
	uint* CEsize,
	const REAL delta, const REAL dt)
{
	__shared__ uint s_bvtts[2][BVTT_SHARED_SIZE];
	__shared__ uint s_bvttNum;
	__shared__ uint s_endNum;
	uint id = blockDim.x * blockIdx.x + threadIdx.x;

	BVHNode node, nodeA, nodeB;
	uint RTriA, RTriB, ind, icomp;
	uint ibvtt, i;
	if (threadIdx.x == 0u) {
		i = blockDim.x * (blockIdx.x + 1u);
		if (i > param._numFaces)
			s_endNum = i - param._numFaces;
		else s_endNum = 0u;
		s_bvttNum = 0u;
	}
	__syncthreads();

	uint ileaf = bvh._size - param._numFaces;
	uint path, level;

	bool isLoop = true;

	if (id >= param._numFaces)
		isLoop = false;
	else {
		path = id;
		level = bvh._maxLevel;
		if (id >= bvh._pivot) {
			path = id - bvh._pivot + (bvh._pivot >> 1u);
			level--;
		}
		i = getBVHIndex(path, level);
		getBVHNode(node, bvh, i);
		path = 0u;
		level = 1u;
	}
	uint num = 0u;
	do {
		__syncthreads();
		if (isLoop) {
			icomp = getBVHIndex(path, level);
			getBVHAABB(nodeA._aabb, bvh, icomp);
			bool isIntersect = intersect(node._aabb, nodeA._aabb);
			if (isIntersect && icomp < ileaf) {
				level++;
				path <<= 1u;
			}
			else {
				if (isIntersect) {
					uint f = bvh._faces[icomp - ileaf];
					if (node._face < f) {
						ibvtt = atomicAdd(&s_bvttNum, 1u);
						s_bvtts[0][ibvtt] = node._face;
						s_bvtts[1][ibvtt] = f;
					}
				}

				if ((path & 1u) == 0u)
					path |= 1u;
				else {
					if (path + 1u >> level)
						isLoop = false;
					else {
						do {
							level--;
							path >>= 1u;
						} while (path & 1u);
						path |= 1u;
					}
				}
			}
			if (!isLoop) atomicAdd(&s_endNum, 1u);
		}
		__syncthreads();
		if (s_endNum >= blockDim.x || s_bvttNum >= BVTT_SHARED_SIZE - blockDim.x) {
			for (i = threadIdx.x; i < s_bvttNum; i += blockDim.x) {
				RTriA = RTri._info[s_bvtts[0][i]];
				RTriB = RTri._info[s_bvtts[1][i]];
				getNumSelfContactElementsCCD_device(
					param, s_bvtts[0][i], s_bvtts[1][i], RTriA, RTriB, dt, num);
			}
			__syncthreads();
			if (threadIdx.x == 0u)
				s_bvttNum = 0u;
		}
	} while (s_endNum < blockDim.x);

	s_bvtts[0][threadIdx.x] = num;
	for (i = blockDim.x >> 1u; i > 32u; i >>= 1u) {
		__syncthreads();
		if (threadIdx.x < i)
			s_bvtts[0][threadIdx.x] += s_bvtts[0][threadIdx.x + i];
	}
	__syncthreads();
	if (threadIdx.x < 32u) {
		warpSum(s_bvtts[0], threadIdx.x);
		if (threadIdx.x == 0u)
			atomicAdd(CEsize, s_bvtts[0][0]);
	}
}
__global__ void getNumSelfCECCD_kernel(
	ObjParam param,
	BVHParam bvh, RTriParam RTri,
	uint* CEsize,
	const REAL delta, const REAL dt)
{
	/*__shared__ uint s_bvtts[2][BVTT_SHARED_SIZE];
	__shared__ uint s_bvttNum;
	__shared__ uint s_endNum;
	uint id = blockDim.x * blockIdx.x + threadIdx.x;

	BVHNode node, nodeA, nodeB;
	uint RTriA, RTriB, ind, icomp;
	uint ibvtt, i;
	if (threadIdx.x == 0u) {
		i = blockDim.x * (blockIdx.x + 1u);
		if (i > param._numFaces)
			s_endNum = i - param._numFaces;
		else s_endNum = 0u;
		s_bvttNum = 0u;
	}
	__syncthreads();

	uint ileaf = bvh._size - param._numFaces;
	uint path, level, lastPath;

	bool isLoop = true;
	bool goback = false;

	if (id >= param._numFaces)
		isLoop = false;
	else {
		path = id;
		level = bvh._maxLevel;
		if (id >= bvh._pivot) {
			path -= bvh._pivot;
			path += (bvh._pivot >> 1u);
			level--;
		}
		i = getBVHIndex(path, level);
		getBVHNode(node, bvh, i);
		{
			uint half = bvh._numFaces >> 1u;
			lastPath = (id + half);
			if (!(bvh._numFaces & 1u) && id >= half)
				lastPath--;
			if (lastPath >= bvh._numFaces)
				lastPath -= bvh._numFaces;

			if (lastPath >= bvh._pivot)
				lastPath += lastPath - bvh._pivot;
			uint test;
			test = path << bvh._maxLevel - level;
		}

		if ((path & 1u) == 0u)
			path |= 1u;
		else {
			if (path + 1u >> level) {
				if (goback) {
					path = 0u;
					level = 1u;
					goback = false;
				}
			}
			else {
				do {
					level--;
					path >>= 1u;
				} while (path & 1u);
				path |= 1u;
			}
		}
	}
	uint num = 0u;
	do {
		__syncthreads();
		if (isLoop) {
			icomp = getBVHIndex(path, level);
			getBVHAABB(nodeA._aabb, bvh, icomp);
			bool isIntersect = intersect(node._aabb, nodeA._aabb);
			if (isIntersect && icomp < ileaf) {
				level++;
				path <<= 1u;
			}
			else {
				if (isIntersect) {
					ibvtt = atomicAdd(&s_bvttNum, 1u);
					s_bvtts[0][ibvtt] = node._face;
					s_bvtts[1][ibvtt] = bvh._faces[icomp - ileaf];
				}

				if ((path & 1u) == 0u)
					path |= 1u;
				else {
					if (path + 1u >> level) {
						if (goback) {
							path = 0u;
							level = 1u;
							goback = false;
						}
						else isLoop = false;
					}
					else {
						do {
							level--;
							path >>= 1u;
						} while (path & 1u);
						path |= 1u;
					}
				}
			}
			if (!goback && path > (lastPath >> bvh._maxLevel - level))
				isLoop = false;
			if (!isLoop) atomicAdd(&s_endNum, 1u);
		}
		__syncthreads();
		if (s_endNum >= blockDim.x || s_bvttNum >= BVTT_SHARED_SIZE - blockDim.x) {
			for (i = threadIdx.x; i < s_bvttNum; i += blockDim.x) {
				RTriA = RTri._info[s_bvtts[0][i]];
				RTriB = RTri._info[s_bvtts[1][i]];
				getNumSelfContactElementsCCD_device(
					param, s_bvtts[0][i], s_bvtts[1][i], RTriA, RTriB, dt, num);
			}
			__syncthreads();
			if (threadIdx.x == 0u)
				s_bvttNum = 0u;
		}
	} while (s_endNum < blockDim.x);

	s_bvtts[0][threadIdx.x] = num;
	for (i = blockDim.x >> 1u; i > 32u; i >>= 1u) {
		__syncthreads();
		if (threadIdx.x < i)
			s_bvtts[0][threadIdx.x] += s_bvtts[0][threadIdx.x + i];
	}
	__syncthreads();
	if (threadIdx.x < 32u) {
		warpSum(s_bvtts[0], threadIdx.x);
		if (threadIdx.x == 0u)
			atomicAdd(CEsize, s_bvtts[0][0]);
	}*/

	__shared__ uint s_bvtts[2][BVTT_SHARED_SIZE];
	__shared__ uint s_bvttNum;
	__shared__ uint s_endNum;
	uint id = blockDim.x * blockIdx.x + threadIdx.x;

	BVHNode node, nodeA, nodeB;
	uint RTriA, RTriB, ind, icomp;
	uint ibvtt, i;
	if (threadIdx.x == 0u) {
		i = blockDim.x * (blockIdx.x + 1u);
		if (i > param._numFaces)
			s_endNum = i - param._numFaces;
		else s_endNum = 0u;
		s_bvttNum = 0u;
	}
	__syncthreads();

	uint ileaf = bvh._size - param._numFaces;
	uint path, level;

	bool isLoop = true;

	if (id >= param._numFaces)
		isLoop = false;
	else {
		path = id;
		level = bvh._maxLevel;
		if (id >= bvh._pivot) {
			path = id - bvh._pivot + (bvh._pivot >> 1u);
			level--;
		}
		i = getBVHIndex(path, level);
		getBVHNode(node, bvh, i);
		path = 0u;
		level = 1u;
	}
	uint num = 0u;
	do {
		__syncthreads();
		if (isLoop) {
			icomp = getBVHIndex(path, level);
			getBVHAABB(nodeA._aabb, bvh, icomp);
			bool isIntersect = intersect(node._aabb, nodeA._aabb);
			if (isIntersect && icomp < ileaf) {
				level++;
				path <<= 1u;
			}
			else {
				if (isIntersect) {
					uint f = bvh._faces[icomp - ileaf];
					if (node._face < f) {
						ibvtt = atomicAdd(&s_bvttNum, 1u);
						s_bvtts[0][ibvtt] = node._face;
						s_bvtts[1][ibvtt] = f;
					}
				}

				if ((path & 1u) == 0u)
					path |= 1u;
				else {
					if (path + 1u >> level)
						isLoop = false;
					else {
						do {
							level--;
							path >>= 1u;
						} while (path & 1u);
						path |= 1u;
					}
				}
			}
			if (!isLoop) atomicAdd(&s_endNum, 1u);
		}
		__syncthreads();
		if (s_endNum >= blockDim.x || s_bvttNum >= BVTT_SHARED_SIZE - blockDim.x) {
			for (i = threadIdx.x; i < s_bvttNum; i += blockDim.x) {
				RTriA = RTri._info[s_bvtts[0][i]];
				RTriB = RTri._info[s_bvtts[1][i]];
				getNumSelfContactElementsCCD_device(
					param, s_bvtts[0][i], s_bvtts[1][i], RTriA, RTriB, dt, num);
			}
			__syncthreads();
			if (threadIdx.x == 0u)
				s_bvttNum = 0u;
		}
	} while (s_endNum < blockDim.x);

	s_bvtts[0][threadIdx.x] = num;
	for (i = blockDim.x >> 1u; i > 32u; i >>= 1u) {
		__syncthreads();
		if (threadIdx.x < i)
			s_bvtts[0][threadIdx.x] += s_bvtts[0][threadIdx.x + i];
	}
	__syncthreads();
	if (threadIdx.x < 32u) {
		warpSum(s_bvtts[0], threadIdx.x);
		if (threadIdx.x == 0u)
			atomicAdd(CEsize, s_bvtts[0][0]);
	}
}
//__global__ void getSelfCEProximity_kernel(
//	ObjParam param,
//	BVHParam bvh, RTriParam RTri,
//	ContactElemParam ceParam,
//	const REAL delta, const REAL dt)
//{
//	__shared__ uint s_bvtts[2][BVTT_SHARED_SIZE];
//	__shared__ uint s_bvttNum;
//	__shared__ uint s_endNum;
//	uint id = blockDim.x * blockIdx.x + threadIdx.x;
//
//	BVHNode node, nodeA, nodeB;
//	uint RTriA, RTriB, ind, icomp;
//	uint ibvtt, i;
//	bool isLoop = true;
//	if (threadIdx.x == 0u)
//		s_bvttNum = s_endNum = 0u;
//	__syncthreads();
//
//	uint ileaf = bvh._size - param._numFaces;
//	uint path, level;
//	if (id >= param._numFaces - 1u) {
//		if (id == param._numFaces - 1u)
//			s_endNum = blockDim.x - threadIdx.x;
//		isLoop = false;
//	}
//	else {
//		if (id < bvh._pivot) {
//			ind = getBVHIndex(id, bvh._maxLevel);
//			path = id;
//		}
//		else {
//			ind = id - bvh._pivot;
//			path = ind + (bvh._pivot >> 1u);
//			ind += ileaf;
//		}
//
//		getBVHNode(node, bvh, ind);
//		level = node._level;
//
//		if ((path & 1u) == 0u)
//			path |= 1u;
//		else {
//			if (path + 1u >> level) {
//				atomicAdd(&s_endNum, 1u);
//				isLoop = false;
//			}
//			else {
//				do {
//					level--;
//					path >>= 1u;
//				} while (path & 1u);
//				path |= 1u;
//			}
//		}
//	}
//
//	do {
//		__syncthreads();
//		if (isLoop) {
//			icomp = getBVHIndex(path, level);
//			getBVHAABB(nodeA._aabb, bvh, icomp);
//			bool isIntersect = intersect(node._aabb, nodeA._aabb);
//			if (isIntersect && icomp < ileaf) {
//				level++;
//				path <<= 1u;
//			}
//			else {
//				if (isIntersect) {
//					ibvtt = atomicAdd(&s_bvttNum, 1u);
//					s_bvtts[0][ibvtt] = node._face;
//					s_bvtts[1][ibvtt] = bvh._faces[icomp - ileaf];
//				}
//
//				if ((path & 1u) == 0u)
//					path |= 1u;
//				else {
//					if (path + 1u >> level) {
//						atomicAdd(&s_endNum, 1u);
//						isLoop = false;
//					}
//					else {
//						do {
//							level--;
//							path >>= 1u;
//						} while (path & 1u);
//						path |= 1u;
//					}
//				}
//			}
//		}
//		__syncthreads();
//		if (s_endNum >= blockDim.x || s_bvttNum >= BVTT_SHARED_SIZE - blockDim.x) {
//			for (i = threadIdx.x; i < s_bvttNum; i += blockDim.x) {
//				RTriA = RTri._info[s_bvtts[0][i]];
//				RTriB = RTri._info[s_bvtts[1][i]];
//				getSelfContactElementsProximity_device(
//					param, ceParam, s_bvtts[0][i], s_bvtts[1][i], RTriA, RTriB, delta);
//			}
//			__syncthreads();
//			if (threadIdx.x == 0u)
//				s_bvttNum = 0u;
//		}
//	} while (s_endNum < blockDim.x);
//}
__global__ void getSelfCEProximity_kernel(
	ObjParam param,
	BVHParam bvh, RTriParam RTri,
	ContactElemParam ceParam,
	const REAL delta, const REAL dt)
{
	/*__shared__ uint s_bvtts[2][BVTT_SHARED_SIZE];
	__shared__ uint s_bvttNum;
	__shared__ uint s_endNum;
	uint id = blockDim.x * blockIdx.x + threadIdx.x;

	BVHNode node, nodeA, nodeB;
	uint RTriA, RTriB, ind, icomp;
	uint ibvtt, i;
	if (threadIdx.x == 0u) {
		i = blockDim.x * (blockIdx.x + 1u);
		if (i > param._numFaces)
			s_endNum = i - param._numFaces;
		else s_endNum = 0u;
		s_bvttNum = 0u;
	}
	__syncthreads();

	uint ileaf = bvh._size - param._numFaces;
	uint path, level, lastPath;

	bool isLoop = true;
	bool goback = false;

	if (id >= param._numFaces)
		isLoop = false;
	else {
		path = id;
		level = bvh._maxLevel;
		if (id >= bvh._pivot) {
			path = id - bvh._pivot + (bvh._pivot >> 1u);
			level--;
		}
		i = getBVHIndex(path, level);
		getBVHNode(node, bvh, i);
		{
			uint half = bvh._numFaces >> 1u;
			lastPath = (id + half);
			if (!(bvh._numFaces & 1u) && id >= half)
				lastPath--;
			if (lastPath >= bvh._numFaces) {
				goback = true;
				lastPath -= bvh._numFaces;
			}

			if (lastPath >= bvh._pivot)
				lastPath += lastPath - bvh._pivot;
		}

		if ((path & 1u) == 0u)
			path |= 1u;
		else {
			if (path + 1u >> level) {
				if (goback) {
					path = 0u;
					level = 1u;
					goback = false;
				}
			}
			else {
				do {
					level--;
					path >>= 1u;
				} while (path & 1u);
				path |= 1u;
			}
		}
		if (!goback) {
			if (path > (lastPath >> bvh._maxLevel - level)) {
				atomicAdd(&s_endNum, 1u);
				isLoop = false;
			}
		}
	}

	do {
		__syncthreads();
		if (isLoop) {
			icomp = getBVHIndex(path, level);
			getBVHAABB(nodeA._aabb, bvh, icomp);
			bool isIntersect = intersect(node._aabb, nodeA._aabb);
			if (isIntersect && icomp < ileaf) {
				level++;
				path <<= 1u;
			}
			else {
				if (isIntersect) {
					ibvtt = atomicAdd(&s_bvttNum, 1u);
					s_bvtts[0][ibvtt] = node._face;
					s_bvtts[1][ibvtt] = bvh._faces[icomp - ileaf];
				}

				if ((path & 1u) == 0u)
					path |= 1u;
				else {
					if (path + 1u >> level) {
						if (goback) {
							path = 0u;
							level = 1u;
							goback = false;
						}
						else isLoop = false;
					}
					else {
						do {
							level--;
							path >>= 1u;
						} while (path & 1u);
						path |= 1u;
					}
					if (!goback) {
						if (path > (lastPath >> bvh._maxLevel - level))
							isLoop = false;
					}
				}
			}
			if (!isLoop) atomicAdd(&s_endNum, 1u);
		}
		__syncthreads();
		if (s_endNum >= blockDim.x || s_bvttNum >= BVTT_SHARED_SIZE - blockDim.x) {
			for (i = threadIdx.x; i < s_bvttNum; i += blockDim.x) {
				RTriA = RTri._info[s_bvtts[0][i]];
				RTriB = RTri._info[s_bvtts[1][i]];
				getSelfContactElementsProximity_device(
					param, ceParam, s_bvtts[0][i], s_bvtts[1][i], RTriA, RTriB, delta);
			}
			__syncthreads();
			if (threadIdx.x == 0u)
				s_bvttNum = 0u;
		}
	} while (s_endNum < blockDim.x);*/
	__shared__ uint s_bvtts[2][BVTT_SHARED_SIZE];
	__shared__ uint s_bvttNum;
	__shared__ uint s_endNum;
	uint id = blockDim.x * blockIdx.x + threadIdx.x;

	BVHNode node, nodeA, nodeB;
	uint RTriA, RTriB, ind, icomp;
	uint ibvtt, i;
	if (threadIdx.x == 0u) {
		i = blockDim.x * (blockIdx.x + 1u);
		if (i > param._numFaces)
			s_endNum = i - param._numFaces;
		else s_endNum = 0u;
		s_bvttNum = 0u;
	}
	__syncthreads();

	uint ileaf = bvh._size - param._numFaces;
	uint path, level;

	bool isLoop = true;

	if (id >= param._numFaces)
		isLoop = false;
	else {
		path = id;
		level = bvh._maxLevel;
		if (id >= bvh._pivot) {
			path = id - bvh._pivot + (bvh._pivot >> 1u);
			level--;
		}
		i = getBVHIndex(path, level);
		getBVHNode(node, bvh, i);
		//getBVHNode(node, bvh, ileaf + id);
		path = 0u;
		level = 1u;
	}
	uint num = 0u;
	do {
		__syncthreads();
		if (isLoop) {
			icomp = getBVHIndex(path, level);
			getBVHAABB(nodeA._aabb, bvh, icomp);
			bool isIntersect = intersect(node._aabb, nodeA._aabb);
			if (isIntersect && icomp < ileaf) {
				level++;
				path <<= 1u;
			}
			else {
				if (isIntersect) {
					uint f = bvh._faces[icomp - ileaf];
					if (node._face < f) {
						ibvtt = atomicAdd(&s_bvttNum, 1u);
						s_bvtts[0][ibvtt] = node._face;
						s_bvtts[1][ibvtt] = f;
					}
				}

				if ((path & 1u) == 0u)
					path |= 1u;
				else {
					if (path + 1u >> level)
						isLoop = false;
					else {
						do {
							level--;
							path >>= 1u;
						} while (path & 1u);
						path |= 1u;
					}
				}
			}
			if (!isLoop) atomicAdd(&s_endNum, 1u);
		}
		__syncthreads();
		if (s_endNum >= blockDim.x || s_bvttNum >= BVTT_SHARED_SIZE - blockDim.x) {
			for (i = threadIdx.x; i < s_bvttNum; i += blockDim.x) {
				RTriA = RTri._info[s_bvtts[0][i]];
				RTriB = RTri._info[s_bvtts[1][i]];
				getSelfContactElementsProximity_device(
					param, ceParam, s_bvtts[0][i], s_bvtts[1][i], RTriA, RTriB, delta);
			}
			__syncthreads();
			if (threadIdx.x == 0u)
				s_bvttNum = 0u;
		}
	} while (s_endNum < blockDim.x);
}
__global__ void getSelfCECCD_kernel(
	ObjParam param,
	BVHParam bvh, RTriParam RTri,
	ContactElemParam ceParam,
	const REAL delta, const REAL dt)
{
	/*__shared__ uint s_bvtts[2][BVTT_SHARED_SIZE];
	__shared__ uint s_bvttNum;
	__shared__ uint s_endNum;
	uint id = blockDim.x * blockIdx.x + threadIdx.x;

	BVHNode node, nodeA, nodeB;
	uint RTriA, RTriB, ind, icomp;
	uint ibvtt, i;
	if (threadIdx.x == 0u) {
		i = blockDim.x * (blockIdx.x + 1u);
		if (i > param._numFaces)
			s_endNum = i - param._numFaces;
		else s_endNum = 0u;
		s_bvttNum = 0u;
	}
	__syncthreads();

	uint ileaf = bvh._size - param._numFaces;
	uint path, level, lastPath;

	bool isLoop = true;
	bool goback = false;

	if (id >= param._numFaces)
		isLoop = false;
	else {
		path = id;
		level = bvh._maxLevel;
		if (id >= bvh._pivot) {
			path = id - bvh._pivot + (bvh._pivot >> 1u);
			level--;
		}
		i = getBVHIndex(path, level);
		getBVHNode(node, bvh, i);
		{
			uint half = bvh._numFaces >> 1u;
			lastPath = (id + half);
			if (!(bvh._numFaces & 1u) && id >= half)
				lastPath--;
			if (lastPath >= bvh._numFaces) {
				goback = true;
				lastPath -= bvh._numFaces;
			}

			if (lastPath >= bvh._pivot)
				lastPath += lastPath - bvh._pivot;
		}

		if ((path & 1u) == 0u)
			path |= 1u;
		else {
			if (path + 1u >> level) {
				if (goback) {
					path = 0u;
					level = 1u;
					goback = false;
				}
			}
			else {
				do {
					level--;
					path >>= 1u;
				} while (path & 1u);
				path |= 1u;
			}
		}
		if (!goback) {
			if (path > (lastPath >> bvh._maxLevel - level)) {
				atomicAdd(&s_endNum, 1u);
				isLoop = false;
			}
		}
	}

	do {
		__syncthreads();
		if (isLoop) {
			icomp = getBVHIndex(path, level);
			getBVHAABB(nodeA._aabb, bvh, icomp);
			bool isIntersect = intersect(node._aabb, nodeA._aabb);
			if (isIntersect && icomp < ileaf) {
				level++;
				path <<= 1u;
			}
			else {
				if (isIntersect) {
					ibvtt = atomicAdd(&s_bvttNum, 1u);
					s_bvtts[0][ibvtt] = node._face;
					s_bvtts[1][ibvtt] = bvh._faces[icomp - ileaf];
				}

				if ((path & 1u) == 0u)
					path |= 1u;
				else {
					if (path + 1u >> level) {
						if (goback) {
							path = 0u;
							level = 1u;
							goback = false;
						}
						else isLoop = false;
					}
					else {
						do {
							level--;
							path >>= 1u;
						} while (path & 1u);
						path |= 1u;
					}
					if (!goback) {
						if (path > (lastPath >> bvh._maxLevel - level))
							isLoop = false;
					}
				}
			}
			if (!isLoop) atomicAdd(&s_endNum, 1u);
		}
		__syncthreads();
		if (s_endNum >= blockDim.x || s_bvttNum >= BVTT_SHARED_SIZE - blockDim.x) {
			for (i = threadIdx.x; i < s_bvttNum; i += blockDim.x) {
				RTriA = RTri._info[s_bvtts[0][i]];
				RTriB = RTri._info[s_bvtts[1][i]];
				getSelfContactElementsCCD_device(
					param, ceParam, s_bvtts[0][i], s_bvtts[1][i], RTriA, RTriB, dt);
			}
			__syncthreads();
			if (threadIdx.x == 0u)
				s_bvttNum = 0u;
		}
	} while (s_endNum < blockDim.x);*/
	__shared__ uint s_bvtts[2][BVTT_SHARED_SIZE];
	__shared__ uint s_bvttNum;
	__shared__ uint s_endNum;
	uint id = blockDim.x * blockIdx.x + threadIdx.x;

	BVHNode node, nodeA, nodeB;
	uint RTriA, RTriB, ind, icomp;
	uint ibvtt, i;
	if (threadIdx.x == 0u) {
		i = blockDim.x * (blockIdx.x + 1u);
		if (i > param._numFaces)
			s_endNum = i - param._numFaces;
		else s_endNum = 0u;
		s_bvttNum = 0u;
	}
	__syncthreads();

	uint ileaf = bvh._size - param._numFaces;
	uint path, level;

	bool isLoop = true;

	if (id >= param._numFaces)
		isLoop = false;
	else {
		path = id;
		level = bvh._maxLevel;
		if (id >= bvh._pivot) {
			path = id - bvh._pivot + (bvh._pivot >> 1u);
			level--;
		}
		i = getBVHIndex(path, level);
		getBVHNode(node, bvh, i);
		//getBVHNode(node, bvh, ileaf + id);
		path = 0u;
		level = 1u;
	}
	uint num = 0u;
	do {
		__syncthreads();
		if (isLoop) {
			icomp = getBVHIndex(path, level);
			getBVHAABB(nodeA._aabb, bvh, icomp);
			bool isIntersect = intersect(node._aabb, nodeA._aabb);
			if (isIntersect && icomp < ileaf) {
				level++;
				path <<= 1u;
			}
			else {
				if (isIntersect) {
					uint f = bvh._faces[icomp - ileaf];
					if (node._face < f) {
						ibvtt = atomicAdd(&s_bvttNum, 1u);
						s_bvtts[0][ibvtt] = node._face;
						s_bvtts[1][ibvtt] = f;
					}
				}

				if ((path & 1u) == 0u)
					path |= 1u;
				else {
					if (path + 1u >> level)
						isLoop = false;
					else {
						do {
							level--;
							path >>= 1u;
						} while (path & 1u);
						path |= 1u;
					}
				}
			}
			if (!isLoop) atomicAdd(&s_endNum, 1u);
		}
		__syncthreads();
		if (s_endNum >= blockDim.x || s_bvttNum >= BVTT_SHARED_SIZE - blockDim.x) {
			for (i = threadIdx.x; i < s_bvttNum; i += blockDim.x) {
				RTriA = RTri._info[s_bvtts[0][i]];
				RTriB = RTri._info[s_bvtts[1][i]];
				getSelfContactElementsCCD_device(
					param, ceParam, s_bvtts[0][i], s_bvtts[1][i], RTriA, RTriB, dt);
			}
			__syncthreads();
			if (threadIdx.x == 0u)
				s_bvttNum = 0u;
		}
	} while (s_endNum < blockDim.x);
}
//-------------------------------------------------------------------------
__global__ void getNodeCEs_kernel(
	ContactElemParam ceParam, uint2* NodeCEs)
{
	uint id = blockDim.x * blockIdx.x + threadIdx.x;
	if (id >= ceParam.h_size)
		return;

	ContactElem ce = ceParam._elems[id];
	uint2 nodeCE;
	uint pos = id;
	nodeCE.y = id;
	for (uint i = 0u; i < 4u; i++) {
		nodeCE.x = ce._i[i];
		NodeCEs[pos] = nodeCE;
		pos += ceParam.h_size;
	}
}
__global__ void reorderNodeCEIds_kernel(
	uint2* NodeCEs, uint* iNodeCEs, uint numNodes, uint numNodeCEs)
{
	extern __shared__ uint s_ids[];
	uint id = blockDim.x * blockIdx.x + threadIdx.x;
	uint curr;
	if (id < numNodeCEs) {
		uint2 tmp = NodeCEs[id];
		curr = tmp.x;
		s_ids[threadIdx.x + 1u] = curr;
		if (id > 0u && threadIdx.x == 0u) {
			tmp = NodeCEs[id - 1u];
			s_ids[0] = tmp.x;
		}
	}
	__syncthreads();

	if (id < numNodeCEs) {
		uint i;
		uint prev = s_ids[threadIdx.x];
		if (id == 0u || prev != curr) {
			if (id == 0u) {
				iNodeCEs[0] = 0u;
				prev = 0u;
			}
			for (i = prev + 1u; i <= curr; i++)
				iNodeCEs[i] = id;
		}
		if (id == numNodeCEs - 1u) {
			iNodeCEs[curr + 1u] = id + 1u;
		}
	}
}
__global__ void initDepth_kernel(
	uint2* NodeCEs, uint* iNodeCEs, uint* icurrs, uint* isEnds, uint numNodes)
{
	uint id = blockDim.x * blockIdx.x + threadIdx.x;
	if (id >= numNodes)
		return;

	uint istart = iNodeCEs[id];
	uint iend = iNodeCEs[id + 1u];
	uint ino;
	uint2 nodeCE;

	ino = istart;
	if (ino < iend) {
		nodeCE = NodeCEs[ino];
		atomicAdd(isEnds + nodeCE.y, 1u);
	}
	else ino = 0xffffffff;
	icurrs[id] = ino;
}
__global__ void nextDepth_kernel(
	uint2* NodeCEs, uint* iNodeCEs, uint* icurrs, uint* isEnds, uint numNodes, uint* isApplied)
{
	uint id = blockDim.x * blockIdx.x + threadIdx.x;

	if (id >= numNodes)
		return;

	uint ino = icurrs[id];
	if (ino == 0xffffffff)
		return;

	uint iend = iNodeCEs[id + 1u];
	uint isEnd;
	uint2 nodeCE;

	nodeCE = NodeCEs[ino];
	isEnd = isEnds[nodeCE.y];
	if (isEnd == 0xffffffff) {
		for (ino++; ino < iend; ino++) {
			nodeCE = NodeCEs[ino];
			isEnd = isEnds[nodeCE.y];
			if (isEnd != 0xffffffff) {
				atomicAdd(isEnds + nodeCE.y, 1u);
				break;
			}
		}
		if (ino >= iend)
			ino = 0xffffffff;
		icurrs[id] = ino;
	}
	if (ino < iend) *isApplied = 1u;
}
__global__ void resolveSelfCollisionProximity_kernel(
	ObjParam param, ContactElemParam ceParam,
	uint* isEnds, const REAL delta, const REAL dt)
{
	uint id = blockDim.x * blockIdx.x + threadIdx.x;
	if (id >= ceParam.h_size)
		return;

	uint isEnd = isEnds[id];
	if (isEnd < 4u || isEnd == 0xffffffff)
		return;
	isEnds[id] = 0xffffffff;

	ContactElem ce = ceParam._elems[id];
	resolveSelfCollisionProximity(
		ce._isfv, ce._i[0], ce._i[1], ce._i[2], ce._i[3],
		param, delta, dt);
}
__global__ void resolveSelfCollisionProximity_kernel(
	ObjParam param,
	ContactElemParam ceParam,
	const REAL delta, const REAL dt)
{
	uint ino;
	ContactElem ce;

	for (ino = 0u; ino < ceParam.h_size; ino++) {
		ce = ceParam._elems[ino];
		resolveSelfCollisionProximity(
			ce._isfv, ce._i[0], ce._i[1], ce._i[2], ce._i[3],
			param, delta, dt);
	}
}
__global__ void resolveSelfCollisionCCD_kernel(
	ObjParam param, ContactElemParam ceParam,
	uint* isEnds, const REAL delta, const REAL dt)
{
	uint id = blockDim.x * blockIdx.x + threadIdx.x;
	if (id >= ceParam.h_size)
		return;

	uint isEnd = isEnds[id];
	if (isEnd < 4u || isEnd == 0xffffffff)
		return;
	isEnds[id] = 0xffffffff;

	ContactElem ce = ceParam._elems[id];
	resolveSelfCollisionCCD(
		ce._isfv, ce._i[0], ce._i[1], ce._i[2], ce._i[3],
		param, delta, dt);
}
__global__ void resolveSelfCollisionCCD_kernel(
	ObjParam param,
	ContactElemParam ceParam,
	const REAL delta, const REAL dt)
{
	uint ino;
	ContactElem ce;

	for (ino = 0u; ino < ceParam.h_size; ino++) {
		ce = ceParam._elems[ino];
		resolveSelfCollisionCCD(
			ce._isfv, ce._i[0], ce._i[1], ce._i[2], ce._i[3],
			param, delta, dt);
	}
}
//-------------------------------------------------------------------------
__global__ void ApplyRigidImpactZone_kernel(
	ObjParam param,
	uint* impactZones,
	uint* sizes,
	uint impactZoneNum)
{
	uint id = blockDim.x * blockIdx.x + threadIdx.x;
	if (id >= impactZoneNum)
		return;


	uint istart = sizes[id];
	uint iend = sizes[id + 1];

	/*for (uint i = istart; i < iend; i++) {
		uint ino = impactZones[i];
		for (uint j = i + 1; j < sizes[impactZoneNum]; j++) {
			uint jno = impactZones[j];
			if (ino == jno)
				printf("Error Impact %d, %d\n", ino, jno);
		}
	}*/

	uint nodeNum = iend - istart;
	REAL3 gc = make_REAL3(0.0);
	REAL3 av = make_REAL3(0.0);
	for (uint i = istart; i < iend; i++) {
		uint ino = impactZones[i] * 3;
		const REAL3 p = make_REAL3(param._ns[ino + 0], param._ns[ino + 1], param._ns[ino + 2]);
		const REAL3 v = make_REAL3(param._vs[ino + 0], param._vs[ino + 1], param._vs[ino + 2]);
		gc += p;
		av += v;
	}
	gc *= 1.0 / (REAL)nodeNum;
	av *= 1.0 / (REAL)nodeNum;
	REAL3 L = make_REAL3(0.0);
	REAL I[9] = { 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0 };
	for (uint i = istart; i < iend; i++) {
		uint ino = impactZones[i] * 3;
		const REAL3 p = make_REAL3(param._ns[ino + 0], param._ns[ino + 1], param._ns[ino + 2]);
		const REAL3 v = make_REAL3(param._vs[ino + 0], param._vs[ino + 1], param._vs[ino + 2]);
		L += Cross(p - gc, v - av);
		const REAL3 q = p - gc;
		REAL vdotv = Dot(v, v);
		I[0] += vdotv - q.x * q.x;	I[1] += -q.x * q.y;			I[2] += -q.x * q.z;
		I[3] += -q.y * q.x;			I[4] += vdotv - q.y * q.y;	I[5] += -q.y * q.z;
		I[6] += -q.z * q.x;			I[7] += -q.z * q.y;			I[8] += vdotv - q.z * q.z;
	}
	REAL Iinv[9];
	CalcInvMat3(Iinv, I);
	REAL3 omg;
	omg.x = Iinv[0] * L.x + Iinv[1] * L.y + Iinv[2] * L.z;
	omg.y = Iinv[3] * L.x + Iinv[4] * L.y + Iinv[5] * L.z;
	omg.z = Iinv[6] * L.x + Iinv[7] * L.y + Iinv[8] * L.z;

	for (uint i = istart; i < iend; i++) {
		uint ino = impactZones[i] * 3;
		const REAL3 p = make_REAL3(param._ns[ino + 0], param._ns[ino + 1], param._ns[ino + 2]);
		const REAL3 rot = Cross(p - gc, omg) * -1.0;
		param._vs[ino + 0] = (av.x + rot.x);
		param._vs[ino + 1] = (av.y + rot.y);
		param._vs[ino + 2] = (av.z + rot.z);
	}
}
//-------------------------------------------------------------------------
void SelfCollisionSolver::ResolveSelfCollisionProximity(
	ContactElems& ceParam,
	const ObjParam& clothParams,
	BVHParam& clothBvh, RTriParam& RTri, 
	Dvector<uint2>& NodeCEs,
	Dvector<uint>& iNodeCEs,
	Dvector<uint>& icurrs,
	Dvector<uint>& isEnds,
	const REAL delta, const REAL dt)
{
#ifdef COLLISION_TESTTIMER
	CUDA_CHECK(cudaDeviceSynchronize());
	ctimer timer = CNOW;
#endif

	ceParam.d_size.memset(0);
	getNumSelfCEProximity_kernel << <divup(clothParams._numFaces, BLOCKSIZE), BLOCKSIZE >> > (
		clothParams, clothBvh, RTri, ceParam.d_size(), delta, dt);
	CUDA_CHECK(cudaPeekAtLastError());
	CUDA_CHECK(cudaMemcpy(&ceParam.h_size, ceParam.d_size(), sizeof(uint), cudaMemcpyDeviceToHost));

#ifdef COLLISION_TESTTIMER
	CUDA_CHECK(cudaDeviceSynchronize());
	printf("getNumSelfCEProximity_kernel: %lf msec\n", (CNOW - timer) / 10000.0);
	timer = CNOW;
#endif

	if (ceParam.h_size > 0u) {
		if (ceParam._elems.size() < ceParam.h_size)
			ceParam.resize();

		ceParam.d_size.memset(0);
		getSelfCEProximity_kernel << <divup(clothParams._numFaces, BLOCKSIZE), BLOCKSIZE >> > (
			clothParams, clothBvh, RTri, ceParam.param(), delta, dt);
		CUDA_CHECK(cudaPeekAtLastError());

#ifdef COLLISION_TESTTIMER
		CUDA_CHECK(cudaDeviceSynchronize());
		printf("getSelfCEProximity_kernel: %lf msec\n", (CNOW - timer) / 10000.0);
		timer = CNOW;
#endif

		thrust::sort(thrust::device_ptr<ContactElem>(ceParam._elems.begin()),
			thrust::device_ptr<ContactElem>(ceParam._elems.begin() + ceParam.h_size), ContactElem_CMP());

#ifdef COLLISION_TESTTIMER
		CUDA_CHECK(cudaDeviceSynchronize());
		printf("thrust::sort: %lf msec\n", (CNOW - timer) / 10000.0);
		timer = CNOW;
#endif

//		uint numNodeCEs = ceParam.h_size << 2u;
//		/*Dvector<uint2> NodeCEs(ceParam.h_size << 2u);
//		Dvector<uint> iNodeCEs(clothParams._numNodes + 1u);
//		Dvector<uint> icurrs(clothParams._numNodes);
//		Dvector<uint> isEnds(ceParam.h_size);*/
//		if (NodeCEs.size() < numNodeCEs)
//			NodeCEs.resize(ceParam.h_size << 2u);
//		iNodeCEs.resize(clothParams._numNodes + 1u);
//		icurrs.resize(clothParams._numNodes);
//		isEnds.resize(ceParam.h_size);
//		iNodeCEs.memset(0);
//		isEnds.memset(0);
//		{
//			getNodeCEs_kernel << <divup(ceParam.h_size, BLOCKSIZE), BLOCKSIZE >> > (
//				ceParam.param(), NodeCEs());
//			CUDA_CHECK(cudaPeekAtLastError());
//
//#ifdef COLLISION_TESTTIMER
//			CUDA_CHECK(cudaDeviceSynchronize());
//			printf("getNodeCEs_kernel: %lf msec\n", (CNOW - timer) / 10000.0);
//			timer = CNOW;
//#endif
//
//			thrust::sort(thrust::device_ptr<uint2>(NodeCEs.begin()),
//				thrust::device_ptr<uint2>(NodeCEs.begin() + numNodeCEs), uint2_CMP());
//
//#ifdef COLLISION_TESTTIMER
//			CUDA_CHECK(cudaDeviceSynchronize());
//			printf("thrust::sort: %lf msec\n", (CNOW - timer) / 10000.0);
//			timer = CNOW;
//#endif
//
//			reorderNodeCEIds_kernel << <divup(numNodeCEs, BLOCKSIZE), BLOCKSIZE, (BLOCKSIZE + 1u) * sizeof(uint) >> > (
//				NodeCEs(), iNodeCEs(), clothParams._numNodes, numNodeCEs);
//			CUDA_CHECK(cudaPeekAtLastError());
//
//#ifdef COLLISION_TESTTIMER
//			CUDA_CHECK(cudaDeviceSynchronize());
//			printf("reorderNodeCEIds_kernel: %lf msec\n", (CNOW - timer) / 10000.0);
//			timer = CNOW;
//#endif
//		}
//
//		uint isApplied;
//		//uint itr = 0u;
//
//		initDepth_kernel << < divup(clothParams._numNodes, BLOCKSIZE), BLOCKSIZE >> > (
//			NodeCEs(), iNodeCEs(), icurrs(), isEnds(), clothParams._numNodes);
//		CUDA_CHECK(cudaPeekAtLastError());
//		do {
//			//itr++;
//			resolveSelfCollisionProximity_kernel << < divup(ceParam.h_size, BLOCKSIZE), BLOCKSIZE >> > (
//				clothParams, ceParam.param(), isEnds(), delta, dt);
//			CUDA_CHECK(cudaPeekAtLastError());
//
//			ceParam.d_size.memset(0);
//			nextDepth_kernel << < divup(clothParams._numNodes, BLOCKSIZE), BLOCKSIZE >> > (
//				NodeCEs(), iNodeCEs(), icurrs(), isEnds(), clothParams._numNodes, ceParam.d_size());
//			CUDA_CHECK(cudaPeekAtLastError());
//			CUDA_CHECK(cudaMemcpy(&isApplied, ceParam.d_size(), sizeof(uint), cudaMemcpyDeviceToHost));
//			/*if (itr > 1000) {
//				vector<uint> test;
//				isEnds.copyToHost(test);
//				for (int i = 0; i < test.size(); i++) {
//					if (test[i] < 4u && test[i] > 0u)
//						printf("%d %d, %d\n", ceParam.h_size, i, test[i]);
//				}
//				vector<ContactElem> test2;
//				ceParam._elems.copyToHost(test2);
//				printf("%d %d %d %d\n", test2[0]._i[0], test2[0]._i[1], test2[0]._i[2], test2[0]._i[3]);
//				vector<uint2> ctest;
//				vector<uint> ictest;
//				NodeCEs.copyToHost(ctest);
//				iNodeCEs.copyToHost(ictest);
//				for (int n = 0; n < 4; n++) {
//					printf("asdfasdf %d %d %d\n", test2[0]._i[n], ictest[test2[0]._i[n]], ictest[test2[0]._i[n] + 1]);
//					for (int i = ictest[test2[0]._i[n]]; i < ictest[test2[0]._i[n] + 1]; i++) {
//						if (ctest[i].y == 0)
//							printf("asdf %d\n", i);
//					}
//				}
//			}*/
//		} while (isApplied);
//		//printf("Iteration %d\n", itr);

		resolveSelfCollisionProximity_kernel << < 1, 1 >> > (
			clothParams, ceParam.param(), delta, dt);
		CUDA_CHECK(cudaPeekAtLastError());

#ifdef COLLISION_TESTTIMER
		CUDA_CHECK(cudaDeviceSynchronize());
		printf("resolveSelfCollisionProximity_kernel: %lf msec\n", (CNOW - timer) / 10000.0);
		timer = CNOW;
#endif
	}

#ifdef COLLISION_TESTTIMER
	CUDA_CHECK(cudaDeviceSynchronize());
	printf("Resolve Self Collision Proximity: %lf msec\n", (CNOW - timer) / 10000.0);
#endif
}
void SelfCollisionSolver::ResolveSelfCollisionCCD(
	ContactElems& ceParam,
	const ObjParam& clothParams,
	BVHParam& clothBvh, RTriParam& RTri,
	Dvector<uint2>& NodeCEs,
	Dvector<uint>& iNodeCEs,
	Dvector<uint>& icurrs,
	Dvector<uint>& isEnds,
	const REAL delta, const REAL dt)
{
#ifdef COLLISION_TESTTIMER
	CUDA_CHECK(cudaDeviceSynchronize());
	ctimer timer = CNOW;
#endif

	ceParam.d_size.memset(0);
	getNumSelfCECCD_kernel << <divup(clothParams._numFaces, BLOCKSIZE), BLOCKSIZE >> > (
		clothParams, clothBvh, RTri, ceParam.d_size(), delta, dt);
	CUDA_CHECK(cudaPeekAtLastError());
	CUDA_CHECK(cudaMemcpy(&ceParam.h_size, ceParam.d_size(), sizeof(uint), cudaMemcpyDeviceToHost));

#ifdef COLLISION_TESTTIMER
	CUDA_CHECK(cudaDeviceSynchronize());
	printf("getNumSelfCECCD_kernel: %lf msec\n", (CNOW - timer) / 10000.0);
	timer = CNOW;
#endif

	if (ceParam.h_size > 0u) {
		if (ceParam._elems.size() < ceParam.h_size)
			ceParam.resize();

		ceParam.d_size.memset(0);
		getSelfCECCD_kernel << <divup(clothParams._numFaces, BLOCKSIZE), BLOCKSIZE >> > (
			clothParams, clothBvh, RTri, ceParam.param(), delta, dt);
		CUDA_CHECK(cudaPeekAtLastError());

#ifdef COLLISION_TESTTIMER
		CUDA_CHECK(cudaDeviceSynchronize());
		printf("getSelfCECCD_kernel: %lf msec\n", (CNOW - timer) / 10000.0);
		timer = CNOW;
#endif

		thrust::sort(thrust::device_ptr<ContactElem>(ceParam._elems.begin()),
			thrust::device_ptr<ContactElem>(ceParam._elems.begin() + ceParam.h_size), ContactElem_CMP());
		/*{
			vector<ContactElem> test;
			ceParam._elems.copyToHost(test);
			for (int i = 0; i < (int)ceParam.h_size - 1; i++) {
				for (int j = i + 1; j < (int)ceParam.h_size; j++) {
					if (test[i]._i[0] == test[j]._i[0] && test[i]._i[1] == test[j]._i[1] &&
						test[i]._i[2] == test[j]._i[2] && test[i]._i[3] == test[j]._i[3])
						printf("asdfasdf\n");
				}
			}
		}*/

#ifdef COLLISION_TESTTIMER
		CUDA_CHECK(cudaDeviceSynchronize());
		printf("thrust::sort: %lf msec\n", (CNOW - timer) / 10000.0);
		timer = CNOW;
#endif

		uint numNodeCEs = ceParam.h_size << 2u;
		/*Dvector<uint2> NodeCEs(ceParam.h_size << 2u);
		Dvector<uint> iNodeCEs(clothParams._numNodes + 1u);
		Dvector<uint> icurrs(clothParams._numNodes);
		Dvector<uint> isEnds(ceParam.h_size);*/
		if (NodeCEs.size() < numNodeCEs)
			NodeCEs.resize(ceParam.h_size << 2u);
		iNodeCEs.resize(clothParams._numNodes + 1u);
		icurrs.resize(clothParams._numNodes);
		if (isEnds.size() < ceParam.h_size)
			isEnds.resize(ceParam.h_size);
		iNodeCEs.memset(0);
		isEnds.memset(0);
		{
			getNodeCEs_kernel << <divup(ceParam.h_size, BLOCKSIZE), BLOCKSIZE >> > (
				ceParam.param(), NodeCEs());
			CUDA_CHECK(cudaPeekAtLastError());

#ifdef COLLISION_TESTTIMER
			CUDA_CHECK(cudaDeviceSynchronize());
			printf("getNodeCEs_kernel: %lf msec\n", (CNOW - timer) / 10000.0);
			timer = CNOW;
#endif

			thrust::sort(thrust::device_ptr<uint2>(NodeCEs.begin()),
				thrust::device_ptr<uint2>(NodeCEs.begin() + numNodeCEs), uint2_CMP());

#ifdef COLLISION_TESTTIMER
			CUDA_CHECK(cudaDeviceSynchronize());
			printf("thrust::sort: %lf msec\n", (CNOW - timer) / 10000.0);
			timer = CNOW;
#endif

			reorderNodeCEIds_kernel << <divup(numNodeCEs, BLOCKSIZE), BLOCKSIZE, (BLOCKSIZE + 1u) * sizeof(uint) >> > (
				NodeCEs(), iNodeCEs(), clothParams._numNodes, numNodeCEs);
			CUDA_CHECK(cudaPeekAtLastError());

#ifdef COLLISION_TESTTIMER
			CUDA_CHECK(cudaDeviceSynchronize());
			printf("reorderNodeCEIds_kernel: %lf msec\n", (CNOW - timer) / 10000.0);
			timer = CNOW;
#endif
		}

		uint isApplied;
		//uint itr = 0u;

		initDepth_kernel << < divup(clothParams._numNodes, BLOCKSIZE), BLOCKSIZE >> > (
			NodeCEs(), iNodeCEs(), icurrs(), isEnds(), clothParams._numNodes);
		CUDA_CHECK(cudaPeekAtLastError());
		do {
			//itr++;
			resolveSelfCollisionCCD_kernel << < divup(ceParam.h_size, BLOCKSIZE), BLOCKSIZE >> > (
				clothParams, ceParam.param(), isEnds(), delta, dt);
			CUDA_CHECK(cudaPeekAtLastError());

			ceParam.d_size.memset(0);
			nextDepth_kernel << < divup(clothParams._numNodes, BLOCKSIZE), BLOCKSIZE >> > (
				NodeCEs(), iNodeCEs(), icurrs(), isEnds(), clothParams._numNodes, ceParam.d_size());
			CUDA_CHECK(cudaPeekAtLastError());
			CUDA_CHECK(cudaMemcpy(&isApplied, ceParam.d_size(), sizeof(uint), cudaMemcpyDeviceToHost));
		} while (isApplied);
		//printf("Iteration %d\n", itr);

		/*resolveSelfCollisionCCD_kernel << < 1, 1 >> > (
			clothParams, ceParam.param(), delta, dt);
		CUDA_CHECK(cudaPeekAtLastError());*/
#ifdef COLLISION_TESTTIMER
		CUDA_CHECK(cudaDeviceSynchronize());
		printf("resolveSelfCollisionCCD_kernel: %lf msec\n", (CNOW - timer) / 10000.0);
		timer = CNOW;
#endif
	}

#ifdef COLLISION_TESTTIMER
	CUDA_CHECK(cudaDeviceSynchronize());
	printf("Resolve Self Collision CCD: %lf msec\n", (CNOW - timer) / 10000.0);
#endif
}
//uint SelfCollisionSolver::ResolveRigidImpactZone(
//	ContactElemParam& ceParam,
//	const ObjParam& clothParams,
//	const DPrefixArray<uint>& clothEdges,
//	BVH* clothBvh,
//	const REAL delta, const REAL dt)
//{
//#ifdef COLLISION_TESTTIMER
//	CUDA_CHECK(cudaDeviceSynchronize());
//	ctimer timer = CNOW;
//#endif
//	CUDA_CHECK(cudaMemset(ceParam.d_size, 0, sizeof(uint)));
//	getSelfCECCD_kernel_NoBuffer << <divup(clothParams._numFaces, BLOCKSIZE), BLOCKSIZE >> >
//		(clothParams, clothBvh->d_tree, ceParam, delta, dt);
//	CUDA_CHECK(cudaPeekAtLastError());
//	CUDA_CHECK(cudaMemcpy(&ceParam.h_size, ceParam.d_size, sizeof(uint), cudaMemcpyDeviceToHost));
//
//	uint ceSize = ceParam.h_size;
//	if (ceSize > 0)
//		ApplyRigidImpactZone(clothParams, ceParam, clothEdges);
//#ifdef COLLISION_TESTTIMER
//	CUDA_CHECK(cudaDeviceSynchronize());
//	printf("Resolve Self Collision RIZ: %lf msec\n", (CNOW - timer) / 10000.0);
//#endif
//	return ceSize;
//}
//-------------------------------------------------------------------------
bool SelfCollisionSolver::ResolveSelfCollision(
	const ObjParam& clothParams,
	const DPrefixArray<uint>& clothEdges,
	BVH* clothBvh, RTriangle* RTri,
	const REAL delta, const REAL dt)
{
	bool result = false;
	ContactElems ceParam;
	Dvector<uint2> NodeCEs;
	Dvector<uint> iNodeCEs;
	Dvector<uint> icurrs;
	Dvector<uint> isEnds;
	uint num = 0u;
	//----------------------------------------------------------------
	printf("=================================\n");
	clothBvh->refit(clothParams, delta, dt, false);
	ResolveSelfCollisionProximity(
		ceParam,clothParams, clothBvh->param(), RTri->param(),
		NodeCEs, iNodeCEs, icurrs, isEnds, delta, dt);
	if (ceParam.h_size > 0u)
		result = true;
	printf("Self Proximity: %d\n", ceParam.h_size);
	//----------------------------------------------------------------
	printf("-------------------\n");
	for (int itr = 0; itr < CCD_ITERATION; itr++) {
		clothBvh->refit(clothParams, COL_CULLING, dt, true);
		ResolveSelfCollisionCCD(
			ceParam, clothParams, clothBvh->param(), RTri->param(),
			NodeCEs, iNodeCEs, icurrs, isEnds, delta, dt);
		if (ceParam.h_size == 0u)
			return result;
		printf("Self CCD %d: %d\n", itr, ceParam.h_size);
		result = true;
	}
	return result;

	////----------------------------------------------------------------
	//printf("-------------------\n");
	//for (int itr = 0; itr < RIGID_ITERATION; itr++) {
	//	clothBvh->refit(clothParams, COL_CULLING, dt, true);
	//	CEnum = ResolveRigidImpactZone(
	//		ceParam, clothParams, clothEdges, clothBvh, delta, dt);
	//	if (CEnum == 0)
	//		break;
	//	printf("Self RIZ %d: %d\n", itr, CEnum);
	//	result = true;
	//	if (itr == 199)
	//		exit(0);
	//}
	//ParamManager::destroyContactElemParamDevice(ceParam);
	//return result;
}