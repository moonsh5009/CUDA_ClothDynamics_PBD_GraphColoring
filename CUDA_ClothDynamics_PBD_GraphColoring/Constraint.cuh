#include "Constraint.h"
#include "DeviceManager.cuh"

__global__ void getEdgeNeisIds_kernel(
	uint* nbNs, uint* inbNs, uint* ineis, 
	uint* ins, uint numSprings)
{
	uint id = blockDim.x * blockIdx.x + threadIdx.x;
	if (id >= numSprings)
		return;

	/*uint2 e = es[id];
	uint istart, iend;
	uint num = 0u;

	istart = inbNs[e.x];
	iend = inbNs[e.x + 1u];
	iend -= istart;
	if (iend)
		num += iend - 1u;

	istart = inbNs[e.y];
	iend = inbNs[e.y + 1u];
	iend -= istart;
	if (iend)
		num += iend - 1u;*/

	uint ino = id << 1u;
	uint ino0 = ins[ino + 0u];
	uint ino1 = ins[ino + 1u];
	uint istart, iend, estart, eend;
	uint eno;
	uint tmp;
	uint num = 0u;

	istart = inbNs[ino0];
	iend = inbNs[ino0 + 1u];
	for (ino = istart; ino < iend; ino++) {
		tmp = nbNs[ino];
		if (tmp >= ino1)
			break;
		num++;
	}

	istart = inbNs[ino1];
	iend = inbNs[ino1 + 1u];
	for (ino = istart; ino < iend; ino++) {
		tmp = nbNs[ino];
		if (tmp >= ino0)
			break;
		num++;
	}

	if (id == 0u)
		ineis[0] = 0u;

	ineis[id + 1u] = num;
}
__global__ void getFaceNeisIds_kernel(
	uint* nbFs, uint* inbFs, uint* ineis,
	uint* ins, uint numSprings)
{
	uint id = blockDim.x * blockIdx.x + threadIdx.x;
	if (id >= numSprings)
		return;

	uint inos[3], jnos[3];
	uint ino = id * 3u;
	inos[0] = ins[ino + 0u];
	inos[1] = ins[ino + 1u];
	inos[2] = ins[ino + 2u];

	uint istart, iend;
	uint tmp, n, i, j;
	uint num = 0u;
	uint flag;
	for (n = 0u; n < 3u; n++) {
		istart = inbFs[inos[n]];
		iend = inbFs[inos[n] + 1u];
		for (ino = istart; ino < iend; ino++) {
			tmp = nbFs[ino];
			if (tmp < id) {
				flag = 1u;
				if (n > 0u) {
					tmp *= 3u;
					jnos[0] = ins[tmp + 0u];
					jnos[1] = ins[tmp + 1u];
					jnos[2] = ins[tmp + 2u];
					for (i = 0u; i < n; i++) {
						for (j = 0u; j < 3u; j++)
							if (inos[i] == jnos[j])
								break;

						if (j < 3u) {
							flag = 0u;
							break;
						}
					}
				}
				num += flag;
			}
		}
	}

	if (id == 0u)
		ineis[0] = 0u;

	ineis[id + 1u] = num;
}

__global__ void getEdgeNeis_kernel(
	uint* es, uint* ies, uint* nbNs, uint* inbNs, 
	uint* neis, uint* ineis, uint*ins, uint numSprings)
{
	uint id = blockDim.x * blockIdx.x + threadIdx.x;
	if (id >= numSprings)
		return;

	/*uint2 e = es[id], tmp;
	uint ns[2];
	ns[0] = e.x; ns[1] = e.y;

	uint enos[2], inos[2], iends[2];
	uint estart, eend;
	uint e0, e1, eno, n;
	uint pos = ineis[id];

	for (n = 0; n < 2; n++) {
		inos[n] = inbNs[ns[n]];
		iends[n] = inbNs[ns[n] + 1u];
		if (inos[n] >= iends[n]) {
			enos[n] = 0xffffffff;
			continue;
		}
		e0 = ns[n];
		e1 = nbNs[inos[n]];
		if (e1 == ns[n ^ 1u]) {
			if (++inos[n] >= iends[n]) {
				enos[n] = 0xffffffff;
				continue;
			}
			e1 = nbNs[inos[n]];
		}
		
		if (e0 > e1) {
			e0 = e1;
			e1 = ns[n];
		}
		estart = ies[e0];
		eend = ies[e0 + 1u];
		for (eno = estart; eno < eend; eno++) {
			tmp = es[eno];
			if (tmp.y == e1) {
				enos[n] = eno;
				break;
			}
		}
	}
	while (enos[0] != 0xffffffff || enos[1] != 0xffffffff) {
		n = enos[0] > enos[1];
		neis[pos++] = enos[n];

		if (++inos[n] >= iends[n]) {
			enos[n] = 0xffffffff;
			continue;
		}
		e0 = ns[n];
		e1 = nbNs[inos[n]];
		if (e1 == ns[n ^ 1u]) {
			if (++inos[n] >= iends[n]) {
				enos[n] = 0xffffffff;
				continue;
			}
			e1 = nbNs[inos[n]];
		}
		if (e0 > e1) {
			e0 = e1;
			e1 = ns[n];
		}
		estart = ies[e0];
		eend = ies[e0 + 1u];
		for (eno = estart; eno < eend; eno++) {
			tmp = es[eno];
			if (tmp.y == e1) {
				enos[n] = eno;
				break;
			}
		}
	}*/
	uint ns[2];
	uint ino = id << 1u;
	ns[0] = ins[ino + 0u];
	ns[1] = ins[ino + 1u];

	uint enos[2], inos[2], iends[2];
	uint estart, eend;
	uint e0, e1, eno, tmp, n;
	uint pos = ineis[id];

	for (n = 0u; n < 2u; n++) {
		inos[n] = inbNs[ns[n]];
		iends[n] = inbNs[ns[n] + 1u];
		if (inos[n] >= iends[n]) {
			enos[n] = 0xffffffff;
			continue;
		}
		e0 = ns[n];
		e1 = nbNs[inos[n]];
		if (e1 >= ns[n ^ 1u]) {
			enos[n] = 0xffffffff;
			continue;
		}

		if (e0 > e1) {
			e0 = e1;
			e1 = ns[n];
		}
		estart = ies[e0];
		eend = ies[e0 + 1u];
		for (eno = estart; eno < eend; eno++) {
			tmp = es[eno];
			if (tmp == e1) {
				enos[n] = eno;
				break;
			}
		}
	}
	while (enos[0] != 0xffffffff || enos[1] != 0xffffffff) {
		n = enos[0] > enos[1];
		neis[pos++] = enos[n];

		if (++inos[n] >= iends[n]) {
			enos[n] = 0xffffffff;
			continue;
		}
		e0 = ns[n];
		e1 = nbNs[inos[n]];
		if (e1 >= ns[n ^ 1u]) {
			enos[n] = 0xffffffff;
			continue;
		}
		if (e0 > e1) {
			e0 = e1;
			e1 = ns[n];
		}
		estart = ies[e0];
		eend = ies[e0 + 1u];
		for (eno = estart; eno < eend; eno++) {
			tmp = es[eno];
			if (tmp == e1) {
				enos[n] = eno;
				break;
			}
		}
	}
}
__global__ void getFaceNeis_kernel(
	uint* nbFs, uint* inbFs, uint* neis, uint* ineis, uint* ins, uint numSprings)
{
	uint id = blockDim.x * blockIdx.x + threadIdx.x;
	if (id >= numSprings)
		return;

	uint ns[3], jnos[3];
	uint ino = id * 3u;
	ns[0] = ins[ino + 0u];
	ns[1] = ins[ino + 1u];
	ns[2] = ins[ino + 2u];

	uint fnos[3], inos[3], iends[3];
	uint tmp, n, i, j;
	uint pos = ineis[id];
	bool flag;

	for (n = 0u; n < 3u; n++) {
		inos[n] = inbFs[ns[n]];
		iends[n] = inbFs[ns[n] + 1u];
		do {
			if (inos[n] >= iends[n]) {
				fnos[n] = 0xffffffff;
				break;
			}
			fnos[n] = nbFs[inos[n]];
			flag = false;
			if (fnos[n] < id) {
				if (n > 0u) {
					tmp = fnos[n] * 3u;
					jnos[0] = ins[tmp + 0u];
					jnos[1] = ins[tmp + 1u];
					jnos[2] = ins[tmp + 2u];
					for (i = 0u; i < n; i++) {
						for (j = 0u; j < 3u; j++)
							if (ns[i] == jnos[j])
								break;

						if (j < 3u) {
							flag = true;
							break;
						}
					}
				}
			}
			else flag = true;

			inos[n]++;
		} while (flag);
	}
	while (fnos[0] != 0xffffffff || fnos[1] != 0xffffffff || fnos[2] != 0xffffffff) {
		n = 0u;
		if (fnos[n] > fnos[1]) n = 1u;
		if (fnos[n] > fnos[2]) n = 2u;

		neis[pos++] = fnos[n];

		do {
			if (inos[n] >= iends[n]) {
				fnos[n] = 0xffffffff;
				break;
			}
			fnos[n] = nbFs[inos[n]];
			flag = false;
			if (fnos[n] < id) {
				if (n > 0u) {
					tmp = fnos[n] * 3u;
					jnos[0] = ins[tmp + 0u];
					jnos[1] = ins[tmp + 1u];
					jnos[2] = ins[tmp + 2u];
					for (i = 0u; i < n; i++) {
						for (j = 0u; j < 3u; j++)
							if (ns[i] == jnos[j])
								break;

						if (j < 3u) {
							flag = true;
							break;
						}
					}
				}
			}
			else flag = true;

			inos[n]++;
		} while (flag);
	}
}

__global__ void getSubTreeParents_kernel(uint* neis, uint* ineis, uint* parents, uint numEdges)
{
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= numEdges)
		return;

	uint ino = ineis[id];
	uint iend = ineis[id + 1u];
	if (ino >= iend)
		return;
	uint nei, jno;
	for (; ino < iend; ino++) {
		nei = neis[ino];

	}
	/*
	uint prevCol, nei, a, b, ap;
	prevCol = neis[ino++];
	for (; ino < iend; ino++) {
		nei = neis[ino];
		a = prevCol;
		b = nei;
		while (a < b) {
			ap = atomicMin(parents + a, b);
			if (ap == 0xffffffff)
				break;
			if (ap > b) {
				a = b;
				b = ap;
			}
			else  a = ap;
		}
		prevCol = nei;
	}
	a = prevCol;
	b = id;
	while (a < b) {
		ap = atomicMin(parents + a, b);
		if (ap == 0xffffffff)
			break;
		if (ap > b) {
			a = b;
			b = ap;
		}
		else  a = ap;
	}*/
}
__global__ void reorderSubTree_kernel(uint* depths, uint* parents, uint* maxDepth, uint numEdges) {
	extern __shared__ uint s_maxDepth[];
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= numEdges) {
		s_maxDepth[threadIdx.x] = 0u;
		return;
	}

	uint parent = parents[id];
	uint depth = 0u;
	for (; parent != 0xffffffff; depth++)
		parent = parents[parent];
	depths[id] = depth;

	s_maxDepth[threadIdx.x] = depth;

	for (uint s = blockDim.x >> 1u; s > 32u; s >>= 1u) {
		__syncthreads();
		if (threadIdx.x < s)
			if (s_maxDepth[threadIdx.x] < s_maxDepth[threadIdx.x + s])
				s_maxDepth[threadIdx.x] = s_maxDepth[threadIdx.x + s];
	}
	__syncthreads();
	if (threadIdx.x < 32) {
		warpMax(s_maxDepth, threadIdx.x);
		if (threadIdx.x == 0)
			atomicMax(maxDepth, s_maxDepth[0]);
	}
}

__global__ void compColoring_kernel(uint* neis, uint* ineis, uint* depths, uint* colors, uint currDepth, uint numEdges) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= numEdges)
		return;

	uint depth = depths[id];
	if (depth != currDepth)
		return;

	uint istart = ineis[id];
	uint iend = ineis[id + 1u];
	uint ino, nei;
	uint color = 0u, icolor;
	bool flag;
	do {
		flag = false;
		for (ino = istart; ino < iend; ino++) {
			nei = neis[ino];
			icolor = colors[nei];
			if (color == icolor) {
				color++;
				flag = true;
				break;
			}
		}
	} while (flag);
	colors[id] = color;
}
__global__ void compColoring_kernel(uint* neis, uint* ineis, uint* icurrs, uint* colors, uint numSprings, uint* isApplied) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= numSprings)
		return;

	uint pos = icurrs[id];
	if (pos == 0xffffffff)
		return;

	uint istart = ineis[id];
	uint iend = ineis[id + 1u];
	uint ino, nei, icolor;

	for (ino = istart + pos; ino < iend; ino++) {
		nei = neis[ino];
		icolor = colors[nei];
		if (icolor == 0xffffffff)
			break;
		pos++;
	}

	if (pos == iend - istart) {
		uint color = 0u;
		bool flag;
		do {
			flag = false;
			for (ino = istart; ino < iend; ino++) {
				nei = neis[ino];
				icolor = colors[nei];
				if (color == icolor) {
					color++;
					flag = true;
					break;
				}
			}
		} while (flag);
		colors[id] = color;
		icurrs[id] = 0xffffffff;
	}
	else {
		icurrs[id] = pos;
		*isApplied = 1u;
	}
}

__global__ void initEdgeConstraintsIds_kernel(uint* es, uint* ies, EdgeSpring springs, uint numVertices) {
	uint ino0 = threadIdx.x + blockDim.x * blockIdx.x, ino1;
	if (ino0 >= numVertices)
		return;

	uint istart = ies[ino0];
	uint iend = ies[ino0 + 1u];
	uint i, ino;
	for (i = istart; i < iend; i++) {
		ino1 = es[i];
		ino = i << 1u;
		springs._inos[ino + 0u] = ino0;
		springs._inos[ino + 1u] = ino1;
	}
}
__global__ void initEdgeConstraints_kernel(REAL* ns, EdgeSpring springs) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= springs._numSprings)
		return;

	uint ino = id << 1u;
	uint ino0 = springs._inos[ino + 0u];
	uint ino1 = springs._inos[ino + 1u];
	REAL3 v0, v1;
	ino0 *= 3u; ino1 *= 3u;
	v0.x = ns[ino0 + 0u]; v0.y = ns[ino0 + 1u]; v0.z = ns[ino0 + 2u];
	v1.x = ns[ino1 + 0u]; v1.y = ns[ino1 + 1u]; v1.z = ns[ino1 + 2u];

	REAL length = Length(v1 - v0);
	springs._rs[id] = length;
}

__device__ void transform_device(REAL* R, REAL3& v0, REAL3& v1, REAL3& v2) {
	REAL c = (REAL)1.0 - v0.z * v0.z;
	if (c > (REAL)1.0e-10) {
		c = sqrt(c);
		REAL invC = (REAL)1.0 / c;

		R[0] = v0.y * invC;
		R[1] = -v0.x * invC;
		R[2] = -R[1] * v0.z;
		R[3] = R[0] * v0.z;
		v1 = make_REAL3(
			v1.x * R[0] + v1.y * R[1],
			v1.x * R[2] + v1.y * R[3] + v1.z * -c, 0.);
		v2 = make_REAL3(
			v2.x * R[0] + v2.y * R[1], 
			v2.x * R[2] + v2.y * R[3] + v2.z * -c, 0.);
	}
}
__global__ void initSBConstraints_kernel(REAL* ns, SBSpring springs, REAL k) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= springs._numSprings)
		return;

	uint ino = id * 3u;
	uint ino0 = springs._inos[ino + 0u];
	uint ino1 = springs._inos[ino + 1u];
	uint ino2 = springs._inos[ino + 2u];

	REAL3 v0, v1, v2;
	ino0 *= 3u; ino1 *= 3u; ino2 *= 3u;
	v0.x = ns[ino0 + 0u]; v0.y = ns[ino0 + 1u]; v0.z = ns[ino0 + 2u];
	v1.x = ns[ino1 + 0u]; v1.y = ns[ino1 + 1u]; v1.z = ns[ino1 + 2u];
	v2.x = ns[ino2 + 0u]; v2.y = ns[ino2 + 1u]; v2.z = ns[ino2 + 2u];
	
	v1 = v1 - v0;
	v2 = v2 - v0;
	v0 = Cross(v1, v2);
	REAL area = LengthSquared(v0);
	v0 *= 1.0 / sqrt(area);

	REAL R[4];
	transform_device(R, v0, v1, v2);

	REAL inv = (REAL)1.0 / (v1.x * v2.y - v1.y * v2.x);
	R[0] = v2.y * inv;
	R[1] = -v1.y * inv;
	R[2] = -v2.x * inv;
	R[3] = v1.x * inv;

	springs._as[id] = area;
	springs._ks[id] = k;
	ino = id << 2u;
	springs._invQs[ino + 0u] = R[0];
	springs._invQs[ino + 1u] = R[1];
	springs._invQs[ino + 2u] = R[2];
	springs._invQs[ino + 3u] = R[3];
}

__global__ void project_kernel(REAL* ns, REAL* invMs, EdgeSpring springs, REAL invdt2, uint currColor) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= springs._numSprings)
		return;

	uint color = springs._colors[id];
	if (color != currColor)
		return;

	uint ino = id << 1u;
	uint ino0 = springs._inos[ino + 0u];
	uint ino1 = springs._inos[ino + 1u];

	REAL w0 = invMs[ino0];
	REAL w1 = invMs[ino1];

	ino0 *= 3u; ino1 *= 3u;
	REAL3 v0, v1;
	v0.x = ns[ino0 + 0u]; v0.y = ns[ino0 + 1u]; v0.z = ns[ino0 + 2u];
	v1.x = ns[ino1 + 0u]; v1.y = ns[ino1 + 1u]; v1.z = ns[ino1 + 2u];

	REAL material = springs._material;
	springs._material *= invdt2;
	REAL restLength = springs._rs[id];
	REAL lambda = springs._lambdas[id];

	REAL3 dir = v1 - v0;
	REAL length = Length(dir);
	REAL constraint = length - restLength;
	REAL dt_lambda = (-constraint - material * lambda) / ((w0 + w1) + material);
	dir = dt_lambda / (length + FLT_EPSILON) * dir;
	lambda += dt_lambda;

	v0 -= w0 * dir;
	v1 += w1 * dir;
	springs._lambdas[id] = lambda;
	ns[ino0 + 0u] = v0.x; ns[ino0 + 1u] = v0.y; ns[ino0 + 2u] = v0.z;
	ns[ino1 + 0u] = v1.x; ns[ino1 + 1u] = v1.y; ns[ino1 + 2u] = v1.z;
}
__global__ void project_kernel(REAL* ns, REAL* invMs, SBSpring springs, uint currColor) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= springs._numSprings)
		return;

	uint color = springs._colors[id];
	if (color != currColor)
		return;

	uint ino0, ino1, ino2;
	uint ino = id * 3u;
	ino0 = springs._inos[ino + 0u];
	ino1 = springs._inos[ino + 1u];
	ino2 = springs._inos[ino + 2u];

	REAL w0 = invMs[ino0];
	REAL w1 = invMs[ino1];
	REAL w2 = invMs[ino2];
	REAL k = springs._ks[id];

	REAL3 v0, v1, v2, force0, force1, force2;
	ino0 *= 3u; ino1 *= 3u; ino2 *= 3u;
	v0.x = ns[ino0 + 0u]; v0.y = ns[ino0 + 1u]; v0.z = ns[ino0 + 2u];
	v1.x = ns[ino1 + 0u]; v1.y = ns[ino1 + 1u]; v1.z = ns[ino1 + 2u];
	v2.x = ns[ino2 + 0u]; v2.y = ns[ino2 + 1u]; v2.z = ns[ino2 + 2u];

	REAL3 p1 = v1 - v0;
	REAL3 p2 = v2 - v0;

	REAL p1Dotp1 = Dot(p1, p1);
	REAL p1Dotp2 = Dot(p1, p2);
	REAL p2Dotp2 = Dot(p2, p2);

	REAL restArea = springs._as[id];
	REAL area = p1Dotp1 * p2Dotp2 - p1Dotp2 * p1Dotp2;
	REAL s = area - restArea;
	REAL3 d1 = p2Dotp2 * p1 - p1Dotp2 * p2; d1 += d1;
	REAL3 d2 = p1Dotp1 * p2 - p1Dotp2 * p1; d2 += d2;
	REAL3 d0 = make_REAL3(-d1.x - d2.x, -d1.y - d2.y, -d1.z - d2.z);
	REAL lambda = -s / (
		LengthSquared(d0) * w0 +
		LengthSquared(d1) * w1 +
		LengthSquared(d2) * w2 + FLT_EPSILON);

	force0 = lambda * d0;
	force1 = lambda * d1;
	force2 = lambda * d2;

	REAL Q[4];
	ino = id << 2u;
	Q[0] = springs._invQs[ino + 0u];
	Q[1] = springs._invQs[ino + 1u];
	Q[2] = springs._invQs[ino + 2u];
	Q[3] = springs._invQs[ino + 3u];

	REAL3 f1 = p1 * Q[0] + p2 * Q[1];
	REAL3 f2 = p1 * Q[2] + p2 * Q[3];

	s = Dot(f1, f2);
	d1 = f1 * Q[2] + f2 * Q[0];
	d2 = f1 * Q[3] + f2 * Q[1];
	d0 = make_REAL3(-d1.x - d2.x, -d1.y - d2.y, -d1.z - d2.z);
	lambda = -s / (
		LengthSquared(d0) * w0 +
		LengthSquared(d1) * w1 +
		LengthSquared(d2) * w2 + FLT_EPSILON);

	force0 += lambda * d0;
	force1 += lambda * d1;
	force2 += lambda * d2;

	s = Dot(f1, f1) - 1.0;
	d1 = f1 * Q[0]; d1 += d1;
	d2 = f1 * Q[1]; d2 += d2;
	d0 = make_REAL3(-d1.x - d2.x, -d1.y - d2.y, -d1.z - d2.z);
	lambda = -s / (
		LengthSquared(d0) * w0 +
		LengthSquared(d1) * w1 +
		LengthSquared(d2) * w2 + FLT_EPSILON);

	force0 += lambda * d0;
	force1 += lambda * d1;
	force2 += lambda * d2;

	s = Dot(f2, f2) - 1.0;
	d1 = f2 * Q[2]; d1 += d1;
	d2 = f2 * Q[3]; d2 += d2;
	d0 = make_REAL3(-d1.x - d2.x, -d1.y - d2.y, -d1.z - d2.z);
	lambda = -s / (
		LengthSquared(d0) * w0 +
		LengthSquared(d1) * w1 +
		LengthSquared(d2) * w2 + FLT_EPSILON);

	force0 += lambda * d0;
	force1 += lambda * d1;
	force2 += lambda * d2;

	v0 += w0 * k * force0;
	v1 += w1 * k * force1;
	v2 += w2 * k * force2;

	ns[ino0 + 0u] = v0.x; ns[ino0 + 1u] = v0.y; ns[ino0 + 2u] = v0.z;
	ns[ino1 + 0u] = v1.x; ns[ino1 + 1u] = v1.y; ns[ino1 + 2u] = v1.z;
	ns[ino2 + 0u] = v2.x; ns[ino2 + 1u] = v2.y; ns[ino2 + 2u] = v2.z;
}
__global__ void project_strain_kernel(REAL* ns, REAL* invMs, SBSpring springs, uint currColor) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= springs._numSprings)
		return;

	uint color = springs._colors[id];
	if (color != currColor)
		return;

	uint ino0, ino1, ino2;
	uint ino = id * 3u;
	ino0 = springs._inos[ino + 0u];
	ino1 = springs._inos[ino + 1u];
	ino2 = springs._inos[ino + 2u];

	REAL w0 = invMs[ino0];
	REAL w1 = invMs[ino1];
	REAL w2 = invMs[ino2];
	REAL k = springs._ks[id];

	REAL3 v0, v1, v2;
	ino0 *= 3u; ino1 *= 3u; ino2 *= 3u;
	v0.x = ns[ino0 + 0u]; v0.y = ns[ino0 + 1u]; v0.z = ns[ino0 + 2u];
	v1.x = ns[ino1 + 0u]; v1.y = ns[ino1 + 1u]; v1.z = ns[ino1 + 2u];
	v2.x = ns[ino2 + 0u]; v2.y = ns[ino2 + 1u]; v2.z = ns[ino2 + 2u];

	REAL3 p1 = v1 - v0;
	REAL3 p2 = v2 - v0;

	REAL3 force0, force1, force2;
	REAL3 d0, d1, d2;
	REAL s, lambda;

	REAL Q[4];
	ino = id << 2u;
	Q[0] = springs._invQs[ino + 0u];
	Q[1] = springs._invQs[ino + 1u];
	Q[2] = springs._invQs[ino + 2u];
	Q[3] = springs._invQs[ino + 3u];

	REAL3 f1 = p1 * Q[0] + p2 * Q[1];
	REAL3 f2 = p1 * Q[2] + p2 * Q[3];

	s = Dot(f1, f2);
	d1 = f1 * Q[2] + f2 * Q[0];
	d2 = f1 * Q[3] + f2 * Q[1];
	d0 = make_REAL3(-d1.x - d2.x, -d1.y - d2.y, -d1.z - d2.z);
	lambda = -s / (
		LengthSquared(d0) * w0 +
		LengthSquared(d1) * w1 +
		LengthSquared(d2) * w2 + FLT_EPSILON);

	force0 = lambda * d0;
	force1 = lambda * d1;
	force2 = lambda * d2;

	s = Dot(f1, f1) - 1.0;
	d1 = f1 * Q[0]; d1 += d1;
	d2 = f1 * Q[1]; d2 += d2;
	d0 = make_REAL3(-d1.x - d2.x, -d1.y - d2.y, -d1.z - d2.z);
	lambda = -s / (
		LengthSquared(d0) * w0 +
		LengthSquared(d1) * w1 +
		LengthSquared(d2) * w2 + FLT_EPSILON);

	force0 += lambda * d0;
	force1 += lambda * d1;
	force2 += lambda * d2;

	s = Dot(f2, f2) - 1.0;
	d1 = f2 * Q[2]; d1 += d1;
	d2 = f2 * Q[3]; d2 += d2;
	d0 = make_REAL3(-d1.x - d2.x, -d1.y - d2.y, -d1.z - d2.z);
	lambda = -s / (
		LengthSquared(d0) * w0 +
		LengthSquared(d1) * w1 +
		LengthSquared(d2) * w2 + FLT_EPSILON);

	force0 += lambda * d0;
	force1 += lambda * d1;
	force2 += lambda * d2;

	v0 += w0 * k * force0;
	v1 += w1 * k * force1;
	v2 += w2 * k * force2;

	ns[ino0 + 0u] = v0.x; ns[ino0 + 1u] = v0.y; ns[ino0 + 2u] = v0.z;
	ns[ino1 + 0u] = v1.x; ns[ino1 + 1u] = v1.y; ns[ino1 + 2u] = v1.z;
	ns[ino2 + 0u] = v2.x; ns[ino2 + 1u] = v2.y; ns[ino2 + 2u] = v2.z;
}
__global__ void project_area_kernel(REAL* ns, REAL* invMs, SBSpring springs, uint currColor) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= springs._numSprings)
		return;

	uint color = springs._colors[id];
	if (color != currColor)
		return;

	uint ino0, ino1, ino2;
	uint ino = id * 3u;
	ino0 = springs._inos[ino + 0u];
	ino1 = springs._inos[ino + 1u];
	ino2 = springs._inos[ino + 2u];

	REAL w0 = invMs[ino0];
	REAL w1 = invMs[ino1];
	REAL w2 = invMs[ino2];
	REAL k = springs._ks[id];

	REAL3 v0, v1, v2;
	ino0 *= 3u; ino1 *= 3u; ino2 *= 3u;
	v0.x = ns[ino0 + 0u]; v0.y = ns[ino0 + 1u]; v0.z = ns[ino0 + 2u];
	v1.x = ns[ino1 + 0u]; v1.y = ns[ino1 + 1u]; v1.z = ns[ino1 + 2u];
	v2.x = ns[ino2 + 0u]; v2.y = ns[ino2 + 1u]; v2.z = ns[ino2 + 2u];

	REAL3 p1 = v1 - v0;
	REAL3 p2 = v2 - v0;

	REAL p1Dotp1 = Dot(p1, p1);
	REAL p1Dotp2 = Dot(p1, p2);
	REAL p2Dotp2 = Dot(p2, p2);

	REAL restArea = springs._as[id];
	REAL area = p1Dotp1 * p2Dotp2 - p1Dotp2 * p1Dotp2;
	REAL s = area - restArea;
	REAL3 d1 = p2Dotp2 * p1 - p1Dotp2 * p2; d1 += d1;
	REAL3 d2 = p1Dotp1 * p2 - p1Dotp2 * p1; d2 += d2;
	REAL3 d0 = make_REAL3(-d1.x - d2.x, -d1.y - d2.y, -d1.z - d2.z);
	REAL lambda = -s / (
		LengthSquared(d0) * w0 +
		LengthSquared(d1) * w1 +
		LengthSquared(d2) * w2 + FLT_EPSILON);

	w0 *= k * lambda;
	w1 *= k * lambda;
	w2 *= k * lambda;

	v0 += w0 * d0;
	v1 += w1 * d1;
	v2 += w2 * d2;

	ns[ino0 + 0u] = v0.x; ns[ino0 + 1u] = v0.y; ns[ino0 + 2u] = v0.z;
	ns[ino1 + 0u] = v1.x; ns[ino1 + 1u] = v1.y; ns[ino1 + 2u] = v1.z;
	ns[ino2 + 0u] = v2.x; ns[ino2 + 1u] = v2.y; ns[ino2 + 2u] = v2.z;
}