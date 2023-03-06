#include "ClothDynamics.h"
#include "DeviceManager.cuh"

__global__ void compExternalForce_kernel(REAL* forces, REAL* ms, REAL3 gravity, uint numVertices) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= numVertices)
		return;

	REAL m = ms[id];
	id *= 3u;

	REAL3 force;
	force.x = forces[id + 0u];
	force.y = forces[id + 1u];
	force.z = forces[id + 2u];

	force += m * gravity;
	forces[id + 0u] = force.x;
	forces[id + 1u] = force.y;
	forces[id + 2u] = force.z;
}

__global__ void compPredictPosition_kernel(REAL* ns, REAL* vels, REAL* forces, REAL* invMs, REAL dt, uint numVertices) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= numVertices)
		return;

	REAL invM = invMs[id];
	id *= 3u;

	REAL3 v, vel, force;
	v.x = ns[id + 0u];
	v.y = ns[id + 1u];
	v.z = ns[id + 2u];
	vel.x = vels[id + 0u];
	vel.y = vels[id + 1u];
	vel.z = vels[id + 2u];
	force.x = forces[id + 0u];
	force.y = forces[id + 1u];
	force.z = forces[id + 2u];

	vel += invM * dt * force; v += dt * vel;
	ns[id + 0u] = v.x;
	ns[id + 1u] = v.y;
	ns[id + 2u] = v.z;
	vels[id + 0u] = vel.x;
	vels[id + 1u] = vel.y;
	vels[id + 2u] = vel.z;
}
__global__ void updateVelocitiy_kernel(REAL* ns, REAL* v0s, REAL* vels, REAL invdt, uint numVertices) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= numVertices)
		return;

	id *= 3u;
	REAL3 v, v0, vel;
	v.x = ns[id + 0u];
	v.y = ns[id + 1u];
	v.z = ns[id + 2u];
	v0.x = v0s[id + 0u];
	v0.y = v0s[id + 1u];
	v0.z = v0s[id + 2u];
	vel.x = vels[id + 0u];
	vel.y = vels[id + 1u];
	vel.z = vels[id + 2u];

	vel = (v - v0) * invdt;
	vels[id + 0u] = vel.x;
	vels[id + 1u] = vel.y;
	vels[id + 2u] = vel.z;
}
__global__ void updatePosition_kernel(REAL* ns, REAL* vs, REAL dt, uint numVertices) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= numVertices)
		return;

	id *= 3u;
	REAL3 n, v;
	n.x = ns[id + 0u];
	n.y = ns[id + 1u];
	n.z = ns[id + 2u];
	v.x = vs[id + 0u];
	v.y = vs[id + 1u];
	v.z = vs[id + 2u];
	n += v * dt;

	ns[id + 0u] = n.x;
	ns[id + 1u] = n.y;
	ns[id + 2u] = n.z;
}

__global__ void compFnorms_kernel(uint* fs, REAL* ns, REAL* fNorms, uint numFaces) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= numFaces)
		return;

	id *= 3u;
	uint iv0 = fs[id + 0u]; iv0 *= 3u;
	uint iv1 = fs[id + 1u]; iv1 *= 3u;
	uint iv2 = fs[id + 2u]; iv2 *= 3u;

	REAL3 v0, v1, v2;
	v0.x = ns[iv0 + 0u]; v0.y = ns[iv0 + 1u]; v0.z = ns[iv0 + 2u];
	v1.x = ns[iv1 + 0u]; v1.y = ns[iv1 + 1u]; v1.z = ns[iv1 + 2u];
	v2.x = ns[iv2 + 0u]; v2.y = ns[iv2 + 1u]; v2.z = ns[iv2 + 2u];

	REAL3 norm = Cross(v1 - v0, v2 - v0);
	Normalize(norm);

	fNorms[id + 0u] = norm.x;
	fNorms[id + 1u] = norm.y;
	fNorms[id + 2u] = norm.z;
}