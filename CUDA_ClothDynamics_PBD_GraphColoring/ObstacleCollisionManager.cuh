//#include "CollisionManager.cuh"
//
//__device__ bool boundaryCollision(
//	int i,
//	const ObjectParams clothParams,
//	const double3 minBoun, const double3 maxBoun,
//	const double delta, const double friction, const double dt)
//{
//	bool apply = false;
//	if (clothParams.masses[i] <= 0.0)
//		return false;
//	double minb[3] = { minBoun.x, minBoun.y, minBoun.z };
//	double maxb[3] = { maxBoun.x, maxBoun.y, maxBoun.z };
//	double p[3] = { clothParams.vertices[i * 3 + 0], clothParams.vertices[i * 3 + 1], clothParams.vertices[i * 3 + 2] };
//	double vel[3] = { clothParams.velocities[i * 3 + 0], clothParams.velocities[i * 3 + 1], clothParams.velocities[i * 3 + 2] };
//	for (int j = 0; j < 3; j++) {
//		double N[3] = { 0.0, 0.0, 0.0 };
//		N[j] = 1.0;
//		double dist = p[j] - minb[j] - delta;
//		double newDist = dist + vel[j] * dt;
//		if (newDist < 0.0) {
//			double A_newVelN = vel[j] - newDist / dt;
//			double3 relVelT = make_double3(vel[0] - N[0] * vel[j], vel[1] - N[1] * vel[j], vel[2] - N[2] * vel[j]);
//			double3 newRelVelT = relVelT *
//				max(0.0, 1.0 - friction * (A_newVelN - vel[j]) / Length(relVelT));
//			vel[0] = (N[0] * A_newVelN + newRelVelT.x);
//			vel[1] = (N[1] * A_newVelN + newRelVelT.y);
//			vel[2] = (N[2] * A_newVelN + newRelVelT.z);
//			apply = true;
//		}
//		N[j] = -1.0;
//		dist = maxb[j] - p[j] - delta;
//		newDist = dist - vel[j] * dt;
//		if (newDist < 0.0) {
//			double A_newVelN = -vel[j] - newDist / dt;
//			double3 relVelT = make_double3(vel[0] + N[0] * vel[j], vel[1] + N[1] * vel[j], vel[2] + N[2] * vel[j]);
//			double3 newRelVelT = relVelT *
//				max(0.0, 1.0 - friction * (A_newVelN + vel[j]) / Length(relVelT));
//			vel[0] = (N[0] * A_newVelN + newRelVelT.x);
//			vel[1] = (N[1] * A_newVelN + newRelVelT.y);
//			vel[2] = (N[2] * A_newVelN + newRelVelT.z);
//			apply = true;
//		}
//	}
//	if (apply) {
//		clothParams.velocities[i * 3 + 0] = vel[0];
//		clothParams.velocities[i * 3 + 1] = vel[1];
//		clothParams.velocities[i * 3 + 2] = vel[2];
//	}
//	return apply;
//}
//
//__device__ __forceinline__ bool isObsContactTV_CCD(
//	const double3 p0, const double3 p1, const double3 p2, const double3 p3,
//	const double3 q0, const double3 q1, const double3 q2, const double3 q3,
//	double delta, double* t, double* w0 = nullptr, double* w1 = nullptr)
//{
//	if (DetectionTV_CCD(p0, p1, p2, p3, q0, q1, q2, q3, COL_TV, t, w0, w1)) {
//#ifdef CCD_PRINT_TRI_TV
//		printf("r tri0 %f %f %f\n", *w0, *w1, *t);
//#endif
//		return true;
//	}
//	if (isContactTV_Proximity(q0, q1, q2, q3, delta, w0, w1)) {
//		*t = -1.0;
//#ifdef CCD_PRINT_TRI_TV
//		printf("r tri0 %f %f %f\n", *w0, *w1, *t);
//#endif
//		return true;
//	}
//	return false;
//}
//__device__ __forceinline__ bool isObsContactEE_CCD(
//	const double3 p0, const double3 p1, const double3 p2, const double3 p3,
//	const double3 q0, const double3 q1, const double3 q2, const double3 q3,
//	double delta, double* t, double* w0 = nullptr, double* w1 = nullptr)
//{
//	if (DetectionEE_CCD(p0, p1, p2, p3, q0, q1, q2, q3, COL_EE, t, w0, w1)) {
//#ifdef CCD_PRINT_EDGE_EE
//		printf("r edge0 %f %f %f\n", *w0, *w1, *t);
//#endif
//		return true;
//	}
//	if (isContactEE_Proximity(q0, q1, q2, q3, delta, w0, w1)) {
//		*t = -1.0;
//#ifdef CCD_PRINT_TRI_TV
//		printf("r edge0 %f %f %f\n", *w0, *w1, *t);
//#endif
//		return true;
//	}
//	return false;
//}
//
//__device__ __forceinline__ void resolveImpulse_obs(
//	int i0,
//	const ImpulseParams impParams,
//	double3 imp, double pene)
//{
//	unsigned long long ret = __double_as_longlong(impParams.penetrations[i0]);
//	unsigned long long old;
//	while (pene > __longlong_as_double(ret))
//	{
//		old = ret;
//		if ((ret = atomicCAS((unsigned long long*)(impParams.penetrations + i0),
//			old, __double_as_longlong(pene))) == old)
//		{
//			ret = __double_as_longlong(impParams.impulses[i0 * 3 + 0]);
//			do {
//				if (pene < impParams.penetrations[i0]) return;
//				old = ret;
//			} while ((ret = atomicCAS((unsigned long long*)(impParams.impulses + i0 * 3 + 0),
//				old, __double_as_longlong(imp.x))) != old);
//
//			ret = __double_as_longlong(impParams.impulses[i0 * 3 + 1]);
//			do {
//				if (pene < impParams.penetrations[i0]) return;
//				old = ret;
//			} while ((ret = atomicCAS((unsigned long long*)(impParams.impulses + i0 * 3 + 1),
//				old, __double_as_longlong(imp.y))) != old);
//
//			ret = __double_as_longlong(impParams.impulses[i0 * 3 + 2]);
//			do {
//				if (pene < impParams.penetrations[i0]) return;
//				old = ret;
//			} while ((ret = atomicCAS((unsigned long long*)(impParams.impulses + i0 * 3 + 2),
//				old, __double_as_longlong(imp.z))) != old);
//
//			break;
//		}
//	}
//}
//
//__device__ bool resolveObstacleResponseTV(
//	int i0, int i1, int i2, int i3,
//	const double3& p0, const double3& p1, const double3& p2, const double3& p3,
//	const double3& v0, const double3& v1, const double3& v2, const double3& v3,
//	const ObjectParams clothParams,
//	const ObjectParams obsParams,
//	const ImpulseParams impParams,
//	const double delta, const double friction, const double dt)
//{
//	double3 q0 = p0 + v0 * dt;
//	double3 q1 = p1 + v1 * dt;
//	double3 q2 = p2 + v2 * dt;
//	double3 q3 = p3 + v3 * dt;
//	double t, w0, w1;
//	if (!isObsContactTV_CCD(p0, p1, p2, p3, q0, q1, q2, q3, delta, &t, &w0, &w1))
//		return false;
//
//	bool isProximity = t < 0.0;
//	if (isProximity)
//		t = 1.0;
//
//	double w2 = 1.0 - w0 - w1;
//	double ht = t * dt * COL_HALFTIME;
//	double3 q0m = p0 + v0 * ht;
//	double3 q1m = p1 + v1 * ht;
//	double3 q2m = p2 + v2 * ht;
//	double3 q3m = p3 + v3 * ht;
//	double3 norm = q0m * w0 + q1m * w1 + q2m * w2 - q3m;
//	Normalize(norm);
//
//	double3 vc = v0 * w0 + v1 * w1 + v2 * w2;
//	double dist = Dot(p0 * w0 + p1 * w1 + p2 * w2 - p3, norm) - delta;
//	double newDist = dist + Dot(vc - v3, norm) * dt;
//	if (newDist >= 0.0)
//		return true;
//	newDist /= -dt;
//
//	double vcn = Dot(vc, norm);
//	double v3n = Dot(v3, norm);
//	double3 vct = vc - norm * vcn;
//	double3 v3t = v3 - norm * v3n;
//	double3 relVT = vct - v3t;
//	double3 n_relVT = make_double3(0.0);
//	double l_relVT = Length(relVT);
//	if (l_relVT != 0.0)
//		n_relVT = relVT * max(1.0 - friction * newDist / l_relVT, 0.0);
//
//	double3 n_vct = v3t + n_relVT;
//	double3 imp = norm * newDist + n_vct - vct;
//	imp *= COL_HALF;
//	if (clothParams.masses[i0] > 0.0)
//		resolveImpulse_obs(i0, impParams, imp, newDist);
//	if (clothParams.masses[i1] > 0.0)
//		resolveImpulse_obs(i1, impParams, imp, newDist);
//	if (clothParams.masses[i2] > 0.0)
//		resolveImpulse_obs(i2, impParams, imp, newDist);
//	return true;
//}
//__device__ bool resolveObstacleResponseVT(
//	int i0, int i1, int i2, int i3,
//	const double3& p0, const double3& p1, const double3& p2, const double3& p3,
//	const double3& v0, const double3& v1, const double3& v2, const double3& v3,
//	const ObjectParams clothParams,
//	const ObjectParams obsParams,
//	const ImpulseParams impParams,
//	const double delta, const double friction, const double dt)
//{
//	double3 q0 = p0 + v0 * dt;
//	double3 q1 = p1 + v1 * dt;
//	double3 q2 = p2 + v2 * dt;
//	double3 q3 = p3 + v3 * dt;
//	double t, w0, w1;
//	if (!isObsContactTV_CCD(p0, p1, p2, p3, q0, q1, q2, q3, delta, &t, &w0, &w1))
//		return true;
//
//	bool isProximity = t < 0.0;
//	if (isProximity)
//		t = 1.0;
//
//	double w2 = 1.0 - w0 - w1;
//	double ht = t * dt * COL_HALFTIME;
//	double3 q0m = p0 + v0 * ht;
//	double3 q1m = p1 + v1 * ht;
//	double3 q2m = p2 + v2 * ht;
//	double3 q3m = p3 + v3 * ht;
//	double3 norm = q3m - q0m * w0 - q1m * w1 - q2m * w2;
//	Normalize(norm);
//
//	double3 vc = v0 * w0 + v1 * w1 + v2 * w2;
//	double dist = Dot(p3 - p0 * w0 - p1 * w1 - p2 * w2, norm) - delta;
//	double newDist = dist + Dot(v3 - vc, norm) * dt;
//	if (newDist >= 0.0)
//		return false;
//	newDist /= -dt;
//
//	double v3n = Dot(v3, norm);
//	double vcn = Dot(vc, norm);
//	double3 v3t = v3 - norm * v3n;
//	double3 vct = vc - norm * vcn;
//	double3 relVT = v3t - vct;
//	double3 n_relVT = make_double3(0.0);
//	double l_relVT = Length(relVT);
//	if (l_relVT != 0.0)
//		n_relVT = relVT * max(1.0 - friction * newDist / l_relVT, 0.0);
//
//	double3 n_v3t = vct + n_relVT;
//	double3 imp = norm * newDist + n_v3t - v3t;
//	imp *= COL_HALF;
//	if (clothParams.masses[i3] > 0.0)
//		resolveImpulse_obs(i3, impParams, imp, newDist);
//
//	return true;
//}
//__device__ bool resolveObstacleResponseEE(
//	int i0, int i1, int i2, int i3,
//	const double3& p0, const double3& p1, const double3& p2, const double3& p3,
//	const double3& v0, const double3& v1, const double3& v2, const double3& v3,
//	const ObjectParams clothParams,
//	const ObjectParams obsParams,
//	const ImpulseParams impParams,
//	const double delta, const double friction, const double dt)
//{
//	double3 q0 = p0 + v0 * dt;
//	double3 q1 = p1 + v1 * dt;
//	double3 q2 = p2 + v2 * dt;
//	double3 q3 = p3 + v3 * dt;
//	double t, w0, w1;
//	if (!isObsContactEE_CCD(p0, p1, p2, p3, q0, q1, q2, q3, delta, &t, &w0, &w1))
//		return true;
//
//	bool isProximity = t < 0.0;
//	if (isProximity)
//		t = 1.0;
//
//	double ht = t * dt * COL_HALFTIME;
//	double3 q0m = p0 + v0 * ht;
//	double3 q1m = p1 + v1 * ht;
//	double3 q2m = p2 + v2 * ht;
//	double3 q3m = p3 + v3 * ht;
//	double3 norm = q0m + (q1m - q0m) * w0 - q2m - (q3m - q2m) * w1;
//	Normalize(norm);
//
//	double3 v01 = v0 + (v1 - v0) * w0;
//	double3 v23 = v2 + (v3 - v2) * w1;
//	double dist = Dot(p0 + (p1 - p0) * w0 - p2 - (p3 - p2) * w1, norm) - delta;
//	double newDist = dist + Dot(v01 - v23, norm) * dt;
//	if (newDist >= 0.0)
//		return false;
//	newDist /= -dt;
//
//	double v01n = Dot(v01, norm);
//	double v23n = Dot(v23, norm);
//	double3 v01t = v01 - norm * v01n;
//	double3 v23t = v23 - norm * v23n;
//	double3 relVT = v01t - v23t;
//	double3 n_relVT = make_double3(0.0);
//	double l_relVT = Length(relVT);
//	if (l_relVT != 0.0)
//		n_relVT = relVT * max(1.0 - friction * newDist / l_relVT, 0.0);
//
//	double3 n_v01t = v23t + n_relVT;
//	double3 imp = norm * newDist + n_v01t - v01t;
//	imp *= COL_HALF;
//	if (clothParams.masses[i0] > 0.0)
//		resolveImpulse_obs(i0, impParams, imp, newDist);
//	if (clothParams.masses[i1] > 0.0)
//		resolveImpulse_obs(i1, impParams, imp, newDist);
//	return true;
//}
//
//
//__global__ void BoundaryCollisionKernel(
//	ObjectParams params,
//	const double3 minBoundary, const double3 maxBoundary,
//	const double delta, const double friction, const double dt,
//	bool* isApplied)
//{
//	uint i = blockDim.x * blockIdx.x + threadIdx.x;
//	if (i >= params.vnum)
//		return;
//
//	if (boundaryCollision(i, params, minBoundary, maxBoundary, delta, friction, dt))
//		*isApplied = true;
//}