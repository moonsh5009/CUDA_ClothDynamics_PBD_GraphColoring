#include "CollisionSolver.h"
#include "DeviceManager.cuh"

//----------------------------------------------
//#define USED_EV_VV

//#define CCD_PRINT_TRI_TV
//#define CCD_PRINT_TRI_EV
//#define CCD_PRINT_TRI_VV
//
//#define CCD_PRINT_EDGE_EE
//#define CCD_PRINT_EDGE_EV
//#define CCD_PRINT_EDGE_VV

//----------------------------------------------
#define RESPONSE_PENETRATION	0
#define RESPONSE_AVERAGE		1
//----------------------------------------------
#define COL_STABILITY		0.2
#define COL_HALFTIME		0.99
#define COL_HALF			0.5

#define COL_CULLING			1.0e-10
#define COL_EPS				1.0e-10

#define COL_VV				1.0e-5
#define COL_EV				1.0e-5
#define COL_EE				1.0e-5
#define COL_TV				1.0e-5
//----------------------------------------------
#define COL_STABILITY_PROXIMITY		0.25
#define COL_STABILITY_CCD			0.2
#define COL_THICKNESS				0.1
#define COL_FRICTION				0.3
//----------------------------------------------

#define CUBIC_ITERATOR		30
inline __host__ __device__ __forceinline__
void PrintCE(const ContactElem& ce) {
	printf("%d, %d, %d, %d, %d, %f\n", (int)ce._isfv, ce._i[0], ce._i[1], ce._i[2], ce._i[3], ce._info);
}
inline __device__ __forceinline__
void makeSelfCE(ContactElem& ce, bool is_fv, REAL info, uint i0, uint i1, uint i2, uint i3) {
	ce._isfv = is_fv;
	ce._info = info;
	//ce._i[0] = i0;  ce._i[1] = i1;  ce._i[2] = i2;  ce._i[3] = i3;
	if (is_fv) {
		if (i0 < i1 && i0 < i2 && i1 < i2) {
			ce._i[0] = i0;  ce._i[1] = i1;  ce._i[2] = i2;  ce._i[3] = i3;
		}
		else if (i0 < i1 && i0 < i2 && i2 < i1) {
			ce._i[0] = i0;  ce._i[1] = i2;  ce._i[2] = i1;  ce._i[3] = i3;
		}
		else if (i1 < i0 && i1 < i2 && i0 < i2) {
			ce._i[0] = i1;  ce._i[1] = i0;  ce._i[2] = i2;  ce._i[3] = i3;
		}
		else if (i1 < i0 && i1 < i2 && i2 < i0) {
			ce._i[0] = i1;  ce._i[1] = i2;  ce._i[2] = i0;  ce._i[3] = i3;
		}
		else if (i2 < i0 && i2 < i1 && i0 < i1) {
			ce._i[0] = i2;  ce._i[1] = i0;  ce._i[2] = i1;  ce._i[3] = i3;
		}
		else if (i2 < i0 && i2 < i1 && i1 < i0) {
			ce._i[0] = i2;  ce._i[1] = i1;  ce._i[2] = i0;  ce._i[3] = i3;
		}
		else {
			printf("Error: tri %d, %d, %d, %d\n", i0, i1, i2, i3);
			//assert(0);
		}
	}
	else {
		if (i0 < i1 && i0 < i2 && i0 < i3 && i2 < i3) {
			ce._i[0] = i0;  ce._i[1] = i1;  ce._i[2] = i2;  ce._i[3] = i3;
		}
		else if (i0 < i1 && i0 < i2 && i0 < i3 && i3 < i2) {
			ce._i[0] = i0;  ce._i[1] = i1;  ce._i[2] = i3;  ce._i[3] = i2;
		}
		else if (i1 < i0 && i1 < i2 && i1 < i3 && i2 < i3) {
			ce._i[0] = i1;  ce._i[1] = i0;  ce._i[2] = i2;  ce._i[3] = i3;
		}
		else if (i1 < i0 && i1 < i2 && i1 < i3 && i3 < i2) {
			ce._i[0] = i1;  ce._i[1] = i0;  ce._i[2] = i3;  ce._i[3] = i2;
		}
		else if (i2 < i0 && i2 < i1 && i2 < i3 && i0 < i1) {
			ce._i[0] = i2;  ce._i[1] = i3;  ce._i[2] = i0;  ce._i[3] = i1;
		}
		else if (i2 < i0 && i2 < i1 && i2 < i3 && i1 < i0) {
			ce._i[0] = i2;  ce._i[1] = i3;  ce._i[2] = i1;  ce._i[3] = i0;
		}
		else if (i3 < i0 && i3 < i1 && i3 < i2 && i0 < i1) {
			ce._i[0] = i3;  ce._i[1] = i2;  ce._i[2] = i0;  ce._i[3] = i1;
		}
		else if (i3 < i0 && i3 < i1 && i3 < i2 && i1 < i0) {
			ce._i[0] = i3;  ce._i[1] = i2;  ce._i[2] = i1;  ce._i[3] = i0;
		}
		else {
			printf("Error: edge %d, %d, %d, %d\n", i0, i1, i2, i3);
			//assert(0);
		}
	}
}
inline __device__ __forceinline__
void addCE(const ContactElemParam& ceParams, const ContactElem* ces, const uint num) {
	uint ind = atomicAdd(ceParams.d_size, num);
	for (uint i = 0u; i < num; i++) {
		ceParams._elems[i + ind] = ces[i];
		//PrintCE(ces[i]);
	}
}

inline __device__ __forceinline__ REAL QuadraticEval(REAL a, REAL b, REAL c, REAL x) {
	return (a * x + b) * x + c;
}
inline __device__ int QuadraticSolver(REAL fa, REAL fb, REAL fc, REAL* times)
{
	int num = 0;
	if (fa == 0.0) {
		if (fc == 0.0)
			times[num++] = 0.0;
		else if (fb != 0.0) {
			REAL t = -fc / fb;
			if (t >= 0.0 && t <= 1.0)
				times[num++] = t;
		}
		return num;
	}
	else if (fa < 0.0) { fa = -fa; fb = -fb; fc = -fc; }

	if (fc > 0.0) {
		if (fb >= 0.0)
			return 0;

		REAL cx = -fb / (fa + fa);
		REAL cy = QuadraticEval(fa, fb, fc, cx);
		if (cy > 0.0)
			return 0;
		else if (cy == 0.0) {
			if (cx >= 0.0 && cx <= 1.0)
				times[num++] = cx;
			return num;
		}
		REAL det = fb * fb - 4.0 * fa * fc;
		fa = 1.0 / (fa + fa);
		det = sqrt(det);
		REAL t = (-fb - det) * fa;
		if (t >= 0.0 && t <= 1.0)
			times[num++] = t;
		t = (-fb + det) * fa;
		if (t >= 0.0 && t <= 1.0)
			times[num++] = t;
		return num;
	}
	else if (fc < 0.0) {
		REAL r1 = QuadraticEval(fa, fb, fc, 1.0);
		if (r1 < 0.0)
			return 0;
		else if (r1 == 0.0)
			times[num++] = 1.0;
		else {
			REAL t = (-fb + sqrt(fb * fb - 4.0 * fa * fc)) / (fa + fa);
			if (t >= 0.0 && t <= 1.0)
				times[num++] = t;
		}
		return num;
	}
	times[num++] = 0.0;
	REAL t = -fb / fa;
	if (t >= 0.0 && t <= 1.0)
		times[num++] = t;
	return num;
}
inline __device__ __forceinline__ REAL CubicEval(REAL a, REAL b, REAL c, REAL d, REAL x) {
	return ((a * x + b) * x + c) * x + d;
}
inline __device__ REAL CubicSolver(REAL a, REAL b, REAL c, REAL d, REAL r0, REAL r1, REAL x0, REAL x1) {
	REAL r2, x2;
	if (x0 == 0.0)	return r0;
	else if (x1 == 0.0) return r1;
	for (int itr = 0; itr < CUBIC_ITERATOR; itr++) {
		r2 = 0.5 * (r0 + r1);
		x2 = CubicEval(a, b, c, d, r2);
		if (x2 == 0.0)
			return r2;

		if (x0 * x2 < 0)	r1 = r2;
		else				r0 = r2;
	}
	return (r0 + r1) * 0.5;
}
inline __device__ int FinePlaneCoTime(
	const REAL3 a0, const REAL3 b0, const REAL3 c0, const REAL3 d0,
	const REAL3 a1, const REAL3 b1, const REAL3 c1, const REAL3 d1,
	REAL* times)
{
	int num = 0;
	REAL3 v01_0 = b0 - a0; REAL3 v01_1 = b1 - a1 - v01_0;
	REAL3 v02_0 = c0 - a0; REAL3 v02_1 = c1 - a1 - v02_0;
	REAL3 v0p_0 = d0 - a0; REAL3 v0p_1 = d1 - a1 - v0p_0;
	auto cross0 = Cross(v01_1, v02_1);
	auto cross1 = Cross(v01_0, v02_1) + Cross(v01_1, v02_0);
	auto cross2 = Cross(v01_0, v02_0);

	REAL a = Dot(cross0, v0p_1);
	REAL b = Dot(cross0, v0p_0) + Dot(cross1, v0p_1);
	REAL c = Dot(cross1, v0p_0) + Dot(cross2, v0p_1);
	REAL d = Dot(cross2, v0p_0);

	if (a == 0.0)
		return QuadraticSolver(b, c, d, times);
	else if (a < 0.0)
	{
		a = -a; b = -b; c = -c; d = -d;
	}

	REAL r0 = 0.0;
	REAL r1 = +1.0;
	REAL f0 = CubicEval(a, b, c, d, r0);
	REAL f1 = CubicEval(a, b, c, d, r1);
	REAL det = b * b - 3.0 * a * c;
	if (det >= 0.0)
	{
		REAL r3 = (-b - sqrt(det)) / (3.0 * a);
		REAL r4 = (-b + sqrt(det)) / (3.0 * a);
		const REAL f3 = CubicEval(a, b, c, d, r3);
		const REAL f4 = CubicEval(a, b, c, d, r4);

		if (r3 > 0.0 && r4 < 1.0) {
			if (f0 * f3 <= 0.0)	times[num++] = CubicSolver(a, b, c, d, r0, r3, f0, f3);
			if (f3 * f4 <= 0.0)	times[num++] = CubicSolver(a, b, c, d, r3, r4, f3, f4);
			if (f4 * f1 <= 0.0)	times[num++] = CubicSolver(a, b, c, d, r4, r1, f4, f1);
			return num;
		}
		else if (r3 > 0.0 && r3 < 1.0) {
			if (f0 * f3 <= 0.0)	times[num++] = CubicSolver(a, b, c, d, r0, r3, f0, f3);
			if (f3 * f1 <= 0.0)	times[num++] = CubicSolver(a, b, c, d, r3, r1, f3, f1);
			return num;
		}
		else if (r4 > 0.0 && r4 < 1.0) {
			if (f0 * f4 <= 0.0)	times[num++] = CubicSolver(a, b, c, d, r0, r4, f0, f4);
			if (f4 * f1 <= 0.0)	times[num++] = CubicSolver(a, b, c, d, r4, r1, f4, f1);
			return num;
		}
	}
	if (f0 * f1 > 0)
		return 0;
	times[num++] = CubicSolver(a, b, c, d, r0, r1, f0, f1);
	return num;
}

inline __device__ __forceinline__ REAL getDistanceTV(
	const REAL3 v0, const REAL3 v1, const REAL3 v2, const REAL3 p,
	REAL error, REAL* wa = nullptr, REAL* wb = nullptr)
{
	/*const REAL3 v01 = v1 - v0;
	const REAL3 v02 = v2 - v0;
	REAL v0pDotv01 = Dot(p - v0, v01);
	REAL v0pDotv02 = Dot(p - v0, v02);
	REAL v01Dotv01 = LengthSquared(v01);
	REAL v02Dotv02 = LengthSquared(v02);
	REAL v01Dotv02 = Dot(v01, v02);
	REAL v1pDotv12 = v0pDotv02 - v0pDotv01 - v01Dotv02 + v01Dotv01;

	REAL3 normal = Cross(v01, v02);
	REAL det = LengthSquared(normal);
	if (fabs(det) <= 1.0e-20)
		return DBL_MAX;
	normal /= sqrt(det);

	REAL3 rv = p;
	REAL result = 0.0;
	REAL dist = Dot(p - v0, normal);
	if (v0pDotv01 < 0.0 && v0pDotv02 < 0.0) {
		rv -= v0;
		w0 = 1.0;
		w1 = 0.0;
	}
	else if (v0pDotv02 + v01Dotv01 < v0pDotv01 + v01Dotv02 && v01Dotv01 < v0pDotv01) {
		rv -= v1;
		w0 = 0.0;
		w1 = 1.0;
	}
	else if (v02Dotv02 < v0pDotv02 && v0pDotv01 + v02Dotv02 < v0pDotv02 + v01Dotv02) {
		rv -= v2;
		w0 = 0.0;
		w1 = 0.0;
	}
	else if (v0pDotv01 * v01Dotv02 > v0pDotv02 * v01Dotv01 && v0pDotv01 >= 0.0 && v01Dotv01 >= v0pDotv01) {
		rv -= v0;
		auto tmp = v0pDotv01 * v0pDotv01 / v01Dotv01;
		result -= tmp;
		w1 = sqrt(tmp / v01Dotv01);
		w0 = 1.0 - w1;
	}
	else if ((v0pDotv01 - v01Dotv01) * (v02Dotv02 - v01Dotv02) > (v0pDotv02 - v01Dotv02) * (v01Dotv02 - v01Dotv01) && v0pDotv01 + v02Dotv02 >= v0pDotv02 + v01Dotv02) {
		rv -= v1;
		auto sqr = v01Dotv01 + v02Dotv02 - v01Dotv02 - v01Dotv02;
		auto tmp = v1pDotv12 * v1pDotv12 / sqr;
		result -= tmp;
		w0 = 0.0;
		w1 = 1.0 - sqrt(tmp / sqr);
	}
	else if (v0pDotv02 * v01Dotv02 > v0pDotv01 * v02Dotv02) {
		rv -= v0;
		auto tmp = v0pDotv02 * v0pDotv02 / v02Dotv02;
		result -= tmp;
		w0 = 1.0 - sqrt(tmp / v02Dotv02);
		w1 = 0.0;
	}
	else {
		REAL invdet = 1.0 / det;
		w1 = (v02Dotv02 * v0pDotv01 - v01Dotv02 * v0pDotv02) * invdet;
		w0 = 1.0 - w1 - (v01Dotv01 * v0pDotv02 - v01Dotv02 * v0pDotv01) * invdet;
		return fabs(dist);
	}
	result += LengthSquared(rv) - dist * dist;
	if (result <= 0.0)
		return fabs(dist);
	else if (sqrt(result) > error)
		return DBL_MAX;
	return fabs(dist);*/
	REAL w0, w1;
	REAL3 v20 = v0 - v2;
	REAL3 v21 = v1 - v2;
	REAL t0 = Dot(v20, v20);
	REAL t1 = Dot(v21, v21);
	REAL t2 = Dot(v20, v21);
	REAL t3 = Dot(v20, p - v2);
	REAL t4 = Dot(v21, p - v2);
	REAL det = t0 * t1 - t2 * t2;
	if (fabs(det) <= 1.0e-20)
		return DBL_MAX;
	REAL invdet = 1.0 / det;
	w0 = (+t1 * t3 - t2 * t4) * invdet;
	if (w0 < 0.0 || w0 > 1.0)
		return DBL_MAX;
	w1 = (-t2 * t3 + t0 * t4) * invdet;
	if (w1 < 0.0 || w1 > 1.0)
		return DBL_MAX;
	const REAL w2 = 1 - w0 - w1;
	if (w2 < 0.0 || w2 > 1.0)
		return DBL_MAX;
	REAL3 pw = v0 * w0 + v1 * w1 + v2 * w2;
	if (wa) { *wa = w0; *wb = w1; }
	return Length(pw - p);
}
inline __device__ __forceinline__ REAL getDistanceEE(
	const REAL3 pa, const REAL3 pb, const REAL3 pc, const REAL3 pd,
	REAL error, REAL* wa = nullptr, REAL* wb = nullptr, bool* isParallel = nullptr)
{
	/*REAL3 vp = pb - pa;
	REAL3 vq = pd - pc;
	REAL t0 = Dot(vp, vp);
	REAL t1 = Dot(vq, vq);
	REAL t2 = Dot(vp, vq);
	REAL det = t0 * t1 - t2 * t2;
	if (fabs(det) < 1.0e-20) {
		isParallel = true;
		REAL lp0 = Dot(pa, vp);
		REAL lp1 = Dot(pb, vp);
		REAL lq0 = Dot(pc, vp);
		REAL lq1 = Dot(pd, vp);
		REAL p_min = (lp0 < lp1) ? lp0 : lp1;
		REAL p_max = (lp0 > lp1) ? lp0 : lp1;
		REAL q_min = (lq0 < lq1) ? lq0 : lq1;
		REAL q_max = (lq0 > lq1) ? lq0 : lq1;
		REAL lm;
		if (p_max < q_min)		lm = (p_max + q_min) * 0.5;
		else if (p_min > q_max) lm = (q_max + p_min) * 0.5;
		else if (p_max < q_max)
			if (p_min < q_min)	lm = (p_max + q_min) * 0.5;
			else				lm = (p_max + p_min) * 0.5;
		else
			if (p_min < q_min)	lm = (q_max + q_min) * 0.5;
			else				lm = (q_max + p_min) * 0.5;
		w0 = (lm - lp0) / (lp1 - lp0);
		w1 = (lm - lq0) / (lq1 - lq0);
		if (w0 < 0.0 || w0 > 1.0 || w1 < 0.0 || w1 > 1.0) {
			REAL dist = 0.0;
			REAL t0 = Dot(vp, vp);
			if (w0 < 0) {
				dist = sqrt(t0) * -w0;
				w0 = 0.0;
			}
			else if (w0 > 1.0) {
				dist = sqrt(t0) * (w0 - 1.0);
				w0 = 1.0;
			}
			if (w1 < 0)			w1 = 0.0;
			else if (w1 > 1.0)	w1 = 1.0;
			if (dist + dist > error)
				return DBL_MAX;
		}
		Normalize(vp);
		auto ppc = pa - pc;
		auto vert = ppc - vp * Dot(ppc, vp);
		return Length(vert);
	}
	isParallel = false;
	REAL t3 = Dot(vp, pc - pa);
	REAL t4 = Dot(vq, pc - pa);
	REAL invdet = 1.0 / det;
	w0 = (+t1 * t3 - t2 * t4) * invdet;
	w1 = (+t2 * t3 - t0 * t4) * invdet;
	if (w0 < 0.0 || w0 > 1.0 || w1 < 0.0 || w1 > 1.0) {
		REAL w0_0 = w0;
		REAL w1_0 = w1;
		REAL inv_t0 = 1.0 / t0;
		REAL inv_t1 = 1.0 / t1;
		if (w1_0 < 0.0)			w0 = t3 * inv_t0;
		else if (w1_0 > 1.0)	w0 = (t3 + t2) * inv_t0;
		if (w0_0 < 0.0)			w1 = -t4 * inv_t1;
		else if (w0_0 > 1.0)	w1 = (-t4 + t2) * inv_t1;
		if (w0 < 0.0)			w0 = 0.0;
		else if (w0 > 1.0)		w0 = 1.0;
		if (w1 < 0.0)			w1 = 0.0;
		else if (w1 > 1.0)		w1 = 1.0;
		REAL3 normal = Cross(vp, vq);
		Normalize(normal);
		REAL3 dir = pa + vp * w0 - pc - vq * w1;
		REAL h = Dot(dir, normal);
		REAL distSqr = LengthSquared(dir) - h * h;
		if (distSqr <= 0.0)
			return fabs(h);
		else if (sqrt(distSqr) > error)
			return DBL_MAX;
		return fabs(h);
	}
	return Length(pa + vp * w0 - pc - vq * w1);*/
	REAL w0, w1;
	REAL3 vp = pb - pa;
	REAL3 vq = pd - pc;
	REAL t0 = Dot(vp, vp);
	REAL t1 = Dot(vq, vq);
	REAL t2 = Dot(vp, vq);
	REAL det = t0 * t1 - t2 * t2;
	if (fabs(det) < 1.0e-20) {
		if (isParallel)
			*isParallel = true;
		REAL lp0 = Dot(pa, vp);
		REAL lp1 = Dot(pb, vp);
		REAL lq0 = Dot(pc, vp);
		REAL lq1 = Dot(pd, vp);
		REAL p_min = (lp0 < lp1) ? lp0 : lp1;
		REAL p_max = (lp0 > lp1) ? lp0 : lp1;
		REAL q_min = (lq0 < lq1) ? lq0 : lq1;
		REAL q_max = (lq0 > lq1) ? lq0 : lq1;
		REAL lm;
		if (p_max < q_min)		lm = (p_max + q_min) * 0.5;
		else if (p_min > q_max) lm = (q_max + p_min) * 0.5;
		else if (p_max < q_max)
			if (p_min < q_min)	lm = (p_max + q_min) * 0.5;
			else				lm = (p_max + p_min) * 0.5;
		else
			if (p_min < q_min)	lm = (q_max + q_min) * 0.5;
			else				lm = (q_max + p_min) * 0.5;
		w0 = (lm - lp0) / (lp1 - lp0);
		if (w0 < 0.0 || w0 > 1.0)
			return DBL_MAX;
		w1 = (lm - lq0) / (lq1 - lq0);
		if (w1 < 0.0 || w1 > 1.0)
			return DBL_MAX;
		Normalize(vp);
		auto ppc = pa - pc;
		auto vert = ppc - vp * Dot(ppc, vp);
		if (wa) { *wa = w0; *wb = w1; }
		return Length(vert);
	}
	if (isParallel)
		*isParallel = false;
	REAL t3 = Dot(vp, pc - pa);
	REAL t4 = Dot(vq, pc - pa);
	REAL invdet = 1.0 / det;
	w0 = (+t1 * t3 - t2 * t4) * invdet;
	if (w0 < 0.0 || w0 > 1.0)
		return DBL_MAX;
	w1 = (+t2 * t3 - t0 * t4) * invdet;
	if (w1 < 0.0 || w1 > 1.0)
		return DBL_MAX;
	if (wa) { *wa = w0; *wb = w1; }
	return Length(pa + vp * w0 - pc - vq * w1);
}
//inline __device__ __forceinline__ REAL getDistanceEV(
//	const REAL3 v0, const REAL3 v1, const REAL3 p,
//	REAL& w, REAL error)
//{
//	REAL3 normal = v1 - v0;
//	Normalize(normal);
//	REAL t1 = Dot(p - v0, normal);
//	REAL t2 = Dot(p - v1, normal);
//	if (isnan(normal.x) || isnan(normal.y) || isnan(normal.z))
//		printf("Error!!!\n");
//	if (t1 < 0.0) {
//		if (-t1 > error)
//			return DBL_MAX;
//		w = 0.0;
//	}
//	else if (t2 > 0.0) {
//		if (t2 > error)
//			return DBL_MAX;
//		w = 1.0;
//	}
//	else w = t1 / (t1 - t2);
//	REAL distsqr = LengthSquared(p - v0) - t1 * t1;
//	if (distsqr < 0.0)	return 0.0;
//	return sqrt(distsqr);
//}
//inline __device__ __forceinline__ bool DetectionVV_CCD(
//	const REAL3 p0, const REAL3 p1, const REAL3 q0, const REAL3 q1,
//	REAL delta, REAL& t)
//{
//	if (min(p0.x, q1.x) - delta > max(p1.x, q1.x) ||
//		min(p0.y, q1.y) - delta > max(p1.y, q1.y) ||
//		min(p0.z, q1.z) - delta > max(p1.z, q1.z) ||
//		max(p0.x, q1.x) + delta < min(p1.x, q1.x) ||
//		max(p0.y, q1.y) + delta < min(p1.y, q1.y) ||
//		max(p0.z, q1.z) + delta < min(p1.z, q1.z))
//		return false;
//
//	REAL deltaSqr = delta * delta;
//	REAL3 x1 = p1 - p0;
//	REAL3 v1 = q1 - q0 - x1;
//	REAL fa = LengthSquared(v1);
//	REAL fb = Dot(x1, v1);
//
//	t = -fb / fa;
//	if (t >= 0.0 && t <= 1.0) {
//		if (LengthSquared(x1 + v1 * t) <= deltaSqr)
//			return true;
//	}
//	t = 1.0;
//	if (LengthSquared(q1 - q0) <= deltaSqr)
//		return true;
//	return false;
//}
//inline __device__ __forceinline__ void addTimesQuadratic(
//	REAL a, REAL b, REAL c, REAL* ts, uint& size)
//{
//	REAL tmp[2];
//	int tn = 0;
//	tn = QuadraticSolver(a, b, c, tmp);
//	if (size == 0) {
//		for (int i = 0; i < tn; i++) {
//			if (tmp[i] < 0.0 || tmp[i] > 1.0)
//				continue;
//			ts[size++] = tmp[i];
//		}
//		return;
//	}
//	int istart = 0;
//	for (int i = 0; i < tn; i++) {
//		if (tmp[i] < 0.0 || tmp[i] > 1.0)
//			continue;
//		if (ts[size - 1] < tmp[i]) {
//			ts[size++] = tmp[i];
//			continue;
//		}
//		for (int j = size - 1; j >= istart; j--) {
//			if (tmp[i] == ts[j])
//				break;
//			if (tmp[i] < ts[j]) {
//				for (int k = size++ - 1; k >= j; k--)
//					ts[k + 1] = ts[k];
//				ts[j] = tmp[i];
//				istart = j + 1;
//				break;
//			}
//			if (j == istart) {
//				for (int k = size++ - 1; k >= j; k--)
//					ts[k + 1] = ts[k];
//				ts[j] = tmp[i];
//			}
//		}
//	}
//}
//inline __device__ __forceinline__ bool DetectionEV_CCD(
//	const REAL3 p0, const REAL3 p1, const REAL3 p2,
//	const REAL3 q0, const REAL3 q1, const REAL3 q2,
//	REAL delta, REAL& t, REAL& w)
//{
//	if (min(p2.x, q2.x) - delta > max(max(max(p0.x, p1.x), q0.x), q1.x) ||
//		min(p2.y, q2.y) - delta > max(max(max(p0.y, p1.y), q0.y), q1.y) ||
//		min(p2.z, q2.z) - delta > max(max(max(p0.z, p1.z), q0.z), q1.z) ||
//		max(p2.x, q2.x) + delta < min(min(min(p0.x, p1.x), q0.x), q1.x) ||
//		max(p2.y, q2.y) + delta < min(min(min(p0.y, p1.y), q0.y), q1.y) ||
//		max(p2.z, q2.z) + delta < min(min(min(p0.z, p1.z), q0.z), q1.z))
//		return false;
//
//	REAL3 x1 = p1 - p0; REAL3 x2 = p2 - p0;
//	REAL3 v1 = q1 - q0 - x1; REAL3 v2 = q2 - q0 - x2;
//	REAL3 va = Cross(v1, v2);
//	REAL3 vb = Cross(x1, v2) + Cross(v1, x2);
//	REAL3 vc = Cross(x1, x2);
//	REAL ts[9];
//	uint tn = 0;
//	{
//		addTimesQuadratic(va.y, vb.y, vc.y, ts, tn);
//		addTimesQuadratic(va.y, vb.y, vc.y, ts, tn);
//		addTimesQuadratic(va.z, vb.z, vc.z, ts, tn);
//
//		addTimesQuadratic(Dot(v1, v2), Dot(x1, v2) + Dot(v1, x2), Dot(x1, x2), ts, tn);
//	}
//	x1 = p0 - p1; x2 = p2 - p1;
//	v1 = q0 - q1 - x1; v2 = q2 - q1 - x2;
//	va = Cross(v1, v2);
//	vb = Cross(x1, v2) + Cross(v1, x2);
//	vc = Cross(x1, x2);
//	{
//		addTimesQuadratic(va.x, vb.x, vc.x, ts, tn);
//		addTimesQuadratic(va.y, vb.y, vc.y, ts, tn);
//		addTimesQuadratic(va.z, vb.z, vc.z, ts, tn);
//
//		addTimesQuadratic(Dot(v1, v2), Dot(x1, v2) + Dot(v1, x2), Dot(x1, x2), ts, tn);
//	}
//	if (ts[tn - 1] < 1.0)
//		ts[tn++] = 1.0;
//	for (uint i = 0; i < tn; i++) {
//		t = ts[i];
//		if (t < 0.0 || t > 1.0)
//			continue;
//		REAL3 p0m = p0 + (q0 - p0) * t;
//		REAL3 p1m = p1 + (q1 - p1) * t;
//		REAL3 p2m = p2 + (q2 - p2) * t;
//		REAL dist = getDistanceEV(p0m, p1m, p2m, w, COL_EPS);
//		if (dist <= delta)
//			return true;
//	}
//	return false;
//}
inline __device__ __forceinline__ bool DetectionTV_CCD(
	const REAL3 p0, const REAL3 p1, const REAL3 p2, const REAL3 p3,
	const REAL3 q0, const REAL3 q1, const REAL3 q2, const REAL3 q3,
	REAL delta, REAL* t, REAL* w0 = nullptr, REAL* w1 = nullptr)
{
	/*if (min(p3.x, q3.x) - delta > max(max(max(max(max(p1.x, p2.x), p0.x), q1.x), q2.x), q0.x) ||
		min(p3.y, q3.y) - delta > max(max(max(max(max(p1.y, p2.y), p0.y), q1.y), q2.y), q0.y) ||
		min(p3.z, q3.z) - delta > max(max(max(max(max(p1.z, p2.z), p0.z), q1.z), q2.z), q0.z) ||
		max(p3.x, q3.x) + delta < min(min(min(min(min(p1.x, p2.x), p0.x), q1.x), q2.x), q0.x) ||
		max(p3.y, q3.y) + delta < min(min(min(min(min(p1.y, p2.y), p0.y), q1.y), q2.y), q0.y) ||
		max(p3.z, q3.z) + delta < min(min(min(min(min(p1.z, p2.z), p0.z), q1.z), q2.z), q0.z))
		return false;*/

	REAL ts[3], time = 0.0;
	int num = FinePlaneCoTime(p0, p1, p2, p3, q0, q1, q2, q3, ts);
	for (int i = 0; i < num; i++) {
		if (ts[i] < 0.0 || ts[i] > 1.0)
			continue;
		time = ts[i];
		REAL3 p0m = p0 + (q0 - p0) * ts[i];
		REAL3 p1m = p1 + (q1 - p1) * ts[i];
		REAL3 p2m = p2 + (q2 - p2) * ts[i];
		REAL3 p3m = p3 + (q3 - p3) * ts[i];
		REAL dist = getDistanceTV(p0m, p1m, p2m, p3m, COL_EPS, w0, w1);
		if (dist == DBL_MAX)
			continue;
		*t = time;
		return true;
	}
	if (time < 1.0) {
		time = 1.0;
		REAL dist = getDistanceTV(q0, q1, q2, q3, COL_EPS, w0, w1);
		if (dist < delta) {
			*t = time;
			return true;
		}
	}
	return false;
}
inline __device__ __forceinline__ bool DetectionEE_CCD(
	const REAL3 p0, const REAL3 p1, const REAL3 p2, const REAL3 p3,
	const REAL3 q0, const REAL3 q1, const REAL3 q2, const REAL3 q3,
	REAL delta, REAL* t, REAL* w0 = nullptr, REAL* w1 = nullptr)
{
	/*if (min(min(min(p0.x, p1.x), q0.x), q1.x) - delta > max(max(max(p2.x, p3.x), q2.x), q3.x) ||
		min(min(min(p0.y, p1.y), q0.y), q1.y) - delta > max(max(max(p2.y, p3.y), q2.y), q3.y) ||
		min(min(min(p0.z, p1.z), q0.z), q1.z) - delta > max(max(max(p2.z, p3.z), q2.z), q3.z) ||
		max(max(max(p0.x, p1.x), q0.x), q1.x) + delta < min(min(min(p2.x, p3.x), q2.x), q3.x) ||
		max(max(max(p0.y, p1.y), q0.y), q1.y) + delta < min(min(min(p2.y, p3.y), q2.y), q3.y) ||
		max(max(max(p0.z, p1.z), q0.z), q1.z) + delta < min(min(min(p2.z, p3.z), q2.z), q3.z))
		return false;*/

	bool isParallel;
	REAL ts[3], time = 0.0;
	int num = FinePlaneCoTime(p0, p1, p2, p3, q0, q1, q2, q3, ts);
	for (int i = 0; i < num; i++) {
		if (ts[i] < 0.0 || ts[i] > 1.0)
			continue;
		time = ts[i];
		REAL3 p0m = p0 + (q0 - p0) * ts[i];
		REAL3 p1m = p1 + (q1 - p1) * ts[i];
		REAL3 p2m = p2 + (q2 - p2) * ts[i];
		REAL3 p3m = p3 + (q3 - p3) * ts[i];
		REAL dist = getDistanceEE(p0m, p1m, p2m, p3m, COL_EPS, w0, w1, &isParallel);
		if (dist == DBL_MAX)
			continue;
		if (isParallel && dist >= delta)
			continue;
		*t = time;
		return true;
	}
	if (time < 1.0) {
		time = 1.0;
		REAL dist = getDistanceEE(q0, q1, q2, q3, COL_EPS, w0, w1);
		if (dist < delta) {
			*t = time;
			return true;
		}
	}
	return false;
}

inline __device__ bool isContactTV_Proximity(
	const REAL3 p0, const REAL3 p1, const REAL3 p2, const REAL3 p3,
	REAL delta, REAL* pene, REAL* w0 = nullptr, REAL* w1 = nullptr)
{
	/*REAL maxX, minX;
	maxX = minX = p0.x;
	if (maxX < p1.x) maxX = p1.x; if (maxX < p2.x) maxX = p2.x;
	if (minX > p1.x) minX = p1.x; if (minX > p2.x) minX = p2.x;
	if (p3.x - delta > maxX || p3.x + delta < minX)
		return false;
	maxX = minX = p0.y;
	if (maxX < p1.y) maxX = p1.y; if (maxX < p2.y) maxX = p2.y;
	if (minX > p1.y) minX = p1.y; if (minX > p2.y) minX = p2.y;
	if (p3.y - delta > maxX || p3.y + delta < minX)
		return false;
	maxX = minX = p0.z;
	if (maxX < p1.z) maxX = p1.z; if (maxX < p2.z) maxX = p2.z;
	if (minX > p1.z) minX = p1.z; if (minX > p2.z) minX = p2.z;
	if (p3.z - delta > maxX || p3.z + delta < minX)
		return false;*/
	
	/*REAL3 minX = p0, maxX = p0;
	if (minX.x > p1.x) minX.x = p1.x; if (minX.x > p2.x) minX.x = p2.x;
	if (minX.y > p1.y) minX.y = p1.y; if (minX.y > p2.y) minX.y = p2.y;
	if (minX.z > p1.z) minX.z = p1.z; if (minX.z > p2.z) minX.z = p2.z;
	if (maxX.x < p1.x) maxX.x = p1.x; if (maxX.x < p2.x) maxX.x = p2.x;
	if (maxX.y < p1.y) maxX.y = p1.y; if (maxX.y < p2.y) maxX.y = p2.y;
	if (maxX.z < p1.z) maxX.z = p1.z; if (maxX.z < p2.z) maxX.z = p2.z;
	if (p3.x + delta < minX.x || p3.y + delta < minX.y || p3.z + delta < minX.z || 
		p3.x - delta > maxX.x || p3.y - delta > maxX.y || p3.z - delta > maxX.z)
		return false;*/

	/*REAL minX[3], maxX[3];
	if (minX[0] > p1.x) minX[0] = p1.x; if (minX[0] > p2.x) minX[0] = p2.x;
	if (minX[1] > p1.y) minX[1] = p1.y; if (minX[1] > p2.y) minX[1] = p2.y;
	if (minX[2] > p1.z) minX[2] = p1.z; if (minX[2] > p2.z) minX[2] = p2.z;
	if (maxX[0] < p1.x) maxX[0] = p1.x; if (maxX[0] < p2.x) maxX[0] = p2.x;
	if (maxX[1] < p1.y) maxX[1] = p1.y; if (maxX[1] < p2.y) maxX[1] = p2.y;
	if (maxX[2] < p1.z) maxX[2] = p1.z; if (maxX[2] < p2.z) maxX[2] = p2.z;
	if (p3.x + delta < minX[0] || p3.y + delta < minX[1] || p3.z + delta < minX[2] ||
		p3.x - delta > maxX[0] || p3.y - delta > maxX[1] || p3.z - delta > maxX[2])
		return false;*/

	/*REAL maxX = p0.x; REAL maxY = p0.y; REAL maxZ = p0.z;
	REAL minX = p0.x; REAL minY = p0.y; REAL minZ = p0.z;
	if (maxX < p1.x) maxX = p1.x; if (maxX < p2.x) maxX = p2.x;
	if (maxY < p1.y) maxY = p1.y; if (maxY < p2.y) maxY = p2.y;
	if (maxZ < p1.z) maxZ = p1.z; if (maxZ < p2.z) maxZ = p2.z;
	if (minX > p1.x) minX = p1.x; if (minX > p2.x) minX = p2.x;
	if (minY > p1.y) minY = p1.y; if (minY > p2.y) minY = p2.y;
	if (minZ > p1.z) minZ = p1.z; if (minZ > p2.z) minZ = p2.z;
	if (p3.x + delta < minX || p3.y + delta < minY || p3.z + delta < minZ ||
		p3.x - delta > maxX || p3.y - delta > maxY || p3.z - delta > maxZ)
		return false;*/

	/*if (p3.x - delta > max(max(p1.x, p2.x), p0.x) ||
		p3.y - delta > max(max(p1.y, p2.y), p0.y) ||
		p3.z - delta > max(max(p1.z, p2.z), p0.z) ||
		p3.x + delta < min(min(p1.x, p2.x), p0.x) ||
		p3.y + delta < min(min(p1.y, p2.y), p0.y) ||
		p3.z + delta < min(min(p1.z, p2.z), p0.z))
		return false;*/

	/*REAL dist = getDistanceTV(p0, p1, p2, p3, COL_EPS, w0, w1);
	if (dist < delta) {
		*pene = delta - dist;
		return true;
	}
	return false;*/

	bool result = false;
	/*REAL maxX = p0.x; REAL maxY = p0.y; REAL maxZ = p0.z;
	REAL minX = p0.x; REAL minY = p0.y; REAL minZ = p0.z;
	if (maxX < p1.x) maxX = p1.x; if (maxX < p2.x) maxX = p2.x;
	if (maxY < p1.y) maxY = p1.y; if (maxY < p2.y) maxY = p2.y;
	if (maxZ < p1.z) maxZ = p1.z; if (maxZ < p2.z) maxZ = p2.z;
	if (minX > p1.x) minX = p1.x; if (minX > p2.x) minX = p2.x;
	if (minY > p1.y) minY = p1.y; if (minY > p2.y) minY = p2.y;
	if (minZ > p1.z) minZ = p1.z; if (minZ > p2.z) minZ = p2.z;
	if (p3.x + delta >= minX && p3.y + delta >= minY && p3.z + delta >= minZ &&
		p3.x - delta <= maxX && p3.y - delta <= maxY && p3.z - delta <= maxZ)
	{*/
	/*REAL minX[3], maxX[3];
	if (minX[0] > p1.x) minX[0] = p1.x; if (minX[0] > p2.x) minX[0] = p2.x;
	if (minX[1] > p1.y) minX[1] = p1.y; if (minX[1] > p2.y) minX[1] = p2.y;
	if (minX[2] > p1.z) minX[2] = p1.z; if (minX[2] > p2.z) minX[2] = p2.z;
	if (maxX[0] < p1.x) maxX[0] = p1.x; if (maxX[0] < p2.x) maxX[0] = p2.x;
	if (maxX[1] < p1.y) maxX[1] = p1.y; if (maxX[1] < p2.y) maxX[1] = p2.y;
	if (maxX[2] < p1.z) maxX[2] = p1.z; if (maxX[2] < p2.z) maxX[2] = p2.z;
	if (p3.x + delta >= minX[0] && p3.y + delta >= minX[1] && p3.z + delta >= minX[2] &&
		p3.x - delta <= maxX[0] && p3.y - delta <= maxX[1] && p3.z - delta <= maxX[2])
	{*/
	/*REAL3 minX = p0, maxX = p0;
	if (minX.x > p1.x) minX.x = p1.x; if (minX.x > p2.x) minX.x = p2.x;
	if (minX.y > p1.y) minX.y = p1.y; if (minX.y > p2.y) minX.y = p2.y;
	if (minX.z > p1.z) minX.z = p1.z; if (minX.z > p2.z) minX.z = p2.z;
	if (maxX.x < p1.x) maxX.x = p1.x; if (maxX.x < p2.x) maxX.x = p2.x;
	if (maxX.y < p1.y) maxX.y = p1.y; if (maxX.y < p2.y) maxX.y = p2.y;
	if (maxX.z < p1.z) maxX.z = p1.z; if (maxX.z < p2.z) maxX.z = p2.z;
	if (p3.x + delta >= minX.x && p3.y + delta >= minX.y && p3.z + delta >= minX.z &&
		p3.x - delta <= maxX.x && p3.y - delta <= maxX.y && p3.z - delta <= maxX.z)
	{*/
		REAL dist = getDistanceTV(p0, p1, p2, p3, COL_EPS, w0, w1);
		if (dist < delta) {
			*pene = delta - dist;
			result = true;
		}
	//}
	return result;
}
inline __device__ bool isContactEE_Proximity(
	const REAL3 p0, const REAL3 p1, const REAL3 p2, const REAL3 p3,
	REAL delta, REAL* pene, REAL* w0 = nullptr, REAL* w1 = nullptr)
{
	/*REAL maxU, minU, maxV, minV;
	maxU = minU = p0.x;
	maxV = minV = p2.x;
	if (p0.x > p1.x) minU = p1.x; else maxU = p1.x;
	if (p2.x > p3.x) minV = p3.x; else maxV = p3.x;
	if (minU - delta > maxV || maxU + delta < minV)
		return false;
	if (p0.y > p1.y) minU = p1.y; else maxU = p1.y;
	if (p2.y > p3.y) minV = p3.y; else maxV = p3.y;
	if (minU - delta > maxV || maxU + delta < minV)
		return false;
	if (p0.z > p1.z) minU = p1.z; else maxU = p1.z;
	if (p2.z > p3.z) minV = p3.z; else maxV = p3.z;
	if (minU - delta > maxV || maxU + delta < minV)
		return false;*/
	/*REAL3 minU = p0, maxU = p0;
		REAL3 minV = p2, maxV = p2;

		if (minU.x > p1.x) minU.x = p1.x; else maxU.x = p1.x;
		if (minU.y > p1.y) minU.y = p1.y; else maxU.y = p1.y;
		if (minU.z > p1.z) minU.z = p1.z; else maxU.z = p1.z;
		if (minV.x > p3.x) minV.x = p3.x; else minV.x = p3.x;
		if (minV.y > p3.y) minV.y = p3.y; else minV.y = p3.y;
		if (minV.z > p3.z) minV.z = p3.z; else minV.z = p3.z;
		if (maxU.x + delta < minV.x || maxU.y + delta < minV.y || maxU.z + delta < minV.z ||
			minU.x - delta > maxV.x || minU.y - delta > maxV.y || minU.z - delta > maxV.z)
			return false;*/
	/*if (min(p0.x, p1.x) - delta > max(p2.x, p3.x) ||
				min(p0.y, p1.y) - delta > max(p2.y, p3.y) ||
				min(p0.z, p1.z) - delta > max(p2.z, p3.z) ||
				max(p0.x, p1.x) + delta < min(p2.x, p3.x) ||
				max(p0.y, p1.y) + delta < min(p2.y, p3.y) ||
				max(p0.z, p1.z) + delta < min(p2.z, p3.z))
				return false;*/

	bool result = false;
	REAL dist = getDistanceEE(p0, p1, p2, p3, COL_EPS, w0, w1);
	if (dist < delta) {
		*pene = delta - dist;
		result = true;
	}
	return result;
}
inline __device__ __forceinline__ bool isSelfContactTV_Proximity(
	int i0, int i1, int i2, int i3,
	const REAL3 p0, const REAL3 p1, const REAL3 p2, const REAL3 p3,
	REAL delta, REAL* pene, REAL* w0 = nullptr, REAL* w1 = nullptr)
{
	if (i3 == i0 || i3 == i1 || i3 == i2)
		return false;

	return isContactTV_Proximity(p0, p1, p2, p3, delta, pene, w0, w1);
}
inline __device__ __forceinline__ bool isSelfContactEE_Proximity(
	int i0, int i1, int i2, int i3,
	const REAL3 p0, const REAL3 p1, const REAL3 p2, const REAL3 p3,
	REAL delta, REAL* pene, REAL* w0 = nullptr, REAL* w1 = nullptr)
{
	if (i0 == i2 || i0 == i3 || i1 == i2 || i1 == i3)
		return false;

	return isContactEE_Proximity(p0, p1, p2, p3, delta, pene, w0, w1);
}

inline __device__ bool isContactTV_CCD(
	const REAL3 p0, const REAL3 p1, const REAL3 p2, const REAL3 p3,
	const REAL3 q0, const REAL3 q1, const REAL3 q2, const REAL3 q3,
	REAL* t, REAL* w0 = nullptr, REAL* w1 = nullptr)
{
	if (DetectionTV_CCD(p0, p1, p2, p3, q0, q1, q2, q3, COL_TV, t, w0, w1)) {
#ifdef CCD_PRINT_TRI_TV
		printf("tri0 %f %f %f\n", *w0, *w1, *t);
#endif
		return true;
	}

#ifdef USED_EV_VV
	if (DetectionEV_CCD(p1, p0, p3, q1, q0, q3, COL_EV, t, w0)) {
		*w1 = 1.0 - *w0;
#ifdef CCD_PRINT_TRI_EV
		printf("r tri1 %f %f %f\n", *w0, *w1, *t);
#endif
		return true;
	}
	if (DetectionEV_CCD(p2, p1, p3, q2, q1, q3, COL_EV, t, w1)) {
		*w0 = 0.0;
#ifdef CCD_PRINT_TRI_EV
		printf("r tri2 %f %f %f\n", *w0, *w1, *t);
#endif
		return true;
	}
	if (DetectionEV_CCD(p2, p0, p3, q2, q0, q3, COL_EV, t, w0)) {
		*w1 = 0.0;
#ifdef CCD_PRINT_TRI_EV
		printf("r tri3 %f %f %f\n", *w0, *w1, *t);
#endif
		return true;
	}
	if (DetectionVV_CCD(p0, p3, q0, q3, COL_VV, t)) {
		*w0 = 1.0; *w1 = 0.0;
#ifdef CCD_PRINT_TRI_VV
		printf("r tri4 %f %f %f\n", *w0, *w1, *t);
#endif
		return true;
	}
	if (DetectionVV_CCD(p1, p3, q1, q3, COL_VV, t)) {
		*w0 = 0.0; *w1 = 1.0;
#ifdef CCD_PRINT_TRI_VV
		printf("r tri5 %f %f %f\n", *w0, *w1, *t);
#endif
		return true;
	}
	if (DetectionVV_CCD(p2, p3, q2, q3, COL_VV, t)) {
		*w0 = 0.0; *w1 = 0.0;
#ifdef CCD_PRINT_TRI_VV
		printf("r tri6 %f %f %f\n", *w0, *w1, *t);
#endif
		return true;
	}
#endif
	return false;
}
inline __device__ bool isContactEE_CCD(
	const REAL3 p0, const REAL3 p1, const REAL3 p2, const REAL3 p3,
	const REAL3 q0, const REAL3 q1, const REAL3 q2, const REAL3 q3,
	REAL* t, REAL* w0 = nullptr, REAL* w1 = nullptr)
{
	if (DetectionEE_CCD(p0, p1, p2, p3, q0, q1, q2, q3, COL_EE, t, w0, w1)) {
#ifdef CCD_PRINT_EDGE_EE
		printf("edge0 %f %f %f\n", *w0, *w1, *t);
#endif
		return true;
	}

#ifdef USED_EV_VV
	if (DetectionEV_CCD(p0, p1, p2, q0, q1, q2, COL_EV, t, w0)) {
		*w1 = 0.0;
#ifdef CCD_PRINT_EDGE_EV
		printf("r edge1 %f %f %f\n", *w0, *w1, *t);
#endif
		return true;
	}
	if (DetectionEV_CCD(p0, p1, p3, q0, q1, q3, COL_EV, t, w0)) {
		*w1 = 1.0;
#ifdef CCD_PRINT_EDGE_EV
		printf("r edge2 %f %f %f\n", *w0, *w1, *t);
#endif
		return true;
	}
	if (DetectionEV_CCD(p2, p3, p0, q2, q3, q0, COL_EV, t, w1)) {
		*w0 = 0.0;
#ifdef CCD_PRINT_EDGE_EV
		printf("r edge3 %f %f %f\n", *w0, *w1, *t);
#endif
		return true;
	}
	if (DetectionEV_CCD(p2, p3, p1, q2, q3, q1, COL_EV, t, w1)) {
		*w0 = 1.0;
#ifdef CCD_PRINT_EDGE_EV
		printf("r edge4 %f %f %f\n", *w0, *w1, *t);
#endif
		return true;
	}
	if (DetectionVV_CCD(p0, p2, q0, q2, COL_VV, t)) {
		*w0 = 0.0; *w1 = 0.0;
#ifdef CCD_PRINT_EDGE_VV
		printf("r edge5 %f %f %f\n", *w0, *w1, *t);
#endif
		return true;
	}
	if (DetectionVV_CCD(p0, p3, q0, q3, COL_VV, t)) {
		*w0 = 0.0; *w1 = 1.0;
#ifdef CCD_PRINT_EDGE_VV
		printf("r edge6 %f %f %f\n", *w0, *w1, *t);
#endif
		return true;
	}
	if (DetectionVV_CCD(p1, p2, q1, q2, COL_VV, t)) {
		*w0 = 1.0; *w1 = 0.0;
#ifdef CCD_PRINT_EDGE_VV
		printf("r edge7 %f %f %f\n", *w0, *w1, *t);
#endif
		return true;
	}
	if (DetectionVV_CCD(p1, p3, q1, q3, COL_VV, t)) {
		*w0 = 1.0; *w1 = 1.0;
#ifdef CCD_PRINT_EDGE_VV
		printf("r edge8 %f %f %f\n", *w0, *w1, *t);
#endif
		return true;
	}
#endif
	return false;
}
inline __device__ __forceinline__ bool isSelfContactTV_CCD(
	int i0, int i1, int i2, int i3,
	const REAL3 p0, const REAL3 p1, const REAL3 p2, const REAL3 p3,
	const REAL3 q0, const REAL3 q1, const REAL3 q2, const REAL3 q3,
	REAL* t, REAL* w0 = nullptr, REAL* w1 = nullptr)
{
	if (i3 == i0 || i3 == i1 || i3 == i2)
		return false;
	return isContactTV_CCD(p0, p1, p2, p3, q0, q1, q2, q3, t, w0, w1);
}
inline __device__ __forceinline__ bool isSelfContactEE_CCD(
	int i0, int i1, int i2, int i3,
	const REAL3 p0, const REAL3 p1, const REAL3 p2, const REAL3 p3,
	const REAL3 q0, const REAL3 q1, const REAL3 q2, const REAL3 q3,
	REAL* t, REAL* w0 = nullptr, REAL* w1 = nullptr)
{
	if (i0 == i2 || i0 == i3 || i1 == i2 || i1 == i3)
		return false;

	return isContactEE_CCD(p0, p1, p2, p3, q0, q1, q2, q3, t, w0, w1);
}

inline __device__ __forceinline__ bool resolveSelfCollisionProximity(
	bool is_fv,
	int i0, int i1, int i2, int i3,
	const ObjParam param,
	REAL delta, REAL dt)
{
	if (is_fv) {
		if (i3 == i0 || i3 == i1 || i3 == i2)
			return false;
	}
	else {
		if (i0 == i2 || i0 == i3 || i1 == i2 || i1 == i3)
			return false;
	}

	REAL stiffness = 10.0;
	REAL invM0 = param._invMs[i0];
	REAL invM1 = param._invMs[i1];
	REAL invM2 = param._invMs[i2];
	REAL invM3 = param._invMs[i3];
	i0 *= 3u; i1 *= 3u; i2 *= 3u; i3 *= 3u;
	REAL3 p0 = make_REAL3(param._ns[i0 + 0u], param._ns[i0 + 1u], param._ns[i0 + 2u]);
	REAL3 p1 = make_REAL3(param._ns[i1 + 0u], param._ns[i1 + 1u], param._ns[i1 + 2u]);
	REAL3 p2 = make_REAL3(param._ns[i2 + 0u], param._ns[i2 + 1u], param._ns[i2 + 2u]);
	REAL3 p3 = make_REAL3(param._ns[i3 + 0u], param._ns[i3 + 1u], param._ns[i3 + 2u]);
	REAL3 v0 = make_REAL3(param._vs[i0 + 0u], param._vs[i0 + 1u], param._vs[i0 + 2u]);
	REAL3 v1 = make_REAL3(param._vs[i1 + 0u], param._vs[i1 + 1u], param._vs[i1 + 2u]);
	REAL3 v2 = make_REAL3(param._vs[i2 + 0u], param._vs[i2 + 1u], param._vs[i2 + 2u]);
	REAL3 v3 = make_REAL3(param._vs[i3 + 0u], param._vs[i3 + 1u], param._vs[i3 + 2u]);
	REAL w0, w1, w2, w3, p_depth;
	if (is_fv) {
		if (!isContactTV_Proximity(p0, p1, p2, p3, delta, &p_depth, &w0, &w1))
			return false;
		w2 = 1.0 - w0 - w1;
		w3 = 1.0;
	}
	else {
		if (!isContactEE_Proximity(p0, p1, p2, p3, delta, &p_depth, &w1, &w3))
			return false;
		w0 = 1.0 - w1;
		w2 = 1.0 - w3;
	}
	REAL iP0 = w0 * w0 * invM0;
	REAL iP1 = w1 * w1 * invM1;
	REAL iP2 = w2 * w2 * invM2;
	REAL iP3 = w3 * w3 * invM3;
	REAL iPt = iP0 + iP1 + iP2 + iP3;
	if (iPt == 0.0)
		return false;

	REAL3 norm, tmp0, tmp1;
	if (is_fv) {
		tmp0 = p0 * w0 + p1 * w1 + p2 * w2;
		tmp1 = p3;
	}
	else {
		tmp0 = p0 + (p1 - p0) * w1;
		tmp1 = p2 + (p3 - p2) * w3;
	}
	norm = tmp1 - tmp0;
	if (!Normalize(norm)) {
		REAL ht = -dt;
		REAL3 q0 = p0 + v0 * ht;
		REAL3 q1 = p1 + v1 * ht;
		REAL3 q2 = p2 + v2 * ht;
		REAL3 q3 = p3 + v3 * ht;
		if (is_fv) {
			tmp0 = q0 * w0 + q1 * w1 + q2 * w2;
			tmp1 = q3;
			printf("Self Proximity VT\n");
		}
		else {
			tmp0 = q0 + (q1 - q0) * w1;
			tmp1 = q2 + (q3 - q2) * w3;
			printf("Self Proximity EE\n");
		}
		norm = tmp1 - tmp0;
		Normalize(norm);
	}
	//REAL p_depth = delta - Dot(tmp1 - tmp0, norm);
	//printf("asdfasdf %f %f\n", pene, p_depth);
	REAL relV;
	if (is_fv)
		relV = Dot(v3 - v0 * w0 - v1 * w1 - v2 * w2, norm);
	else
		relV = Dot(v2 + (v3 - v2) * w3 - v0 - (v1 - v0) * w1, norm);
	
	if (relV >= COL_STABILITY * p_depth / dt)
		return false;

	REAL imp_el = dt * stiffness * p_depth;
	REAL imp_ie = COL_STABILITY * p_depth / dt - relV;
	REAL imp = (imp_el < imp_ie) ? imp_el : imp_ie;
	REAL imp_mod = 2. * imp / iPt;
	imp_mod *= 0.25;

	if (invM0 > 0.0) {
		imp = imp_mod * w0 * invM0;
		v0 -= norm * imp;
		param._vs[i0 + 0u] = v0.x;
		param._vs[i0 + 1u] = v0.y;
		param._vs[i0 + 2u] = v0.z;
	}
	if (invM1 > 0.0) {
		imp = imp_mod * w1 * invM1;
		v1 -= norm * imp;
		param._vs[i1 + 0u] = v1.x;
		param._vs[i1 + 1u] = v1.y;
		param._vs[i1 + 2u] = v1.z;
	}
	if (invM2 > 0.0) {
		imp = imp_mod * w2 * invM2;
		if (is_fv)
			imp = -imp;
		v2 += norm * imp;
		param._vs[i2 + 0u] = v2.x;
		param._vs[i2 + 1u] = v2.y;
		param._vs[i2 + 2u] = v2.z;
	}
	if (invM3 > 0.0) {
		imp = imp_mod * w3 * invM3;
		v3 += norm * imp;
		param._vs[i3 + 0u] = v3.x;
		param._vs[i3 + 1u] = v3.y;
		param._vs[i3 + 2u] = v3.z;
	}
	return true;
}
inline __device__ __forceinline__ bool resolveSelfCollisionCCD(
	bool is_fv,
	int i0, int i1, int i2, int i3,
	const ObjParam param,
	REAL delta, REAL dt)
{
	if (is_fv) {
		if (i3 == i0 || i3 == i1 || i3 == i2)
			return false;
	}
	else {
		if (i0 == i2 || i0 == i3 || i1 == i2 || i1 == i3)
			return false;
	}

	REAL invM0 = param._invMs[i0];
	REAL invM1 = param._invMs[i1];
	REAL invM2 = param._invMs[i2];
	REAL invM3 = param._invMs[i3];
	i0 *= 3u; i1 *= 3u; i2 *= 3u; i3 *= 3u;
	REAL3 p0 = make_REAL3(param._ns[i0 + 0u], param._ns[i0 + 1u], param._ns[i0 + 2u]);
	REAL3 p1 = make_REAL3(param._ns[i1 + 0u], param._ns[i1 + 1u], param._ns[i1 + 2u]);
	REAL3 p2 = make_REAL3(param._ns[i2 + 0u], param._ns[i2 + 1u], param._ns[i2 + 2u]);
	REAL3 p3 = make_REAL3(param._ns[i3 + 0u], param._ns[i3 + 1u], param._ns[i3 + 2u]);
	REAL3 v0 = make_REAL3(param._vs[i0 + 0u], param._vs[i0 + 1u], param._vs[i0 + 2u]);
	REAL3 v1 = make_REAL3(param._vs[i1 + 0u], param._vs[i1 + 1u], param._vs[i1 + 2u]);
	REAL3 v2 = make_REAL3(param._vs[i2 + 0u], param._vs[i2 + 1u], param._vs[i2 + 2u]);
	REAL3 v3 = make_REAL3(param._vs[i3 + 0u], param._vs[i3 + 1u], param._vs[i3 + 2u]);
	REAL3 q0 = p0 + v0 * dt;
	REAL3 q1 = p1 + v1 * dt;
	REAL3 q2 = p2 + v2 * dt;
	REAL3 q3 = p3 + v3 * dt;
	REAL t, w0, w1, w2, w3;
	if (is_fv) {
		if (!isContactTV_CCD(p0, p1, p2, p3, q0, q1, q2, q3, &t, &w0, &w1))
			return false;
		w2 = 1.0 - w0 - w1;
		w3 = 1.0;
	}
	else {
		if (!isContactEE_CCD(p0, p1, p2, p3, q0, q1, q2, q3, &t, &w1, &w3))
			return false;
		w0 = 1.0 - w1;
		w2 = 1.0 - w3;
	}

	REAL iP0 = w0 * w0 * invM0;
	REAL iP1 = w1 * w1 * invM1;
	REAL iP2 = w2 * w2 * invM2;
	REAL iP3 = w3 * w3 * invM3;
	REAL iPt = iP0 + iP1 + iP2 + iP3;
	if (iPt == 0.0)
		return false;

	REAL3 norm, tmp0, tmp1;
	/*REAL ht = t * dt * COL_HALFTIME;
	q0 = p0 + v0 * ht;
	q1 = p1 + v1 * ht;
	q2 = p2 + v2 * ht;
	q3 = p3 + v3 * ht;
	if (is_fv) {
		tmp0 = q0 * w0 + q1 * w1 + q2 * w2;
		norm = q3 - tmp0;
	}
	else {
		tmp0 = q0 + (q1 - q0) * w1;
		tmp1 = q2 + (q3 - q2) * w3;
		norm = tmp1 - tmp0;
	}*/
	if (is_fv) {
		tmp0 = p0 * w0 + p1 * w1 + p2 * w2;
		tmp1 = p3;
	}
	else {
		tmp0 = p0 + (p1 - p0) * w1;
		tmp1 = p2 + (p3 - p2) * w3;
	}
	norm = tmp1 - tmp0;
	if (!Normalize(norm)) {
		//REAL ht = -dt * (1.0 - COL_HALFTIME);
		REAL ht = -dt;
		q0 = p0 + v0 * ht;
		q1 = p1 + v1 * ht;
		q2 = p2 + v2 * ht;
		q3 = p3 + v3 * ht;
		if (is_fv) {
			tmp0 = q0 * w0 + q1 * w1 + q2 * w2;
			tmp1 = q3;
			printf("Self CCD VT\n");
		}
		else {
			tmp0 = q0 + (q1 - q0) * w1;
			tmp1 = q2 + (q3 - q2) * w3;
			printf("Self CCD EE\n");
		}
		norm = tmp1 - tmp0;
		Normalize(norm);
	}
	REAL relV;
	if (is_fv)
		relV = Dot(v3 - v0 * w0 - v1 * w1 - v2 * w2, norm);
	else
		relV = Dot(v2 + (v3 - v2) * w3 - v0 - (v1 - v0) * w1, norm);

	REAL imp = COL_STABILITY * delta / dt - relV;
	if (imp <= 0.0)
		return false;

	REAL imp_mod = 2. * imp / iPt;
	imp_mod *= 0.1;
	if (invM0 > 0.0) {
		imp = imp_mod * w0 * invM0;
		v0 -= norm * imp;
		param._vs[i0 + 0u] = v0.x;
		param._vs[i0 + 1u] = v0.y;
		param._vs[i0 + 2u] = v0.z;
	}
	if (invM1 > 0.0) {
		imp = imp_mod * w1 * invM1;
		v1 -= norm * imp;
		param._vs[i1 + 0u] = v1.x;
		param._vs[i1 + 1u] = v1.y;
		param._vs[i1 + 2u] = v1.z;
	}
	if (invM2 > 0.0) {
		imp = imp_mod * w2 * invM2;
		if (is_fv)
			imp = -imp;
		v2 += norm * imp;
		param._vs[i2 + 0u] = v2.x;
		param._vs[i2 + 1u] = v2.y;
		param._vs[i2 + 2u] = v2.z;
	}
	if (invM3 > 0.0) {
		imp = imp_mod * w3 * invM3;
		v3 += norm * imp;
		param._vs[i3 + 0u] = v3.x;
		param._vs[i3 + 1u] = v3.y;
		param._vs[i3 + 2u] = v3.z;
	}
	return true;
}

inline __device__ __forceinline__ void CalcInvMat3(REAL* ainv, const REAL* a)
{
	const REAL det =
		+a[0] * a[4] * a[8] + a[3] * a[7] * a[2] + a[6] * a[1] * a[5]
		- a[0] * a[7] * a[5] - a[6] * a[4] * a[2] - a[3] * a[1] * a[8];
	const REAL inv_det = 1.0 / det;

	ainv[0] = inv_det * (a[4] * a[8] - a[5] * a[7]);
	ainv[1] = inv_det * (a[2] * a[7] - a[1] * a[8]);
	ainv[2] = inv_det * (a[1] * a[5] - a[2] * a[4]);

	ainv[3] = inv_det * (a[5] * a[6] - a[3] * a[8]);
	ainv[4] = inv_det * (a[0] * a[8] - a[2] * a[6]);
	ainv[5] = inv_det * (a[2] * a[3] - a[0] * a[5]);

	ainv[6] = inv_det * (a[3] * a[7] - a[4] * a[6]);
	ainv[7] = inv_det * (a[1] * a[6] - a[0] * a[7]);
	ainv[8] = inv_det * (a[0] * a[4] - a[1] * a[3]);
}