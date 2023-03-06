#ifndef __DEVICE_MANAGER_H__
#define __DEVICE_MANAGER_H__

#pragma once
#include "device_launch_parameters.h"
#include "cuda_runtime.h"
#include "thrust/device_ptr.h"
#include "thrust/scan.h"
#include "thrust/sort.h"
#include <chrono>
#include <vector>

using namespace std;
typedef unsigned int uint;
typedef unsigned char uchar;
typedef std::chrono::system_clock::time_point ctimer;

#define M_PI				3.14159265359

#if 1
#define REAL				double
#define REAL2				double2
#define REAL3				double3
#define REAL4				double4
#define REAL_AS_INT			__double_as_longlong
#define INT_AS_REAL			__longlong_as_double
#define REAL_INT			unsigned long long
#define REAL_MAX			DBL_MAX
#define REAL_INT_MAX		ULLONG_MAX
#define make_REAL3			make_double3
#else
#define REAL				float
#define REAL2				float2
#define REAL3				float3
#define REAL4				float4
#define REAL_AS_INT			__float_as_int
#define INT_AS_REAL			__int_as_float
#define REAL_INT			unsigned int
#define REAL_MAX			FLT_MAX
#define REAL_INT_MAX		UINT_MAX
#define make_REAL3			make_float3
#endif

#define MAX_DBLOCKSIZE		2048
#define MAX_BLOCKSIZE		1024
#define DBLOCKSIZE			256
#define BLOCKSIZE			128
#define HBLOCKSIZE			64
#define WARPSIZE			32
#define EPS					1.0e-20

#define INV3				0.33333333333333333333333333333

//#define TESTTIMER
#define CUDA_DEBUG

#ifndef CUDA_DEBUG
#define CUDA_CHECK(x)	(x)
#else
#define CUDA_CHECK(x)	do {\
		(x); \
		cudaError_t e = cudaGetLastError(); \
		if (e != cudaSuccess) { \
			printf("cuda failure %s:%d: '%s'\n", \
				__FILE__, __LINE__, cudaGetErrorString(e)); \
			/*exit(1);*/ \
		}\
	} while(0)
#endif

#define CNOW		std::chrono::system_clock::now()

struct uint2_CMP
{
	__host__ __device__
		bool operator()(const uint2& a, const uint2& b) {
		if (a.x != b.x)
			return a.x < b.x;
		return a.y < b.y;
	}
};
//-------------------------------------------------------------------------------------------------------------
class StreamParam {
public:
	vector<cudaStream_t>	_streams;
public:
	StreamParam() {}
	virtual ~StreamParam() {
		freeStream();
	}
public:
	void initStream(uint num) {
		if (_streams.size() > 0)
			freeStream();
		_streams.resize(num);

		for (int i = 0; i < num; i++)
			CUDA_CHECK(cudaStreamCreate(&_streams[i]));
	}
	void freeStream(void) {
		for (int i = 0; i < _streams.size(); i++)
			CUDA_CHECK(cudaStreamDestroy(_streams[i]));
	}
public:
	inline cudaStream_t* begin(void) {
		return &_streams[0];
	}
	inline cudaStream_t* end(void) {
		return (&_streams[0]) + _streams.size();
	}
	inline cudaStream_t& operator[](size_t i) {
		if (i >= _streams.size()) {
			printf("Error : StreamParam_[] : index out\n");
			exit(1);
		}
		return _streams[i];
	}
};
//-------------------------------------------------------------------------------------------------------------
static __inline__ __host__ __device__ uint divup(uint x, uint y) {
	return (x + y - 1) / y;
}
static __inline__ __host__ __device__ void setV(REAL3 & a, int i, REAL x)
{
	if (i == 0) a.x = x;
	else if (i == 1) a.y = x;
	else a.z = x;
}
static __inline__ __host__ __device__ REAL& getV(const REAL3 & a, int i)
{
	/*if (i == 0) return a.x;
	else if (i == 1) return a.y;
	return a.z;*/
	return ((REAL*)&a)[i];
}
static __inline__ __host__ __device__ REAL3 operator+(const REAL3 a, const REAL3 b)
{
	return make_REAL3(a.x + b.x, a.y + b.y, a.z + b.z);
}
static __inline__ __host__ __device__ REAL3 operator-(const REAL3 a, const REAL3 b)
{
	return make_REAL3(a.x - b.x, a.y - b.y, a.z - b.z);
}
static __inline__ __host__ __device__ REAL operator*(const REAL3 a, const REAL3 b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}
static __inline__ __host__ __device__ REAL3 operator+(const REAL3 a, const REAL b)
{
	return make_REAL3(a.x + b, a.y + b, a.z + b);
}
static __inline__ __host__ __device__ REAL3 operator-(const REAL3 a, const REAL b)
{
	return make_REAL3(a.x - b, a.y - b, a.z - b);
}
static __inline__ __host__ __device__ REAL3 operator*(const REAL3 a, const REAL b)
{
	return make_REAL3(a.x * b, a.y * b, a.z * b);
}
static __inline__ __host__ __device__ REAL3 operator/(const REAL3 a, const REAL b)
{
	return make_REAL3(a.x / b, a.y / b, a.z / b);
}
static __inline__ __host__ __device__ REAL3 operator*(const REAL b, const REAL3 a)
{
	return make_REAL3(a.x * b, a.y * b, a.z * b);
}
static __inline__ __host__ __device__ REAL3 operator/(const REAL b, const REAL3 a)
{
	return make_REAL3(b / a.x, b / a.y, b / a.z);
}
static __inline__ __host__ __device__ void operator+=(REAL3 & a, const REAL3 b)
{
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
}
static __inline__ __host__ __device__ void operator-=(REAL3 & a, const REAL3 b)
{
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
}
static __inline__ __host__ __device__ void operator*=(REAL3 & a, const REAL3 b)
{
	a.x *= b.x;
	a.y *= b.y;
	a.z *= b.z;
}
static __inline__ __host__ __device__ void operator/=(REAL3 & a, const REAL3 b)
{
	a.x /= b.x;
	a.y /= b.y;
	a.z /= b.z;
}
static __inline__ __host__ __device__ void operator+=(REAL3 & a, const REAL b)
{
	a.x += b;
	a.y += b;
	a.z += b;
}
static __inline__ __host__ __device__ void operator-=(REAL3 & a, const REAL b)
{
	a.x -= b;
	a.y -= b;
	a.z -= b;
}
static __inline__ __host__ __device__ void operator*=(REAL3 & a, const REAL b)
{
	a.x *= b;
	a.y *= b;
	a.z *= b;
}
static __inline__ __host__ __device__ void operator/=(REAL3 & a, const REAL b)
{
	a.x /= b;
	a.y /= b;
	a.z /= b;
}
static __inline__ __host__ __device__ bool operator==(const REAL3 a, const REAL3 b)
{
	return a.x == b.x && a.y == b.y && a.z == b.z;
}
static __inline__ __host__ __device__ bool operator!=(const REAL3 a, const REAL3 b)
{
	return a.x != b.x && a.y != b.y && a.z != b.z;
}
static __inline__ __host__ __device__ bool operator==(const REAL3 a, const REAL b)
{
	return a.x == b && a.y == b && a.z == b;
}
static __inline__ __host__ __device__ bool operator!=(const REAL3 a, const REAL b)
{
	return a.x != b && a.y != b && a.z != b;
}
static __inline__ __host__ __device__ bool operator<(const REAL3 a, const REAL3 b)
{
	return a.x < b.x&& a.y < b.y&& a.z < b.z;
}
static __inline__ __host__ __device__ bool operator>(const REAL3 a, const REAL3 b)
{
	return a.x > b.x && a.y > b.y && a.z > b.z;
}
static __inline__ __host__ __device__ bool operator<=(const REAL3 a, const REAL3 b)
{
	return a.x <= b.x && a.y <= b.y && a.z <= b.z;
}
static __inline__ __host__ __device__ bool operator>=(const REAL3 a, const REAL3 b)
{
	return a.x >= b.x && a.y >= b.y && a.z >= b.z;
}
static __inline__ __host__ __device__ bool operator<(const REAL3 a, const REAL b)
{
	return a.x < b&& a.y < b&& a.z < b;
}
static __inline__ __host__ __device__ bool operator>(const REAL3 a, const REAL b)
{
	return a.x > b && a.y > b && a.z > b;
}
static __inline__ __host__ __device__ bool operator<=(const REAL3 a, const REAL b)
{
	return a.x <= b && a.y <= b && a.z <= b;
}
static __inline__ __host__ __device__ bool operator>=(const REAL3 a, const REAL b)
{
	return a.x >= b && a.y >= b && a.z >= b;
}
static __inline__ __host__ __device__ REAL3 make_REAL3(REAL s)
{
	return make_REAL3(s, s, s);
}
static __inline__ __host__ __device__ REAL3 make_REAL3(const REAL * a)
{
	return make_REAL3(a[0], a[1], a[2]);
}
static __inline__ __host__ __device__ REAL3 make_REAL3(const REAL2 & a)
{
	return make_REAL3(a.x, a.y, 0.0f);
}
static __inline__ __host__ __device__ REAL3 make_REAL3(const REAL2 & a, REAL s)
{
	return make_REAL3(a.x, a.y, s);
}
static __inline__ __host__ __device__ REAL3 make_REAL3(const REAL3 a)
{
	return make_REAL3(a.x, a.y, a.z);
}
static __inline__ __host__ __device__ REAL3 make_REAL3(const REAL4 & a)
{
	return make_REAL3(a.x, a.y, a.z);
}

static __inline__ __host__ __device__ REAL3 Add(REAL* a, const REAL* b)
{
	a[0] += b[0]; a[1] += b[1]; a[2] += b[2];
}
static __inline__ __host__ __device__ REAL3 Add(REAL* a, const REAL b)
{
	a[0] += b; a[1] += b; a[2] += b;
}
static __inline__ __host__ __device__ REAL3 Add(REAL* a, const REAL* b, const REAL* c)
{
	a[0] = b[0] + c[0]; a[1] = b[1] + c[1]; a[2] = b[2] + c[2];
}
static __inline__ __host__ __device__ REAL3 Add(REAL* a, const REAL* b, const REAL c)
{
	a[0] = b[0] + c; a[1] = b[1] + c; a[2] = b[2] + c;
}
static __inline__ __host__ __device__ REAL3 Sub(REAL* a, const REAL* b)
{
	a[0] -= b[0]; a[1] -= b[1]; a[2] -= b[2];
}
static __inline__ __host__ __device__ REAL3 Sub(REAL* a, const REAL b)
{
	a[0] -= b; a[1] -= b; a[2] -= b;
}
static __inline__ __host__ __device__ REAL3 Sub(REAL* a, const REAL* b, const REAL* c)
{
	a[0] = b[0] - c[0]; a[1] = b[1] - c[1]; a[2] = b[2] - c[2];
}
static __inline__ __host__ __device__ REAL3 Sub(REAL* a, const REAL* b, const REAL c)
{
	a[0] = b[0] - c; a[1] = b[1] - c; a[2] = b[2] - c;
}
static __inline__ __host__ __device__ void Mult(REAL* a, const REAL* b)
{
	a[0] *= b[0]; a[1] *= b[1]; a[2] *= b[2];
}
static __inline__ __host__ __device__ void Mult(REAL* a, const REAL b)
{
	a[0] *= b; a[1] *= b; a[2] *= b;
}
static __inline__ __host__ __device__ REAL3 Mult(REAL* a, const REAL* b, const REAL* c)
{
	a[0] = b[0] * c[0]; a[1] = b[1] * c[1]; a[2] = b[2] * c[2];
}
static __inline__ __host__ __device__ REAL3 Mult(REAL* a, const REAL* b, const REAL c)
{
	a[0] = b[0] * c; a[1] = b[1] * c; a[2] = b[2] * c;
}
static __inline__ __host__ __device__ REAL3 Mult(REAL* a, const REAL* b, const REAL c, const REAL* d, const REAL e)
{
	a[0] = b[0] * c + d[0] * e; 
	a[1] = b[1] * c + d[1] * e; 
	a[2] = b[2] * c + d[2] * e;
}
static __inline__ __host__ __device__ void Div(REAL* a, const REAL* b)
{
	a[0] /= b[0]; a[1] /= b[1]; a[2] /= b[2];
}
static __inline__ __host__ __device__ void Div(REAL* a, const REAL b)
{
	a[0] /= b; a[1] /= b; a[2] /= b;
}
static __inline__ __host__ __device__ REAL3 Div(REAL* a, const REAL* b, const REAL* c)
{
	a[0] = b[0] / c[0]; a[1] = b[1] / c[1]; a[2] = b[2] / c[2];
}
static __inline__ __host__ __device__ REAL3 Div(REAL* a, const REAL* b, const REAL c)
{
	a[0] = b[0] / c; a[1] = b[1] / c; a[2] = b[2] / c;
}

static __inline__ __host__ __device__ REAL3 minVec(const REAL3 a, const REAL3 b)
{
	REAL3 x;
	if (a.x <= b.x) x.x = a.x;
	else x.x = b.x;
	if (a.y <= b.y) x.y = a.y;
	else x.y = b.y;
	if (a.z <= b.z) x.z = a.z;
	else x.z = b.z;
	return x;
}
static __inline__ __host__ __device__ REAL3 maxVec(const REAL3 a, const REAL3 b)
{
	REAL3 x;
	if (a.x >= b.x) x.x = a.x;
	else x.x = b.x;
	if (a.y >= b.y) x.y = a.y;
	else x.y = b.y;
	if (a.z >= b.z) x.z = a.z;
	else x.z = b.z;
	return x;
}
static __inline__ __host__ __device__ void Print(const REAL3 a)
{
	printf("%f, %f, %f\n", a.x, a.y, a.z);
}
static __inline__ __host__ __device__ REAL LengthSquared(const REAL3 a)
{
	return a * a;
}
static __inline__ __host__ __device__ REAL LengthSquared(REAL x, REAL y, REAL z)
{
	return x * x + y * y + z * z;
}
static __inline__ __host__ __device__ REAL Length(const REAL3 a)
{
	return sqrt(a * a);
}
static __inline__ __host__ __device__ bool Normalize(REAL3 & a)
{
	REAL norm = Length(a);
	if (norm == 0) {
		//printf("Error Normalize Length 0\n"); 
		return false;
	}
	a *= 1.0 / norm;
	return true;
}
static __inline__ __host__ __device__ REAL Dot(REAL3 a, REAL3 b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}
static __inline__ __host__ __device__ REAL3 Cross(REAL3 a, REAL3 b)
{
	return make_REAL3(a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x);
}
static __inline__ __host__ __device__ REAL Dot(REAL* a, REAL* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
static __inline__ __host__ __device__ REAL3 Cross(REAL* a, REAL* b)
{
	return make_REAL3(a[1] * b[2] - a[2] * b[1],
		a[2] * b[0] - a[0] * b[2],
		a[0] * b[1] - a[1] * b[0]);
}
static __inline__ __host__ __device__ void Rotate(REAL3 & v, float degreeX, float degreeY) {
	float cosTheta = cosf(degreeX * M_PI / 180.0f);
	float sinTheta = sinf(degreeX * M_PI / 180.0f);
	v.y = v.y * cosTheta - v.z * sinTheta;
	v.z = v.y * sinTheta + v.z * cosTheta;

	cosTheta = cosf(degreeY * M_PI / 180.0f);
	sinTheta = sinf(degreeY * M_PI / 180.0f);
	v.x = v.x * cosTheta + v.z * sinTheta;
	v.z = -v.x * sinTheta + v.z * cosTheta;
}
static __inline__ __host__ __device__ REAL Squared(REAL num)
{
	return num * num;
}
static __inline__ __host__ __device__ uint Log2(uint num)
{
	uint k = 2, n = 0;
	while (k << n <= num) n++;
	return n;
}
static __inline__ __host__ __device__ uint MaxBinary(uint num)
{
	uint n = 1;
	while (n < num)
		n = n << 1;
	return n;
}
//-------------------------------------------------------------------------------------------------------------
struct AABB {
	REAL3 _min;
	REAL3 _max;
};
static __inline__ __host__ __device__ void printAABB(const AABB & a) {
	printf("Min: (%f, %f, %f)\nMax: (%f, %f, %f)\n", a._min.x, a._min.y, a._min.z, a._max.x, a._max.y, a._max.z);
}
inline __host__ __device__ __forceinline__ void resetAABB(AABB& aabb) {
	aabb._min.x = aabb._min.y = aabb._min.z = DBL_MAX;
	aabb._max.x = aabb._max.y = aabb._max.z = -DBL_MAX;
}
static __inline__ __host__ __device__ void setAABB(AABB& a, const REAL3& min, const REAL3& max) {
	a._min = min;
	a._max = max;
}
static __inline__ __host__ __device__ void setAABB(AABB& a, const REAL3& cen, const REAL delta) {
	a._min = cen - delta;
	a._max = cen + delta;
}
static __inline__ __host__ __device__ void setAABB(AABB& a, const AABB& b) {
	a._min = b._min;
	a._max = b._max;
}
static __inline__ __host__ __device__ void addAABB(AABB& a, const REAL3& x) {
	if (a._min.x > x.x)
		a._min.x = x.x;
	else if (a._max.x < x.x)
		a._max.x = x.x;
	if (a._min.y > x.y)
		a._min.y = x.y;
	else if (a._max.y < x.y)
		a._max.y = x.y;
	if (a._min.z > x.z)
		a._min.z = x.z;
	else if (a._max.z < x.z)
		a._max.z = x.z;
}
static __inline__ __host__ __device__ void addAABB(AABB& a, const REAL3& x, REAL delta) {
	if (a._min.x > x.x - delta)
		a._min.x = x.x - delta;
	else if (a._max.x < x.x + delta)
		a._max.x = x.x + delta;
	if (a._min.y > x.y - delta)
		a._min.y = x.y - delta;
	else if (a._max.y < x.y + delta)
		a._max.y = x.y + delta;
	if (a._min.z > x.z - delta)
		a._min.z = x.z - delta;
	else if (a._max.z < x.z + delta)
		a._max.z = x.z + delta;
}
static __inline__ __host__ __device__ void addAABB(AABB& a, const AABB& x) {
	if (a._min.x > x._min.x)
		a._min.x = x._min.x;
	else if (a._max.x < x._max.x)
		a._max.x = x._max.x;
	if (a._min.y > x._min.y)
		a._min.y = x._min.y;
	else if (a._max.y < x._max.y)
		a._max.y = x._max.y;
	if (a._min.z > x._min.z)
		a._min.z = x._min.z;
	else if (a._max.z < x._max.z)
		a._max.z = x._max.z;
}
inline  __host__ __device__ __forceinline__ bool intersect(const AABB& a, const AABB& b) {
	return a._min.x <= b._max.x
		&& a._min.y <= b._max.y
		&& a._min.z <= b._max.z
		&& a._max.x >= b._min.x
		&& a._max.y >= b._min.y
		&& a._max.z >= b._min.z;
}

#endif