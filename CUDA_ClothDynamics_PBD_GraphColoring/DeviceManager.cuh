#ifndef __DEVICE_MANAGER_CUH__
#define __DEVICE_MANAGER_CUH__

#pragma once
#include "DeviceManager.h"

//#define MIN(x, y)			x < y ? x : y
//#define MAX(x, y)			x > y ? x : y
#define MIN(x, y)			min(x, y)
#define MAX(x, y)			max(x, y)

static __inline__ __device__ __forceinline__ REAL atomicMax_REAL(REAL* address, REAL val)
{
	REAL_INT ret = REAL_AS_INT(*address);
	while (val > INT_AS_REAL(ret))
	{
		REAL_INT old = ret;
		if ((ret = atomicCAS((REAL_INT*)address, old, REAL_AS_INT(val))) == old)
			break;
	}
	return INT_AS_REAL(ret);
}
static __inline__ __device__ __forceinline__ REAL atomicMin_REAL(REAL* address, REAL val)
{
	REAL_INT ret = REAL_AS_INT(*address);
	while (val < INT_AS_REAL(ret))
	{
		REAL_INT old = ret;
		if ((ret = atomicCAS((REAL_INT*)address, old, REAL_AS_INT(val))) == old)
			break;
	}
	return INT_AS_REAL(ret);
}
static __inline__ __device__ __forceinline__ REAL atomicAdd_REAL(REAL* address, REAL val)
{
	REAL_INT* address_as_ull = (REAL_INT*)address;
	REAL_INT old = *address_as_ull, assumed;

	do {
		assumed = old;
		old = atomicCAS(address_as_ull, assumed,
			REAL_AS_INT(val + INT_AS_REAL(assumed)));
	} while (assumed != old);

	return INT_AS_REAL(old);
}
static __inline__ __device__ __forceinline__ void atomicPairAdd_REAL(
	REAL_INT* address0, REAL_INT* address1, REAL val0, REAL val1)
{
	REAL_INT old0 = *address0, assumed0;
	REAL_INT old1 = *address0, assumed1;
	bool flag0 = true; bool flag1 = true;
	do {
		if (flag0) {
			assumed0 = old0;
			old0 = atomicCAS(address0, assumed0,
				REAL_AS_INT(val0 + INT_AS_REAL(assumed0)));

			if (assumed0 == old0)
				flag0 = false;
		}
		if (flag1) {
			assumed1 = old1;
			old1 = atomicCAS(address1, assumed1,
				REAL_AS_INT(val1 + INT_AS_REAL(assumed1)));

			if (assumed1 == old1)
				flag1 = false;
		}
	} while (flag0 || flag1);
}
static __inline__ __device__ __forceinline__ REAL atomicExch_REAL(REAL* address, REAL val)
{
	REAL_INT* address_as_ull = (REAL_INT*)address;
	REAL_INT old = *address_as_ull, assumed;

	do {
		assumed = old;
		old = atomicCAS(address_as_ull, assumed, REAL_AS_INT(val));
	} while (assumed != old);

	return INT_AS_REAL(old);
}
static __inline__ __device__ __forceinline__ void atomicAdd_REAL3(REAL* address, REAL* val)
{
	/*REAL* raddress = (REAL*)address;
	REAL3 old;
	old.x = atomicAdd_REAL(raddress + 0, val.x);
	old.y = atomicAdd_REAL(raddress + 1, val.y);
	old.z = atomicAdd_REAL(raddress + 2, val.z);
	return old;*/
	REAL_INT* address_as_ull = (REAL_INT*)address;
	REAL_INT old, assumed;

	uint i, n = 7u;
	do {
		for (i = 0u; i < 3u; i++) {
			if ((n >> i)) {
				assumed = *address_as_ull;
				old = atomicCAS(address_as_ull, assumed,
					REAL_AS_INT(val[i] + INT_AS_REAL(assumed)));
				n &= ((uint)(assumed != old)) << i;
			}
		}
	} while (n);
}
static __inline__ __device__ __forceinline__ uint atomicMax_flag(uint* address, uint val)
{
	uint ret = *address;
	uint old;
	while (val > ret || ret == 0xffffffff)
	{
		old = ret;
		if ((ret = atomicCAS(address, old, val)) == old)
			break;
	}
	return ret;
}
static __inline__ __device__ __forceinline__ REAL MyAtomicOr(uint* address, uint val)
{
	uint ret = *address;
	while ((ret & val) != val)
	{
		uint old = ret;
		if ((ret = atomicCAS(address, old, old | val)) == old)
			break;
	}
	return ret;
}
static __inline__ __device__ __forceinline__ void cudaLock(uint* mutex, uint i) {
	while (atomicCAS(mutex + i, 0, 1) != 0);
}
static __inline__ __device__ __forceinline__ void cudaUnLock(uint* mutex, uint i) {
	atomicExch(mutex + i, 0);
}

static __inline__ __device__ __forceinline__ void warpSum(volatile uint* s_data, uint tid) {
	s_data[tid] += s_data[tid + 32];
	s_data[tid] += s_data[tid + 16];
	s_data[tid] += s_data[tid + 8];
	s_data[tid] += s_data[tid + 4];
	s_data[tid] += s_data[tid + 2];
	s_data[tid] += s_data[tid + 1];
}
static __inline__ __device__ __forceinline__ void warpSum(volatile REAL* s_data, uint tid) {
	s_data[tid] += s_data[tid + 32];
	s_data[tid] += s_data[tid + 16];
	s_data[tid] += s_data[tid + 8];
	s_data[tid] += s_data[tid + 4];
	s_data[tid] += s_data[tid + 2];
	s_data[tid] += s_data[tid + 1];
}
static __inline__ __device__ __forceinline__ void warpMin(volatile uint* s_data, uint tid) {
	if (s_data[tid] > s_data[tid + 32])
		s_data[tid] = s_data[tid + 32];
	if (s_data[tid] > s_data[tid + 16])
		s_data[tid] = s_data[tid + 16];
	if (s_data[tid] > s_data[tid + 8])
		s_data[tid] = s_data[tid + 8];
	if (s_data[tid] > s_data[tid + 4])
		s_data[tid] = s_data[tid + 4];
	if (s_data[tid] > s_data[tid + 2])
		s_data[tid] = s_data[tid + 2];
	if (s_data[tid] > s_data[tid + 1])
		s_data[tid] = s_data[tid + 1];
}
static __inline__ __device__ __forceinline__ void warpMin(volatile REAL* s_data, uint tid) {
	if (s_data[tid] > s_data[tid + 32])
		s_data[tid] = s_data[tid + 32];
	if (s_data[tid] > s_data[tid + 16])
		s_data[tid] = s_data[tid + 16];
	if (s_data[tid] > s_data[tid + 8])
		s_data[tid] = s_data[tid + 8];
	if (s_data[tid] > s_data[tid + 4])
		s_data[tid] = s_data[tid + 4];
	if (s_data[tid] > s_data[tid + 2])
		s_data[tid] = s_data[tid + 2];
	if (s_data[tid] > s_data[tid + 1])
		s_data[tid] = s_data[tid + 1];
}
static __inline__ __device__ __forceinline__ void warpMax(volatile uint* s_data, uint tid) {
	if (s_data[tid] < s_data[tid + 32])
		s_data[tid] = s_data[tid + 32];
	if (s_data[tid] < s_data[tid + 16])
		s_data[tid] = s_data[tid + 16];
	if (s_data[tid] < s_data[tid + 8])
		s_data[tid] = s_data[tid + 8];
	if (s_data[tid] < s_data[tid + 4])
		s_data[tid] = s_data[tid + 4];
	if (s_data[tid] < s_data[tid + 2])
		s_data[tid] = s_data[tid + 2];
	if (s_data[tid] < s_data[tid + 1])
		s_data[tid] = s_data[tid + 1];
}
static __inline__ __device__ __forceinline__ void warpMax(volatile REAL* s_data, uint tid) {
	if (s_data[tid] < s_data[tid + 32])
		s_data[tid] = s_data[tid + 32];
	if (s_data[tid] < s_data[tid + 16])
		s_data[tid] = s_data[tid + 16];
	if (s_data[tid] < s_data[tid + 8])
		s_data[tid] = s_data[tid + 8];
	if (s_data[tid] < s_data[tid + 4])
		s_data[tid] = s_data[tid + 4];
	if (s_data[tid] < s_data[tid + 2])
		s_data[tid] = s_data[tid + 2];
	if (s_data[tid] < s_data[tid + 1])
		s_data[tid] = s_data[tid + 1];
}

#endif