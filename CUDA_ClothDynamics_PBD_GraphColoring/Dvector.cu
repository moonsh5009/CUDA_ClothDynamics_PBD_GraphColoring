#include "Dvector.h"
#include "DeviceManager.cuh"

__global__ void rotateX_kernel(REAL* x, size_t size, REAL radian, REAL c, REAL s) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= size)
		return;
	uint id0 = id * 3u + 1u;
	uint id1 = id * 3u + 2u;
	REAL a = x[id0];
	REAL b = x[id1];
	a = a * c - b * s;
	b = a * s + b * c;
	x[id0] = a;
	x[id1] = b;
}

__global__ void vectorAdd_kernel(REAL* x, size_t size, REAL delta) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= size)
		return;
	REAL a = x[id];
	a += delta;
	x[id] = a;
}
__global__ void vectorSub_kernel(REAL* x, size_t size, REAL delta) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= size)
		return;
	REAL a = x[id];
	a -= delta;
	x[id] = a;
}
__global__ void vectorMulti_kernel(REAL* x, size_t size, REAL delta) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= size)
		return;
	REAL a = x[id];
	a *= delta;
	x[id] = a;
}
__global__ void vectorDivide_kernel(REAL* x, size_t size, REAL delta) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= size)
		return;
	REAL a = x[id];
	a /= delta;
	x[id] = a;
}

__global__ void vectorAdd_kernel(REAL* x, size_t size, REAL3 delta) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= size)
		return;
	REAL a = x[id];
	REAL d = (&delta.x)[id % 3];
	a += d;
	x[id] = a;
}
__global__ void vectorSub_kernel(REAL* x, size_t size, REAL3 delta) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= size)
		return;
	REAL a = x[id];
	REAL d = (&delta.x)[id % 3];
	a -= d;
	x[id] = a;
}
__global__ void vectorMulti_kernel(REAL* x, size_t size, REAL3 delta) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= size)
		return;
	REAL a = x[id];
	REAL d = (&delta.x)[id % 3];
	a *= d;
	x[id] = a;
}
__global__ void vectorDivide_kernel(REAL* x, size_t size, REAL3 delta) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= size)
		return;
	REAL a = x[id];
	REAL d = (&delta.x)[id % 3];
	a /= d;
	x[id] = a;
}

__global__ void vectorAdd_kernel(REAL* x, size_t size, REAL* y) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= size)
		return;
	REAL a = x[id];
	REAL b = y[id];
	a += b;
	x[id] = a;
}
__global__ void vectorSub_kernel(REAL* x, size_t size, REAL* y) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= size)
		return;
	REAL a = x[id];
	REAL b = y[id];
	a -= b;
	x[id] = a;
}
__global__ void vectorMulti_kernel(REAL* x, size_t size, REAL* y) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= size)
		return;
	REAL a = x[id];
	REAL b = y[id];
	a *= b;
	x[id] = a;
}
__global__ void vectorDivide_kernel(REAL* x, size_t size, REAL* y) {
	uint id = threadIdx.x + blockDim.x * blockIdx.x;
	if (id >= size)
		return;
	REAL a = x[id];
	REAL b = y[id];
	a /= b;
	x[id] = a;
}
__global__ void vectorMax_kernel(uint* X, size_t size, uint* result) {
	extern __shared__ uint s_maxs[];
	uint id = threadIdx.x + (blockDim.x * blockIdx.x << 1u);
	if (id < size) {
		uint x = X[id];
		if (id + blockDim.x < size) {
			uint tmp = X[id + blockDim.x];
			if (x < tmp) x = tmp;
		}
		s_maxs[threadIdx.x] = x;
	}
	else s_maxs[threadIdx.x] = 0.0;

	for (uint s = blockDim.x >> 1u; s > 32u; s >>= 1u) {
		__syncthreads();
		if (threadIdx.x < s)
			if (s_maxs[threadIdx.x] < s_maxs[threadIdx.x + s])
				s_maxs[threadIdx.x] = s_maxs[threadIdx.x + s];
	}
	__syncthreads();
	if (threadIdx.x < 32) {
		warpMax(s_maxs, threadIdx.x);
		if (threadIdx.x == 0)
			atomicMax(result, s_maxs[0]);
	}
}

void rotateX(Dvector<REAL> v, REAL degree) {
	REAL radian = degree * 3.141592 / 180.0;
	REAL c = cos(radian);
	REAL s = sin(radian);
	size_t length = v.size();
	rotateX_kernel << <divup(length, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		v.begin(), length, radian, c, s);
	CUDA_CHECK(cudaPeekAtLastError());
}

Dvector<REAL> operator+(const Dvector<REAL>& a, const REAL b) {
	Dvector<REAL> result;
	result = a;
	uint length = result.size();
	vectorAdd_kernel << <divup(length, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		result.begin(), length, b);
	CUDA_CHECK(cudaPeekAtLastError());
	return result;
}
Dvector<REAL> operator-(const Dvector<REAL>& a, const REAL b) {
	Dvector<REAL> result;
	result = a;
	uint length = result.size();
	vectorSub_kernel << <divup(length, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		result.begin(), length, b);
	CUDA_CHECK(cudaPeekAtLastError());
	return result;
}
Dvector<REAL> operator*(const Dvector<REAL>& a, const REAL b) {
	Dvector<REAL> result;
	result = a;
	uint length = result.size();
	vectorMulti_kernel << <divup(length, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		result.begin(), length, b);
	CUDA_CHECK(cudaPeekAtLastError());
	return result;
}
Dvector<REAL> operator/(const Dvector<REAL>& a, const REAL b) {
	Dvector<REAL> result;
	result = a;
	uint length = result.size();
	vectorDivide_kernel << <divup(length, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		result.begin(), length, b);
	CUDA_CHECK(cudaPeekAtLastError());
	return result;
}

Dvector<REAL> operator+(const Dvector<REAL>& a, const REAL3 b) {
	Dvector<REAL> result;
	result = a;
	uint length = result.size();
	vectorAdd_kernel << <divup(length, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		result.begin(), length, b);
	CUDA_CHECK(cudaPeekAtLastError());
	return result;
}
Dvector<REAL> operator-(const Dvector<REAL>& a, const REAL3 b) {
	Dvector<REAL> result;
	result = a;
	uint length = result.size();
	vectorSub_kernel << <divup(length, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		result.begin(), length, b);
	CUDA_CHECK(cudaPeekAtLastError());
	return result;
}
Dvector<REAL> operator*(const Dvector<REAL>& a, const REAL3 b) {
	Dvector<REAL> result;
	result = a;
	uint length = result.size();
	vectorMulti_kernel << <divup(length, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		result.begin(), length, b);
	CUDA_CHECK(cudaPeekAtLastError());
	return result;
}
Dvector<REAL> operator/(const Dvector<REAL>& a, const REAL3 b) {
	Dvector<REAL> result;
	result = a;
	uint length = result.size();
	vectorDivide_kernel << <divup(length, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		result.begin(), length, b);
	CUDA_CHECK(cudaPeekAtLastError());
	return result;
}

Dvector<REAL> operator+(const Dvector<REAL>& a, const Dvector<REAL>& b) {
	Dvector<REAL> result;
	result = a;
	uint length = result.size();
	vectorAdd_kernel << <divup(length, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		result.begin(), length, b.begin());
	CUDA_CHECK(cudaPeekAtLastError());
	return result;
}
Dvector<REAL> operator-(const Dvector<REAL>& a, const Dvector<REAL>& b) {
	Dvector<REAL> result;
	result = a;
	uint length = result.size();
	vectorSub_kernel << <divup(length, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		result.begin(), length, b.begin());
	CUDA_CHECK(cudaPeekAtLastError());
	return result;
}
Dvector<REAL> operator*(const Dvector<REAL>& a, const Dvector<REAL>& b) {
	Dvector<REAL> result;
	result = a;
	uint length = result.size();
	vectorMulti_kernel << <divup(length, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		result.begin(), length, b.begin());
	CUDA_CHECK(cudaPeekAtLastError());
	return result;

}
Dvector<REAL> operator/(const Dvector<REAL>& a, const Dvector<REAL>& b) {
	Dvector<REAL> result;
	result = a;
	uint length = result.size();
	vectorDivide_kernel << <divup(length, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		result.begin(), length, b.begin());
	CUDA_CHECK(cudaPeekAtLastError());
	return result;
}
void operator+=(Dvector<REAL>& a, const REAL b) {
	uint length = a.size();
	vectorAdd_kernel << <divup(length, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		a.begin(), length, b);
	CUDA_CHECK(cudaPeekAtLastError());
}
void operator*=(Dvector<REAL>& a, const REAL b) {
	uint length = a.size();
	vectorMulti_kernel << <divup(length, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		a.begin(), length, b);
	CUDA_CHECK(cudaPeekAtLastError());
}
void operator+=(Dvector<REAL>& a, const REAL3 b) {
	uint length = a.size();
	vectorAdd_kernel << <divup(length, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		a.begin(), length, b);
	CUDA_CHECK(cudaPeekAtLastError());
}
void operator*=(Dvector<REAL>& a, const REAL3 b) {
	uint length = a.size();
	vectorMulti_kernel << <divup(length, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		a.begin(), length, b);
	CUDA_CHECK(cudaPeekAtLastError());
}
void operator+=(Dvector<REAL>& a, const Dvector<REAL>& b) {
	uint length = a.size();
	vectorAdd_kernel << <divup(length, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		a.begin(), length, b.begin());
	CUDA_CHECK(cudaPeekAtLastError());
}
void operator-=(Dvector<REAL>& a, const Dvector<REAL>& b) {
	uint length = a.size();
	vectorSub_kernel << <divup(length, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		a.begin(), length, b.begin());
	CUDA_CHECK(cudaPeekAtLastError());
}
void operator*=(Dvector<REAL>& a, const Dvector<REAL>& b) {
	uint length = a.size();
	vectorMulti_kernel << <divup(length, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		a.begin(), length, b.begin());
	CUDA_CHECK(cudaPeekAtLastError());
}
void operator/=(Dvector<REAL>& a, const Dvector<REAL>& b) {
	uint length = a.size();
	vectorDivide_kernel << <divup(length, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		a.begin(), length, b.begin());
	CUDA_CHECK(cudaPeekAtLastError());
}
void getDvectorMax(const Dvector<uint>& X, uint* deviceResult) {
	uint length = X.size();
	vectorMax_kernel << <divup(length, MAX_BLOCKSIZE << 1u), MAX_BLOCKSIZE, MAX_BLOCKSIZE * sizeof(uint) >> > (
		X._list, length, deviceResult);
	CUDA_CHECK(cudaPeekAtLastError());
}