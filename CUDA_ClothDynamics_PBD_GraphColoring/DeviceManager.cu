//#include "DeviceManager.cuh"
//
//__global__ void Sum_kernel(REAL* X, size_t size, REAL* result) {
//	extern __shared__ REAL s_sums[];
//	uint id = threadIdx.x + (blockDim.x * blockIdx.x << 1u);
//	if (id < size) {
//		REAL x = X[id];
//		if (id + blockDim.x < size) {
//			REAL tmp = X[id + blockDim.x];
//			x += tmp;
//		}
//		s_sums[threadIdx.x] = x;
//	}
//	else s_sums[threadIdx.x] = (REAL)0.0;
//
//	for (uint s = blockDim.x >> 1u; s > 32u; s >>= 1u) {
//		__syncthreads();
//		if (threadIdx.x < s)
//			s_sums[threadIdx.x] += s_sums[threadIdx.x + s];
//	}
//	__syncthreads();
//	if (threadIdx.x < 32) {
//		warpSum(s_sums, threadIdx.x);
//		if (threadIdx.x == 0)
//			atomicAdd_REAL(result, s_sums[0]);
//	}
//}
//
//__global__ void Min_kernel(REAL* X, size_t size, REAL* result) {
//	extern __shared__ REAL s_mins[];
//	uint id = threadIdx.x + (blockDim.x * blockIdx.x << 1u);
//	if (id < size) {
//		REAL x = X[id];
//		if (id + blockDim.x < size) {
//			REAL tmp = X[id + blockDim.x];
//			if (x > tmp) x = tmp;
//		}
//		s_mins[threadIdx.x] = x;
//	}
//	else s_mins[threadIdx.x] = REAL_MAX;
//
//	for (uint s = blockDim.x >> 1u; s > 32u; s >>= 1u) {
//		__syncthreads();
//		if (threadIdx.x < s)
//			if (s_mins[threadIdx.x] > s_mins[threadIdx.x + s])
//				s_mins[threadIdx.x] = s_mins[threadIdx.x + s];
//	}
//	__syncthreads();
//	if (threadIdx.x < 32) {
//		warpMin(s_mins, threadIdx.x);
//		if (threadIdx.x == 0)
//			atomicMin_REAL(result, s_mins[0]);
//	}
//}
//__global__ void AbsMin_kernel(REAL* X, size_t size, REAL* result) {
//	extern __shared__ REAL s_mins[];
//	uint id = threadIdx.x + (blockDim.x * blockIdx.x << 1u);
//	if (id < size) {
//		REAL x = X[id];
//		x = abs(x);
//		if (id + blockDim.x < size) {
//			REAL tmp = X[id + blockDim.x];
//			tmp = abs(x);
//			if (x > tmp) x = tmp;
//		}
//		s_mins[threadIdx.x] = x;
//	}
//	else s_mins[threadIdx.x] = REAL_MAX;
//
//	for (uint s = blockDim.x >> 1u; s > 32u; s >>= 1u) {
//		__syncthreads();
//		if (threadIdx.x < s)
//			if (s_mins[threadIdx.x] > s_mins[threadIdx.x + s])
//				s_mins[threadIdx.x] = s_mins[threadIdx.x + s];
//	}
//	__syncthreads();
//	if (threadIdx.x < 32) {
//		warpMin(s_mins, threadIdx.x);
//		if (threadIdx.x == 0)
//			atomicMin_REAL(result, s_mins[0]);
//	}
//}
//__global__ void Max_REAL_kernel(REAL* X, size_t size, REAL* result) {
//	extern __shared__ REAL s_maxs[];
//	uint id = threadIdx.x + (blockDim.x * blockIdx.x << 1u);
//	if (id < size) {
//		REAL x = X[id];
//		if (id + blockDim.x < size) {
//			REAL tmp = X[id + blockDim.x];
//			if (x < tmp) x = tmp;
//		}
//		s_maxs[threadIdx.x] = x;
//	}
//	else s_maxs[threadIdx.x] = 0.0;
//
//	for (uint s = blockDim.x >> 1u; s > 32u; s >>= 1u) {
//		__syncthreads();
//		if (threadIdx.x < s)
//			if (s_maxs[threadIdx.x] < s_maxs[threadIdx.x + s])
//				s_maxs[threadIdx.x] = s_maxs[threadIdx.x + s];
//	}
//	__syncthreads();
//	if (threadIdx.x < 32) {
//		warpMax(s_maxs, threadIdx.x);
//		if (threadIdx.x == 0)
//			atomicMax_REAL(result, s_maxs[0]);
//	}
//}
//__global__ void Max_uint_kernel(uint* X, size_t size, uint* result) {
//	extern __shared__ uint s_maxs[];
//	uint id = threadIdx.x + (blockDim.x * blockIdx.x << 1u);
//	if (id < size) {
//		uint x = X[id];
//		if (id + blockDim.x < size) {
//			uint tmp = X[id + blockDim.x];
//			if (x < tmp) x = tmp;
//		}
//		s_maxs[threadIdx.x] = x;
//	}
//	else s_maxs[threadIdx.x] = 0.0;
//
//	for (uint s = blockDim.x >> 1u; s > 32u; s >>= 1u) {
//		__syncthreads();
//		if (threadIdx.x < s)
//			if (s_maxs[threadIdx.x] < s_maxs[threadIdx.x + s])
//				s_maxs[threadIdx.x] = s_maxs[threadIdx.x + s];
//	}
//	__syncthreads();
//	if (threadIdx.x < 32) {
//		warpMax(s_maxs, threadIdx.x);
//		if (threadIdx.x == 0)
//			atomicMax(result, s_maxs[0]);
//	}
//}