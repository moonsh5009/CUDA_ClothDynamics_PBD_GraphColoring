#ifndef __REPRESENTIVE_TRIANGLE_H__
#define __REPRESENTIVE_TRIANGLE_H__

#pragma once
#include "Mesh.h"

inline __host__ __device__ __forceinline__ void setRTriVertex(uint& info, uint n) {
	info |= (1u << n);
}
inline __host__ __device__ __forceinline__ void setRTriEdge(uint& info, uint n) {
	info |= (1u << (n + 3u));
}
inline __device__ __forceinline__ bool RTriVertex(const uint info, uint n) {
	return (info >> n) & 1u;
}
inline __device__ __forceinline__ bool RTriEdge(const uint info, uint n) {
	return (info >> (n + 3u)) & 1u;
}

struct RTriParam {
	uint* _info;
	uint _size;
};

class RTriangle {
public:
	Dvector<uint> _rtris;
public:
	RTriangle() {};
	~RTriangle() {};
public:
	inline RTriParam param(void) const {
		RTriParam p;
		p._info = _rtris._list;
		p._size = _rtris.size();
		return p;
	}
public:
	void init(Dvector<uint>& fs, DPrefixArray<uint> nbFs);
};

#endif