#ifndef __CUDA_PREFIX_ARRAY_H__
#define __CUDA_PREFIX_ARRAY_H__

#pragma once
#include <set>
#include "Dvector.h"

template <typename T>
class DPrefixArray;

template <typename T>
class PrefixArray {
public:
	vector<T>		_array;
	vector<uint>	_index;
public:
	PrefixArray() {}
	virtual ~PrefixArray() {}
public:
	inline void init(vector<vector<T>>& v) {
		uint indexSize = v.size();
		_index.resize(indexSize + 1u);
		_index[0] = 0u;
		uint num = 0u;
		for (uint i = 0u; i < indexSize; i++) {
			uint unitSize = v[i].size();
			num += unitSize;
			_index[i + 1u] = num;
		}
		_array.resize(num);
		num = 0u;
		for (uint i = 0u; i < indexSize; i++) {
			uint ino = _index[i];
			for (auto x : v[i])
				_array[num++] = x;
		}
	}
	inline void init(vector<set<T>>& v) {
		uint indexSize = v.size();
		_index.resize(indexSize + 1u);
		_index[0] = 0u;
		uint num = 0u;
		for (uint i = 0u; i < indexSize; i++) {
			uint unitSize = v[i].size();
			num += unitSize;
			_index[i + 1u] = num;
		}
		_array.resize(num);
		num = 0u;
		for (uint i = 0u; i < indexSize; i++) {
			uint ino = _index[i];
			for (auto x : v[i])
				_array[num++] = x;
		}
	}
	inline size_t arraySize(void) {
		return _array.size();
	}
	inline size_t indexSize(void) {
		return _index.size();
	}
public:
	inline void copyToHost(PrefixArray<T>& host) {
		host._array = _array;
		host._index = _index;
	}
	inline void copyToDevice(DPrefixArray<T>& device) {
		device._array = _array;
		device._index = _index;
	}
	inline void copyFromHost(const PrefixArray<T>& host) {
		_array = host._array;
		_index = host._index;
	}
	inline void copyFromDevice(const DPrefixArray<T>& device) {
		device._array.copyToHost(_array);
		device._index.copyToHost(_index);
	}
public:
	inline PrefixArray<T>& operator=(const PrefixArray<T>& v) {
		copyFromHost(v);
		return *this;
	}
	inline PrefixArray<T>& operator=(const DPrefixArray<T>& v) {
		copyFromDevice(v);
		return *this;
	}
	/*__inline__ __host__ __device__ */
	inline T* operator[](size_t i) {
		if (i >= _index.size() - 1u) {
			printf("Error : DPrefixArray[] : index out\n");
			exit(1);
		}
		return &_array[_index[i]];
	}
};

template <typename T>
class DPrefixArray {
public:
	Dvector<T>		_array;
	Dvector<uint>	_index;
public:
	DPrefixArray() {}
	virtual ~DPrefixArray() {
		_array.freeChk();
		_index.freeChk();
	}
public:
	inline size_t arraySize(void) const {
		return _array.size();
	}
	inline size_t indexSize(void) const {
		return _index.size();
	}
public:
	inline void copyToHost(PrefixArray<T>& host) const {
		_array.copyToHost(host._array);
		_index.copyToHost(host._index);
	}
	inline void copyToDevice(DPrefixArray<T>& device) const {
		device._array = _array;
		device._index = _index;
	}
	inline void copyFromHost(const PrefixArray<T>& host) {
		_array = host._array;
		_index = host._index;
	}
	inline void copyFromDevice(const DPrefixArray<T>& device) {
		_array = device._array;
		_index = device._index;
	}
	inline void copyToHost(PrefixArray<T>& host, cudaStream_t* s) const {
		_array.copyToHost(host._array, s[0]);
		_index.copyToHost(host._index, s[1]);
	}
	inline void copyToDevice(DPrefixArray<T>& device, cudaStream_t* s) const {
		_array.copyToDevice(device._array, s[0]);
		_index.copyToDevice(device._index, s[1]);
	}
	inline void copyFromHost(const PrefixArray<T>& host, cudaStream_t* s) {
		_array.copyFromHost(host._array, s[0]);
		_index.copyFromHost(host._index, s[1]);
	}
	inline void copyFromDevice(const DPrefixArray<T>& device, cudaStream_t* s) {
		_array.copyFromDevice(device._array, s[0]);
		_index.copyFromDevice(device._index, s[1]);
	}
public:
	inline DPrefixArray<T>& operator=(const PrefixArray<T>& v) {
		copyFromHost(v);
		return *this;
	}
	inline DPrefixArray<T>& operator=(const DPrefixArray<T>& v) {
		copyFromDevice(v);
		return *this;
	}
	/*__inline__ __host__ __device__ */
	/*inline T* operator[](size_t i) {
		if (i >= _index._size - 1u) {
			printf("Error : DPrefixArray[] : index out\n");
			exit(1);
		}
		return &_array[_list[i]];
	}*/
};

#endif