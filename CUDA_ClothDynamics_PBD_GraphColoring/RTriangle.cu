#include "RTriangle.cuh"

void RTriangle::init(Dvector<uint>& fs, DPrefixArray<uint> nbFs) {
	printf("CUDA RTri Build: ");
	CUDA_CHECK(cudaDeviceSynchronize());
	ctimer timer = CNOW;

	uint numFaces = fs.size() / 3u;
	_rtris.resize(numFaces);
	RTriBuild_Kernel << <divup(numFaces, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		fs(), nbFs._array(), nbFs._index(), param());
	CUDA_CHECK(cudaPeekAtLastError());

	CUDA_CHECK(cudaDeviceSynchronize());
	printf("%lf msec\n", (CNOW - timer) / 10000.0);
}