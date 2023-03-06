#include "BVH.cuh"

struct TriInfo_CMP
{
	__host__ __device__
		bool operator()(const TriInfo& a, const TriInfo& b) {
		if (a._id != b._id)
			return a._id < b._id;
		return a._pos < b._pos;
	}
};

void BVH::build(Dvector<uint>& fs, Dvector<REAL>& ns, const DPrefixArray<uint>& nbFs) {
	CUDA_CHECK(cudaDeviceSynchronize());
	ctimer timer = CNOW;

	uint numFaces = fs.size() / 3u;
	Dvector<TriInfo> infos(numFaces);

	initBVHTreeDevice(numFaces);
	initBVHInfo_kernel << < divup(numFaces, BLOCKSIZE), BLOCKSIZE >> > (
		fs(), ns(), infos(), param());
	CUDA_CHECK(cudaPeekAtLastError());
	InitMinMaxKernel << < divup(numFaces, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		param());
	CUDA_CHECK(cudaPeekAtLastError());
	updateBVHInfo_kernel << < divup(numFaces, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		infos(), param());
	CUDA_CHECK(cudaPeekAtLastError());
	thrust::sort(thrust::device_ptr<TriInfo>(infos.begin()),
		thrust::device_ptr<TriInfo>(infos.end()), TriInfo_CMP());

	for (uint level = 0u; level < _maxLevel - 2u; level++) {
		uint size = 1u << level;
		updateMinMax_kernel << < divup(size, BLOCKSIZE), BLOCKSIZE >> > (
			param(), level, size);
		CUDA_CHECK(cudaPeekAtLastError());
		subdivBVH_kernel << < divup(numFaces, BLOCKSIZE), BLOCKSIZE >> > (
			infos(), param());
		CUDA_CHECK(cudaPeekAtLastError());
		updateBVHInfo_kernel << < divup(numFaces, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
			infos(), param());
		CUDA_CHECK(cudaPeekAtLastError());
		thrust::sort(thrust::device_ptr<TriInfo>(infos.begin()),
			thrust::device_ptr<TriInfo>(infos.end()), TriInfo_CMP());
		/*{
			vector<TriInfo> hinfo;
			infos.copyToHost(hinfo);
			printf("%d: ", level);
			for (auto t : hinfo) {
				printf("%d ", t._id);
			}
			printf("\n");
		}*/
	}
	buildBVH_kernel << < divup(_size, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		infos(), param());
	CUDA_CHECK(cudaPeekAtLastError());

	{
		/*vector<TriInfo> hinfo;
		infos.copyToHost(hinfo);
		for (auto t : hinfo) {
			printf("%d\n", t._id);
		}
		printf("\n");*/
	}

	_test = 0u;

	CUDA_CHECK(cudaDeviceSynchronize());
	printf("BVH Build: %lf msec\n", (CNOW - timer) / 10000.0);
}
void BVH::refit(ObjParam obj, REAL delta, const REAL dt, bool isCCD) {
	//CUDA_CHECK(cudaDeviceSynchronize());
	//ctimer timer = CNOW;

	if (!isCCD)	delta *= 0.5;

	uint currLevel = _maxLevel - 1u;
	int num = 1u << (currLevel--);
	RefitLeafBVHKernel << < divup(num, REFIT_BLOCKSIZE), REFIT_BLOCKSIZE >> >
		(obj, param(), num, delta, dt, isCCD);
	CUDA_CHECK(cudaPeekAtLastError());

	while (currLevel > 10u) {
		currLevel--;
		num >>= 1;
		RefitBVHKernel << < divup(num, REFIT_BLOCKSIZE), REFIT_BLOCKSIZE >> >
			(obj, param(), num, delta, dt, isCCD);
		CUDA_CHECK(cudaPeekAtLastError());
	}
	RefitNodeBVHKernel << < 1u, MAX_BLOCKSIZE >> >
		(obj, param(), currLevel, delta, dt, isCCD);
	CUDA_CHECK(cudaPeekAtLastError());

	//CUDA_CHECK(cudaDeviceSynchronize());
	//printf("BVH Refit: %lf msec\n", (CNOW - timer) / 10000.0);
}
void BVH::draw(const AABB& aabb) {
	glDisable(GL_LIGHTING);
	glPushMatrix();
	glLineWidth(1.0f);
	//glColor3d(0.6, 0.6, 0.6);
	glColor3d(1.0, 0.0, 0.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glBegin(GL_LINES);
	glVertex3d(aabb._min.x, aabb._min.y, aabb._min.z);
	glVertex3d(aabb._min.x, aabb._min.y, aabb._max.z);
	glVertex3d(aabb._min.x, aabb._max.y, aabb._min.z);
	glVertex3d(aabb._min.x, aabb._max.y, aabb._max.z);
	glVertex3d(aabb._max.x, aabb._min.y, aabb._min.z);
	glVertex3d(aabb._max.x, aabb._min.y, aabb._max.z);
	glVertex3d(aabb._max.x, aabb._max.y, aabb._min.z);
	glVertex3d(aabb._max.x, aabb._max.y, aabb._max.z);
	glEnd();
	glTranslated(0, 0, aabb._min.z);
	glRectd(aabb._min.x, aabb._min.y, aabb._max.x, aabb._max.y);
	glTranslated(0, 0, aabb._max.z - aabb._min.z);
	glRectd(aabb._min.x, aabb._min.y, aabb._max.x, aabb._max.y);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glPopMatrix();
	glEnable(GL_LIGHTING);
}
void BVH::draw(void) {
	vector<REAL> mins[3];
	vector<REAL> maxs[3];
	vector<uint> levels;
	for (uint i = 0u; i < 3u; i++) {
		_mins[i].copyToHost(mins[i]);
		_maxs[i].copyToHost(maxs[i]);
	}
	_levels.copyToHost(levels);
	//uint testlevel = levels[_test];
	//uint testpath = _test - (1u << levels[_test]) + 1;
	for (uint i = 0u; i < _size; i++) {
		AABB aabb;
		aabb._min = make_REAL3(mins[0][i], mins[1][i], mins[2][i]);
		aabb._max = make_REAL3(maxs[0][i], maxs[1][i], maxs[2][i]);
		draw(aabb);
		/*if (levels[i] > testlevel && i > _size - _numFaces) {
			uint path = i - (1u << levels[i]) + 1;
			path >>= (levels[i] - testlevel);
			if (path == testpath) {
				draw(aabb);
			}
		}
		if (i == _test) {
			draw(aabb);
		}*/
	}
}