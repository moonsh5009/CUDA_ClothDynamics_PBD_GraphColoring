#include "Constraint.cuh"

void Constraint::compColoring(const Dvector<uint>& neis, const Dvector<uint>& ineis)
{
	Dvector<uint> icurrs(_numSprings);
	_colors.resize(_numSprings);
	_colors.memset(0xffffffff);
	icurrs.memset(0);

	uint* d_tmp;
	CUDA_CHECK(cudaMalloc((void**)&d_tmp, sizeof(uint)));

	uint flag = 1u;
	//uint num = 0u;
	while (flag) {
		CUDA_CHECK(cudaMemset(d_tmp, 0, sizeof(uint)));
		compColoring_kernel << <divup(_numSprings, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
			neis(), ineis(), icurrs(), _colors(), _numSprings, d_tmp);
		CUDA_CHECK(cudaPeekAtLastError());
		CUDA_CHECK(cudaMemcpy(&flag, d_tmp, sizeof(uint), cudaMemcpyDeviceToHost));
		//num++;
	}
	//printf("%d\n", num);

	{
		CUDA_CHECK(cudaMemset(d_tmp, 0, sizeof(uint)));
		getDvectorMax(_colors, d_tmp);
		CUDA_CHECK(cudaMemcpy(&_colorSize, d_tmp, sizeof(uint), cudaMemcpyDeviceToHost));
		_colorSize++;
		printf("ColorSize %d\n", _colorSize);
	}

	CUDA_CHECK(cudaFree(d_tmp));
}

void Constraint::getEdgeNeighbors(
	const DPrefixArray<uint>& es, const DPrefixArray<uint>& nbNs,
	Dvector<uint>& neis, Dvector<uint>& ineis, uint numVertices) 
{
	ineis.resize(_numSprings + 1u);
	
	getEdgeNeisIds_kernel << <divup(_numSprings, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		nbNs._array(), nbNs._index(), ineis(), _inos(), _numSprings);
	CUDA_CHECK(cudaPeekAtLastError());
	cudaDeviceSynchronize();

	uint numNeis;
	thrust::inclusive_scan(thrust::device_ptr<uint>(ineis.begin()),
		thrust::device_ptr<uint>(ineis.end()), thrust::device_ptr<uint>(ineis.begin()));
	CUDA_CHECK(cudaMemcpy(&numNeis, ineis() + _numSprings, sizeof(uint), cudaMemcpyDeviceToHost));
	cudaDeviceSynchronize();

	neis.resize(numNeis);
	getEdgeNeis_kernel << <divup(_numSprings, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		es._array(), es._index(), nbNs._array(), nbNs._index(), neis(), ineis(), _inos(), _numSprings);
	CUDA_CHECK(cudaPeekAtLastError());

	/*{
		vector<uint> h_neis;
		vector<uint> h_ineis;
		neis.copyToHost(h_neis);
		ineis.copyToHost(h_ineis);
		for (int i = 0; i < _numSprings; i++) {
			printf("%d: ", i);
			for (int j = h_ineis[i]; j < h_ineis[i + 1]; j++) {
				printf("%d ", h_neis[j]);
			}
			printf("\n");
		}
	}*/

	//printf("%d\n", numNeis);
}
void Constraint::getFaceNeighbors(
	const DPrefixArray<uint>& nbFs,
	Dvector<uint>& neis, Dvector<uint>& ineis, uint numVertices)
{
	ineis.resize(_numSprings + 1u);

	getFaceNeisIds_kernel << <divup(_numSprings, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		nbFs._array(), nbFs._index(), ineis(), _inos(), _numSprings);
	CUDA_CHECK(cudaPeekAtLastError());
	cudaDeviceSynchronize();

	uint numNeis;
	thrust::inclusive_scan(thrust::device_ptr<uint>(ineis.begin()),
		thrust::device_ptr<uint>(ineis.end()), thrust::device_ptr<uint>(ineis.begin()));
	CUDA_CHECK(cudaMemcpy(&numNeis, ineis() + _numSprings, sizeof(uint), cudaMemcpyDeviceToHost));
	cudaDeviceSynchronize();

	neis.resize(numNeis);
	getFaceNeis_kernel << <divup(_numSprings, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		nbFs._array(), nbFs._index(), neis(), ineis(), _inos(), _numSprings);
	CUDA_CHECK(cudaPeekAtLastError());

	/*{
		vector<uint> h_neis;
		vector<uint> h_ineis;
		neis.copyToHost(h_neis);
		ineis.copyToHost(h_ineis);
		for (int i = 0; i < _numSprings; i++) {
			printf("%d: ", i);
			for (int j = h_ineis[i]; j < h_ineis[i + 1]; j++) {
				printf("%d ", h_neis[j]);
			}
			printf("\n");
		}
	}*/

	//printf("%d\n", numNeis);
}

void EdgeConstraint::init(
	const DPrefixArray<uint>& es, const DPrefixArray<uint>& nbNs, 
	const Dvector<REAL>& ns, REAL material) 
{
	cudaDeviceSynchronize();
	ctimer timer = CNOW;

	uint numVertices = ns.size() / 3u;
	uint numEdges = es.arraySize();

	printf("numSpring: %d\n", numEdges);

	initSprings(material, numEdges);
	initEdgeConstraintsIds_kernel << <divup(numVertices, CONSTRAINT_BLOCKSIZE), CONSTRAINT_BLOCKSIZE >> > (
		es._array(), es._index(), springs(), numVertices);
	CUDA_CHECK(cudaPeekAtLastError());

	initEdgeConstraints_kernel << <divup(numEdges, CONSTRAINT_BLOCKSIZE), CONSTRAINT_BLOCKSIZE >> > (
		ns(), springs());
	CUDA_CHECK(cudaPeekAtLastError());

	Dvector<uint> neis;
	Dvector<uint> ineis;
	getEdgeNeighbors(es, nbNs, neis, ineis, numVertices);
	compColoring(neis, ineis);

	cudaDeviceSynchronize();
	printf("init constraints %f msec\n", (CNOW - timer) / 10000.0);
}
void EdgeConstraint::project(const Dvector<REAL>& ns, const Dvector<REAL>& invMs, REAL invdt2) {
	_lambdas.memset(0);
	
	/*vector<uint> colors;
	vector<uint2> edges;
	_colors.copyToHost(colors);
	es.copyToHost(edges);
	set<uint> test;
	set<uint> test2;
	for (uint c = 0; c < _colorSize; c++) {
		test.clear();
		for (int i = 0; i < numEdges; i++) {
			if (colors[i] == c) {
				test2.insert(i);
				if (test.find(edges[i].x) != test.end()) {
					printf("asdfasdfasdf\n");
				}
				test.insert(edges[i].x);
				if (test.find(edges[i].y) != test.end()) {
					printf("asdfasdfasdf\n");
				}
				test.insert(edges[i].y);
			}
		}
	}
	for (int i = 0; i < numEdges; i++) {
		if (test2.find(i) == test2.end() || test2.size() != numEdges) {
			printf("asdfasdfasdf\n");
		}
	}*/
	for (uint i = 0; i < _colorSize; i++) {
		project_kernel << <divup(_numSprings, CONSTRAINT_BLOCKSIZE), CONSTRAINT_BLOCKSIZE >> > (
			ns(), invMs(), springs(), invdt2, i);
		CUDA_CHECK(cudaPeekAtLastError());
	}
}

void SBConstraint::init(
	const Dvector<uint>& fs, const DPrefixArray<uint>& nbFs,
	const Dvector<REAL>& ns, REAL k)
{
	cudaDeviceSynchronize();
	ctimer timer = CNOW;

	uint numVertices = ns.size() / 3u;
	uint numFaces = fs.size() / 3u;

	printf("numSpring: %d\n", numFaces);

	initSprings(numFaces);
	_inos = fs;
	initSBConstraints_kernel << <divup(numFaces, CONSTRAINT_BLOCKSIZE), CONSTRAINT_BLOCKSIZE >> > (
		ns(), springs(), k);
	CUDA_CHECK(cudaPeekAtLastError());

	Dvector<uint> neis;
	Dvector<uint> ineis;
	getFaceNeighbors(nbFs, neis, ineis, numVertices);
	compColoring(neis, ineis);

	cudaDeviceSynchronize();
	printf("init constraints %f msec\n", (CNOW - timer) / 10000.0);
}
void SBConstraint::project(const Dvector<REAL>& ns, const Dvector<REAL>& invMs) {
	uint blockSize = divup(_numSprings, CONSTRAINT_BLOCKSIZE);

	for (uint i = 0; i < _colorSize; i++) {
		project_kernel << <blockSize, CONSTRAINT_BLOCKSIZE >> > (
			ns(), invMs(), springs(), i);
		CUDA_CHECK(cudaPeekAtLastError());
	}

	/*for (uint i = 0; i < _colorSize; i++) {
		project_strain_kernel << <blockSize, CONSTRAINT_BLOCKSIZE >> > (
			ns(), invMs(), springs(), i);
		CUDA_CHECK(cudaPeekAtLastError());
	}
	for (uint i = 0; i < _colorSize; i++) {
		project_area_kernel << <blockSize, CONSTRAINT_BLOCKSIZE >> > (
			ns(), invMs(), springs(), i);
		CUDA_CHECK(cudaPeekAtLastError());
	}*/
}
void SBConstraint::project_strain(const Dvector<REAL>& ns, const Dvector<REAL>& invMs) {
	uint blockSize = divup(_numSprings, CONSTRAINT_BLOCKSIZE);
	for (uint i = 0; i < _colorSize; i++) {
		project_strain_kernel << <blockSize, CONSTRAINT_BLOCKSIZE >> > (
			ns(), invMs(), springs(), i);
		CUDA_CHECK(cudaPeekAtLastError());
	}
}
void SBConstraint::project_area(const Dvector<REAL>& ns, const Dvector<REAL>& invMs) {
	uint blockSize = divup(_numSprings, CONSTRAINT_BLOCKSIZE);
	for (uint i = 0; i < _colorSize; i++) {
		project_area_kernel << <blockSize, CONSTRAINT_BLOCKSIZE >> > (
			ns(), invMs(), springs(), i);
		CUDA_CHECK(cudaPeekAtLastError());
	}
}