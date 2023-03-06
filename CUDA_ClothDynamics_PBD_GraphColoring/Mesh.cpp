#include "Mesh.h"

void Mesh::loadObj(const char* filename, REAL3 center, REAL scale) {
	CUDA_CHECK(cudaDeviceSynchronize());
	ctimer timer = CNOW;

	_fs.clear();
	_ns.clear();

	bool flag = true;
	ifstream fin;
	fin.open(filename);
	if (fin.is_open()) {
		while (!fin.eof()) {
			string head;
			fin >> head;
			if (head.length() > 1)
				continue;
			if (head[0] == 'v') {
				REAL3 x;
				fin >> x.x >> x.y >> x.z;
				_ns.push_back(x.x);
				_ns.push_back(x.y);
				_ns.push_back(x.z);
				if (flag) {
					_aabb._min = _aabb._max = x;
					flag = false;
				}
				else addAABB(_aabb, x);
			}
			else if(head[0] == 'f') {
				uint3 x;
				fin >> x.x >> x.y >> x.z;
				_fs.push_back(x.x - 1u);
				_fs.push_back(x.y - 1u);
				_fs.push_back(x.z - 1u);
			}
		}
		fin.close();
	}
	if (_ns.empty() || _fs.empty()) {
		printf("Error : Mesh_init : Object Load Error\n");
		exit(1);
		return;
	}
	moveCenter(center, scale);
	buildAdjacency();

	_numFaces = _fs.size() / 3u;
	_numVertices = _ns.size() / 3u;
	CUDA_CHECK(cudaDeviceSynchronize());
	printf("Num of Faces: %d, Num of Vertices: %d, %f ms\n", _numFaces, _numVertices, (CNOW - timer) / 10000.);
}

void Mesh::moveCenter(REAL3 center, REAL scale) {
	REAL3 size = _aabb._max - _aabb._min;
	REAL max_length = size.x;
	if (max_length < size.y)
		max_length = size.y;
	if (max_length < size.z)
		max_length = size.z;
	max_length = 2.0 * scale / max_length;

	REAL3 prevCenter = (_aabb._min + _aabb._max) * (REAL)0.5;

	bool flag = true;
	uint vlen = _ns.size();
	for (uint i = 0u; i < vlen; i += 3u) {
		REAL3 pos = make_REAL3(_ns[i], _ns[i + 1u], _ns[i + 2u]);
		REAL3 grad = pos - prevCenter;
		grad *= max_length;
		pos = center + grad;
		_ns[i] = pos.x;
		_ns[i + 1u] = pos.y;
		_ns[i + 2u] = pos.z;
		if (flag) {
			_aabb._min = _aabb._max = pos;
			flag = false;
		}
		else addAABB(_aabb, pos);
	}
}
void Mesh::buildAdjacency(void)
{
	uint numFaces = _fs.size() / 3u;
	uint numVertices = _ns.size() / 3u;
	vector<set<uint>> nbFs(numVertices);
	vector<set<uint>> nbNs(numVertices);

	for (uint i = 0u; i < numFaces; i++) {
		uint ino = i * 3u;
		for (uint j = 0u; j < 3u; j++)
			nbFs[_fs[ino + j]].insert(i);
	}

	for (uint i = 0u; i < numVertices; i++) {
		for (auto inbf : nbFs[i]) {
			uint ino = inbf * 3u;
			uint ino0, ino1;
			uint pivot;
			if (_fs[ino + 0u] == i) {
				pivot = 0u;
				ino0 = _fs[ino + 1u];
				ino1 = _fs[ino + 2u];
			}
			else if (_fs[ino + 1u] == i) {
				pivot = 1u;
				ino0 = _fs[ino + 2u];
				ino1 = _fs[ino + 0u];
			}
			else {
				pivot = 2u;
				ino0 = _fs[ino + 0u];
				ino1 = _fs[ino + 1u];
			}
			nbNs[i].insert(ino0);
			nbNs[i].insert(ino1);
		}
	}
	vector<set<uint>> ses(numVertices);
	vector<set<uint>> bes(numVertices);
	for (uint i = 0u; i < numVertices; i++) {
		for (auto inbv : nbNs[i])
			if (i < inbv)
				ses[i].insert(inbv);
	}

	/*uint numEdges = _ses.size();
	for (uint i = 0u; i < numEdges; i++) {
		uint ns[2];
		uint ino0 = _ses[i].x;
		uint ino1 = _ses[i].y;
		int num = 0;
		for (auto inbf0 : nbFs[ino0]) {
			for (auto inbf1 : nbFs[ino1]) {
				if (inbf0 == inbf1) {
					uint iface = inbf0 * 3u;
					for (uint j = 0u; j < 3u; j++) {
						uint ino = _fs[iface + j];
						if (ino != ino0 && ino != ino1)
							ns[num++] = ino;
					}
					break;
				}
			}
			if (num == 2) break;
		}
		if (num == 2)
			_bes.push_back(make_uint2(ns[0], ns[1]));
		else if (num == 1)
			_bes.push_back(make_uint2(ns[0], 0xffffffff));
		else {
			printf("Error : Mesh_buildAdjacency : Edge Error\n");
			exit(1);
		}
	}*/
	uint ns[2];
	for (uint ino0 = 0u; ino0 < numVertices; ino0++) {
		for (auto ino1 : ses[ino0]) {
			uint num = 0u;
			for (auto inbf0 : nbFs[ino0]) {
				for (auto inbf1 : nbFs[ino1]) {
					if (inbf0 == inbf1) {
						uint iface = inbf0 * 3u;
						for (uint j = 0u; j < 3u; j++) {
							uint ino = _fs[iface + j];
							if (ino != ino0 && ino != ino1)
								ns[num++] = ino;
						}
						break;
					}
				}
				if (num == 2u) break;
			}
			if (num == 2u) {
				num = ns[0] > ns[1];
				bes[ns[num]].insert(ns[num ^ 1u]);
			}
		}
	}
	_nbFs.init(nbFs);
	_nbNs.init(nbNs);
	_ses.init(ses);
	_bes.init(bes);
}