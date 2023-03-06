//#ifndef __CLOTH_DYNAMICS_H__
//#define __CLOTH_DYNAMICS_H__
//
//#pragma once
//#include "CollisionSolver.h"
//
//class Obstacle
//{
//public:
//	Mesh					*_mesh;
//	BVH						*_bvh;
//	RTri					*_Rtri;
//public:
//	Dvector<uint>			d_fs;
//	Dvector<REAL>			d_ns;
//	Dvector<REAL>			d_n0s;
//	Dvector<REAL>			d_vs;
//	Dvector<REAL>			d_ms;
//	Dvector<REAL>			d_invMs;
//	DPrefixArray<uint>		d_ses;
//	DPrefixArray<uint>		d_nbFs;
//	DPrefixArray<uint>		d_nbNs;
//	Dvector<REAL>			d_fNorms;
//	Dvector<REAL>			d_nNorms;
//public:
//	vector<uint>			h_fs;
//	vector<REAL>			h_ns;
//	vector<REAL>			h_ms;
//	vector<REAL>			h_invMs;
//	PrefixArray<uint>		h_ses;
//	PrefixArray<uint>		h_nbFs;
//	PrefixArray<uint>		h_nbNs;
//	vector<REAL>			h_fNorms;
//	vector<REAL>			h_nNorms;
//public:
//	vector<uint>			h_fs0;
//	vector<REAL>			h_ns0;
//	vector<REAL>			h_ms0;
//	vector<REAL>			h_invMs0;
//	PrefixArray<uint>		h_ses0;
//	PrefixArray<uint>		h_nbFs0;
//	PrefixArray<uint>		h_nbNs0;
//public:
//	uint					_numFaces;
//	uint					_numVertices;
//public:
//	StreamParam				_streams;
//public:
//	REAL					_mass;
//	REAL					_dt;
//	REAL					_dt2;
//	REAL					_invdt;
//	REAL					_invdt2;
//	uint					_maxIter;
//	AABB					_boundary;
//public:
//	Obstacle() {}
//	Obstacle(Mesh* mesh, REAL mass) {
//		init(mesh, mass);
//	}
//	virtual ~Obstacle() {}
//public:
//	inline ObjParam param(void) {
//		ObjParam p;
//		p._fs = d_fs._list;
//		p._ns = d_ns._list;
//		p._vs = d_vs._list;
//		p._invMs = d_invMs._list;
//		p._ms = d_ms._list;
//		p._numFaces = _numFaces;
//		return p;
//	}
//public:
//	void	init(Mesh* mesh, REAL mass);
//	void	initVelocities(void);
//	void	initMasses(void);
//	void	initNoramls(void);
//	void	initElements(void);
//	void	initConstraints(void);
//public:
//	void	computeExternalForce(void);
//	void	update(void);
//	void	simulation(void);
//	void	computeNormal(void);
//public:
//	void	draw(void);
//	void	drawWire(void);
//	void	drawSurface(void);
//	void	drawBoundary(void);
//public:
//	void	reset(void);
//	void	copyToDevice(void);
//	void	copyToHost(void);
//	void	copyNbToDevice(void);
//	void	copyNbToHost(void);
//	void	copyMassToDevice(void);
//	void	copyMassToHost(void);
//	void	copyNormToDevice(void);
//	void	copyNormToHost(void);
//};
//#endif