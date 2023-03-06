#ifndef __CONSTRAINT_H__
#define __CONSTRAINT_H__

#pragma once
#include "Mesh.h"

#define CONSTRAINT_BLOCKSIZE		128u

struct BaseSpring {
	uint	*_inos;
	uint	*_colors;
	uint	_numSprings;
};
struct EdgeSpring : public BaseSpring {
	REAL	*_rs;
	REAL	*_lambdas;
	REAL	_material;
};
struct SBSpring : public BaseSpring {
	REAL	*_invQs;
	REAL	*_as;
	REAL	*_ks;
};

class Constraint
{
public:
	Dvector<uint>	_colors;
	uint			_colorSize;
public:
	Dvector<uint>	_inos;
	uint			_numSprings;
public:
	Constraint() {}
	virtual ~Constraint() {}
public:
	inline size_t numSprings(void) {
		return _numSprings;
	}
public:
	void compColoring(const Dvector<uint>& neis, const Dvector<uint>& ineis);
public:
	void getEdgeNeighbors(
		const DPrefixArray<uint>& es, const DPrefixArray<uint>& nbNs,
		Dvector<uint>& neis, Dvector<uint>& ineis, uint numVertices);
	void getFaceNeighbors(
		const DPrefixArray<uint>& nbFs,
		Dvector<uint>& neis, Dvector<uint>& ineis, uint numEdges);
};

class EdgeConstraint : public Constraint
{
public:
	Dvector<REAL>	_lambdas;
	Dvector<REAL>	_rs;
	REAL			_material;
public:
	EdgeConstraint() {}
	virtual ~EdgeConstraint() {}
public:
	inline void initSprings(REAL material, uint numSprings) {
		_inos.resize(numSprings << 1u);
		_rs.resize(numSprings);
		_lambdas.resize(numSprings);
		_material = material;
		_numSprings = numSprings;
	}
	inline EdgeSpring springs(void) {
		EdgeSpring s;
		s._inos = _inos._list;
		s._rs = _rs._list;
		s._lambdas = _lambdas._list;
		s._colors = _colors._list;
		s._material = _material;
		s._numSprings = _numSprings;
		return s;
	}
public:
	void init(
		const DPrefixArray<uint>& es, const DPrefixArray<uint>& nbNs, 
		const Dvector<REAL>& ns, REAL material);
	void project(const Dvector<REAL>& ns, const Dvector<REAL>& invMs, REAL invdt2);
};
class SBConstraint : public Constraint
{
public:
	Dvector<REAL> _invQs;
	Dvector<REAL> _as;
	Dvector<REAL> _ks;
public:
	SBConstraint() {}
	virtual ~SBConstraint() {}
public:
	inline void initSprings(uint numSprings) {
		_inos.resize(numSprings * 3u);
		_invQs.resize(numSprings << 2u);
		_as.resize(numSprings);
		_ks.resize(numSprings);
		_numSprings = numSprings;
	}
	inline SBSpring springs(void) {
		SBSpring s;
		s._inos = _inos._list;
		s._invQs = _invQs._list;
		s._as = _as._list;
		s._ks = _ks._list;
		s._colors = _colors._list;
		s._numSprings = _numSprings;
		return s;
	}
public:
	void init(
		const Dvector<uint>& fs, const DPrefixArray<uint>& nbFs,
		const Dvector<REAL>& ns, REAL k);
	void project(const Dvector<REAL>& ns, const Dvector<REAL>& invMs);
	void project_strain(const Dvector<REAL>& ns, const Dvector<REAL>& invMs);
	void project_area(const Dvector<REAL>& ns, const Dvector<REAL>& invMs);
};

#endif