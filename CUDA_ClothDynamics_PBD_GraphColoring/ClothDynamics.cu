#include "ClothDynamics.cuh"

void ClothDynamics::init(Mesh* mesh, REAL mass) {
	_maxIter = 30u;
	_dt = (REAL)0.01;
	{
		_dt2 = _dt * _dt;
		_invdt = 1.0 / _dt;
		_invdt2 = _invdt * _invdt;
	}

	_mesh = mesh;
	_mass = mass;

	_boundary._min = make_REAL3(-1.5);
	_boundary._max = make_REAL3(1.5);

	_streams.initStream(10);
	initElements();
	initConstraints();
	_RTri = new RTriangle();
	_bvh = new BVH();
	_RTri->init(d_fs, d_nbFs);
	_bvh->build(d_fs, d_ns, d_nbFs);
}

void ClothDynamics::initVelocities(void) {
	d_vs.resize(_numNodes * 3u);
	d_vs.memset((REAL)0.0, _streams[2]);
}
void ClothDynamics::initMasses(void) {
	h_ms.clear();
	h_invMs.clear();
	h_ms.resize(_numNodes, _mass);
	h_invMs.resize(_numNodes, 1.0 / _mass);

	double offset = 0.000001;
	for (int i = 0; i < _numNodes; i++) {
		// hanging cloth1
		REAL3 v = make_REAL3(h_ns[i * 3 + 0], h_ns[i * 3 + 1], h_ns[i * 3 + 2]);
		if (v.x <= _mesh->_aabb._min.x + offset && v.y >= _mesh->_aabb._max.y - offset) {
			h_ms[i] = h_invMs[i] = 0.;
		}
		/*if (v.y >= _mesh->_aabb._max.y - offset) {
			h_ms[i] = h_invMs[i] = 0.;
		}*/
	}

	copyMassToDevice();
}
void ClothDynamics::initNoramls(void) {
	h_fNorms.resize(_numFaces * 3u);
	h_nNorms.resize(_numNodes * 3u);
	d_fNorms.resize(_numFaces * 3u);
	d_nNorms.resize(_numNodes * 3u);
	computeNormal();
}
void ClothDynamics::initElements(void) {
	_numFaces = _mesh->_numFaces;
	_numNodes = _mesh->_numVertices;

	initVelocities();

	h_fs = _mesh->_fs;
	h_ns = _mesh->_ns;
	h_ses = _mesh->_ses;
	h_bes = _mesh->_bes;
	h_nbFs = _mesh->_nbFs;
	h_nbNs = _mesh->_nbNs;

	initMasses();

	copyToDevice();
	copyNbToDevice();
	initNoramls();

	h_fs0 = h_fs;
	h_ns0 = h_ns;
	h_ses0 = h_ses;
	h_bes0 = h_bes;
	h_nbFs0 = h_nbFs;
	h_nbNs0 = h_nbFs;
	h_ms0 = h_ms;
	h_invMs0 = h_invMs;
}
void ClothDynamics::initConstraints(void) {
	_constraints.init(d_ses, d_nbNs, d_ns, 0.00000000000016);
	_fconstraints.init(d_fs, d_nbFs, d_ns, 0.9);
}
void ClothDynamics::computeExternalForce(void)
{
	d_forces.resize(_numNodes * 3u);
	d_forces.memset(0);

	REAL3 gravity = make_REAL3(0.0, -2.8, 0.0);
	compExternalForce_kernel << <divup(_numNodes, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		d_forces(), d_ms(), gravity, _numNodes);
	CUDA_CHECK(cudaPeekAtLastError());

	d_n0s = d_ns;
}

void ClothDynamics::update(void) {

	//d_vs += d_forces * d_invMs * _dt;
	//d_ns += d_vs * _dt;
	compPredictPosition_kernel << <divup(_numNodes, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		d_ns(), d_vs(), d_forces(), d_invMs(), _dt, _numNodes);
	CUDA_CHECK(cudaPeekAtLastError());

	for (int k = 0; k < _maxIter; k++) {
		_constraints.project(d_ns, d_invMs, _invdt);
		//_fconstraints.project(d_ns, d_invMs);
		//_fconstraints.project_strain(d_ns, d_invMs);
		_fconstraints.project_area(d_ns, d_invMs);
	}

	//d_vs = (d_ns - d_n0s) * _invdt;
	updateVelocitiy_kernel << <divup(_numNodes, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		d_ns(), d_n0s(), d_vs(), _invdt, _numNodes);
	CUDA_CHECK(cudaPeekAtLastError());
}
void ClothDynamics::collision(void) {
	bool isApplied = false;
	REAL contactClearance = 0.2;
	REAL thickness = 0.01;
	REAL boundaryFriction = 0.2;
	REAL obstacleFriction = 0.1;

	d_ns = d_n0s;
	//_bvh->build(d_fs, d_ns, d_nbFs);
	isApplied |= CollisionSolver::ResolveSelfCollision(param(), d_ses, _bvh, _RTri, contactClearance, _dt);

	updatePosition_kernel << <divup(_numNodes, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		d_ns(), d_vs(), _dt, _numNodes);
	CUDA_CHECK(cudaPeekAtLastError());
}
void ClothDynamics::simulation(void)
{
	cudaDeviceSynchronize();
	ctimer timer = CNOW;

	computeNormal();
	computeExternalForce();
	update();
	//collision();

	copyToHost();

	cudaDeviceSynchronize();
	printf("Simulation %f msec\n", (CNOW - timer) / 10000.0);
}

void ClothDynamics::computeNormal(void) {
	compFnorms_kernel << <divup(_numFaces, MAX_BLOCKSIZE), MAX_BLOCKSIZE >> > (
		d_fs(), d_ns(), d_fNorms(), _numFaces);
	CUDA_CHECK(cudaPeekAtLastError());
	copyNormToHost();
}
void ClothDynamics::draw(void)
{
	drawWire();
	drawSurface();
	//drawBoundary();
	//_bvh->draw();
}
void ClothDynamics::drawBoundary(void)
{
	glPushMatrix();
	glDisable(GL_LIGHTING);
	glColor3d(1, 1, 1);
	glLineWidth(3.0f);

	glBegin(GL_LINES);
	glVertex3f(_boundary._min.x, _boundary._min.y, _boundary._min.z);
	glVertex3f(_boundary._max.x, _boundary._min.y, _boundary._min.z);
	glVertex3f(_boundary._min.x, _boundary._min.y, _boundary._min.z);
	glVertex3f(_boundary._min.x, _boundary._max.y, _boundary._min.z);
	glVertex3f(_boundary._min.x, _boundary._min.y, _boundary._min.z);
	glVertex3f(_boundary._min.x, _boundary._min.y, _boundary._max.z);
	glVertex3f(_boundary._max.x, _boundary._max.y, _boundary._max.z);
	glVertex3f(_boundary._min.x, _boundary._max.y, _boundary._max.z);
	glVertex3f(_boundary._max.x, _boundary._max.y, _boundary._max.z);
	glVertex3f(_boundary._max.x, _boundary._min.y, _boundary._max.z);
	glVertex3f(_boundary._max.x, _boundary._max.y, _boundary._max.z);
	glVertex3f(_boundary._max.x, _boundary._max.y, _boundary._min.z);
	glVertex3f(_boundary._max.x, _boundary._min.y, _boundary._min.z);
	glVertex3f(_boundary._max.x, _boundary._min.y, _boundary._max.z);
	glVertex3f(_boundary._max.x, _boundary._min.y, _boundary._min.z);
	glVertex3f(_boundary._max.x, _boundary._max.y, _boundary._min.z);
	glVertex3f(_boundary._min.x, _boundary._max.y, _boundary._min.z);
	glVertex3f(_boundary._min.x, _boundary._max.y, _boundary._max.z);
	glVertex3f(_boundary._min.x, _boundary._max.y, _boundary._min.z);
	glVertex3f(_boundary._max.x, _boundary._max.y, _boundary._min.z);
	glVertex3f(_boundary._min.x, _boundary._min.y, _boundary._max.z);
	glVertex3f(_boundary._min.x, _boundary._max.y, _boundary._max.z);
	glVertex3f(_boundary._min.x, _boundary._min.y, _boundary._max.z);
	glVertex3f(_boundary._max.x, _boundary._min.y, _boundary._max.z);
	glEnd();

	glLineWidth(1.0f);
	glEnable(GL_LIGHTING);
	glShadeModel(GL_FLAT);
	glPopMatrix();
}
void ClothDynamics::drawWire(void)
{
	glPushMatrix();
	glDisable(GL_LIGHTING);
	glColor3d(0, 0, 0);

	for (int i = 0; i < _numFaces; i++) {
		glBegin(GL_LINE_LOOP);
		for (int j = 0; j < 3; j++) {
			auto x = h_ns[h_fs[i * 3 + j] * 3 + 0];
			auto y = h_ns[h_fs[i * 3 + j] * 3 + 1];
			auto z = h_ns[h_fs[i * 3 + j] * 3 + 2];
			glVertex3f(x, y, z);
		}
		glEnd();
	}

	glEnable(GL_LIGHTING);
	glShadeModel(GL_FLAT);
	glPopMatrix();
}
void ClothDynamics::drawSurface(void)
{
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1); // turn on two-sided lighting.
	float purple[] = { 0.0f, 0.44705882352941176470588235294118f, 0.66666666666666666666666666666667f, 1.0f };
	float yellow[] = { 0.6f, 0.6f, 0.0f, 1.0f };
	float white[] = { 0.4f, 0.4f, 0.4f, 1.0f };
	float black[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, purple);
	glMaterialfv(GL_FRONT, GL_SPECULAR, white);
	glMaterialf(GL_FRONT, GL_SHININESS, 64);
	glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, yellow); // back material
	glMaterialfv(GL_BACK, GL_SPECULAR, black); // no specular highlights

	for (uint i = 0u; i < _numFaces; i++) {
		uint ino0 = h_fs[i * 3u + 0u];
		uint ino1 = h_fs[i * 3u + 1u];
		uint ino2 = h_fs[i * 3u + 2u];
		REAL3 a = make_REAL3(h_ns[ino0 * 3u + 0u], h_ns[ino0 * 3u + 1u], h_ns[ino0 * 3u + 2u]);
		REAL3 b = make_REAL3(h_ns[ino1 * 3u + 0u], h_ns[ino1 * 3u + 1u], h_ns[ino1 * 3u + 2u]);
		REAL3 c = make_REAL3(h_ns[ino2 * 3u + 0u], h_ns[ino2 * 3u + 1u], h_ns[ino2 * 3u + 2u]);
		glBegin(GL_TRIANGLES);
		glNormal3f(h_fNorms[i * 3u + 0u], h_fNorms[i * 3u + 1u], h_fNorms[i * 3u + 2u]);
		glVertex3f(a.x, a.y, a.z);
		glVertex3f(b.x, b.y, b.z);
		glVertex3f(c.x, c.y, c.z);
		glEnd();
	}
}

void ClothDynamics::reset(void) {
	_numFaces = _mesh->_numFaces;
	_numNodes = _mesh->_numVertices;

	initVelocities();

	h_fs = h_fs0;
	h_ns = h_ns0;
	h_ses = h_ses0;
	h_bes = h_bes0;
	h_nbFs = h_nbFs0;
	h_nbNs = h_nbFs0;

	copyToDevice();
	copyNbToDevice();
	copyMassToDevice();
	initNoramls();

	h_ms = h_ms0;
	h_invMs = h_invMs0;
}
void ClothDynamics::copyToDevice(void) {
	d_fs.copyFromHost(h_fs, _streams[0]);
	d_ns.copyFromHost(h_ns, _streams[1]);
}
void ClothDynamics::copyToHost(void) {
	d_fs.copyToHost(h_fs, _streams[0]);
	d_ns.copyToHost(h_ns, _streams[1]);
}
void ClothDynamics::copyNbToDevice(void) {
	d_ses.copyFromHost(h_ses, &_streams[2]);
	d_nbFs.copyFromHost(h_nbFs, &_streams[4]);
	d_nbNs.copyFromHost(h_nbNs, &_streams[6]);
}
void ClothDynamics::copyNbToHost(void) {
	d_ses.copyToHost(h_ses, &_streams[2]);
	d_nbFs.copyToHost(h_nbFs, &_streams[4]);
	d_nbNs.copyToHost(h_nbNs, &_streams[6]);
}
void ClothDynamics::copyMassToDevice(void) {
	d_ms.copyFromHost(h_ms, _streams[8]);
	d_invMs.copyFromHost(h_invMs, _streams[9]);
}
void ClothDynamics::copyMassToHost(void) {
	d_ms.copyToHost(h_ms, _streams[8]);
	d_invMs.copyToHost(h_invMs, _streams[9]);
}
void ClothDynamics::copyNormToDevice(void) {
	d_fNorms.copyFromHost(h_fNorms, _streams[0]);
	d_nNorms.copyFromHost(h_nNorms, _streams[1]);
}
void ClothDynamics::copyNormToHost(void) {
	d_fNorms.copyToHost(h_fNorms, _streams[0]);
	d_nNorms.copyToHost(h_nNorms, _streams[1]);
}