#include "RTriangle.h"
#include "DeviceManager.cuh"

__global__ void RTriBuild_Kernel(
	uint* fs, uint* nbFs, uint* inbFs, RTriParam rtri)
{
	uint id = blockDim.x * blockIdx.x + threadIdx.x;
	if (id >= rtri._size)
		return;

	uint inos[3];
	uint ino = id * 3u, jno;
	inos[0] = fs[ino + 0u];
	inos[1] = fs[ino + 1u];
	inos[2] = fs[ino + 2u];

	uint istart, iend, jstart, jend;
	uint i, j, itri, jtri;
	uint info = 0u;

	jstart = inbFs[inos[0]];
	jend = inbFs[inos[0] + 1u];
	for (i = 0u; i < 3u; i++) {
		ino = jstart;
		iend = jend;
		itri = nbFs[ino];
		if (itri == id)
			setRTriVertex(info, i);

		j = (i + 1u) % 3u;
		jstart = inbFs[inos[j]];
		jend = inbFs[inos[j] + 1u];
		jno = jstart;
		jtri = nbFs[jno];

		while (ino < iend && jno < jend) {
			if (itri < jtri)
				itri = nbFs[++ino];
			else if (itri > jtri)
				jtri = nbFs[++jno];
			else {
				if (itri == id) {
					itri = nbFs[++ino];
					jtri = nbFs[++jno];
				}
				else {
					if (id > itri)
						setRTriEdge(info, i);
					break;
				}
			}
		}
		/*flag = true;
		for (; ino < iend;ino++) {
			itri = nbFs[ino];
			for (jno = jstart; jno < jend; jno++) {
				jtri = nbFs[jno];
				if (itri == jtri)
					break;
			}
			if (jno < jend && itri != id) {
				setRTriEdge(info, i);
				break;
			}
		}*/
	}
	rtri._info[id] = info;
	//rtri._info[id] = 0xffffffff;
}