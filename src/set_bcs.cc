#include <iostream>
#include "grid.h"
#include "inputs.h"
#include "bc.h"

extern int rank;

bool within_box(Vec3D centroid, Vec3D box_1, Vec3D box_2);

void set_bcs(Grid& grid, InputFile &input, BC &bc) {
	// Loop through each boundary condition region and apply sequentially
	numberedSubsection bcSection=input.section["boundaryConditions"].numberedSubsections["BC"];
	for (unsigned int b=0;b<bcSection.count;++b) {
		BCregion bcRegion;
		bcRegion.type=bcSection.strings[b]["type"];
		bcRegion.kind=bcSection.strings[b]["kind"];
		bcRegion.rho=bcSection.doubles[b]["rho"];
		bcRegion.p=bcSection.doubles[b]["p"];
		bcRegion.v=bcSection.Vec3Ds[b]["v"];
		bc.region.push_back(bcRegion);
		if (bcSection.strings[b]["region"]=="box") {
			Vec3D box_1=bcSection.Vec3Ds[b]["box_1"];
			Vec3D box_2=bcSection.Vec3Ds[b]["box_2"];
			for (unsigned int f = 0;f < grid.faceCount;++f) {
				// if the face is not already marked as internal or partition boundary
				// And if the face centroid falls within the defined box
				if ((grid.face[f].bc>=0 || grid.face[f].bc==-2 ) && within_box(grid.face[f].centroid,box_1,box_2)) {				
					if (bcSection.strings[b]["pick"]=="override") {
						grid.face[f].bc=b; // real boundary conditions are marked as positive
					} else if (bcSection.strings[b]["pick"]=="unassigned" && grid.face[f].bc==-2) {
						grid.face[f].bc=b; // real boundary conditions are marked as positive
					}
				}
			}
		}
	}

	// Mark nodes that touch boundaries
	for (unsigned int f=0;f<grid.faceCount;++f) {
		if (grid.face[f].bc>=0) { // if a boundary face
			for (unsigned int fn=0;fn<grid.face[f].nodeCount;++fn) {
				grid.face[f].node(fn).bcs.insert(grid.face[f].bc);
			}
		}
	}

	return;
}
