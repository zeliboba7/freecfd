#include <iostream>
#include "grid.h"
#include "inputs.h"
#include "bc.h"

bool within_box(Vec3D centroid, Vec3D box_1, Vec3D box_2);
	
void set_bcs(Grid& grid, InputFile &input, BC &bc) {
	// Loop through each boundary condition region and apply sequentially
	numberedSubsection bcSection=input.section["boundaryConditions"].numberedSubsections["BC"];
	for (unsigned int b=0;b<bcSection.count;++b) {
		BCregion bcRegion;
		bcRegion.type=bcSection.strings[b]["type"];
		bcRegion.rho=bcSection.doubles[b]["rho"];
		bcRegion.p=bcSection.doubles[b]["p"];
		bcRegion.v=bcSection.Vec3Ds[b]["v"];		
		bc.region.push_back(bcRegion);	
		if (bcSection.strings[b]["region"]=="box") {
			for (unsigned int f = 0;f < grid.faceCount;++f) {
				Vec3D box_1=bcSection.Vec3Ds[b]["box_1"];
				Vec3D box_2=bcSection.Vec3Ds[b]["box_2"];
				if (within_box(grid.face[f].centroid,box_1,box_2)) {
					// The face centroid is inside the box region
					grid.face[f].bc=b;
				}
				if (grid.face[f].neighbor>=0) grid.face[f].bc=-1; // internal face
			}
		}
	}
	
}
