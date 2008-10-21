#include <iostream>
#include "grid.h"
#include "inputs.h"
#include "bc.h"

extern int Rank;

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
		bcRegion.k=bcSection.doubles[b]["k"];
		bcRegion.omega=bcSection.doubles[b]["omega"];
		bcRegion.v=bcSection.Vec3Ds[b]["v"];
		bcRegion.area=0.;
		bcRegion.areaVec=0.;
		bcRegion.momentum=0.;
		if (bcRegion.kind=="") {
			cout << "[I Rank=" << Rank << "] BC_" << b+1 << " assigned as " << bcRegion.type << endl;
		} else {
			cout << "[I Rank=" << Rank << "] BC_" << b+1 << " assigned as " << bcRegion.type << " of " << bcRegion.kind << " kind" << endl;
		}	
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

	// Integrate boundary areas
	for (unsigned int f=0;f<grid.faceCount;++f) {
		int bcIndex=grid.face[f].bc;
		if (bcIndex>=0) {
			bc.region[bcIndex].area+=grid.face[f].area;
			bc.region[bcIndex].areaVec+=grid.face[f].area*grid.face[f].normal;
		}
	}
	
	return;
}
