
#include "grid.h"
#include "inputs.h"

bool within_box(Vec3D centroid, Vec3D box_1, Vec3D box_2);

void initialize(Grid &grid, InputFile &input) {
	
	// Loop through each initial condition region and apply sequentially
	numberedSubsection region=input.section["initialConditions"].numberedSubsections["region"];
	for (unsigned int r=0;r<region.count;++r) {
		if (region.strings[r]["type"]=="box") {
			for (unsigned int c = 0;c < grid.cellCount;++c) {
				Vec3D box_1=region.Vec3Ds[r]["box_1"];
				Vec3D box_2=region.Vec3Ds[r]["box_2"];
				if (within_box(grid.cell[c].centroid,box_1,box_2)) {
					// The cell centroid is inside the box region
					grid.cell[c].rho = region.doubles[r]["rho"];
					grid.cell[c].v = region.Vec3Ds[r]["v"];
					grid.cell[c].p = region.doubles[r]["p"];
				}
			}
		} else if (region.strings[r]["type"]=="circle") {
			double radius=region.doubles[r]["radius"];
			Vec3D center=region.Vec3Ds[r]["center"];
			Vec3D zComp=center;
			Vec3D center2cell;
			zComp.comp[0]=0.;zComp.comp[1]=0.;
			center=center-zComp;
			for (unsigned int c = 0;c < grid.cellCount;++c) {
			zComp=grid.cell[c].centroid;  zComp.comp[0]=0.;zComp.comp[1]=0.;
			center2cell=(grid.cell[c].centroid-zComp)-center;
				if ( fabs(center2cell)<=radius ) {
					// The cell centroid is inside the box region
					grid.cell[c].rho = region.doubles[r]["rho"];
					grid.cell[c].v = region.Vec3Ds[r]["v"];
					grid.cell[c].v.comp[0]=(region.Vec3Ds[r]["v"].comp[0]*center2cell).comp[0];
					grid.cell[c].v.comp[1]=(region.Vec3Ds[r]["v"].comp[0]*center2cell).comp[1];					
					grid.cell[c].p = region.doubles[r]["p"];
				}
			}
		}
	}

	for (unsigned int c=0;c<grid.cellCount;++c) {
		for (unsigned int i=0;i<5;++i) {
			grid.cell[c].flux[i]=0.;
		}
	}
	
return;
}
