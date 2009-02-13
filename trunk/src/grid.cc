/************************************************************************
	
	Copyright 2007-2008 Emre Sozer & Patrick Clark Trizila

	Contact: emresozer@freecfd.com , ptrizila@freecfd.com

	This file is a part of Free CFD

	Free CFD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    Free CFD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    For a copy of the GNU General Public License,
    see <http://www.gnu.org/licenses/>.

*************************************************************************/
#include "commons.h"
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <iomanip>
#include<cmath>
using namespace std;
#include <cgnslib.h>
#include <parmetis.h>

#include "inputs.h"
#include "bc.h"

extern InputFile input;
extern BC bc;

string int2str(int number) ;
double minmod(double a, double b);
double doubleMinmod(double a, double b);
double harmonic(double a, double b);
double superbee(double a, double b);


extern GridRawData raw;

Grid::Grid() {
	;
}

int Grid::read(string fname) {
	fstream file;
	fileName=fname;
	file.open(fileName.c_str());
	if (file.is_open()) {
		if (Rank==0) cout << "[I] Found grid file " << fileName  << endl;
		file.close();
		readCGNS();
		//readTEC();
		//reorderRCM();
		scale();
		rotate();
		partition();
		create_nodes_cells();
		mesh2dual();
		create_faces();
		create_ghosts();
		areas_volumes();
		return 1;
	} else {
		if (Rank==0) cerr << "[E] Grid file "<< fileName << " could not be found." << endl;
		return 0;
	}
}

// int Grid::reorderRCM() {
// 	//void genrcmi (const int n, const int flags, const int *xadj, const int * adj, int * perm, signed char * mask, int * deg)
// 	
// 	vector <int> adjIndex;
// 	adjIndex.assign(raw.cellConnIndex.begin(),raw.cellConnIndex.end());
// 	adjIndex.push_back(raw.cellConnectivity.size());
// 	vector<int> RCMpermutation;
// 	vector<signed char> mask;
// 	vector<int> deg;
// 	RCMpermutation.resize(raw.cellConnIndex.size());
// 	
// 	genrcmi(raw.cellConnIndex.size(),0,&adjIndex[0],&raw.cellConnectivity[0],&RCMpermutation[0],&mask[0],&deg[0]);
// 	return 1;
// }

int Grid::scale() {
	if (input.section("grid").get_Vec3D("scaleBy").is_found) {
		Vec3D scaleFactor=input.section("grid").get_Vec3D("scaleBy");
		for (unsigned int n=0;n<globalNodeCount;++n) {
			raw.x[n]*=scaleFactor[0];
			raw.y[n]*=scaleFactor[1];
			raw.z[n]*=scaleFactor[2];
		}
		if (Rank==0) cout << "[I] Grid is scaled by " << scaleFactor << endl;
	}
	return 1;
}

int Grid::rotate() {
	if (input.section("grid").get_Vec3D("rotationCenter").is_found) {
		Vec3D center=input.section("grid").get_Vec3D("rotationCenter");
		Vec3D angles=input.section("grid").get_Vec3D("rotationAngles");
		// Convert angles to radian
		angles*=4.*atan(1.)/180.;
		double x,y,z;
		for (unsigned int n=0;n<globalNodeCount;++n) {
			// Translate to the new center
			raw.x[n]-=center[0];
			raw.y[n]-=center[1];
			raw.z[n]-=center[2];			
			x=raw.x[n]; y=raw.y[n]; z=raw.z[n];
			// Rotate around x axis
			raw.y[n]=y*cos(angles[0])-z*sin(angles[0]);
			raw.z[n]=y*sin(angles[0])+z*cos(angles[0]);
			// Rotate around y axis
			raw.x[n]=x*cos(angles[1])+z*sin(angles[1]);
			raw.z[n]=-x*sin(angles[1])+z*cos(angles[1]);
			// Rotate around z axis
			raw.x[n]=x*cos(angles[2])-y*sin(angles[2]);
			raw.y[n]=x*sin(angles[2])+y*cos(angles[2]);
			// Translate back to the original center
			raw.x[n]+=center[0];
			raw.y[n]+=center[1];
			raw.z[n]+=center[2];
		}
		if (Rank==0) cout << "[I] Grid rotation applied" << endl;
	}
	return 1;
}

int Grid::areas_volumes() {
	// NOTE The methodology here is generic hence slower than it could be:
	// i.e Even though volume/centroid of tetrahedra is easy enough, we still break it down to
	// smaller tetrahedras

	// Now loop through faces and calculate centroids and areas
	for (unsigned int f=0;f<faceCount;++f) {
		Vec3D centroid=0.;
		Vec3D areaVec=0.;
		Vec3D patchCentroid,patchArea;
		// Find an approxiamate centroid (true centroid for triangle)
		for (unsigned int n=0;n<face[f].nodeCount;++n) {
			centroid+=face[f].node(n);
		}
		centroid/=double(face[f].nodeCount);
		
		// Sum the area as a patch of triangles formed by connecting two nodes and an interior point
		face[f].centroid=0.;
		areaVec=0.;
		int next;
		// First calculate face normal
		face[f].normal=(face[f].node(2)-face[f].node(1)).cross(face[f].node(0)-face[f].node(1)).norm();
		for (unsigned int n=0;n<face[f].nodeCount;++n) {
			next=n+1;
			if (next==face[f].nodeCount) next=0;
			patchArea=0.5*(face[f].node(n)-centroid).cross(face[f].node(next)-centroid);
			patchCentroid=1./3.*(face[f].node(n)+face[f].node(next)+centroid);
			face[f].centroid+=patchCentroid*patchArea.dot(face[f].normal);
			areaVec+=patchArea;
		}
		face[f].area=fabs(areaVec);
		face[f].centroid/=face[f].area;

	}
			
	// Loop through the cells and calculate the volumes
	double totalVolume=0.;
	for (unsigned int c=0;c<cellCount;++c) {
		double volume=0.;
		double patchVolume=0.;
		Vec3D height,basePatchArea,patchCentroid;
		unsigned int f,next;
		// Calculate cell centroid
		// First calculate an approximate one
		Vec3D centroid=0.;
		for (unsigned int cn=0;cn<cell[c].nodeCount;++cn) {
			centroid+=cell[c].node(cn);
		}
		centroid/=double(cell[c].nodeCount);
		
		// Break-up the volume into tetrahedras and add the volumes
		// Calculate the centroid of the cell by taking moments of each tetra
		cell[c].centroid=0.;
		double sign;
		for (unsigned int cf=0;cf<cell[c].faceCount;++cf) {
			f=cell[c].faces[cf];
			// Every cell face is broken to triangles
			for (unsigned int n=0;n<face[f].nodeCount;++n) {
				next=n+1;
				if (next==face[f].nodeCount) next=0;
				// Triangle area
				basePatchArea=0.5*(face[f].node(n)-face[f].centroid).cross(face[f].node(next)-face[f].centroid);
				// Height of the tetrahedra
				height=(face[f].centroid-centroid).dot(face[f].normal)*face[f].normal;
				// Fix face orientation issue
				sign=-1.;
				if (face[f].parent==c) sign=1.;
				patchVolume=sign*basePatchArea.dot(height)/3.;
				patchCentroid=0.25*(face[f].centroid+face[f].node(n)+face[f].node(next)+centroid);
				cell[c].centroid+=patchVolume*patchCentroid;
				// TODO Keep an eye on this
				//if (patchVolume<0.) cout << "[E Rank=" << Rank << "] Encountered error when calculating volume of cell " << c << endl;
				volume+=patchVolume;
			}
		}
		cell[c].volume=volume;
		cell[c].centroid/=volume;
		totalVolume+=volume;
	}
	
	cout << "[I Rank=" << Rank << "] Total Volume= " << setw(16) << setprecision(8) << scientific << totalVolume << endl;
	double globalTotalVolume=0.;
	MPI_Allreduce (&totalVolume,&globalTotalVolume,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	if (np>1 && Rank==0) cout << "[I] Global Total Volume= " << globalTotalVolume << endl;
	
	for (unsigned int f=0;f<faceCount;++f) {
		if (face[f].normal.dot(face[f].centroid-cell[face[f].parent].centroid)<=0.) {
			// [TBM] Need to swap the face and reflect the area vector
			cout << "[W Rank=" << Rank << "] Face " << f << " normal is pointing in to its parent ... fixing " << endl;
			cout << face[f].parent << endl;
			face[f].normal*=-1.;
			vector<int>::reverse_iterator rit;
			face[f].nodes.assign(face[f].nodes.rbegin(),face[f].nodes.rend());
		}
	}
	
	return 0;

}


Node::Node(double x, double y, double z) {
	comp[0]=x;
	comp[1]=y;
	comp[2]=z;
}

Cell::Cell(void) {
	;
}


Node& Cell::node(int n) {
	return grid.node[nodes[n]];
};

Face& Cell::face(int f) {
	return grid.face[faces[f]];
};

Node& Face::node(int n) {
	return grid.node[nodes[n]];
};


void Grid::gradMaps() {
	
	unsigned int f;
	map<int,double>::iterator it;
	Vec3D areaVec;
	for (unsigned int c=0;c<cellCount;++c) {
		cell[c].gradMap.clear();
		for (unsigned int cf=0;cf<cell[c].faceCount;++cf) {
			f=cell[c].faces[cf];
			if (face[f].bc==INTERNAL || face[f].bc==GHOST ) { // if internal or interpartition face
				areaVec=face[f].normal*face[f].area/cell[c].volume;
				if (face[f].parent!=c) areaVec*=-1.;
				for (it=face[f].average.begin();it!=face[f].average.end();it++) {
					if (cell[c].gradMap.find((*it).first)!=cell[c].gradMap.end()) { // if the cell contributing to the face average is already contained in the cell gradient map
						cell[c].gradMap[(*it).first]+=(*it).second*areaVec;
					} else {
						cell[c].gradMap.insert(pair<int,Vec3D>((*it).first,(*it).second*areaVec));
					}
				}
			} // end if internal face
		} // end cell face loop
	} // end cell loop
	
} // end Grid::gradMaps()
	

void Grid::gradients(void) {

	// Calculate cell gradients
	map<int,Vec3D>::iterator it;
	map<int,double>::iterator fit;
	unsigned int f;
	Vec3D faceVel,areaVec;
	double faceP,faceT,faceRho,faceK,faceOmega,faceMach;
	
	for (unsigned int c=0;c<cellCount;++c) {
		// Initialize all gradients to zero
		for (unsigned int i=0;i<7;++i) cell[c].grad[i]=0.;
		// Add internal and interpartition face contributions
		for (it=cell[c].gradMap.begin();it!=cell[c].gradMap.end(); it++ ) {
			if ((*it).first>=0) { // if contribution is coming from a real cell
				cell[c].grad[0]+=(*it).second*cell[(*it).first].p;
				for (unsigned int i=1;i<4;++i) cell[c].grad[i]+=(*it).second*cell[(*it).first].v.comp[i-1];
				cell[c].grad[4]+=(*it).second*cell[(*it).first].T;
				cell[c].grad[5]+=(*it).second*cell[(*it).first].k;
				cell[c].grad[6]+=(*it).second*cell[(*it).first].omega;
			} else { // if contribution is coming from a ghost cell
				cell[c].grad[0]+=(*it).second*ghost[-1*((*it).first+1)].p;
				for (unsigned int i=1;i<4;++i) cell[c].grad[i]+=(*it).second*ghost[-1*((*it).first+1)].v.comp[i-1];
				cell[c].grad[4]+=(*it).second*ghost[-1*((*it).first+1)].T;
				cell[c].grad[5]+=(*it).second*ghost[-1*((*it).first+1)].k;
				cell[c].grad[6]+=(*it).second*ghost[-1*((*it).first+1)].omega;
			}
		} // end gradMap loop

 		// Add boundary face contributions
		for (unsigned int cf=0;cf<cell[c].faceCount;++cf) {
			f=cell[c].faces[cf];
			if (face[f].bc>=0) { // if a boundary face
				areaVec=face[f].normal*face[f].area/cell[c].volume;
				
				faceP=0.; faceVel=0.;faceT=0.;;faceK=0.;faceOmega=0.;
				for (fit=face[f].average.begin();fit!=face[f].average.end();fit++) {
					if ((*fit).first>=0) { // if contribution is coming from a real cell
						faceP+=(*fit).second*cell[(*fit).first].p;
						faceVel+=(*fit).second*cell[(*fit).first].v;
						faceT+=(*fit).second*cell[(*fit).first].T;
						faceK+=(*fit).second*cell[(*fit).first].k;
						faceOmega+=(*fit).second*cell[(*fit).first].omega;
					} else { // if contribution is coming from a ghost cell
						faceP+=(*fit).second*ghost[-1*((*fit).first+1)].p;
						faceVel+=(*fit).second*ghost[-1*((*fit).first+1)].v;
						faceT+=(*fit).second*ghost[-1*((*fit).first+1)].T;
						faceK+=(*fit).second*ghost[-1*((*fit).first+1)].k;
						faceOmega+=(*fit).second*ghost[-1*((*fit).first+1)].omega;
					}
				}

 				faceRho=eos.rho(faceP,faceT);
				
				if (bc.region[face[f].bc].thermalType==ADIABATIC) {
					faceT=cell[face[f].parent].T;
					faceRho=eos.rho(faceP,faceT);
				}
				
				if (bc.region[face[f].bc].specified==BC_STATE) {
					faceP=bc.region[face[f].bc].p;
					faceT=bc.region[face[f].bc].T;
					faceRho=bc.region[face[f].bc].rho;
				} else if (bc.region[face[f].bc].specified==BC_P) {
					faceP=bc.region[face[f].bc].p;
					faceRho=eos.rho(faceP,faceT); // temperature is extrapolated
				} else if (bc.region[face[f].bc].specified==BC_T) {
					faceT=bc.region[face[f].bc].T;
					faceRho=eos.rho(faceP,faceT); // pressure is extrapolated
				} else if (bc.region[face[f].bc].specified==BC_RHO) {
					faceRho=bc.region[face[f].bc].rho;
					faceT=eos.T(faceP,faceRho); // pressure is extrapolated
				} // If nothing is specified, everything is extrapolated
				
				if (bc.region[face[f].bc].type==INLET) {
					faceVel=0.5*(faceVel+bc.region[face[f].bc].v);
					faceK=bc.region[face[f].bc].k;
					faceOmega=bc.region[face[f].bc].omega;
				} else if (bc.region[face[f].bc].type==SLIP) { 
					faceVel-=faceVel.dot(face[f].normal)*face[f].normal;
					faceK=0.;
					faceOmega=60.*viscosity/(faceRho*0.075*pow(fabs((cell[face[f].parent].centroid-face[f].centroid).dot(face[f].normal)),2.));	
				} else if (bc.region[face[f].bc].type==SYMMETRY) {
					// Symmetry mirrors everything
					faceP=cell[face[f].parent].p;
					faceVel=cell[face[f].parent].v;
					faceVel-=faceVel.dot(face[f].normal)*face[f].normal;
					faceK=cell[face[f].parent].k;
					faceOmega=cell[face[f].parent].omega;
				} else if (bc.region[face[f].bc].type==NOSLIP) {
					faceVel=0.;
					faceK=0.;
					faceOmega=60.*viscosity/(faceRho*0.075*pow(fabs((cell[face[f].parent].centroid-face[f].centroid).dot(face[f].normal)),2.));
				}
				
				cell[c].grad[0]+=faceP*areaVec;
				cell[c].grad[1]+=faceVel[0]*areaVec;
				cell[c].grad[2]+=faceVel[1]*areaVec;
				cell[c].grad[3]+=faceVel[2]*areaVec;
				cell[c].grad[4]+=faceT*areaVec;
				cell[c].grad[5]+=faceK*areaVec;
				cell[c].grad[6]+=faceOmega*areaVec;

			} // end if a boundary face
		} // end cell face loop
	} // end cell loop
	
 	
} // end Grid::gradients(void)

void Grid::limit_gradients(void) {
	
	unsigned int neighbor,g;
	Vec3D maxGrad[7],minGrad[7];
	if(LIMITER==NONE) {
		for (unsigned int c=0;c<cellCount;++c) for (unsigned int var=0;var<7;++var) for (unsigned int comp=0;comp<3;++comp) cell[c].limited_grad[var].comp[comp]=cell[c].grad[var].comp[comp];
	} else {
		for (unsigned int c=0;c<cellCount;++c) {
	
			// Initialize min and max to current cells values
			for (unsigned int i=0;i<7;++i) maxGrad[i]=minGrad[i]=cell[c].grad[i];
			// Find extremes in neighboring real cells
			for (unsigned int cc=0;cc<cell[c].neighborCellCount;++cc) {
				neighbor=cell[c].neighborCells[cc];
				for (unsigned int var=0;var<7;++var) {
					for (unsigned int comp=0;comp<3;++comp) {
						maxGrad[var].comp[comp]=max(
								maxGrad[var].comp[comp],
								(1.-limiter_sharpening)*cell[neighbor].grad[var].comp[comp]
									+limiter_sharpening*cell[c].grad[var].comp[comp]);
						minGrad[var].comp[comp]=min(minGrad[var].comp[comp],
								(1.-limiter_sharpening)*cell[neighbor].grad[var].comp[comp]
									+limiter_sharpening*cell[c].grad[var].comp[comp]);
					}
				}
			}
			// Find extremes in neighboring ghost cells
			for (unsigned int cg=0;cg<cell[c].ghostCount;++cg) {
				g=cell[c].ghosts[cg];
				for (unsigned int var=0;var<7;++var) {
					for (unsigned int comp=0;comp<3;++comp) {
						maxGrad[var].comp[comp]=max(
								maxGrad[var].comp[comp],
								(1.-limiter_sharpening)*ghost[g].grad[var].comp[comp]
									+limiter_sharpening*cell[c].grad[var].comp[comp]);
						minGrad[var].comp[comp]=min(
								minGrad[var].comp[comp],
								(1.-limiter_sharpening)*ghost[g].grad[var].comp[comp]
									+limiter_sharpening*cell[c].grad[var].comp[comp]);
					}
				}
			}
			
			if(LIMITER==MINMOD) for (unsigned int var=0;var<7;++var) for (unsigned int comp=0;comp<3;++comp) cell[c].limited_grad[var].comp[comp]=minmod(maxGrad[var].comp[comp],minGrad[var].comp[comp]);
			if(LIMITER==DOUBLEMINMOD) for (unsigned int var=0;var<7;++var) for (unsigned int comp=0;comp<3;++comp) cell[c].limited_grad[var].comp[comp]=doubleMinmod(maxGrad[var].comp[comp],minGrad[var].comp[comp]);
			if(LIMITER==HARMONIC) for (unsigned int var=0;var<7;++var) for (unsigned int comp=0;comp<3;++comp) cell[c].limited_grad[var].comp[comp]=harmonic(maxGrad[var].comp[comp],minGrad[var].comp[comp]);
			if(LIMITER==SUPERBEE) for (unsigned int var=0;var<7;++var) for (unsigned int comp=0;comp<3;++comp) cell[c].limited_grad[var].comp[comp]=superbee(maxGrad[var].comp[comp],minGrad[var].comp[comp]);
			
			
		}
	}
	
	return;

} // end Grid::limit_gradients()

void Grid::lengthScales(void) {
	// Loop through the cells and calculate length scales
	for (unsigned int c=0;c<cellCount;++c) {
		cell[c].lengthScale=1.e20;
		double height;
		unsigned int f;
		for (unsigned int cf=0;cf<cell[c].faceCount;++cf) {
			f=cell[c].faces[cf];
			height=fabs(face[f].normal.dot(face[f].centroid-cell[c].centroid));
			bool skipScale=false;
			if (face[f].bc>=0) {
				if (bc.region[face[f].bc].type==SYMMETRY) {
					skipScale=true;
				}
			}
			if (!skipScale) cell[c].lengthScale=min(cell[c].lengthScale,height);
		}
	}
	return;
} // end Grid::lengthScales
