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
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <iomanip>
#include<cmath>
using namespace std;
#include <cgnslib.h>
#include <parmetis.h>


#include "grid.h"
#include "bc.h"

#define EPS 1e-10

				 
extern Grid grid;
extern BC bc;
extern int np, Rank;
extern double Gamma;
extern double Pref;

double superbee(double a, double b);
double minmod(double a, double b);
int gelimd(double **a,double *b,double *x, int n);

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

int Grid::areas_volumes() {

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
		
		// Sum the area as a patch of triangles formed by connected two nodes and an interior point
		face[f].centroid=0.;
		areaVec=0.;
		int next;
		for (unsigned int n=0;n<face[f].nodeCount;++n) {
			next=n+1;
			if (next==face[f].nodeCount) next=0;
			patchArea=0.5*(face[f].node(n)-centroid).cross(face[f].node(next)-centroid);
			patchCentroid=1./3.*(face[f].node(n)+face[f].node(next)+centroid);
			face[f].centroid+=patchCentroid*fabs(patchArea);
			areaVec+=patchArea;
			
		}
		if (areaVec.dot(centroid-cell[face[f].parent].centroid) <=0.) {
			// [TBM] Need to swap the face and reflect the area vector
			cout << "face " << f << " should be swapped" << endl;
		}
		face[f].area=fabs(areaVec);
		face[f].centroid/=face[f].area;
		face[f].normal=areaVec/face[f].area;
	}
			
	// Loop through the cells and calculate the volumes
	double totalVolume=0.;
	for (unsigned int c=0;c<cellCount;++c) {
		double volume=0.,height;
		unsigned int f;
		for (unsigned int cf=0;cf<cell[c].faceCount;++cf) {
			f=cell[c].faces[cf];
			height=fabs(face[f].normal.dot(face[f].centroid-cell[c].centroid));
			volume+=1./3.*face[f].area*height;
		}
		cell[c].volume=volume;
		totalVolume+=volume;
	}
	
	cout << "[I Rank=" << Rank << "] Total Volume= " << setw(16) << setprecision(8) << scientific << totalVolume << endl;
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

class InterpolationTriangle {
	public:
		unsigned int c1,c2,c3;
		Vec3D cv1,cv2,cv3,centroid;
		double edge1,edge2,edge3;
		double perimeter;
		double area;
		double minEdge,maxEdge;
		double distance;
		double quality;
		double cw1,cw2,cw3;
		double weight;
};


void Grid::nodeAverages() {
	
	double **a;
	double *b;
	double *weights;
	
	for (unsigned int n=0;n<nodeCount;++n) {
		set<int> stencil; //interpolation stencil
		set<int>::iterator sit,sit1,sit2;
		int c1,c2,c3,c4;
		// Initialize stencil to nearest neighbor cells
		for (int nc=0;nc<node[n].cells.size();++nc) stencil.insert(node[n].cells[nc]);
		string method;
		Vec3D planeNormal;
		// if stencil doesn't have at least 4 points, expand it to include 2nd nearest neighbor cells
		if (stencil.size()<4) {
			for (int nc=0;nc<node[n].cells.size();++nc) {
				int ncell=node[n].cells[nc];
				for (int cc=0;cc<cell[ncell].neighborCells.size();++cc) {
					stencil.insert(cell[ncell].neighborCells[cc]);
				}
			}
		}
		if (stencil.size()>=4) { // Candidate for tetra interpolation
			// Pick first 3 points in the stencil
			sit=stencil.begin();
			c1=*sit; sit++;
			c2=*sit; sit++;
			c3=*sit;
			// Calculate the normal vector of the plane they form
			planeNormal=((cell[c2].centroid-cell[c1].centroid).cross(cell[c3].centroid-cell[c2].centroid)).norm();
			// Check the remaining stencil cell centroids to see if they are on the same plane too
			method="tri";
			for (sit=sit;sit!=stencil.end();sit++) {
				// If the centroid lies on the plane fomed by first three
				if (fabs(planeNormal.dot((cell[*sit].centroid-cell[*(stencil.begin())].centroid).norm()))>1.e-7) {
					method="tetra";
					// Finding on off-plane point is enough
					// That means tetra method can be used
					break;
				}
			}
		}
		else if (stencil.size()==3) { method="tri";}
		else if (stencil.size()==2) { method="line";}
		else if (stencil.size()==1) { method="point";}

		// TODO implement tetra method
		if (method=="tri") {
			vector<InterpolationTriangle> triangles;
			// Find out best combination neigboring cells for triangulation
			double edge1,edge2,edge3,maxEdge,minEdge,perimeter,area;
			double areaEquil,quality,moment;
			double minMoment=1.e20,bestQuality=0.;
			Vec3D p1,p2,p3,centroid;
			unsigned int cq1,cq2,cq3;
			bool foundOne=false;
			// Loop through the stencil and construct triangle combinations
			double weightSum=0.;
			for (sit=stencil.begin();sit!=stencil.end();sit++) {
				sit1=sit; sit1++;
				for (sit1=sit1;sit1!=stencil.end();sit1++) {
					sit2=sit; sit2++;
					for (sit2=sit2;sit2!=stencil.end();sit2++) {
						InterpolationTriangle temp;
						temp.c1=*sit;
						temp.c2=*sit1;
						temp.c3=*sit2;
						temp.cv1=cell[*sit].centroid;
						temp.cv2=cell[*sit1].centroid;
						temp.cv3=cell[*sit2].centroid;
						temp.centroid=(temp.cv1+temp.cv2+temp.cv3)/3.;
						temp.area=0.5*fabs((temp.cv2-temp.cv1).cross(temp.cv3-temp.cv1));
						temp.edge1=fabs(temp.cv1-temp.cv2);
						temp.edge2=fabs(temp.cv1-temp.cv3);
						temp.edge3=fabs(temp.cv2-temp.cv3);
						temp.perimeter=temp.edge1+temp.edge2+temp.edge3;
						temp.maxEdge=max(temp.edge1,temp.edge2); temp.maxEdge=max(temp.maxEdge,temp.edge3);
						temp.minEdge=min(temp.edge1,temp.edge2); temp.minEdge=min(temp.minEdge,temp.edge3);
						temp.quality=temp.area/((temp.perimeter*temp.perimeter/9.)*sqrt(3.)/4.);
						temp.distance=fabs(node[n]-temp.centroid);
						temp.weight=temp.quality/(temp.distance*temp.distance);
						if (temp.quality>1.e-7) {
							triangles.push_back(temp);
							weightSum+=temp.weight;
						}
					}
				}
			} // Loop through stencil
			if (triangles.size()==0) cerr << "[E] An interpolation triangle couldn't be found for node " << n << endl;
			// Now loop over the triangles and store interpolation weights
			for (int t=0;t<triangles.size();++t) {
				triangles[t].weight/=weightSum;
				
				Vec3D pp,origin,basis1,basis2;
				origin=cell[triangles[t].c1].centroid;
				basis1=(cell[triangles[t].c2].centroid-origin).norm();
				planeNormal=(basis1.cross((cell[triangles[t].c3].centroid-origin).norm()));
				// Normalize the plane normal vector
				planeNormal/=fabs(planeNormal);
				basis2=-1.*(basis1.cross(planeNormal)).norm();
				// Project the node point to the plane
				pp=node[n]-origin;
				pp-=pp.dot(planeNormal)*planeNormal;
				pp+=origin;
				// Form the linear system
				a = new double* [3];
				for (int i=0;i<3;++i) a[i]=new double[3];
				b = new double[3];
				weights= new double [3];
				a[0][0]=(cell[triangles[t].c1].centroid).dot(basis1);
				a[0][1]=(cell[triangles[t].c2].centroid).dot(basis1);
				a[0][2]=(cell[triangles[t].c3].centroid).dot(basis1);
				a[1][0]=(cell[triangles[t].c1].centroid).dot(basis2);
				a[1][1]=(cell[triangles[t].c2].centroid).dot(basis2);
				a[1][2]=(cell[triangles[t].c3].centroid).dot(basis2);
				a[2][0]=1.;
				a[2][1]=1.;
				a[2][2]=1.;
				b[0]=pp.dot(basis1);
				b[1]=pp.dot(basis2);
				b[2]=1.;
				// Solve the 3x3 linear system by Gaussion Elimination
				gelimd(a,b,weights,3);
				triangles[t].cw1=weights[0];
				triangles[t].cw2=weights[1];
				triangles[t].cw3=weights[2];
			}

			node[n].average.clear();
			for (int t=0;t<triangles.size();++t) {
				if (node[n].average.find(triangles[t].c1)!=node[n].average.end()) { // if the cell contributing to the node average is already contained in the map
					node[n].average[triangles[t].c1]+=triangles[t].weight*triangles[t].cw1;
				} else {
					node[n].average.insert(pair<int,double>(triangles[t].c1,triangles[t].weight*triangles[t].cw1));
				}
				
				if (node[n].average.find(triangles[t].c2)!=node[n].average.end()) { // if the cell contributing to the node average is already contained in the map
					node[n].average[triangles[t].c2]+=triangles[t].weight*triangles[t].cw2;
				} else {
					node[n].average.insert(pair<int,double>(triangles[t].c2,triangles[t].weight*triangles[t].cw2));
				}
				
				if (node[n].average.find(triangles[t].c3)!=node[n].average.end()) { // if the cell contributing to the node average is already contained in the map
					node[n].average[triangles[t].c3]+=triangles[t].weight*triangles[t].cw3;
				} else {
					node[n].average.insert(pair<int,double>(triangles[t].c3,triangles[t].weight*triangles[t].cw3));
				}
				
			}

// 				// Insert the weights
// 			node[n].average.clear();
// 			node[n].average.insert(pair<int,double>(triangles[0].c1,triangles[0].cw1));
// 			node[n].average.insert(pair<int,double>(triangles[0].c2,triangles[0].cw2));
// 			node[n].average.insert(pair<int,double>(triangles[0].c3,triangles[0].cw3));
			

		} // end if method is tri
		if (method=="line") {
			// Chose the closest 2 points in stencil
			double distanceMin=1.e20;
			double distance;
			for (sit=stencil.begin();sit!=stencil.end();sit++) {
				distance=fabs(node[n]-cell[*sit].centroid);
				if (distance<distanceMin) {
					distanceMin=distance;
					sit1=sit;
				}
			}
			distanceMin=1.e20;
			for (sit=stencil.begin();sit!=stencil.end();sit++) {
				distance=fabs(node[n]-cell[*sit].centroid);
				if (distance<distanceMin && sit!=sit1) {
					distanceMin=distance;
					sit2=sit;
				}
			}
			Vec3D pp,p1,p2,direction;
			double a,weight1,weight2;
			p1=cell[*sit1].centroid;
			p2=cell[*sit2].centroid;
			direction=(p2-p1).norm();
			a=(node[n]-p1).dot(direction)/fabs(p2-p1);
			weight1=a;
			weight2=1.-a;
			node[n].average.clear();
			node[n].average.insert(pair<int,double>(*sit1,weight1));
			node[n].average.insert(pair<int,double>(*sit2,weight2));
		} // end if method is line
		if (method=="point") {
			node[n].average.clear();
			node[n].average.insert(pair<int,double>(*(stencil.begin()),1.));
		} // end if method is point
		
		if (method!="tri") cout << method << endl;
		std::map<int,double>::iterator it;
		double avg=0.;
		for ( it=node[n].average.begin() ; it != node[n].average.end(); it++ ) {
			avg+=(*it).second*cell[(*it).first].rho;
		}
// 		if (fabs(avg-2.*node[n].comp[0]-2.)>1.e-7) cout << 2.*node[n].comp[0]+2. << "\t" << avg << endl;

	} // node loop
	
	return;
} // end Grid::nodeAverages()

void Grid::faceAverages() {

	unsigned int n;
	map<int,double>::iterator it;
	map<int,double>::iterator fit;

 	for (unsigned int f=0;f<faceCount;++f) {
		double weightSumCheck;
		vector<double> nodeWeights;
		// Mid point of the face
		Vec3D mid=0.,patchArea;
		double patchRatio;
		for (unsigned int fn=0;fn<face[f].nodeCount;++fn) {
			nodeWeights.push_back(0.);
			mid+=face[f].node(fn);
		}
		mid/=double(face[f].nodeCount);
		// Get node weights
		int next;
		for (unsigned int fn=0;fn<face[f].nodeCount;++fn) {
			next=fn+1;
			if (next==face[f].nodeCount) next=0;
			patchArea=0.5*(face[f].node(fn)-mid).cross(face[f].node(next)-mid);
			patchRatio=fabs(patchArea)/face[f].area;
			for (unsigned int fnn=0;fnn<face[f].nodeCount;++fnn) {
				nodeWeights[fnn]+=(1./double(3*face[f].nodeCount))*patchRatio;
			}
			nodeWeights[fn]+=patchRatio/3.;
			nodeWeights[next]+=patchRatio/3.;
		}

		face[f].average.clear();
		for (unsigned int fn=0;fn<face[f].nodeCount;++fn) {
			n=face[f].nodes[fn];
			for ( it=node[n].average.begin() ; it != node[n].average.end(); it++ ) {
				if (face[f].average.find((*it).first)!=face[f].average.end()) { // if the cell contributing to the node average is already contained in the face average map
					face[f].average[(*it).first]+=nodeWeights[fn]*(*it).second;
				} else {
					face[f].average.insert(pair<int,double>((*it).first,nodeWeights[fn]*(*it).second));
				}
			}
		}
// 		weightSumCheck=0.;
// 		for ( fit=face[f].average.begin() ; fit != face[f].average.end(); fit++ ) {
// 			//(*fit).second/=double(face[f].nodeCount);
// 			weightSumCheck+=(*fit).second;
// 		}
// 
// 		cout << nodeWeights.size() << "\t" << weightSumCheck << endl;
		//cout << setw(16) << setprecision(8) << scientific << weightSumCheck << endl;
// 		std::map<int,double>::iterator it;
// 		double avg=0.;
// 		for ( it=face[f].average.begin() ; it != face[f].average.end(); it++ ) {
// 			avg+=(*it).second*cell[(*it).first].rho;
// 		}
// 		if (fabs(avg-2.*face[f].centroid.comp[0]-2.)>1.e-8) cout << 2.*face[f].centroid.comp[0]+2. << "\t" << avg << endl;
	} // end face loop

} // end Grid::faceAverages()

void Grid::gradMaps() {
	
	unsigned int f;
	map<int,double>::iterator it;
	Vec3D areaVec;
	for (unsigned int c=0;c<cellCount;++c) {
		cell[c].gradMap.clear();
		for (unsigned int cf=0;cf<cell[c].faceCount;++cf) {
			f=cell[c].faces[cf];
			if (face[f].bc<0 ) { // if internal or interpartition face
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
	double faceRho,faceP,faceK,faceOmega;
	double Mach;
	
	for (unsigned int c=0;c<cellCount;++c) {
		// Initialize all gradients to zero
		for (unsigned int i=0;i<7;++i) cell[c].grad[i]=0.;
		// Add internal and interpartition face contributions
		for (it=cell[c].gradMap.begin();it!=cell[c].gradMap.end(); it++ ) {
			if ((*it).first>=0) { // if contribution is coming from a real cell
				cell[c].grad[0]+=(*it).second*cell[(*it).first].rho;
				for (unsigned int i=1;i<4;++i) cell[c].grad[i]+=(*it).second*cell[(*it).first].v.comp[i-1];
				cell[c].grad[4]+=(*it).second*cell[(*it).first].p;
				cell[c].grad[5]+=(*it).second*cell[(*it).first].k;
				cell[c].grad[6]+=(*it).second*cell[(*it).first].omega;
			} else { // if contribution is coming from a ghost cell
				cell[c].grad[0]+=(*it).second*ghost[-1*((*it).first+1)].rho;
				for (unsigned int i=1;i<4;++i) cell[c].grad[i]+=(*it).second*ghost[-1*((*it).first+1)].v.comp[i-1];
				cell[c].grad[4]+=(*it).second*ghost[-1*((*it).first+1)].p;
				cell[c].grad[5]+=(*it).second*ghost[-1*((*it).first+1)].k;
				cell[c].grad[6]+=(*it).second*ghost[-1*((*it).first+1)].omega;
			}
		} // end gradMap loop

 		// Add boundary face contributions

		for (unsigned int cf=0;cf<cell[c].faceCount;++cf) {
			f=cell[c].faces[cf];
			if (face[f].bc>=0) { // if a boundary face
				areaVec=face[f].normal*face[f].area/cell[c].volume;
				if (face[f].parent!=c) areaVec*=-1.;

				faceVel=0.;faceRho=0.;faceP=0.;faceK=0.;faceOmega=0.;
				for (fit=face[f].average.begin();fit!=face[f].average.end();fit++) {
					if ((*fit).first>=0) { // if contribution is coming from a real cell
						faceRho+=(*fit).second*cell[(*fit).first].rho;
						faceVel+=(*fit).second*cell[(*fit).first].v;
						faceP+=(*fit).second*cell[(*fit).first].p;
						faceK+=(*fit).second*cell[(*fit).first].k;
						faceOmega+=(*fit).second*cell[(*fit).first].omega;
					} else { // if contribution is coming from a ghost cell
						faceRho+=(*fit).second*ghost[-1*((*fit).first+1)].rho;
						faceVel+=(*fit).second*ghost[-1*((*fit).first+1)].v;
						faceP+=(*fit).second*ghost[-1*((*fit).first+1)].p;
						faceK+=(*fit).second*ghost[-1*((*fit).first+1)].k;
						faceOmega+=(*fit).second*ghost[-1*((*fit).first+1)].omega;
					}
				}

				//cout << 2.*face[f].centroid.comp[0]+2. << "\t" << faceRho << endl;

				if (bc.region[face[f].bc].type=="inlet") {
					faceRho=bc.region[face[f].bc].rho;
					faceVel=bc.region[face[f].bc].v;
					faceK=bc.region[face[f].bc].k;
					faceOmega=bc.region[face[f].bc].omega;
					//faceP=cell[face[f].parent].p;
					//faceP=bc.region[face[f].bc].p;
				}

				if (bc.region[grid.face[f].bc].type=="outlet" &&
					bc.region[grid.face[f].bc].kind=="fixedPressure") {
					// find Mach number
// 					Mach=(cell[c].v.dot(face[f].normal))/sqrt(Gamma*(cell[c].p+Pref)/cell[c].rho);
// 					if (Mach<1.) faceP=bc.region[face[f].bc].p;
					faceP=bc.region[face[f].bc].p;//-0.5*faceRho*faceVel.dot(faceVel);
				}
				// Kill the wall normal component for slip or symmetry, pressure and rho is extrapolated
				if (bc.region[face[f].bc].type=="slip") faceVel-=faceVel.dot(face[f].normal)*face[f].normal;
				if (bc.region[face[f].bc].type=="symmetry") {
					faceRho=cell[face[f].parent].rho;
					faceVel=cell[face[f].parent].v;
					faceP=cell[face[f].parent].p;
					faceVel-=faceVel.dot(face[f].normal)*face[f].normal;
				}
				// Kill the velocity for no-slip, pressure and rho is extrapolated
				if (bc.region[face[f].bc].type=="noslip") faceVel=0.;

				cell[c].grad[0]+=faceRho*areaVec;
				cell[c].grad[1]+=faceVel.comp[0]*areaVec;
				cell[c].grad[2]+=faceVel.comp[1]*areaVec;
				cell[c].grad[3]+=faceVel.comp[2]*areaVec;
				cell[c].grad[4]+=faceP*areaVec;

			} // end if a boundary face
		} // end cell face loop
		//cout << cell[c].grad[0].comp[0] << endl;
		//cout << c << "\t" << areaVecSum/cell[c].volume << endl;
	} // end cell loop
	
 	
} // end Grid::gradients(void)

void Grid::limit_gradients(string limiter, double sharpeningFactor) {
	
	unsigned int neighbor,g;
	Vec3D maxGrad[7],minGrad[7];
	
	for (unsigned int c=0;c<cellCount;++c) {

		// Initialize min and max to current cells values
		for (unsigned int i=0;i<7;++i) maxGrad[i]=minGrad[i]=cell[c].grad[i];
		// Find extremes in neighboring real cells
		for (unsigned int cc=0;cc<cell[c].neighborCellCount;++cc) {
			neighbor=cell[c].neighborCells[cc];
			for (unsigned int var=0;var<7;++var) {
				for (unsigned int comp=0;comp<3;++comp) {
					maxGrad[var].comp[comp]=max(maxGrad[var].comp[comp],(1.-sharpeningFactor)*cell[neighbor].grad[var].comp[comp]+sharpeningFactor*cell[c].grad[var].comp[comp]);
					minGrad[var].comp[comp]=min(minGrad[var].comp[comp],(1.-sharpeningFactor)*cell[neighbor].grad[var].comp[comp]+sharpeningFactor*cell[c].grad[var].comp[comp]);
				}
			}
		}
		// Find extremes in neighboring ghost cells
		for (unsigned int cg=0;cg<cell[c].ghostCount;++cg) {
			g=cell[c].ghosts[cg];
			for (unsigned int var=0;var<7;++var) {
				for (unsigned int comp=0;comp<3;++comp) {
					maxGrad[var].comp[comp]=max(maxGrad[var].comp[comp],(1.-sharpeningFactor)*ghost[g].grad[var].comp[comp]+sharpeningFactor*cell[c].grad[var].comp[comp]);
					minGrad[var].comp[comp]=min(minGrad[var].comp[comp],(1.-sharpeningFactor)*ghost[g].grad[var].comp[comp]+sharpeningFactor*cell[c].grad[var].comp[comp]);
				}
			}
		}
		if(limiter=="superbee") for (unsigned int var=0;var<7;++var) for (unsigned int comp=0;comp<3;++comp) cell[c].limited_grad[var].comp[comp]=superbee(maxGrad[var].comp[comp],minGrad[var].comp[comp]);
		if(limiter=="minmod") for (unsigned int var=0;var<7;++var) for (unsigned int comp=0;comp<3;++comp) cell[c].limited_grad[var].comp[comp]=minmod(maxGrad[var].comp[comp],minGrad[var].comp[comp]);

	}

} // end Grid::limit_gradients(string limiter, double sharpeningFactor)

int gelimd(double **a,double *b,double *x, int n)
{
	double tmp,pvt,*t;
	int i,j,k,itmp;

	for (i=0;i<n;i++) {             // outer loop on rows
		pvt = a[i][i];              // get pivot value
		if (fabs(pvt) < EPS) {
			for (j=i+1;j<n;j++) {
				if(fabs(pvt = a[j][i]) >= EPS) break;
			}
			if (fabs(pvt) < EPS) return 1;     // nowhere to run!
			t=a[j];                 // swap matrix rows...
			a[j]=a[i];
			a[i]=t;
			tmp=b[j];               // ...and result vector
			b[j]=b[i];
			b[i]=tmp;
		}
// (virtual) Gaussian elimination of column
		for (k=i+1;k<n;k++) {       // alt: for (k=n-1;k>i;k--)
			tmp = a[k][i]/pvt;
			for (j=i+1;j<n;j++) {   // alt: for (j=n-1;j>i;j--)
				a[k][j] -= tmp*a[i][j];
			}
			b[k] -= tmp*b[i];
		}
	}
// Do back substitution
	for (i=n-1;i>=0;i--) {
		x[i]=b[i];
		for (j=n-1;j>i;j--) {
			x[i] -= a[i][j]*x[j];
		}
		x[i] /= a[i][i];
	}
	return 0;
}

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
				if (bc.region[face[f].bc].type=="symmetry") {
					skipScale=true;
				}
			}
			if (skipScale==false) cell[c].lengthScale=min(cell[c].lengthScale,height);
		}
	}
	return;
} // end Grid::lengthScales
