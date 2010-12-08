/************************************************************************
	
	Copyright 2007-2010 Emre Sozer

	Contact: emresozer@freecfd.com

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
#include "grid.h"
#include "inputs.h"
#include <iomanip>
using namespace std;

extern vector<Grid> grid;
extern InputFile input;

void gradient_maps(int gid) {

	
	bool DEBUG=false;
	
	set<int> stencil;
	set<int>::iterator sit,sit1,sit2,sit3;
	Vec3D my_normal,other_normal,cell2face;
	int opposite_face;
	double product,min_product,max_product;
	vector<int> face_pairs (6,-1);
	vector<int> cell_pairs (6,-1);
	vector<Vec3D> Jac; // Changes in othogonal system and curvilinear system
	Jac.resize(3);
	vector<Vec3D> Jac_invT; // Inverse Jacobian Transposed
	Jac_invT.resize(3);
	double det; // Determinant of the Jacobian
	int cell_plus,cell_minus;
	/*
	Jac:
	| dx/dxsi  dy/dxsi  dz/dxsi  | 
	| dx/deta  dy/deta  dz/deta  | 
 	| dx/dzeta dy/dzeta dz/dzeta | 
	*/
	
	// Loop all the cells
	for (int c=0;c<grid[gid].cellCount;++c) {
			
		for (int cf=0;cf<grid[gid].cell[c].faceCount;++cf) stencil.insert(cf); // Note that these are not actual face indices
		
		if (grid[gid].cell[c].nodeCount==8) { //If hexa cell
			 // Put the faces of the cell into a set (except the first one)
			 for (int cf=0;cf<grid[gid].cell[c].faceCount;++cf) stencil.insert(cf); // Note that these are not actual face indices	 
			 // Loop the stencil
			 int counter=0;
			 for (sit1=stencil.begin();sit1!=stencil.end();sit1++) {
				 // For each face, find the corresponding opposite
				 min_product=1.e20;
				 my_normal=grid[gid].cellFace(c,*sit1).normal;
				 if (grid[gid].cellFace(c,*sit1).parent!=c) my_normal*=-1.;
				 face_pairs[2*counter+1]=*sit1;
				 sit2=sit1; sit2++;
				 for (sit2=sit2;sit2!=stencil.end();sit2++) {
				 other_normal=grid[gid].cellFace(c,*sit2).normal;
				 // Make sure these normals are pointing outward
				 if (grid[gid].cellFace(c,*sit2).parent!=c) other_normal*=-1.;
				 // Pick the face with the normal in the most opposite direction to the current normal
				 product=my_normal.dot(other_normal);
				 if (product<min_product) {
				 min_product=product;
				 face_pairs[2*counter]=*sit2;
				 }
				 }
				 stencil.erase(face_pairs[2*counter]); 
				 stencil.erase(face_pairs[2*counter+1]);
				 counter++;
			 }
		} else {
			//First internal face will define the first direction
			sit1=stencil.begin();
			for (sit=stencil.begin();sit!=stencil.end();sit++) {
				if (grid[gid].cellFace(c,*sit).bc==INTERNAL_FACE) {
					sit1=sit;
					break;
				}
			}
			my_normal=grid[gid].cellFace(c,*sit1).normal;
			// Find the most perpendicular face to the first one
			// First restrict the search to non-boundary faces
			sit2=sit1;
			min_product=1.e20;
			for (sit=stencil.begin();sit!=stencil.end();sit++) {
				if (sit!=sit1) {
					if (grid[gid].cellFace(c,*sit).bc==INTERNAL_FACE) {
						other_normal=grid[gid].cellFace(c,*sit).normal;
						product=fabs(my_normal.dot(other_normal));
						if (product<min_product) {
							min_product=product;
							sit2=sit;
						}
					}
				}
			}
			// Check if an internal face is picked, if not remove the restriction
			min_product=1.e20;
			if (sit2==sit1) { // Means no selection made in the previous loop
				for (sit=stencil.begin();sit!=stencil.end();sit++) {
					if (sit!=sit1) {
						other_normal=grid[gid].cellFace(c,*sit).normal;
						product=fabs(my_normal.dot(other_normal));
						if (product<min_product) {
							min_product=product;
							sit2=sit;
						}
					}
				}
			}
			// Find the face with the normal most parallel to cross product of the first two
			max_product=0.;
			my_normal=grid[gid].cellFace(c,*sit1).normal.cross(grid[gid].cellFace(c,*sit2).normal);
			for (sit=stencil.begin();sit!=stencil.end();sit++) {
				if (sit!=sit1 && sit!=sit2) {
					other_normal=grid[gid].cellFace(c,*sit).normal;
					product=fabs(my_normal.dot(other_normal));
					if (product>max_product) {
						max_product=product;
						sit3=sit;
					}
				}
			}		
			face_pairs[1]=*sit1;
			face_pairs[3]=*sit2;
			face_pairs[5]=*sit3;

			// Just initialize these
			face_pairs[0]=*sit2;
			face_pairs[2]=*sit3;
			face_pairs[4]=*sit1;
			
			// Now find the pairs of the 3 faces (most opposing faces)
			min_product=1.e20;
			my_normal=grid[gid].cellFace(c,*sit1).normal;
			if (grid[gid].cellFace(c,*sit1).parent!=c) my_normal*=-1.;
			for (sit=stencil.begin();sit!=stencil.end();sit++) {
				other_normal=grid[gid].cellFace(c,*sit).normal;
				// Make sure these normals are pointing outward
				if (grid[gid].cellFace(c,*sit).parent!=c) other_normal*=-1.;
				// Pick the face with the normal in the most opposite direction to the current normal
				product=my_normal.dot(other_normal);
				if (product<min_product) {
					min_product=product;
					face_pairs[0]=*sit;
				}
			}
			min_product=1.e20;
			my_normal=grid[gid].cellFace(c,*sit2).normal;
			if (grid[gid].cellFace(c,*sit2).parent!=c) my_normal*=-1.;
			for (sit=stencil.begin();sit!=stencil.end();sit++) {
				other_normal=grid[gid].cellFace(c,*sit).normal;
				// Make sure these normals are pointing outward
				if (grid[gid].cellFace(c,*sit).parent!=c) other_normal*=-1.;
				// Pick the face with the normal in the most opposite direction to the current normal
				product=my_normal.dot(other_normal);
				if (product<min_product) {
					if (!(*sit==face_pairs[1] && face_pairs[3]==face_pairs[0])) { // Avoid the same but reverse direction pick
						min_product=product;
						face_pairs[2]=*sit;
					}
				}
			}		
			min_product=1.e20;
			my_normal=grid[gid].cellFace(c,*sit3).normal;
			if (grid[gid].cellFace(c,*sit3).parent!=c) my_normal*=-1.;
			for (sit=stencil.begin();sit!=stencil.end();sit++) {
				other_normal=grid[gid].cellFace(c,*sit).normal;
				// Make sure these normals are pointing outward
				if (grid[gid].cellFace(c,*sit).parent!=c) other_normal*=-1.;
				// Pick the face with the normal in the most opposite direction to the current normal
				product=my_normal.dot(other_normal);
				if (product<min_product) {
					if (!(*sit==face_pairs[1] && face_pairs[5]==face_pairs[0])) { // Avoid the same but reverse direction pick
						if (!(*sit==face_pairs[3] && face_pairs[5]==face_pairs[2])) {
							min_product=product;
							face_pairs[4]=*sit;
						}
					}
				}
			}
		}
		
		stencil.clear();
		
		// Now we know the three directions
		for (int i=0;i<3;++i) {
			
			cell_plus=grid[gid].cellFace(c,face_pairs[2*i+1]).neighbor;
			if (cell_plus==c) cell_plus=grid[gid].cellFace(c,face_pairs[2*i+1]).parent;
			cell_pairs[2*i+1]=cell_plus;
			if (grid[gid].cellFace(c,face_pairs[2*i+1]).bc>=0) { // Face is at a boundary
				cell_pairs[2*i+1]=c;
				Jac[i]=grid[gid].cellFace(c,face_pairs[2*i+1]).centroid;
			} else if (cell_plus<0) { // An inter-partition ghost cell
				Jac[i]=grid[gid].ghost[-cell_plus-1].centroid;
			} else {
				Jac[i]=grid[gid].cell[cell_plus].centroid;
			}
			
			cell_minus=grid[gid].cellFace(c,face_pairs[2*i]).neighbor;
			if (cell_minus==c) cell_minus=grid[gid].cellFace(c,face_pairs[2*i]).parent;
			cell_pairs[2*i]=cell_minus;
			if (grid[gid].cellFace(c,face_pairs[2*i]).bc>=0) { // Face is at a boundary
				cell_pairs[2*i]=c;
				Jac[i]-=grid[gid].cellFace(c,face_pairs[2*i]).centroid;
			} else if (cell_minus<0) { // An inter-partition ghost cell
				Jac[i]-=grid[gid].ghost[-cell_minus-1].centroid;
			} else {
				Jac[i]-=grid[gid].cell[cell_minus].centroid;
			}

		}
		
		// The Jacobian is calculated, now need to invert it
		
		// Determinant
		det=Jac[0][0]*(Jac[1][1]*Jac[2][2]-Jac[1][2]*Jac[2][1])
		   -Jac[0][1]*(Jac[1][0]*Jac[2][2]-Jac[1][2]*Jac[2][0])
		   +Jac[0][2]*(Jac[1][0]*Jac[2][1]-Jac[1][1]*Jac[2][0]);

		Jac_invT[0][0]=(Jac[1][1]*Jac[2][2]-Jac[1][2]*Jac[2][1])/det;
		Jac_invT[1][0]=(Jac[0][2]*Jac[2][1]-Jac[0][1]*Jac[2][2])/det;
		Jac_invT[2][0]=(Jac[0][1]*Jac[1][2]-Jac[0][2]*Jac[1][1])/det;

		Jac_invT[0][1]=(Jac[1][2]*Jac[2][0]-Jac[1][0]*Jac[2][2])/det;
		Jac_invT[1][1]=(Jac[0][0]*Jac[2][2]-Jac[0][2]*Jac[2][0])/det;
		Jac_invT[2][1]=(Jac[0][2]*Jac[1][0]-Jac[0][0]*Jac[1][2])/det;

		Jac_invT[0][2]=(Jac[1][0]*Jac[2][1]-Jac[1][1]*Jac[2][0])/det;
		Jac_invT[1][2]=(Jac[0][1]*Jac[2][0]-Jac[0][0]*Jac[2][1])/det;
		Jac_invT[2][2]=(Jac[0][0]*Jac[1][1]-Jac[0][1]*Jac[1][0])/det;
	
		/*
		 
		 Jac:
		 | dx/dxsi  dy/dxsi  dz/dxsi  | 
		 | dx/deta  dy/deta  dz/deta  | 
		 | dx/dzeta dy/dzeta dz/dzeta | 
		 
		 Jac_invT
		 | dxsi/dx  dxsi/dy  dxsi/dz  | 
		 | deta/dx  deta/dy  deta/dz  | 
		 | dzeta/dx dzeta/dy dzeta/dz | 
		
		 
		dphi/dx=dphi/dxsi*dxsi/dx + dphi/deta*deta/dx + dphi/dzeta*dzeta/dx
		 
		dphi/dx= (cell_pairs[1]-cell_pairs[0])*Jac_invT[0][0]
		        +(cell_pairs[3]-cell_pairs[2])*Jac_invT[1][0]
				+(cell_pairs[5]-cell_pairs[4])*Jac_invT[2][0]		

		dphi/dy= (cell_pairs[1]-cell_pairs[0])*Jac_invT[0][1]
		        +(cell_pairs[3]-cell_pairs[2])*Jac_invT[1][1]
         		+(cell_pairs[5]-cell_pairs[4])*Jac_invT[2][1]	
		
		dphi/dy= (cell_pairs[1]-cell_pairs[0])*Jac_invT[0][2]
				+(cell_pairs[3]-cell_pairs[2])*Jac_invT[1][2]
				+(cell_pairs[5]-cell_pairs[4])*Jac_invT[2][2]
		 
		*/
		
		map<int,Vec3D> gradMap;
		
		gradMap.insert(pair<int,Vec3D>(cell_pairs[0],-1.*Jac_invT[0]));
		
		if (gradMap.find(cell_pairs[1])==gradMap.end())	gradMap.insert(pair<int,Vec3D>(cell_pairs[1],Jac_invT[0]));		
		else gradMap[cell_pairs[1]]+=Jac_invT[0];
		
		if (gradMap.find(cell_pairs[2])==gradMap.end())	gradMap.insert(pair<int,Vec3D>(cell_pairs[2],-1.*Jac_invT[1]));		
		else gradMap[cell_pairs[2]]-=Jac_invT[1];
		
		if (gradMap.find(cell_pairs[3])==gradMap.end())	gradMap.insert(pair<int,Vec3D>(cell_pairs[3],Jac_invT[1]));		
		else gradMap[cell_pairs[3]]+=Jac_invT[1];
		
		if (gradMap.find(cell_pairs[4])==gradMap.end())	gradMap.insert(pair<int,Vec3D>(cell_pairs[4],-1.*Jac_invT[2]));		
		else gradMap[cell_pairs[4]]-=Jac_invT[2];
		
		if (gradMap.find(cell_pairs[5])==gradMap.end())	gradMap.insert(pair<int,Vec3D>(cell_pairs[5],Jac_invT[2]));		
		else gradMap[cell_pairs[5]]+=Jac_invT[2];
		
		grid[gid].cell[c].gradMap=gradMap;
		
		gradMap.clear();
		
	} // end of cell loop
		
	stencil.clear();
	
	return;
}
