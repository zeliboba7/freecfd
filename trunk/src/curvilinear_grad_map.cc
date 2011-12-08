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
using namespace std;

extern vector<Grid> grid;
	
void curvilinear_grad_map(int gid,int c) {

	set<int> stencil;
	set<int>::iterator sit,sit1,sit2,sit3,sit4,sit5,sit6;
	Vec3D my_normal,other_normal,cell2face;
	int opposite_face;
	double product,min_product,max_product;
	vector<int> face_pairs (6,-1);
	vector<int> cell_pairs (6,-1);
	vector<int> cell_pairs_temp (6,-1);
	vector<Vec3D> Jac,Jac_temp; // Changes in othogonal system and curvilinear system
	Jac.resize(3),Jac_temp.resize(3);
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
	int counter=0;
			
	for (int cf=0;cf<grid[gid].cell[c].faces.size();++cf) stencil.insert(cf); // Note that these are not actual face indices
	
	if (grid[gid].cell[c].nodes.size()==8) { //If hexa cell
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
		
		// Now we know the three directions
		for (int i=0;i<3;++i) {
			cell_plus=grid[gid].cellFace(c,face_pairs[2*i+1]).neighbor;
			if (cell_plus==c) cell_plus=grid[gid].cellFace(c,face_pairs[2*i+1]).parent;
			cell_pairs[2*i+1]=cell_plus;
			Jac[i]=grid[gid].cell[cell_plus].centroid;
					
			cell_minus=grid[gid].cellFace(c,face_pairs[2*i]).neighbor;
			if (cell_minus==c) cell_minus=grid[gid].cellFace(c,face_pairs[2*i]).parent;
			cell_pairs[2*i]=cell_minus;
			Jac[i]-=grid[gid].cell[cell_minus].centroid;
		}
		
	} else {

		int counter=0;
		int max_counter;
		double max_det=0.;
		
		stencil.clear();
		for (int cf=0;cf<grid[gid].cell[c].faces.size();++cf) stencil.insert(grid[gid].cell[c].faces[cf]); 
		
		for (sit1=stencil.begin();sit1!=stencil.end();sit1++) {
			for (sit2=sit1;sit2!=stencil.end();sit2++) {

				cell_plus=grid[gid].face[*sit2].neighbor;
				if (cell_plus==c) cell_plus=grid[gid].face[*sit2].parent;
				cell_pairs_temp[1]=cell_plus;
				Jac_temp[0]=grid[gid].cell[cell_plus].centroid;
				
				cell_minus=grid[gid].face[*sit1].neighbor;
				if (cell_minus==c) cell_minus=grid[gid].face[*sit1].parent;
				cell_pairs_temp[0]=cell_minus;
				Jac_temp[0]-=grid[gid].cell[cell_minus].centroid;
				
				for (sit3=sit1;sit3!=stencil.end();sit3++) {
					for (sit4=sit3;sit4!=stencil.end();sit4++) {
						
						cell_plus=grid[gid].face[*sit4].neighbor;
						if (cell_plus==c) cell_plus=grid[gid].face[*sit4].parent;
						cell_pairs_temp[3]=cell_plus;
						Jac_temp[1]=grid[gid].cell[cell_plus].centroid;
						
						cell_minus=grid[gid].face[*sit3].neighbor;
						if (cell_minus==c) cell_minus=grid[gid].face[*sit3].parent;
						cell_pairs_temp[2]=cell_minus;
						Jac_temp[1]-=grid[gid].cell[cell_minus].centroid;
						
						for (sit5=sit2;sit5!=stencil.end();sit5++) {
							for (sit6=sit5;sit6!=stencil.end();sit6++) {
								
								cell_plus=grid[gid].face[*sit6].neighbor;
								if (cell_plus==c) cell_plus=grid[gid].face[*sit6].parent;
								cell_pairs_temp[5]=cell_plus;
								Jac_temp[2]=grid[gid].cell[cell_plus].centroid;
								
								cell_minus=grid[gid].face[*sit5].neighbor;
								if (cell_minus==c) cell_minus=grid[gid].face[*sit5].parent;
								cell_pairs_temp[4]=cell_minus;
								Jac_temp[2]-=grid[gid].cell[cell_minus].centroid;
								
								counter++;
								det=Jac_temp[0].dot(Jac_temp[1].cross(Jac_temp[2]));
								if (det>max_det) {
									max_counter=counter;
									max_det=det;
									cell_pairs=cell_pairs_temp;
									Jac=Jac_temp;
								}
							}
						}
					}
				}
			}
		}
	}
		
	
	stencil.clear();
	
	det=Jac[0].dot(Jac[1].cross(Jac[2]));

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
	stencil.clear();
	
	return;
}
