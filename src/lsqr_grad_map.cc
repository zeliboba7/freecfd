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
#include "bc.h"
using namespace std;

extern vector<Grid> grid;
extern InputFile input;
extern vector<vector<BCregion> > bc;

void lsqr_grad_map(int gid,int ci) {

	vector<int> stencil;
	vector<double> w;
	vector<Vec3D> distance;
	vector<int>::iterator sit;
	double a,b,c,d,e,f,det;
	Vec3D point;
	int neighbor;
	Vec3D weight;
	
	// Populate cell neighbors into the stencil
	/*
	for (int cf=0;cf<grid[gid].cell[ci].faceCount;++cf) {
		if (grid[gid].cellFace(ci,cf).bc<0) { // if not a boundary face
			neighbor=grid[gid].cellFace(ci,cf).neighbor;
			if (neighbor==ci) neighbor=grid[gid].cellFace(ci,cf).parent;
			stencil.push_back(neighbor);
		}
	}
	*/
	
	for (int i=0;i<grid[gid].cell[ci].neighborCells.size();++i) {
		neighbor=grid[gid].cell[ci].neighborCells[i];
		if (neighbor!=ci) stencil.push_back(neighbor);
	}
	
	// Loop the stencil and calculate the metrics
	a=0.; b=0.; c=0.; d=0.; e=0; f=0.;
	
	for (int k=0;k<stencil.size();k++) {
		if (stencil[k]>=0) point=grid[gid].cell[stencil[k]].centroid;
		else point=grid[gid].ghost[stencil[k]].centroid;
		distance.push_back(point-grid[gid].cell[ci].centroid);
 	    w.push_back(distance[k].dot(distance[k]));
		a+=w[k]*distance[k][0]*distance[k][0];
		b+=w[k]*distance[k][0]*distance[k][1];
		c+=w[k]*distance[k][0]*distance[k][2];
		d+=w[k]*distance[k][1]*distance[k][1];
		e+=w[k]*distance[k][1]*distance[k][2];
		f+=w[k]*distance[k][2]*distance[k][2];
	}
	
	// Handle the coplanar stencil happening in 2D or 1D runs
	// Add mirror cells to the symmetry faces of the target cell (ci)
	//if (grid[gid].dimension<3) {
	for (int cf=0;cf<grid[gid].cell[ci].faceCount;++cf) {
		int bcno=grid[gid].cellFace(ci,cf).bc;
		if (grid[gid].face[grid[gid].cell[ci].faces[cf]].symmetry) {
			stencil.push_back(ci);
			point=grid[gid].cell[ci].centroid+(grid[gid].cellFace(ci,cf).centroid-grid[gid].cell[ci].centroid).dot(grid[gid].cellFace(ci,cf).normal)*grid[gid].cellFace(ci,cf).normal;
			distance.push_back(point-grid[gid].cell[ci].centroid);
			int k=distance.size()-1;
			w.push_back(distance[k].dot(distance[k]));
			a+=w[k]*distance[k][0]*distance[k][0];
			b+=w[k]*distance[k][0]*distance[k][1];
			c+=w[k]*distance[k][0]*distance[k][2];
			d+=w[k]*distance[k][1]*distance[k][1];
			e+=w[k]*distance[k][1]*distance[k][2];
			f+=w[k]*distance[k][2]*distance[k][2];
		}
	}
	
	det=a*d*f-a*e*e-b*b*f-c*c*d+2.*b*c*e;
	
	weight=0.;
	grid[gid].cell[ci].gradMap.clear();
	for (int k=0;k<stencil.size();k++) {
		weight[0]=w[k]*(distance[k][0]*(d*f-e*e)+distance[k][1]*(c*e-b*f)+distance[k][2]*(b*e-c*d));
		weight[1]=w[k]*(distance[k][0]*(c*e-b*f)+distance[k][1]*(a*f-c*c)+distance[k][2]*(b*c-a*e));
		weight[2]=w[k]*(distance[k][0]*(b*e-c*d)+distance[k][1]*(b*c-a*e)+distance[k][2]*(a*d-b*b));		
		weight/=det;
		grid[gid].cell[ci].gradMap.insert(pair<int,Vec3D>(stencil[k],weight));
		if (grid[gid].cell[ci].gradMap.find(ci)==grid[gid].cell[ci].gradMap.end()) grid[gid].cell[ci].gradMap.insert(pair<int,Vec3D>(ci,-1.*weight));	
		else grid[gid].cell[ci].gradMap[ci]-=weight;
		
	}

	return;
}
