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
#include "ns.h"

template<typename T>
inline T signum(T n)
{
	if (n < 0) return -1;
	if (n > 0) return 1;
	return 0;
}

void NavierStokes::calc_limiter(void) {
	if(limiter_function==NONE || order==FIRST) {
		// Do nothing
	} else if (limiter_function==VK) {
		venkatakrishnan_limiter();
	} else if (limiter_function==BJ) {
		barth_jespersen_limiter();
	}
	return;
}

double limit(double phi,double delta,double delta_p) {
	if (delta_p==0.) return phi;
	if (delta*delta_p<0.) return 0.;
	
	return min(phi,delta*delta/(delta_p*delta_p));
	/*
	 else if (delta_p>0) {
	 if (delta_p>delta) return min(phi,delta/delta_p);
	 else return phi;
	 } else {
	 if (delta_p<delta) return min(phi,delta/delta_p);
	 else return phi;
	 }
	 */
}

void NavierStokes::barth_jespersen_limiter(void) {
	
	int neighbor,g;
	double phi[5];
	double delta,delta_p;
	
	for (int c=0;c<grid[gid].cellCount;++c) {
		
		for (int i=0;i<5;++i) {
			phi[i]=1.;
		}
		
		// First loop through face neighbors to find max values
		for (int cf=0;cf<grid[gid].cell[c].faceCount;++cf) {
			if (grid[gid].cellFace(c,cf).bc<0) { // If not a bouddary face
				c==grid[gid].cellFace(c,cf).parent ? neighbor=grid[gid].cellFace(c,cf).neighbor : neighbor=grid[gid].cellFace(c,cf).parent;
				if (neighbor>=0) { // real cell
					delta=p.cell(neighbor)-p.cell(c);
					delta_p=gradp.cell(c).dot(grid[gid].cellFace(c,cf).centroid-grid[gid].cell[c].centroid);
					phi[0]=limit(phi[0],delta,delta_p);
					
					delta=V.cell(neighbor)[0]-V.cell(c)[0];
					delta_p=gradu.cell(c).dot(grid[gid].cellFace(c,cf).centroid-grid[gid].cell[c].centroid);
					phi[1]=limit(phi[1],delta,delta_p);
					
					delta=V.cell(neighbor)[1]-V.cell(c)[1];
					delta_p=gradv.cell(c).dot(grid[gid].cellFace(c,cf).centroid-grid[gid].cell[c].centroid);
					phi[2]=limit(phi[2],delta,delta_p);
					
					delta=V.cell(neighbor)[2]-V.cell(c)[2];
					delta_p=gradw.cell(c).dot(grid[gid].cellFace(c,cf).centroid-grid[gid].cell[c].centroid);
					phi[3]=limit(phi[3],delta,delta_p);

					delta=T.cell(neighbor)-T.cell(c);
					delta_p=gradT.cell(c).dot(grid[gid].cellFace(c,cf).centroid-grid[gid].cell[c].centroid);
					phi[4]=limit(phi[4],delta,delta_p);
				} else { // neighbor cell is a ghost (partition interface)
					neighbor=-1*neighbor-1;
					delta=p.ghost(neighbor)-p.cell(c);
					delta_p=gradp.cell(c).dot(grid[gid].cellFace(c,cf).centroid-grid[gid].cell[c].centroid);
					phi[0]=limit(phi[0],delta,delta_p);
					
					delta=V.ghost(neighbor)[0]-V.cell(c)[0];
					delta_p=gradu.cell(c).dot(grid[gid].cellFace(c,cf).centroid-grid[gid].cell[c].centroid);
					phi[1]=limit(phi[1],delta,delta_p);
					
					delta=V.ghost(neighbor)[1]-V.cell(c)[1];
					delta_p=gradv.cell(c).dot(grid[gid].cellFace(c,cf).centroid-grid[gid].cell[c].centroid);
					phi[2]=limit(phi[2],delta,delta_p);
					
					delta=V.ghost(neighbor)[2]-V.cell(c)[2];
					delta_p=gradw.cell(c).dot(grid[gid].cellFace(c,cf).centroid-grid[gid].cell[c].centroid);
					phi[3]=limit(phi[3],delta,delta_p);
					
					delta=T.ghost(neighbor)-T.cell(c);
					delta_p=gradT.cell(c).dot(grid[gid].cellFace(c,cf).centroid-grid[gid].cell[c].centroid);
					phi[4]=limit(phi[4],delta,delta_p);
				}
			}
		} // end face loop
		for (int var=0;var<5;++var) limiter[var].cell(c)=phi[var];
		
	} // end cell loop
	
	for (int var=0;var<5;++var) limiter[var].mpi_update();	
	
	return;
}

void NavierStokes::venkatakrishnan_limiter(void) {
	
	int neighbor,g;
	double phi;
	double deltaP,deltaP2,deltaM,deltaM2,eps,eps2;
	double K=limiter_threshold;
	double max_delta[3], min_delta[3];
	
	double uref,pref,tref;
	
	Vec3D normal;
	Vec3D maxGrad[3],minGrad[3];
	for (int c=0;c<grid[gid].cellCount;++c) {
		
		for (int i=0;i<3;++i) {
			phi=1.;
			max_delta[i]=0.;
			min_delta[i]=0.;
		}
		
		// First loop through face neighbors to find max values
		for (int cf=0;cf<grid[gid].cell[c].faceCount;++cf) {
			if (grid[gid].cellFace(c,cf).bc<0) { // If not a boundary face
				if (c==grid[gid].cellFace(c,cf).parent) {
					neighbor=grid[gid].cellFace(c,cf).neighbor;
					normal=grid[gid].cellFace(c,cf).normal;
				} else {
					neighbor=grid[gid].cellFace(c,cf).parent;
					normal=-1.*grid[gid].cellFace(c,cf).normal;
				}
				if (neighbor>=0) { // real cell
					max_delta[0]=max(max_delta[0],p.cell(neighbor)-p.cell(c));
					max_delta[1]=max(max_delta[1],V.cell(neighbor).dot(normal)-V.cell(c).dot(normal));
					max_delta[2]=max(max_delta[2],T.cell(neighbor)-T.cell(c));
					
					min_delta[0]=min(min_delta[0],p.cell(neighbor)-p.cell(c));
					max_delta[1]=min(min_delta[1],V.cell(neighbor).dot(normal)-V.cell(c).dot(normal));
					min_delta[2]=min(min_delta[2],T.cell(neighbor)-T.cell(c));
					
				} else { // neighbor cell is a ghost (partition interface)
					neighbor=-1*neighbor-1;
					max_delta[0]=max(max_delta[0],p.ghost(neighbor)-p.cell(c));
					max_delta[1]=max(max_delta[1],V.ghost(neighbor).dot(normal)-V.cell(c).dot(normal));
					max_delta[2]=max(max_delta[2],T.ghost(neighbor)-T.cell(c));
					
					min_delta[0]=min(min_delta[0],p.ghost(neighbor)-p.cell(c));
					min_delta[1]=min(min_delta[1],V.ghost(neighbor).dot(normal)-V.cell(c).dot(normal));
					min_delta[2]=min(min_delta[2],T.ghost(neighbor)-T.cell(c));
					
				}
			}
		} // end face loop		
		
		pref=p.cell(c)+material.Pref;
		tref=T.cell(c)+material.Tref;
		uref=max(fabs(V.cell(c)),material.a(p.cell(c),T.cell(c)));
		// Second loop through face neigbors to calculate min limiter		
		for (int cf=0;cf<grid[gid].cell[c].faceCount;++cf) {
			if (grid[gid].cellFace(c,cf).bc<0) { // If not a boundary face
				Vec3D cell2face=grid[gid].cellFace(c,cf).centroid-grid[gid].cell[c].centroid;
				eps=pow(K*fabs(cell2face),3);
				//eps=pow(1.e-3/fabs(grid[gid].cellFace(c,cf).centroid-grid[gid].cell[c].centroid),2);
				for (int var=0;var<3;++var) {
					if (var==0) {
						eps2=eps*pref*pref;
						deltaM=gradp.cell(c).dot(cell2face);
					} else if (var==1) {
						eps2=eps*uref*uref;
						deltaM=pow(gradu.cell(c).dot(cell2face),2);
						deltaM+=pow(gradv.cell(c).dot(cell2face),2);
						deltaM+=pow(gradw.cell(c).dot(cell2face),2);
					} else if (var==2) {
						eps2=eps*tref*tref;
						deltaM=gradT.cell(c).dot(cell2face);
					}
					
					if (deltaM>=0.) deltaP=max_delta[var];
					else if (deltaM<0.) deltaP=min_delta[var];

					deltaM2=deltaM*deltaM;
					deltaP2=deltaP*deltaP;
					
					phi=min(phi,((deltaP2+eps2)+2.*deltaM*deltaP)/(deltaP2+2.*deltaM2+deltaM*deltaP+eps2)); 
					
				} // var loop
			} // if not a boundary face
		} // end face loop
		
		for (int var=0;var<5;++var) limiter[var].cell(c)=phi;
		
	} // end cell loop
	for (int var=0;var<5;++var) limiter[var].mpi_update();
	
	return;
} 
