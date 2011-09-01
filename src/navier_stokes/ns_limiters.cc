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
	} else if (limiter_function==MINMOD) {
		minmod_limiter();
	} 
	return;
}

void NavierStokes::barth_jespersen_limiter(void) {
	
	int neighbor,g;
	double phi[5];
	double delta;
	double max_val[5], min_val[5];
	
	for (int c=0;c<grid[gid].cellCount;++c) {
		
		for (int i=0;i<5;++i) phi[i]=1.;
		
		max_val[0]=min_val[0]=p.cell(c);
		max_val[1]=min_val[1]=V.cell(c)[0];
		max_val[2]=min_val[2]=V.cell(c)[1];
		max_val[3]=min_val[3]=V.cell(c)[2];
		max_val[4]=min_val[4]=T.cell(c);
		
		// First loop through face neighbors to find max values
		for (int cf=0;cf<grid[gid].cell[c].faceCount;++cf) {
			if (grid[gid].cellFace(c,cf).bc<0) { // If not a bouddary face
				c==grid[gid].cellFace(c,cf).parent ? neighbor=grid[gid].cellFace(c,cf).neighbor : neighbor=grid[gid].cellFace(c,cf).parent;
				if (neighbor>=0) { // real cell
					max_val[0]=max(max_val[0],p.cell(neighbor));
					max_val[1]=max(max_val[1],V.cell(neighbor)[0]);
					max_val[2]=max(max_val[2],V.cell(neighbor)[1]);
					max_val[3]=max(max_val[3],V.cell(neighbor)[2]);
					max_val[4]=max(max_val[4],T.cell(neighbor));
					
					min_val[0]=min(min_val[0],p.cell(neighbor));
					min_val[1]=min(min_val[1],V.cell(neighbor)[0]);
					min_val[2]=min(min_val[2],V.cell(neighbor)[1]);
					min_val[3]=min(min_val[3],V.cell(neighbor)[2]);
					min_val[4]=min(min_val[4],T.cell(neighbor));
				} else { // neighbor cell is a ghost (partition interface)
					neighbor=-1*neighbor-1;
					max_val[0]=max(max_val[0],p.ghost(neighbor));
					max_val[1]=max(max_val[1],V.ghost(neighbor)[0]);
					max_val[2]=max(max_val[2],V.ghost(neighbor)[1]);
					max_val[3]=max(max_val[3],V.ghost(neighbor)[2]);
					max_val[4]=max(max_val[4],T.ghost(neighbor));
					
					min_val[0]=min(min_val[0],p.ghost(neighbor));
					min_val[1]=min(min_val[1],V.ghost(neighbor)[0]);
					min_val[2]=min(min_val[2],V.ghost(neighbor)[1]);
					min_val[3]=min(min_val[3],V.ghost(neighbor)[2]);
					min_val[4]=min(min_val[4],T.ghost(neighbor));
				}
			}
		} // end face loop
				
		// Second loop through face neigbors to calculate min limiter
		for (int cf=0;cf<grid[gid].cell[c].faceCount;++cf) {
			if (grid[gid].cellFace(c,cf).bc<0) { // If not a boundary face
				for (int var=0;var<5;++var) {
					if (var==0) {
						delta=gradp.cell(c).dot(grid[gid].cellFace(c,cf).centroid-grid[gid].cell[c].centroid);
						if (delta>0.) phi[var]=min(phi[var],min(1.,(max_val[var]-p.cell(c))/delta));
						else if (delta<0.) phi[var]=min(phi[var],min(1.,(min_val[var]-p.cell(c))/delta));
					} else if (var==1) {
						delta=gradu.cell(c).dot(grid[gid].cellFace(c,cf).centroid-grid[gid].cell[c].centroid);
						if (delta>0.) phi[var]=min(phi[var],min(1.,(max_val[var]-V.cell(c)[0])/delta));
						else if (delta<0.) phi[var]=min(phi[var],min(1.,(min_val[var]-V.cell(c)[0])/delta));
					} else if (var==2) {
						delta=gradv.cell(c).dot(grid[gid].cellFace(c,cf).centroid-grid[gid].cell[c].centroid);
						if (delta>0.) phi[var]=min(phi[var],min(1.,(max_val[var]-V.cell(c)[1])/delta));
						else if (delta<0.) phi[var]=min(phi[var],min(1.,(min_val[var]-V.cell(c)[1])/delta));
					} else if (var==3) {
						delta=gradw.cell(c).dot(grid[gid].cellFace(c,cf).centroid-grid[gid].cell[c].centroid);
						if (delta>0.) phi[var]=min(phi[var],min(1.,(max_val[var]-V.cell(c)[2])/delta));
						else if (delta<0.) phi[var]=min(phi[var],min(1.,(min_val[var]-V.cell(c)[2])/delta));
					} else if (var==4) {
						delta=gradT.cell(c).dot(grid[gid].cellFace(c,cf).centroid-grid[gid].cell[c].centroid);
						if (delta>0.) phi[var]=min(phi[var],min(1.,(max_val[var]-T.cell(c))/delta));
						else if (delta<0.) phi[var]=min(phi[var],min(1.,(min_val[var]-T.cell(c))/delta));
					}
				}
			}
		}

		for (int var=0;var<5;++var) limiter[var].cell(c)=phi[var];
		
	} // end cell loop
	
	for (int var=0;var<5;++var) limiter[var].mpi_update();	
	
	return;
}

void NavierStokes::venkatakrishnan_limiter(void) {
	
	int neighbor,g;
	double phi[5];
	double deltaP,deltaP2,deltaM,deltaM2,eps,eps2;
	double K=limiter_threshold;
	double max_delta[5], min_delta[5];
	
	double uref,pref,tref;
	
	Vec3D maxGrad[5],minGrad[5];
	for (int c=0;c<grid[gid].cellCount;++c) {
		
		for (int i=0;i<5;++i) phi[i]=1.;
		
		for (int i=0;i<5;++i) {
			max_delta[i]=0.;
			min_delta[i]=0.;
		}
		
		// First loop through face neighbors to find max values
		for (int cf=0;cf<grid[gid].cell[c].faceCount;++cf) {
			if (grid[gid].cellFace(c,cf).bc<0) { // If not a boundary face
				c==grid[gid].cellFace(c,cf).parent ? neighbor=grid[gid].cellFace(c,cf).neighbor : neighbor=grid[gid].cellFace(c,cf).parent;
				Vec3D neighbor2face=grid[gid].cellFace(c,cf).centroid-grid[gid].cell[neighbor].centroid;
				if (neighbor>=0) { // real cell
					max_delta[0]=max(max_delta[0],p.cell(neighbor)+gradp.cell(neighbor).dot(neighbor2face)-p.cell(c));
					max_delta[1]=max(max_delta[1],V.cell(neighbor)[0]+gradu.cell(neighbor).dot(neighbor2face)-V.cell(c)[0]);
					max_delta[2]=max(max_delta[2],V.cell(neighbor)[1]+gradv.cell(neighbor).dot(neighbor2face)-V.cell(c)[1]);	
					max_delta[3]=max(max_delta[3],V.cell(neighbor)[2]+gradw.cell(neighbor).dot(neighbor2face)-V.cell(c)[2]);
					max_delta[4]=max(max_delta[4],T.cell(neighbor)-T.cell(c));
					
					min_delta[0]=min(min_delta[0],p.cell(neighbor)+gradp.cell(neighbor).dot(neighbor2face)-p.cell(c));
					min_delta[1]=min(min_delta[1],V.cell(neighbor)[0]+gradu.cell(neighbor).dot(neighbor2face)-V.cell(c)[0]);
					min_delta[2]=min(min_delta[2],V.cell(neighbor)[1]+gradv.cell(neighbor).dot(neighbor2face)-V.cell(c)[1]);	
					min_delta[3]=min(min_delta[3],V.cell(neighbor)[2]+gradw.cell(neighbor).dot(neighbor2face)-V.cell(c)[2]);
					min_delta[4]=min(min_delta[4],T.cell(neighbor)-T.cell(c));
					
				} else { // neighbor cell is a ghost (partition interface)
					neighbor=-1*neighbor-1;
					max_delta[0]=max(max_delta[0],p.ghost(neighbor)-p.cell(c));
					max_delta[1]=max(max_delta[1],V.ghost(neighbor)[0]-V.cell(c)[0]);
					max_delta[2]=max(max_delta[2],V.ghost(neighbor)[1]-V.cell(c)[1]);	
					max_delta[3]=max(max_delta[3],V.ghost(neighbor)[2]-V.cell(c)[2]);
					max_delta[4]=max(max_delta[4],T.ghost(neighbor)-T.cell(c));
					
					min_delta[0]=min(min_delta[0],p.ghost(neighbor)-p.cell(c));
					min_delta[1]=min(min_delta[1],V.ghost(neighbor)[0]-V.cell(c)[0]);
					min_delta[2]=min(min_delta[2],V.ghost(neighbor)[1]-V.cell(c)[1]);	
					min_delta[3]=min(min_delta[3],V.ghost(neighbor)[2]-V.cell(c)[2]);
					min_delta[4]=min(min_delta[4],T.ghost(neighbor)-T.cell(c));
					
				}
			}
		} // end face loop		
		
		pref=p.cell(c)+material.Pref;
		tref=T.cell(c)+material.Tref;
		uref=material.a(p.cell(c),T.cell(c));
		
		// Second loop through face neigbors to calculate min limiter		
		for (int cf=0;cf<grid[gid].cell[c].faceCount;++cf) {
			if (grid[gid].cellFace(c,cf).bc<0) { // If not a boundary face
				Vec3D cell2face=grid[gid].cellFace(c,cf).centroid-grid[gid].cell[c].centroid;
				eps=pow(K*fabs(cell2face),3);
				//eps=pow(1.e-3/fabs(grid[gid].cellFace(c,cf).centroid-grid[gid].cell[c].centroid),2);
				for (int var=0;var<5;++var) {
					deltaP=1.;
					if (var==0) {
						eps2=eps*pref*pref;
						deltaM=gradp.cell(c).dot(cell2face);
					} else if (var==1) {
						eps2=eps*uref*uref;
						deltaM=gradu.cell(c).dot(cell2face);
					} else if (var==2) {
						eps2=eps*uref*uref;
						deltaM=gradv.cell(c).dot(cell2face);
					} else if (var==3) {
						eps2=eps*uref*uref;
						deltaM=gradw.cell(c).dot(cell2face);
					} else if (var==4) {
						eps2=eps*tref*tref;
						deltaM=gradT.cell(c).dot(cell2face);
					}
					
					if (deltaM>0.) deltaP=max_delta[var];
					else if (deltaM<0.) deltaP=min_delta[var];

					//deltaM>0. ? deltaM+=1.e-12 : deltaM-=1.e-12;
					deltaM2=deltaM*deltaM;
					deltaP2=deltaP*deltaP;
					
					//phi[var]=min(phi[var],(1./deltaM)*((deltaP2+eps2)*deltaM+2.*deltaM2*deltaP)/(deltaP2+2.*deltaM2+deltaM*deltaP+eps2)); 
					phi[var]=min(phi[var],((deltaP2+eps2)+2.*deltaM*deltaP)/(deltaP2+2.*deltaM2+deltaM*deltaP+eps2)); 
					
				} // var loop
			} // if not a boundary face
		} // end face loop
		
		for (int var=0;var<5;++var) limiter[var].cell(c)=phi[var];
		
	} // end cell loop
	for (int var=0;var<5;++var) limiter[var].mpi_update();
	
	return;
} 

void NavierStokes::minmod_limiter(void) {
	
	int neighbor,g;
	double delta1[5],delta2[5];
	double phi[5];
	double eps=1.e-8;
	
	for (int c=0;c<grid[gid].cellCount;++c) {
		
		for (int i=0;i<5;++i) limiter[i].cell(c)=1.;
		
		// First loop through face neighbors to find max values
		for (int cf=0;cf<grid[gid].cell[c].faceCount;++cf) {
			for (int i=0;i<5;++i) {
				delta1[i]=0.;
				delta2[i]=0.;
				phi[i]=0.;
			}
			if (grid[gid].cellFace(c,cf).bc<0) { // If not a boundary face
				c==grid[gid].cellFace(c,cf).parent ? neighbor=grid[gid].cellFace(c,cf).neighbor : neighbor=grid[gid].cellFace(c,cf).parent;
				Vec3D cell2face=grid[gid].cellFace(c,cf).centroid-grid[gid].cell[c].centroid;
				Vec3D neighbor2face=grid[gid].cellFace(c,cf).centroid-grid[gid].cell[neighbor].centroid;

				delta1[0]=gradp.cell(c).dot(cell2face);
				delta1[1]=gradu.cell(c).dot(cell2face);
				delta1[2]=gradv.cell(c).dot(cell2face);
				delta1[3]=gradw.cell(c).dot(cell2face);
				delta1[4]=gradT.cell(c).dot(cell2face);
				if (neighbor>=0) { // real cell
					delta2[0]=p.cell(neighbor)+gradp.cell(neighbor).dot(neighbor2face)-p.cell(c);
					delta2[1]=V.cell(neighbor)[0]+gradu.cell(neighbor).dot(neighbor2face)-V.cell(c)[0];
					delta2[2]=V.cell(neighbor)[1]+gradv.cell(neighbor).dot(neighbor2face)-V.cell(c)[1];
					delta2[3]=V.cell(neighbor)[2]+gradw.cell(neighbor).dot(neighbor2face)-V.cell(c)[2];
					delta2[4]=T.cell(neighbor)+gradT.cell(neighbor).dot(neighbor2face)-T.cell(c);
				} else { // neighbor cell is a ghost (partition interface)
					neighbor=-1*neighbor-1;
					delta2[0]=p.ghost(neighbor)+gradp.ghost(neighbor).dot(neighbor2face)-p.cell(c);
					delta2[1]=V.ghost(neighbor)[0]+gradu.ghost(neighbor).dot(neighbor2face)-V.cell(c)[0];
					delta2[2]=V.ghost(neighbor)[1]+gradv.ghost(neighbor).dot(neighbor2face)-V.cell(c)[1];
					delta2[3]=V.ghost(neighbor)[2]+gradw.ghost(neighbor).dot(neighbor2face)-V.cell(c)[2];
					delta2[4]=T.ghost(neighbor)+gradT.ghost(neighbor).dot(neighbor2face)-T.cell(c);
				}
			} else { // If a boundary face
				
			}
			
			for (int i=0;i<5;++i) {
				if (delta1[i]*delta2[i]<0.) { // They are of opposite signs
					phi[i]=0.;
				} else {
	//				Put minmod here;
	//				phi[i]=(min(fabs(delta1[i]),fabs(delta2[i]))+eps)/(fabs(delta1[i])+eps);
				}
				limiter[i].cell(c)=min(limiter[i].cell(c),phi[i]);
			}

		} // end face loop		

	} // end cell loop
	for (int var=0;var<5;++var) limiter[var].mpi_update();
	
	return;
} 
