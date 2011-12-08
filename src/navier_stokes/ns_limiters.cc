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
#define sign(x) (( x > 0 ) - ( x < 0 ))

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

double minmod(double a1,double a2) {
	if (a1*a2<0.) return 0.;
	else if (a1>0.) return min(a1,a2);
	else return max(a1,a2);
}

void NavierStokes::minmod_limiter(void) {
	
	int neighbor,g;
	double phi[5];
	double deltaU,deltaUg;
	double small=1.e-12;
	
	for (int c=0;c<grid[gid].cellCount;++c) {
		
		for (int i=0;i<5;++i) phi[i]=1.;
		
		// Repeat the loop to calculate the limiter for each face
		for (int cf=0;cf<grid[gid].cell[c].faces.size();++cf) {
			Vec3D cell2face=grid[gid].cellFace(c,cf).centroid-grid[gid].cell[c].centroid;
			if (c==grid[gid].cellFace(c,cf).parent) {
				neighbor=grid[gid].cellFace(c,cf).neighbor;
			} else {
				neighbor=grid[gid].cellFace(c,cf).parent;
			}
			for (int var=0;var<5;++var) {
				if (var==0) { 
					deltaU=p.cell(neighbor)-p.cell(c);
					deltaUg=gradp.cell(c).dot(cell2face);
				} else if (var==1) { 
					deltaU=V.cell(neighbor)[0]-V.cell(c)[0];
					deltaUg=gradu.cell(c).dot(cell2face);
				} else if (var==2) { 
					deltaU=V.cell(neighbor)[1]-V.cell(c)[1];
					deltaUg=gradv.cell(c).dot(cell2face);
				} else if (var==3) { 
					deltaU=V.cell(neighbor)[2]-V.cell(c)[2];
					deltaUg=gradw.cell(c).dot(cell2face);
				} else if (var==4) { 
					deltaU=T.cell(neighbor)-T.cell(c);
					deltaUg=gradT.cell(c).dot(cell2face);
				}

				deltaU=minmod(deltaU,deltaUg);
				phi[var]=min(phi[var],(deltaU+sign(deltaU)*small)/(deltaUg+sign(deltaUg)*small));
			}
		} // end face loop		

		for (int var=0;var<5;++var) limiter[var].cell(c)=phi[var];
		
	} // end cell loop
	
	return;
}

void NavierStokes::barth_jespersen_limiter(void) {
	
	int neighbor,g;
	double phi[5];
	double deltaP,deltaM;
	double umax[5],umin[5];

	
	for (int c=0;c<grid[gid].cellCount;++c) {
		
		for (int i=0;i<5;++i) phi[i]=1.;
		umax[0]=umin[0]=p.cell(c);
		umax[1]=umin[1]=V.cell(c)[0];
		umax[2]=umin[2]=V.cell(c)[1];
		umax[3]=umin[3]=V.cell(c)[2];
		umax[4]=umin[4]=T.cell(c);
		
		// First loop through face neighbors to find the max and min values
		for (int cf=0;cf<grid[gid].cell[c].faces.size();++cf) {
			if (c==grid[gid].cellFace(c,cf).parent) {
				neighbor=grid[gid].cellFace(c,cf).neighbor;
			} else {
				neighbor=grid[gid].cellFace(c,cf).parent;
			}
			umax[0]=max(umax[0],p.cell(neighbor));
			umax[1]=max(umax[1],V.cell(neighbor)[0]);
			umax[2]=max(umax[2],V.cell(neighbor)[1]);
			umax[3]=max(umax[3],V.cell(neighbor)[2]);
			umax[4]=max(umax[4],T.cell(neighbor));

			umin[0]=min(umin[0],p.cell(neighbor));
			umin[1]=min(umin[1],V.cell(neighbor)[0]);
			umin[2]=min(umin[2],V.cell(neighbor)[1]);
			umin[3]=min(umin[3],V.cell(neighbor)[2]);
			umin[4]=min(umin[4],T.cell(neighbor));
		}

		// Repeat the loop to calculate the limiter for each face
		for (int cf=0;cf<grid[gid].cell[c].faces.size();++cf) {
			Vec3D cell2face=grid[gid].cellFace(c,cf).centroid-grid[gid].cell[c].centroid;
			for (int var=0;var<5;++var) {
				if (var==0) { deltaM=gradp.cell(c).dot(cell2face); deltaP=p.cell(c); } 
				if (var==1) { deltaM=gradu.cell(c).dot(cell2face); deltaP=V.cell(c)[0]; }
				if (var==2) { deltaM=gradv.cell(c).dot(cell2face); deltaP=V.cell(c)[1]; }
				if (var==3) { deltaM=gradw.cell(c).dot(cell2face); deltaP=V.cell(c)[2]; }
				if (var==4) { deltaM=gradT.cell(c).dot(cell2face); deltaP=T.cell(c); }

				if (deltaM>0.) {
					deltaP=umax[var]-deltaP;
				} else if (deltaM<0.) {
					deltaP=umin[var]-deltaP;
				} else {
					deltaP=1.; deltaM=1.;
				}

				phi[var]=min(phi[var],deltaP/deltaM);
			}	

		} // end face loop		
		
		for (int var=0;var<5;++var) limiter[var].cell(c)=phi[var];
	
	} // end cell loop
	
	return;
}

void NavierStokes::venkatakrishnan_limiter(void) {
	
	int neighbor,g;
	double phi[5];
	double deltaP,deltaP2,deltaM,deltaM2,eps,eps2;
	double K=limiter_threshold;
	double umax[5],umin[5];
	//double pref,uref,tref
	double Lref;
	double deltaRef[5];

	find_min_max();

	Lref=grid[gid].lengthScale;
	for (int i=0;i<5;++i) deltaRef[i]=fabs(qmax[i]-qmin[i]);

	for (int c=0;c<grid[gid].cellCount;++c) {
		
		for (int i=0;i<5;++i) {
			phi[i]=1.;
			umax[0]=umin[0]=p.cell(c);
			umax[1]=umin[1]=V.cell(c)[0];
			umax[2]=umin[2]=V.cell(c)[1];
			umax[3]=umin[3]=V.cell(c)[2];
			umax[4]=umin[4]=T.cell(c);
		}
		
		// First loop through face neighbors to find the max and min values
		for (int cf=0;cf<grid[gid].cell[c].faces.size();++cf) {
			if (c==grid[gid].cellFace(c,cf).parent) {
				neighbor=grid[gid].cellFace(c,cf).neighbor;
			} else {
				neighbor=grid[gid].cellFace(c,cf).parent;
			}
			umax[0]=max(umax[0],p.cell(neighbor));
			umax[1]=max(umax[1],V.cell(neighbor)[0]);
			umax[2]=max(umax[2],V.cell(neighbor)[1]);
			umax[3]=max(umax[3],V.cell(neighbor)[2]);
			umax[4]=max(umax[4],T.cell(neighbor));

			umin[0]=min(umin[0],p.cell(neighbor));
			umin[1]=min(umin[1],V.cell(neighbor)[0]);
			umin[2]=min(umin[2],V.cell(neighbor)[1]);
			umin[3]=min(umin[3],V.cell(neighbor)[2]);
			umin[4]=min(umin[4],T.cell(neighbor));
		}

		// Repeat the loop to calculate the limiter for each face
		for (int cf=0;cf<grid[gid].cell[c].faces.size();++cf) {
			Vec3D cell2face=grid[gid].cellFace(c,cf).centroid-grid[gid].cell[c].centroid;
			for (int var=0;var<5;++var) {
				if (var==0) { deltaM=gradp.cell(c).dot(cell2face); deltaP=p.cell(c); } 
				if (var==1) { deltaM=gradu.cell(c).dot(cell2face); deltaP=V.cell(c)[0]; }
				if (var==2) { deltaM=gradv.cell(c).dot(cell2face); deltaP=V.cell(c)[1]; }
				if (var==3) { deltaM=gradw.cell(c).dot(cell2face); deltaP=V.cell(c)[2]; }
				if (var==4) { deltaM=gradT.cell(c).dot(cell2face); deltaP=T.cell(c); }
				eps=deltaRef[var]/Lref;

				if (deltaM>0.) {
					deltaP=umax[var]-deltaP;
				//	deltaM+=1.e-12;
				} else if (deltaM<0.) {
					deltaP=umin[var]-deltaP;
				//	deltaM-=1.e-12;
				} 
		
				deltaP2=deltaP*deltaP;
				deltaM2=deltaM*deltaM;
				eps2=pow(1.e-6*K*eps,2);

				if (deltaM!=0.) {
					phi[var]=min(phi[var],((deltaP2+eps2)+2.*deltaM*deltaP)/(deltaP2+2.*deltaM2+deltaM*deltaP+eps2));
				//	phi[var]=min(
				//		phi[var],
				//		(1./deltaM)*(((deltaP2+eps2)*deltaM+2.*deltaM2*deltaP)/(deltaP2+2.*deltaM2+deltaM*deltaP+eps2))
				//		);
				}

			}	

		} // end face loop		

		for (int var=0;var<5;++var) limiter[var].cell(c)=phi[var];
	
	} // end cell loop
	
	return;
} 
