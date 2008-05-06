#include <cmath>
#include<string>
#include<mpi.h>
#include "grid.h"
#include "bc.h"

extern Grid grid;
extern BC bc;
extern int np, rank;

void RoeFlux(const double gamma, const double qL[], const double qR[], double flux[]);
double superbee(double a, double b);
double minmod(double a, double b);

void hancock_predictor(double gamma, double dt, string limiter) {

// 	double delta[5],deltaMin[5], deltaMax[5]; // These arrays has rho,u,v,w,p sequentially
// 	Vec3D faceVel, faceNormal;
// 	double a,rho, u, v, w, p, uN;
// 	double flux[5];
// 
// 	unsigned int f;
// 	// Loop through cells
// 	for (unsigned int c = 0;c < grid.cellCount;++c) {
// 		// Loop through each face of the cell to get the limited face values
// 		for (unsigned int cf=0; cf<grid.cell[c].faceCount; ++cf) {
// 			f=grid.cell[c].faces[cf];
// 			if (grid.face[f].parent==c) {
// 				faceNormal=grid.face[f].normal;
// 			} else {
// 				faceNormal=-1.*grid.face[f].normal;
// 			}
// 
// 			find_extremas(grid,c,f,deltaMin,deltaMax);
// 			if (limiter=="superbee") {
// 				for (unsigned int i=0;i<5;++i) delta[i]=superbee(deltaMax[i],deltaMin[i]);
// 			} else if (limiter=="minmod") {
// 				for (unsigned int i=0;i<5;++i) delta[i]=minmod(deltaMax[i],deltaMin[i]);
// 			}
// 
// 
// 			
// 			// Values at the face
// 			rho=grid.cell[c].rho+delta[0];
// 			u=grid.cell[c].v.comp[0]+delta[1];
// 			v=grid.cell[c].v.comp[1]+delta[2];
// 			w=grid.cell[c].v.comp[2]+delta[3];
// 			p=grid.cell[c].p+delta[4];
// 
// 			faceVel.comp[0]=u;faceVel.comp[1]=v;faceVel.comp[2]=w;
// 
// 			uN=faceVel.dot(faceNormal);
// 
// 			if (grid.face[f].bc!=-1) { // means a boundary face
// 				if (bc.region[grid.face[f].bc].type=="slip") {
// 					uN=0.;
// 					faceVel-=faceVel.dot(faceNormal) *faceNormal;
// 					u=faceVel.comp[0]; v=faceVel.comp[1]; w=faceVel.comp[2];
// 				} else if (bc.region[grid.face[f].bc].type=="noslip") {
// 					uN=0.;u=0.;v=0.;w=0.;
// 				} else if (bc.region[grid.face[f].bc].type=="outlet") {
// 					/*u=grid.cell[c].v.comp[0];
// 					v=grid.cell[c].v.comp[1];
// 					w=grid.cell[c].v.comp[2];
// 					uN=grid.cell[c].v.dot(faceNormal);*/
// 				} else if (bc.region[grid.face[f].bc].type=="inlet") {
// 					rho=bc.region[grid.face[f].bc].p;
// 					u=bc.region[grid.face[f].bc].v.comp[0];
// 					v=bc.region[grid.face[f].bc].v.comp[1];
// 					w=bc.region[grid.face[f].bc].v.comp[2];
// 					//p=bc.region[grid.face[f].bc].p;
// 					uN=bc.region[grid.face[f].bc].v.dot(faceNormal);
// 				}
// 			}
// 
// 			a=sqrt(gamma*p/rho);
// 
// 			flux[0]= (rho*uN) *grid.face[f].area;
// 			flux[1]= (uN*u+p/rho*faceNormal.comp[0]) *grid.face[f].area;
// 			flux[2]= (uN*v+p/rho*faceNormal.comp[1]) *grid.face[f].area;
// 			flux[3]= (uN*w+p/rho*faceNormal.comp[2]) *grid.face[f].area;
// 			flux[4]= (rho* (uN) *a*a+uN*p) *grid.face[f].area;
// 
// 			for (int i = 0;i < 5;++i) {
// 				grid.cell[c].flux[i] += flux[i];
// 			}
// 
// 		} // cell faces loop
// 
// 		grid.cell[c].rho -= 0.5*dt / grid.cell[c].volume * grid.cell[c].flux[0];
// 		for (int i = 0;i <= 3;++i) grid.cell[c].v.comp[i]-= 0.5*dt / grid.cell[c].volume * grid.cell[c].flux[i+1];
// 		grid.cell[c].p-= 0.5*dt / grid.cell[c].volume * grid.cell[c].flux[4];
// 		for (int i = 0;i < 5;++i) grid.cell[c].flux[i] = 0.;
// 
// 	} // cell loop

} // end function


void hancock_corrector(double gamma, string limiter) {

	double delta[5],deltaMin[5], deltaMax[5]; // These arrays has rho,u,v,w,p sequentially
	double rhoL,pL,uL,vL,wL,rhoR,pR,uR,vR,wR,uNL,uNR,vTL,vTR,wTL,wTR;
	double qL[5], qR[5], fluxNormal[5], flux[5];
	Vec3D faceVel,parent2face,neighbor2face,faceTangent1,faceTangent2;
	unsigned int parent,neighbor;

	// Calculate flux through each face
	for (unsigned int f = 0;f < grid.faceCount;++f) {
		parent=grid.face[f].parent; neighbor=grid.face[f].neighbor;
		// Take 2 nodes of the face to find a vector in the face plane
		faceTangent1= (grid.face[f].node(1)-grid.face[f].node(0));
		faceTangent1/=fabs(faceTangent1);
		// Cross the tangent vector with the normal vector to get the second tangent
		faceTangent2= (grid.face[f].normal).cross(faceTangent1);

		for (unsigned int i=0;i<5;++i) delta[i]=(grid.face[f].centroid-grid.cell[parent].centroid).dot(grid.cell[parent].limited_grad[i]);

		rhoL=grid.cell[parent].rho+delta[0];
		uL=grid.cell[parent].v.comp[0]+delta[1];
		vL=grid.cell[parent].v.comp[1]+delta[2];
		wL=grid.cell[parent].v.comp[2]+delta[3];
		pL=grid.cell[parent].p+delta[4];

		faceVel.comp[0]=uL;faceVel.comp[1]=vL;faceVel.comp[2]=wL;
		uNL=faceVel.dot(grid.face[f].normal);
		vTL=faceVel.dot(faceTangent1);
		wTL=faceVel.dot(faceTangent2);
		
		if (grid.face[f].bc>=0) { // means a real boundary face
			rhoR=rhoL; uNR=uNL; vTR=vTL; wTR=wTL; pR=pL; // outlet condition
			if (bc.region[grid.face[f].bc].type=="slip") uNR=-uNL;
			if (bc.region[grid.face[f].bc].type=="noslip") {uNR=-uNL; vTR=vTL=0.; wTR=wTL=0.;}
			if (bc.region[grid.face[f].bc].type=="inlet") {
				uNR=bc.region[grid.face[f].bc].v.dot(grid.face[f].normal);
				vTR=bc.region[grid.face[f].bc].v.dot(faceTangent1);
				wTR=bc.region[grid.face[f].bc].v.dot(faceTangent2);
				rhoR=bc.region[grid.face[f].bc].rho;
				//pR=bc.region[grid.face[f].bc].p;
			}
		} else if (grid.face[f].bc==-1) { // internal face
			for (unsigned int i=0;i<5;++i) delta[i]=(grid.face[f].centroid-grid.cell[neighbor].centroid).dot(grid.cell[neighbor].limited_grad[i]);

			rhoR=grid.cell[neighbor].rho+delta[0];
			uR=grid.cell[neighbor].v.comp[0]+delta[1];
			vR=grid.cell[neighbor].v.comp[1]+delta[2];
			wR=grid.cell[neighbor].v.comp[2]+delta[3];
			pR=grid.cell[neighbor].p+delta[4];

			faceVel.comp[0]=uR;faceVel.comp[1]=vR;faceVel.comp[2]=wR;
			uNR=faceVel.dot(grid.face[f].normal);
			vTR=faceVel.dot(faceTangent1);
			wTR=faceVel.dot(faceTangent2);
		
		} else { // partition boundary 
			int g=-1*grid.face[f].bc-2;

			for (unsigned int i=0;i<5;++i) {
				delta[i]=(grid.face[f].centroid-grid.ghost[g].centroid).dot(grid.ghost[g].limited_grad[i]);
				//cout << grid.ghost[g].limited_grad[i] << endl;
			}
			
			rhoR=grid.ghost[g].rho+delta[0];
			uR=grid.ghost[g].v.comp[0]+delta[1];
			vR=grid.ghost[g].v.comp[1]+delta[2];
			wR=grid.ghost[g].v.comp[2]+delta[3];
			pR=grid.ghost[g].p+delta[4];
			
			faceVel.comp[0]=uR;faceVel.comp[1]=vR;faceVel.comp[2]=wR;
			uNR=faceVel.dot(grid.face[f].normal);
			vTR=faceVel.dot(faceTangent1);
			wTR=faceVel.dot(faceTangent2);
		}

		qL[0]=rhoL;
		qL[1]=rhoL * uNL;
		qL[2]=rhoL * vTL;
		qL[3]=rhoL * wTL;
		qL[4]=0.5*rhoL* (uNL*uNL+vTL*vTL+wTL*wTL) + pL/ (gamma - 1.);

		qR[0] = rhoR;
		qR[1] = rhoR * uNR;
		qR[2] = rhoR * vTR;
		qR[3] = rhoR * wTR;
		qR[4] = 0.5*rhoR* (uNR*uNR+vTR*vTR+wTR*wTR) + pR/ (gamma - 1.);

		RoeFlux(gamma, qL, qR, fluxNormal);

		flux[0] = fluxNormal[0]*grid.face[f].area;
		flux[1] = (fluxNormal[1]*grid.face[f].normal.comp[0]+fluxNormal[2]*faceTangent1.comp[0]+fluxNormal[3]*faceTangent2.comp[0]) * grid.face[f].area;
		flux[2] = (fluxNormal[1]*grid.face[f].normal.comp[1]+fluxNormal[2]*faceTangent1.comp[1]+fluxNormal[3]*faceTangent2.comp[1]) * grid.face[f].area;
		flux[3] = (fluxNormal[1]*grid.face[f].normal.comp[2]+fluxNormal[2]*faceTangent1.comp[2]+fluxNormal[3]*faceTangent2.comp[2]) * grid.face[f].area;
		flux[4] = fluxNormal[4]*grid.face[f].area;

		for (int i = 0;i < 5;++i) {
			grid.cell[parent].flux[i] += flux[i];
			if (grid.face[f].bc==-1) {  // internal face
				grid.cell[neighbor].flux[i] -= flux[i];
			}
		}

	} // face loop

} // end function


