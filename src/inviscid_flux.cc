#include <cmath>
#include "grid.h"
#include "bc.h"

extern Grid grid;
extern BC bc;
extern int rank;
extern double Gamma;

void update(double dt);

void RoeFlux(const double qL[], const double qR[], double flux[]);

void inviscid_flux(string order, string limiter) {

	double rhoL,rhoR,uNL,uNR,vTL,vTR,wTL,wTR, pL, pR;
	Vec3D faceTangent1, faceTangent2, deltaV;
	double qL[5], qR[5], fluxNormal[5], flux[5];
	double delta[5];

	unsigned int parent, neighbor;

	deltaV=0.;
	for (int i=0; i<5; ++i) delta[i]=0.;
	
	// Calculate flux through each face
	for (unsigned int f = 0;f < grid.faceCount;++f) {

		parent=grid.face[f].parent; neighbor=grid.face[f].neighbor;
										
		// Take 2 nodes of the face to find a vector in the face plane
		faceTangent1= (grid.face[f].node(1)-grid.face[f].node(0));
		faceTangent1/=fabs(faceTangent1);
		// Cross the tangent vector with the normal vector to get the second tangent
		faceTangent2= (grid.face[f].normal).cross(faceTangent1);

		if (order!="first") {
			for (unsigned int i=0;i<5;++i) delta[i]=(grid.face[f].centroid-grid.cell[parent].centroid).dot(grid.cell[parent].limited_grad[i]);
			deltaV.comp[0]=delta[1]; deltaV.comp[1]=delta[2]; deltaV.comp[2]=delta[3];
		}

		rhoL=grid.cell[parent].rho+delta[0];
		pL=grid.cell[parent].p+delta[4];
		uNL=(grid.cell[parent].v+deltaV).dot(grid.face[f].normal);
		vTL=(grid.cell[parent].v+deltaV).dot(faceTangent1);
		wTL=(grid.cell[parent].v+deltaV).dot(faceTangent2);
		
		if (grid.face[f].bc>=0) { // means a real boundary face
			rhoR=rhoL; uNR=uNL; vTR=vTL; wTR=wTL; pR=pL; // outlet condition
			if (bc.region[grid.face[f].bc].type=="slip") uNR=-uNL;
			if (bc.region[grid.face[f].bc].type=="noslip") {uNR=-uNL; vTR=0.; wTR=0.;}
			if (bc.region[grid.face[f].bc].type=="inlet") {
				uNR=bc.region[grid.face[f].bc].v.dot(grid.face[f].normal);
				vTR=bc.region[grid.face[f].bc].v.dot(faceTangent1);
				wTR=bc.region[grid.face[f].bc].v.dot(faceTangent2);
				rhoR=bc.region[grid.face[f].bc].rho;
				//pR=bc.region[grid.face[f].bc].p;
			}
		} else if (grid.face[f].bc==-1) { // internal face

			if (order!="first") {
				for (unsigned int i=0;i<5;++i) delta[i]=(grid.face[f].centroid-grid.cell[neighbor].centroid).dot(grid.cell[neighbor].limited_grad[i]);
				deltaV.comp[0]=delta[1]; deltaV.comp[1]=delta[2]; deltaV.comp[2]=delta[3];
			}
			
			rhoR=grid.cell[neighbor].rho+delta[0];
			pR=grid.cell[neighbor].p+delta[4];
			uNR=(grid.cell[neighbor].v+deltaV).dot(grid.face[f].normal);
			vTR=(grid.cell[neighbor].v+deltaV).dot(faceTangent1);
			wTR=(grid.cell[neighbor].v+deltaV).dot(faceTangent2);
		} else { // partition boundary
			int g=-1*grid.face[f].bc-3;

			if (order!="first") {
				for (unsigned int i=0;i<5;++i) 				delta[i]=(grid.face[f].centroid-grid.ghost[g].centroid).dot(grid.ghost[g].limited_grad[i]);
				deltaV.comp[0]=delta[1]; deltaV.comp[1]=delta[2]; deltaV.comp[2]=delta[3];
			}
			rhoR=grid.ghost[g].rho+delta[0];
			pR=grid.ghost[g].p+delta[4];
			uNR=(grid.ghost[g].v+deltaV).dot(grid.face[f].normal);
			vTR=(grid.ghost[g].v+deltaV).dot(faceTangent1);
			wTR=(grid.ghost[g].v+deltaV).dot(faceTangent2);
		}
		
// 		if (grid.face[f].bc>=0 && (bc.region[grid.face[f].bc].type=="outlet" | bc.region[grid.face[f].bc].type=="inlet")) {
// 			fluxNormal[0]=rhoR*uNR;
// 			fluxNormal[1]=rhoR*uNR*uNR+pR;
// 			fluxNormal[2]=rhoR*uNR*vTR;
// 			fluxNormal[3]=rhoR*uNR*wTR;
// 			fluxNormal[4]=uNR* (0.5*rhoR* (uNR*uNR+vTR*vTR+wTR*wTR) +pR/ (Gamma - 1.) +pR);
// 		} else {
			qL[0]=rhoL;
			qL[1]=qL[0] * uNL;
			qL[2]=qL[0] * vTL;
			qL[3]=qL[0] * wTL;
			qL[4]=0.5*qL[0]* (uNL*uNL+vTL*vTL+wTL*wTL) +pL/ (Gamma - 1.);

			qR[0] = rhoR;
			qR[1] = qR[0] * uNR;
			qR[2] = qR[0] * vTR;
			qR[3] = qR[0] * wTR;
			qR[4] = 0.5*qR[0]* (uNR*uNR+vTR*vTR+wTR*wTR) +pR/ (Gamma - 1.);

			RoeFlux(qL, qR, fluxNormal);
		//} 


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

	return;
} // end function

