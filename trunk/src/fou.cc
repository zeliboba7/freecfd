#include <cmath>
#include "grid.h"
#include "bc.h"

extern Grid grid;
extern BC bc;
extern int rank;

void update(double dt,double gamma);

void RoeFlux(const double gamma, const double qL[], const double qR[], double flux[]);

void fou(double gamma) {

	double rhoL,rhoR,uNL,uNR,vTL,vTR,wTL,wTR, pL, pR;
	Vec3D faceTangent1, faceTangent2;
	double qL[5], qR[5], fluxNormal[5], flux[5];

	unsigned int parent, neighbor;

	// Calculate flux through each face
	for (unsigned int f = 0;f < grid.faceCount;++f) {

		parent=grid.face[f].parent; neighbor=grid.face[f].neighbor;
										
		// Take 2 nodes of the face to find a vector in the face plane
		faceTangent1= (grid.face[f].node(1)-grid.face[f].node(0));
		faceTangent1/=fabs(faceTangent1);
		// Cross the tangent vector with the normal vector to get the second tangent
		faceTangent2= (grid.face[f].normal).cross(faceTangent1);

		uNL=grid.cell[parent].v.dot(grid.face[f].normal);
		vTL=grid.cell[parent].v.dot(faceTangent1);
		wTL=grid.cell[parent].v.dot(faceTangent2);
		rhoL=grid.cell[parent].rho;
		pL=grid.cell[parent].p;
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
		} else { // internal face
			uNR=grid.cell[neighbor].v.dot(grid.face[f].normal);
			vTR=grid.cell[neighbor].v.dot(faceTangent1);
			wTR=grid.cell[neighbor].v.dot(faceTangent2);
			rhoR=grid.cell[neighbor].rho;
			pR=grid.cell[neighbor].p;
		}
		if (grid.face[f].bc>=0 && (bc.region[grid.face[f].bc].type=="outlet" | bc.region[grid.face[f].bc].type=="inlet")) {
			fluxNormal[0]=rhoR*uNR;
			fluxNormal[1]=rhoR*uNR*uNR+pR;
			fluxNormal[2]=rhoR*uNR*vTR;
			fluxNormal[3]=rhoR*uNR*wTR;
			fluxNormal[4]=uNR* (0.5*rhoR* (uNR*uNR+vTR*vTR+wTR*wTR) +pR/ (gamma - 1.) +pR);
		} else {
			qL[0]=rhoL;
			qL[1]=qL[0] * uNL;
			qL[2]=qL[0] * vTL;
			qL[3]=qL[0] * wTL;
			qL[4]=0.5*qL[0]* (uNL*uNL+vTL*vTL+wTL*wTL) +pL/ (gamma - 1.);

			qR[0] = rhoR;
			qR[1] = qR[0] * uNR;
			qR[2] = qR[0] * vTR;
			qR[3] = qR[0] * wTR;
			qR[4] = 0.5*qR[0]* (uNR*uNR+vTR*vTR+wTR*wTR) +pR/ (gamma - 1.);

			RoeFlux(gamma, qL, qR, fluxNormal);
		}


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

