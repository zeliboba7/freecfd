#include <cmath>
#include "grid.h"
#include "bc.h"
#include "petscksp.h"

extern Grid grid;
extern BC bc;
extern int rank;
extern double Gamma;

void RoeFlux(const double qL[], const double qR[], double flux[]);

void jac(Mat impOP) {

	double rhoL,rhoR,pL,pR,uL,vL,wL,uR,vR,wR;
	Vec3D faceTangent1, faceTangent2, qVL, qVR, deltaV;
	double qNL[5], qNR[5], fluxNormal[5], flux[5], fluxPlus[5], deltaL[5], deltaR[5];
	
	unsigned int parent,neighbor,f;

	int row,col;
	PetscScalar value;
	
	for (int i=0;i<5;++i) {
		deltaL[i]=sqrt(std::numeric_limits<double>::epsilon());
		deltaR[i]=sqrt(std::numeric_limits<double>::epsilon());
	}
	deltaV=0.;

	for (unsigned int f=0;f<grid.faceCount;++f) {
			parent=grid.face[f].parent; neighbor=grid.face[f].neighbor;

			if (grid.face[f].bc==-1) { // means a real internal face

				rhoL=grid.cell[parent].rho;
				pL=grid.cell[parent].p;

				rhoR=grid.cell[neighbor].rho;
				pR=grid.cell[neighbor].p;
				
				// Take 2 nodes of the face to find a vector in the face plane
				faceTangent1= (grid.face[f].node(1)-grid.face[f].node(0));
				faceTangent1/=fabs(faceTangent1);	
				// Cross the tangent vector with the normal vector to get the second tangent
				faceTangent2= (grid.face[f].normal).cross(faceTangent1);
		
				qVL.comp[0]=rhoL * uL;
				qVL.comp[1]=rhoL * vL;
				qVL.comp[2]=rhoL * wL;
				
				qNL[0]=rhoL;
				qNL[1]=qVL.dot(grid.face[f].normal);
				qNL[2]=qVL.dot(faceTangent1);
				qNL[3]=qVL.dot(faceTangent2);
				qNL[4]=0.5*rhoL* (uL*uL+vL*vL+wL*wL) +pL/ (Gamma - 1.);
				
				qVR.comp[0]=rhoR * uR;
				qVR.comp[1]=rhoR * vR;
				qVR.comp[2]=rhoR * wR;
				
				qNR[0]=rhoR;
				qNR[1]=qVR.dot(grid.face[f].normal);
				qNR[2]=qVR.dot(faceTangent1);
				qNR[3]=qVR.dot(faceTangent2);
				qNR[4]=0.5*rhoR* (uR*uR+vR*vR+wR*wR) +pR/ (Gamma - 1.);

				
				deltaV=0.;deltaV.comp[0]=deltaL[0];
				deltaL[1]=(qVL+deltaV).dot(grid.face[f].normal)-qNL[1];
				deltaR[1]=(qVR+deltaV).dot(grid.face[f].normal)-qNR[1];
								
				deltaV=0.; deltaV.comp[1]=deltaL[0];
				deltaL[2]=(qVL+deltaV).dot(faceTangent1)-qNL[2];
				deltaR[2]=(qVR+deltaV).dot(faceTangent1)-qNR[2];
				
				deltaV=0.; deltaV.comp[2]=deltaL[0];
				deltaL[3]=(qVL+deltaV).dot(faceTangent2)-qNL[3];
				deltaR[3]=(qVR+deltaV).dot(faceTangent2)-qNR[3];
				
				RoeFlux(qNL, qNR, fluxNormal);

				flux[0] = fluxNormal[0]*grid.face[f].area;
				flux[1] = (fluxNormal[1]*grid.face[f].normal.comp[0]+fluxNormal[2]*faceTangent1.comp[0]+fluxNormal[3]*faceTangent2.comp[0]) * grid.face[f].area;
				flux[2] = (fluxNormal[1]*grid.face[f].normal.comp[1]+fluxNormal[2]*faceTangent1.comp[1]+fluxNormal[3]*faceTangent2.comp[1]) * grid.face[f].area;
				flux[3] = (fluxNormal[1]*grid.face[f].normal.comp[2]+fluxNormal[2]*faceTangent1.comp[2]+fluxNormal[3]*faceTangent2.comp[2]) * grid.face[f].area;
				flux[4] = fluxNormal[4]*grid.face[f].area;

				for (int i=0;i<5;++i) {

					qNL[i]+=deltaL[i];

					RoeFlux(qNL, qNR, fluxNormal);

					qNL[i]-=deltaL[i];

					fluxPlus[0] = fluxNormal[0]*grid.face[f].area;
					fluxPlus[1] = (fluxNormal[1]*grid.face[f].normal.comp[0]+fluxNormal[2]*faceTangent1.comp[0]+fluxNormal[3]*faceTangent2.comp[0]) * grid.face[f].area;
					fluxPlus[2] = (fluxNormal[1]*grid.face[f].normal.comp[1]+fluxNormal[2]*faceTangent1.comp[1]+fluxNormal[3]*faceTangent2.comp[1]) * grid.face[f].area;
					fluxPlus[3] = (fluxNormal[1]*grid.face[f].normal.comp[2]+fluxNormal[2]*faceTangent1.comp[2]+fluxNormal[3]*faceTangent2.comp[2]) * grid.face[f].area;
					fluxPlus[4] = fluxNormal[4]*grid.face[f].area;

					// i is the variable, j is the flux component
					for (int j=0;j<5;++j) {
						row=grid.cell[parent].globalId*5+j;
						col=grid.cell[parent].globalId*5+i;
						value=(fluxPlus[j]-flux[j])/deltaL[i];
						
						//cout << fluxPlus[i] << "\t" << flux[i] << "\t" << delta[i] << "\t" << value << endl;
						MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
						row=grid.cell[neighbor].globalId*5+j;
						//if (col>row) cout << row << "\t" << col << endl;
						value*=-1.;
						//cout << row  << "\t" << col << endl;
						MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
					} // for j

				} // for i
								
				for (int i=0;i<5;++i) {

					qNR[i]+=deltaR[i];

					RoeFlux(qNL, qNR, fluxNormal);

					qNR[i]-=deltaR[i];

					fluxPlus[0] = fluxNormal[0]*grid.face[f].area;
					fluxPlus[1] = (fluxNormal[1]*grid.face[f].normal.comp[0]+fluxNormal[2]*faceTangent1.comp[0]+fluxNormal[3]*faceTangent2.comp[0]) * grid.face[f].area;
					fluxPlus[2] = (fluxNormal[1]*grid.face[f].normal.comp[1]+fluxNormal[2]*faceTangent1.comp[1]+fluxNormal[3]*faceTangent2.comp[1]) * grid.face[f].area;
					fluxPlus[3] = (fluxNormal[1]*grid.face[f].normal.comp[2]+fluxNormal[2]*faceTangent1.comp[2]+fluxNormal[3]*faceTangent2.comp[2]) * grid.face[f].area;
					fluxPlus[4] = fluxNormal[4]*grid.face[f].area;

					for (int j=0;j<5;++j) {
						row=grid.cell[neighbor].globalId*5+j;
						col=grid.cell[neighbor].globalId*5+i;
						value=(flux[j]-fluxPlus[j])/deltaR[i];
						//cout << fluxPlus[i] << "\t" << flux[i] << "\t" << delta[i] << "\t" << value << endl;
						MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
						row=grid.cell[parent].globalId*5+j;
						//if (col>row) cout << row << "\t" << col << endl;
						value*=-1.;
						//cout << row  << "\t" << col << endl;
						MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
					} // for j

				} // for i



				
			} // if internal face
		} // for faces
		//cout << c << endl;

	

	return;
} // end function


