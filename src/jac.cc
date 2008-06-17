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

	double rhoL,rhoR,uNL,uNR,vTL,vTR,wTL,wTR, pL, pR;
	Vec3D faceTangent1, faceTangent2;
	double qL[5], qR[5], fluxNormal[5], flux[5], fluxPlus[5], delta[5];

	unsigned int parent,neighbor,f;

	int row,col;
	PetscScalar value;

	for (int i=0;i<5;++i) delta[i]=0.;
		
	for (unsigned int f=0;f<grid.faceCount;++f) {
			parent=grid.face[f].parent; neighbor=grid.face[f].neighbor;

			if (grid.face[f].bc==-1) { // means a real internal face
				
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

				uNR=grid.cell[neighbor].v.dot(grid.face[f].normal);
				vTR=grid.cell[neighbor].v.dot(faceTangent1);
				wTR=grid.cell[neighbor].v.dot(faceTangent2);
				rhoR=grid.cell[neighbor].rho;
				pR=grid.cell[neighbor].p;
				
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

				flux[0] = fluxNormal[0]*grid.face[f].area;
				flux[1] = (fluxNormal[1]*grid.face[f].normal.comp[0]+fluxNormal[2]*faceTangent1.comp[0]+fluxNormal[3]*faceTangent2.comp[0]) * grid.face[f].area;
				flux[2] = (fluxNormal[1]*grid.face[f].normal.comp[1]+fluxNormal[2]*faceTangent1.comp[1]+fluxNormal[3]*faceTangent2.comp[1]) * grid.face[f].area;
				flux[3] = (fluxNormal[1]*grid.face[f].normal.comp[2]+fluxNormal[2]*faceTangent1.comp[2]+fluxNormal[3]*faceTangent2.comp[2]) * grid.face[f].area;
				flux[4] = fluxNormal[4]*grid.face[f].area;

				//for (int i=0;i<5;++i) delta[i]=max(fabs(qL[i]*1.e-4),1.e-4);
				//for (int i=0;i<5;++i) delta[i]=sqrt(1.e-12)*(1.+fabs(qL[i]));
				for (int i=0;i<5;++i) delta[i]=sqrt(std::numeric_limits<double>::epsilon());
								
				for (int i=0;i<5;++i) {

					qL[0]=rhoL;
					qL[1]=qL[0] * uNL;
					qL[2]=qL[0] * vTL;
					qL[3]=qL[0] * wTL;
					qL[4]=0.5*qL[0]* (uNL*uNL+vTL*vTL+wTL*wTL) +pL/ (Gamma - 1.);

					qL[i]+=delta[i];

					RoeFlux(qL, qR, fluxNormal);

					fluxPlus[0] = fluxNormal[0]*grid.face[f].area;
					fluxPlus[1] = (fluxNormal[1]*grid.face[f].normal.comp[0]+fluxNormal[2]*faceTangent1.comp[0]+fluxNormal[3]*faceTangent2.comp[0]) * grid.face[f].area;
					fluxPlus[2] = (fluxNormal[1]*grid.face[f].normal.comp[1]+fluxNormal[2]*faceTangent1.comp[1]+fluxNormal[3]*faceTangent2.comp[1]) * grid.face[f].area;
					fluxPlus[3] = (fluxNormal[1]*grid.face[f].normal.comp[2]+fluxNormal[2]*faceTangent1.comp[2]+fluxNormal[3]*faceTangent2.comp[2]) * grid.face[f].area;
					fluxPlus[4] = fluxNormal[4]*grid.face[f].area;

					//int j=i;
					for (int j=0;j<5;++j) {
						row=grid.cell[parent].globalId*5+j;
						col=grid.cell[parent].globalId*5+i;
						value=(fluxPlus[j]-flux[j])/delta[j];
						
						//cout << fluxPlus[i] << "\t" << flux[i] << "\t" << delta[i] << "\t" << value << endl;
						MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
						row=grid.cell[neighbor].globalId*5+j;
						//if (col>row) cout << row << "\t" << col << endl;
						value*=-1.;
						//cout << row  << "\t" << col << endl;
						MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
					} // for j

				} // for i

				qL[0]=rhoL;
				qL[1]=qL[0] * uNL;
				qL[2]=qL[0] * vTL;
				qL[3]=qL[0] * wTL;
				qL[4]=0.5*qL[0]* (uNL*uNL+vTL*vTL+wTL*wTL) +pL/ (Gamma - 1.);

				//for (int i=0;i<5;++i) delta[i]=max(fabs(qR[i]*1.e-4),1.e-4);
				//for (int i=0;i<5;++i) delta[i]=sqrt(1.e-12)*(1.+fabs(qR[i]));
				for (int i=0;i<5;++i) delta[i]=sqrt(std::numeric_limits<double>::epsilon());
								
				for (int i=0;i<5;++i) {

					qR[0] = rhoR;
					qR[1] = qR[0] * uNR;
					qR[2] = qR[0] * vTR;
					qR[3] = qR[0] * wTR;
					qR[4] = 0.5*qR[0]* (uNR*uNR+vTR*vTR+wTR*wTR) +pR/ (Gamma - 1.);

					qR[i]+=delta[i];

					RoeFlux(qL, qR, fluxNormal);

					fluxPlus[0] = fluxNormal[0]*grid.face[f].area;
					fluxPlus[1] = (fluxNormal[1]*grid.face[f].normal.comp[0]+fluxNormal[2]*faceTangent1.comp[0]+fluxNormal[3]*faceTangent2.comp[0]) * grid.face[f].area;
					fluxPlus[2] = (fluxNormal[1]*grid.face[f].normal.comp[1]+fluxNormal[2]*faceTangent1.comp[1]+fluxNormal[3]*faceTangent2.comp[1]) * grid.face[f].area;
					fluxPlus[3] = (fluxNormal[1]*grid.face[f].normal.comp[2]+fluxNormal[2]*faceTangent1.comp[2]+fluxNormal[3]*faceTangent2.comp[2]) * grid.face[f].area;
					fluxPlus[4] = fluxNormal[4]*grid.face[f].area;

					//int j=i;
					for (int j=0;j<5;++j) {
						row=grid.cell[neighbor].globalId*5+j;
						col=grid.cell[neighbor].globalId*5+i;
						value=(flux[j]-fluxPlus[j])/delta[j];
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


