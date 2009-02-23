/************************************************************************
	
	Copyright 2007-2009 Emre Sozer & Patrick Clark Trizila

	Contact: emresozer@freecfd.com , ptrizila@freecfd.com

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
#include "commons.h"
#include <cmath>
#include "grid.h"
#include "bc.h"
#include "inputs.h"
#include "petsc_functions.h"
#include "turbulence.h"

extern BC bc;
extern InputFile input;

double EPS=1.e-10;
TurbulenceModel kepsilon,komega;


void get_kOmega(void);
double get_blending(double &k,double &omega,double &rho,double &y,Vec3D &gradK,Vec3D &gradOmega);
void set_turbulence_model_constants(void);

unsigned int parent,neighbor,f;
double leftK,leftOmega, rightK, rightOmega, faceK, faceOmega;
double leftK_center, rightK_center, leftOmega_center, rightOmega_center;
Vec3D faceGradK, faceGradOmega, left2right;
int row,col;
double mu_t=0.;
double value;
double dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz;
double strainRate;

double convectiveFlux[2],diffusiveFlux[2],source[2];
double jacL[2],jacR[2];
double weightL,weightR;
double faceRho;
double closest_wall_distance;
// For blending function
double cross_diffusion,arg,arg1,arg2,arg3,CD;

void set_turbulence_model_constants(void) {
	
	kepsilon.sigma_k=1.;
	kepsilon.sigma_omega=0.856;
	kepsilon.beta=0.0828;
	kepsilon.beta_star=0.09;
	kepsilon.kappa=0.41;
	kepsilon.alpha=kepsilon.beta/kepsilon.beta_star
			-kepsilon.sigma_omega*kepsilon.kappa*kepsilon.kappa/sqrt(kepsilon.beta_star);
	
	komega.sigma_k=0.5;
	komega.sigma_omega=0.5;
	komega.beta=0.075;
	komega.beta_star=0.09;
	komega.kappa=0.41;
	komega.alpha=komega.beta/komega.beta_star
			-komega.sigma_omega*komega.kappa*komega.kappa/sqrt(komega.beta_star);

	return;
}	

void linear_system_turb(void) {

	double blending;
	double sigma_k,sigma_omega,beta,beta_star,kappa,alpha;


	MatZeroEntries(impOP_turb);

	// Loop through faces
	for (f=0;f<grid.faceCount;++f) {
		if (TURBULENCE_MODEL==KOMEGA) blending=1.;
		if (TURBULENCE_MODEL==KEPSILON) blending=0.;
		
		parent=grid.face[f].parent; neighbor=grid.face[f].neighbor;
		
		// Get left, right and face values of k and omega as well as the face normal gradients
		get_kOmega();
		
		// Convective flux is based on mdot calculated through the Riemann solver right after
		// main flow was updated
		
		// These weights are coming from the Riemann solver
		weightL=grid.face[f].weightL; weightR=grid.face[f].weightR;
				
// 		double weightMax=max(weightL,weightR);
// 		weightL+=0.05*weightMax;
// 		weightR+=0.05*weightMax;

		convectiveFlux[0]=grid.face[f].mdot*(weightL*leftK+weightR*rightK)*grid.face[f].area;
		convectiveFlux[1]=grid.face[f].mdot*(weightL*leftOmega+weightR*rightOmega)*grid.face[f].area;

		// Diffusive k and omega fluxes	
		if (TURBULENCE_MODEL==BSL) {
			closest_wall_distance=grid.cell[parent].closest_wall_distance;
			blending=get_blending(faceK,faceOmega,faceRho,closest_wall_distance,faceGradK,faceGradOmega);
		}
		
		sigma_omega=blending*komega.sigma_omega+(1.-blending)*kepsilon.sigma_omega;
		sigma_k=blending*komega.sigma_k+(1.-blending)*kepsilon.sigma_k;
		alpha=blending*komega.alpha+(1.-blending)*kepsilon.alpha;
		beta=blending*komega.beta+(1.-blending)*kepsilon.beta;
		beta_star=blending*komega.beta_star+(1.-blending)*kepsilon.beta_star;
		
		mu_t=fabs(faceRho*faceK/faceOmega);
		mu_t=min(mu_t,viscosityRatioLimit*viscosity);
		
		diffusiveFlux[0]=(viscosity+mu_t*sigma_k)*faceGradK.dot(grid.face[f].normal)*grid.face[f].area;
		diffusiveFlux[1]=(viscosity+mu_t*sigma_omega)*faceGradOmega.dot(grid.face[f].normal)*grid.face[f].area;
		
		// Fill in rhs vector for turbulence scalars
		for (int i=0;i<2;++i) {
			row=(grid.myOffset+parent)*2+i;
			value=diffusiveFlux[i]-convectiveFlux[i];
			VecSetValues(rhs_turb,1,&row,&value,ADD_VALUES);
			if (grid.face[f].bc==INTERNAL) { 
				row=(grid.myOffset+neighbor)*2+i;
				value*=-1.;
				VecSetValues(rhs_turb,1,&row,&value,ADD_VALUES);
			}
		}
		
		// Calculate flux jacobians
		
		// Assumes k flux doesn't change with omega and vice versa 
		// This is true for convective flux (effect of mu_t in diffusive flux ignored)
		// TODO What is density doing in diffusive flux jacobians?
		double leftRho=grid.cell[parent].rho;
		double rightRho=1.; // if right cell is not internal, jacobian is not inserted into the matrix, so this doesn't matter
		if (grid.face[f].bc==INTERNAL) rightRho=grid.cell[neighbor].rho;
		// dF_k/dk_left
		jacL[0]=weightL*grid.face[f].mdot*grid.face[f].area; // convective
		jacL[0]+=(viscosity+mu_t*sigma_k)/(leftRho*left2right.dot(grid.face[f].normal))*grid.face[f].area; // diffusive
		// dF_k/dk_right
		jacR[0]=weightR*grid.face[f].mdot*grid.face[f].area; // convective
		jacR[0]-=(viscosity+mu_t*sigma_k)/(rightRho*left2right.dot(grid.face[f].normal))*grid.face[f].area; // diffusive
		// dF_omega/dOmega_left
		jacL[1]=weightL*grid.face[f].mdot*grid.face[f].area; // convective
		jacL[1]+=(viscosity+mu_t*sigma_omega)/(leftRho*left2right.dot(grid.face[f].normal))*grid.face[f].area; // diffusive
		// dF_omega/dOmega_right
		jacR[1]=weightR*grid.face[f].mdot*grid.face[f].area; // convective
		jacR[1]-=(viscosity+mu_t*sigma_omega)/(rightRho*left2right.dot(grid.face[f].normal))*grid.face[f].area; // diffusive
		
		// Insert flux jacobians for the parent cell
		// left_k/left_k
		row=(grid.myOffset+parent)*2; col=row; value=jacL[0];
		MatSetValues(impOP_turb,1,&row,1,&col,&value,ADD_VALUES);
		// left_omega/left_omega
		row++; col++; value=jacL[1];
		MatSetValues(impOP_turb,1,&row,1,&col,&value,ADD_VALUES);
		if (grid.face[f].bc==INTERNAL) { 
			// left_k/right_k
			row=(grid.myOffset+parent)*2; col=(grid.myOffset+neighbor)*2; value=jacR[0];
			MatSetValues(impOP_turb,1,&row,1,&col,&value,ADD_VALUES);
			// left_omega/right_omega
			row++; col++; value=jacR[1];
			MatSetValues(impOP_turb,1,&row,1,&col,&value,ADD_VALUES);
			
			// Insert flux jacobians for the neighbor cell
			// right_k/right_k
			row=(grid.myOffset+neighbor)*2; col=row; value=-jacR[0];
			MatSetValues(impOP_turb,1,&row,1,&col,&value,ADD_VALUES);
			// right_omega/right_omega
			row++; col++; value=-jacR[1];
			MatSetValues(impOP_turb,1,&row,1,&col,&value,ADD_VALUES);
			// right_k/left_k
			row=(grid.myOffset+neighbor)*2; col=(grid.myOffset+parent)*2; value=-jacL[0];
			MatSetValues(impOP_turb,1,&row,1,&col,&value,ADD_VALUES);
			// right_omega/left_omega
			row++; col++; value=-jacL[1];
			MatSetValues(impOP_turb,1,&row,1,&col,&value,ADD_VALUES);
		} else if (grid.face[f].bc==GHOST) { 
			// left_k/right_k
			row=(grid.myOffset+parent)*2; col=(grid.ghost[-1*neighbor-1].matrix_id)*2; value=jacR[0];
			MatSetValues(impOP_turb,1,&row,1,&col,&value,ADD_VALUES);
			// left_omega/right_omega
			row++; col++; value=jacR[1];
			MatSetValues(impOP_turb,1,&row,1,&col,&value,ADD_VALUES);
		}

	} // for faces
	
	// Now to a cell loop to add the unsteady and source terms
	double divU,Prod_k;
	for (unsigned int c=0;c<grid.cellCount;++c) {
		double d=0.;
		
		if (TIME_STEP_TYPE==CFL_LOCAL) {
			// Determine time step with CFL condition
			double dtLocal=1.E20;
			double a=sqrt(Gamma*(grid.cell[c].p+Pref)/grid.cell[c].rho);
			double lengthScale=grid.cell[c].lengthScale;
			dtLocal=min(dtLocal,CFLlocal*lengthScale/(fabs(grid.cell[c].v[0])+a));
			dtLocal=min(dtLocal,CFLlocal*lengthScale/(fabs(grid.cell[c].v[1])+a));
			dtLocal=min(dtLocal,CFLlocal*lengthScale/(fabs(grid.cell[c].v[2])+a));
			d=grid.cell[c].volume/dtLocal;
		} else {
			d=grid.cell[c].volume/dt;
		}

		// Insert unsteady term
		row=(grid.myOffset+c)*2;
		value=grid.cell[c].rho*d;
		MatSetValues(impOP_turb,1,&row,1,&row,&value,ADD_VALUES);
		row++;
		MatSetValues(impOP_turb,1,&row,1,&row,&value,ADD_VALUES);

		if (TURBULENCE_MODEL==BSL) {
			blending=get_blending(grid.cell[c].k,grid.cell[c].omega,grid.cell[c].rho,grid.cell[c].closest_wall_distance,grid.cell[c].grad_turb[0],grid.cell[c].grad_turb[1]);
		}
		
		sigma_omega=blending*komega.sigma_omega+(1.-blending)*kepsilon.sigma_omega;
		sigma_k=blending*komega.sigma_k+(1.-blending)*kepsilon.sigma_k;
		alpha=blending*komega.alpha+(1.-blending)*kepsilon.alpha;
		beta=blending*komega.beta+(1.-blending)*kepsilon.beta;
		beta_star=blending*komega.beta_star+(1.-blending)*kepsilon.beta_star;
		
		// Calculate the source terms
		mu_t=fabs(grid.cell[c].rho*grid.cell[c].k/grid.cell[c].omega);
		mu_t=min(mu_t,viscosityRatioLimit*viscosity);
		
		dudx=grid.cell[c].grad[1][0];
		dudy=grid.cell[c].grad[1][1];
		dudz=grid.cell[c].grad[1][2];
		
		dvdx=grid.cell[c].grad[2][0];
		dvdy=grid.cell[c].grad[2][1];
		dvdz=grid.cell[c].grad[2][2];
		
		dwdx=grid.cell[c].grad[3][0];
		dwdy=grid.cell[c].grad[3][1];
		dwdz=grid.cell[c].grad[3][2];
		
		divU=dudx+dvdy+dwdz;
		strainRate=sqrt(2.*(dudx*dudx + dvdy*dvdy + dwdz*dwdz)
				+ (dudy+dvdx)*(dudy+dvdx) + (dudz+dwdx)*(dudz+dwdx)
				+ (dvdz+dwdy)*(dvdz+dwdy)
			       );
		
		Prod_k=mu_t*(strainRate*strainRate-2./3.*divU*divU)-2./3.*grid.cell[c].rho*grid.cell[c].k*divU;

		//if (Prod_k>10./grid.cell[c].volume) cout << "[I] k-production term clipped" << endl;
		//Prod_k=min(Prod_k,10./grid.cell[c].volume);
		
		cross_diffusion=0.;
		if (TURBULENCE_MODEL!=KOMEGA) cross_diffusion=grid.cell[c].grad_turb[0].dot(grid.cell[c].grad_turb[1]);
			
		source[0]=Prod_k; // k production
		source[0]-=beta_star*grid.cell[c].rho*grid.cell[c].omega*grid.cell[c].k; // k destruction
		
		source[1]=alpha/(mu_t/grid.cell[c].rho+EPS)*Prod_k; // omega production
		source[1]-=beta*grid.cell[c].rho*grid.cell[c].omega*grid.cell[c].omega; // omega destruction
				
		source[1]+=2.*(1.-blending)*grid.cell[c].rho*komega.sigma_omega*cross_diffusion/grid.cell[c].omega;
		
		// Add source terms to rhs
		for (int i=0;i<2;++i) {
			row=(grid.myOffset+c)*2+i;
			value=source[i]*grid.cell[c].volume;
			VecSetValues(rhs_turb,1,&row,&value,ADD_VALUES);
		}
		
		// Add source jacobians
		// Only include destruction terms
		
		// dS_k/dk
		row=(grid.myOffset+c)*2; col=row;
		value=beta_star*grid.cell[c].rho*grid.cell[c].omega*grid.cell[c].volume; // approximate 
		MatSetValues(impOP_turb,1,&row,1,&col,&value,ADD_VALUES);
		
		// dS_k/dOmega
		col++; 
		value=beta_star*grid.cell[c].rho*grid.cell[c].k*grid.cell[c].volume; 
		MatSetValues(impOP_turb,1,&row,1,&col,&value,ADD_VALUES);
		
		// dS_omega/dOmega
		row++; 
		value=2.*beta*grid.cell[c].rho*grid.cell[c].omega*grid.cell[c].volume; // approximate
		// Add cross-diffusion term jacobian
		value+=2.*(1.-blending)*komega.sigma_omega*komega.sigma_omega*cross_diffusion/(grid.cell[c].omega*grid.cell[c].omega)*grid.cell[c].volume;
		MatSetValues(impOP_turb,1,&row,1,&col,&value,ADD_VALUES);
	
		
	} // end cell loop
	
	return;
} // end function

void get_kOmega() {
	
	// Get left k and omega
	double delta[2];

	if (order==SECOND) {
		for (unsigned int i=0;i<2;++i) {
			delta[i]=(grid.face[f].centroid-grid.cell[parent].centroid).dot(grid.cell[parent].grad_turb[i]);
		}
	} else {
		for (unsigned int i=0;i<2;++i) delta[i]=0.;
	}
	
	leftK_center=grid.cell[parent].k;
	leftK=leftK_center+delta[0];
	leftOmega_center=grid.cell[parent].omega;
	leftOmega=leftOmega_center+delta[1];
	
	leftK=max(leftK,kLowLimit);
	leftOmega=max(leftOmega,omegaLowLimit);

	// Find face averaged quantities
	faceRho=0.; faceGradK=0.; faceGradOmega=0.;
 	map<int,double>::iterator fit;
	for (fit=grid.face[f].average.begin();fit!=grid.face[f].average.end();fit++) {
		if ((*fit).first>=0) { // if contribution is coming from a real cell
			faceRho+=(*fit).second*grid.cell[(*fit).first].rho;
			faceGradK+=(*fit).second*grid.cell[(*fit).first].grad_turb[0];
			faceGradOmega+=(*fit).second*grid.cell[(*fit).first].grad_turb[1];
		} else { // if contribution is coming from a ghost cell
			faceRho+=(*fit).second*grid.ghost[-1*((*fit).first+1)].rho;
			faceGradK+=(*fit).second*grid.ghost[-1*((*fit).first+1)].grad_turb[0];
			faceGradOmega+=(*fit).second*grid.ghost[-1*((*fit).first+1)].grad_turb[1];
		}
	}
	
	// Find distance between left and right centroids 
	Vec3D leftCentroid,rightCentroid;
	leftCentroid=grid.cell[parent].centroid;
	if (grid.face[f].bc==INTERNAL) { rightCentroid=grid.cell[neighbor].centroid;}
	else if (grid.face[f].bc>=0) {
		rightCentroid=leftCentroid+2.*(grid.face[f].centroid-leftCentroid).dot(grid.face[f].normal)*grid.face[f].normal;
	} else { rightCentroid=grid.ghost[-1*grid.face[f].neighbor-1].centroid;}
	
	left2right=rightCentroid-leftCentroid;
		
	// Get right k and omega
	if (grid.face[f].bc==INTERNAL) {// internal face

		if (order==SECOND) {
			for (unsigned int i=0;i<2;++i) {
				delta[i]=(grid.face[f].centroid-grid.cell[neighbor].centroid).dot(grid.cell[neighbor].grad_turb[i]);
			}
		} else {
			for (unsigned int i=0;i<2;++i) delta[i]=0.;
		}

		rightK_center=grid.cell[neighbor].k;
		rightK=rightK_center+delta[0];
		rightOmega_center=grid.cell[neighbor].omega;
		rightOmega=rightOmega_center+delta[1];
	
		rightK=max(rightK,kLowLimit);
		rightOmega=max(rightOmega,omegaLowLimit);
		
	} else if (grid.face[f].bc>=0) { // boundary face
		
		rightK=leftK; rightOmega=leftOmega;

		if (bc.region[grid.face[f].bc].type==NOSLIP) {
			rightK=0.;
			rightOmega=60.*viscosity/(faceRho*0.075*pow(0.5*left2right.dot(grid.face[f].normal),2.));
			rightK_center=2.*rightK-leftK_center; 
			rightOmega_center=2.*rightOmega-leftOmega_center;
		} else if (bc.region[grid.face[f].bc].type==SYMMETRY) {
			rightK_center=leftK_center; 
			rightOmega_center=leftOmega_center;
		} else if (bc.region[grid.face[f].bc].type==SLIP) {
			rightK_center=leftK_center; 
			rightOmega_center=leftOmega_center;
		} else if (bc.region[grid.face[f].bc].type==INLET) {
			rightK=bc.region[grid.face[f].bc].k;
			rightOmega=bc.region[grid.face[f].bc].omega;
			rightK_center=rightK;
			rightOmega_center=rightOmega;
			//cout << bc.region[grid.face[f].bc].omega << "\t" << leftOmega << "\t" << faceOmega << "\t" << leftOmega_center << endl;
		} else if (bc.region[grid.face[f].bc].type==OUTLET) {
			rightK_center=rightK;
			rightOmega_center=rightOmega;
		}
		
	} else { // partition boundary

		int g=-1*grid.face[f].neighbor-1; // ghost cell index

		if (order==SECOND) {
			for (unsigned int i=0;i<2;++i) {
				delta[i]=(grid.face[f].centroid-grid.ghost[g].centroid).dot(grid.ghost[g].grad_turb[i]);
			}
		}
		
		rightK_center=grid.ghost[g].k;
		rightK=rightK_center+delta[0];
		rightOmega_center=grid.ghost[g].omega;
		rightOmega=rightOmega_center+delta[1];
		
		rightK=max(rightK,kLowLimit);
		rightOmega=max(rightOmega,omegaLowLimit);
	}
	

	faceGradK-=faceGradK.dot(grid.face[f].normal)*grid.face[f].normal;
	faceGradK+=(rightK_center-leftK_center)/(left2right.dot(grid.face[f].normal))*grid.face[f].normal;
	
	faceGradOmega-=faceGradOmega.dot(grid.face[f].normal)*grid.face[f].normal;
	faceGradOmega+=(rightOmega_center-leftOmega_center)/(left2right.dot(grid.face[f].normal))*grid.face[f].normal;

	faceK=0.5*(leftK+rightK);
	faceOmega=0.5*(leftOmega+rightOmega);
	
	return;	
}

double get_blending(double &k,double &omega,double &rho,double &y,Vec3D &gradK,Vec3D &gradOmega) {
	double F;
	
	cross_diffusion=gradK[0]*gradOmega[0]+gradK[1]*gradOmega[1]+gradK[2]*gradOmega[2];
	CD=max(2.*rho*komega.sigma_omega*cross_diffusion/omega,1.e-20);
	arg1=sqrt(k)/(0.09*omega*y);
	arg2=500.*viscosity/(rho*y*y*omega);
	arg3=4.*rho*komega.sigma_omega*k/(CD*y*y);
	arg=min(max(arg1,arg2),arg3);
	F=tanh(pow(arg,4.));

	return F;
}

void update_eddy_viscosity(void) {
	
	map<int,double>::iterator fit;
	
	for (unsigned int f=0;f<grid.faceCount;++f) {
		parent=grid.face[f].parent;
		// Find face averaged variables
		faceRho=0.; faceK=0.; faceOmega=0.;
		for (fit=grid.face[f].average.begin();fit!=grid.face[f].average.end();fit++) {
			if ((*fit).first>=0) { // if contribution is coming from a real cell
				faceK+=(*fit).second*grid.cell[(*fit).first].k;
				faceOmega+=(*fit).second*grid.cell[(*fit).first].omega;
				faceRho+=(*fit).second*grid.cell[(*fit).first].rho;
			} else { // if contribution is coming from a ghost cell
				faceK+=(*fit).second*grid.ghost[-1*((*fit).first+1)].k;
				faceOmega+=(*fit).second*grid.ghost[-1*((*fit).first+1)].omega;
				faceRho+=(*fit).second*grid.ghost[-1*((*fit).first+1)].rho;
			}
		}
		
		faceK=max(faceK,kLowLimit);
		faceOmega=max(faceOmega,omegaLowLimit);
		
		if (grid.face[f].bc>=0) { // boundary face
			if (bc.region[grid.face[f].bc].type==NOSLIP) {
				mu_t=0.;
			} else if (bc.region[grid.face[f].bc].type==SYMMETRY) {
				mu_t=grid.cell[parent].rho*grid.cell[parent].k/grid.cell[parent].omega;
			} else if (bc.region[grid.face[f].bc].type==SLIP) {
				mu_t=faceRho*faceK/faceOmega;
			} else if (bc.region[grid.face[f].bc].type==INLET) {
				mu_t=faceRho*bc.region[grid.face[f].bc].k/bc.region[grid.face[f].bc].omega;
			} else if (bc.region[grid.face[f].bc].type==OUTLET) {
				mu_t=faceRho*faceK/faceOmega;
			}
		}
		
		grid.face[f].mu_t=min(fabs(mu_t),viscosityRatioLimit*viscosity);
		
	}
	
} // update_eddy_viscosity




