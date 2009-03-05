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
#include "rans.h"
extern RANS rans;
double EPS=1.e-10;

void get_kOmega(void);
double get_blending(double &k,double &omega,double &rho,double &y,Vec3D &gradK,Vec3D &gradOmega);

unsigned int parent,neighbor,f;
double leftK,leftOmega, rightK, rightOmega, faceK, faceOmega, faceRho;
double leftK_center, rightK_center, leftOmega_center, rightOmega_center;
Vec3D faceGradK, faceGradOmega, left2right;
double mu_t=0.;
double weightL,weightR;
bool extrapolated;

void RANS::terms(void) {

	double blending;
	double sigma_k,sigma_omega,beta,beta_star,kappa,alpha;
	double value;
	double dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz;
	double strainRate;
	int row,col;
	double convectiveFlux[2],diffusiveFlux[2],source[2];
	double jacL[2],jacR[2];
	double cross_diffusion;
	double closest_wall_distance;
	
	MatZeroEntries(impOP); // Flush the implicit operator

	// Loop through faces
	for (f=0;f<grid.faceCount;++f) {
		if (TURBULENCE_MODEL==KOMEGA) blending=1.;
		if (TURBULENCE_MODEL==KEPSILON) blending=0.;
		
		parent=grid.face[f].parent; neighbor=grid.face[f].neighbor;
				
		// Convective flux is based on mdot calculated through the Riemann solver right after
		// main flow was updated
		
		// These weights are coming from the Riemann solver
		weightL=grid.face[f].weightL;
		extrapolated=false;
		// Get left, right and face values of k and omega as well as the face normal gradients
		get_kOmega();
		weightR=1.-weightL;

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
		
		// Fill in rhs vector for rans scalars
		for (int i=0;i<2;++i) {
			row=(grid.myOffset+parent)*2+i;
			value=diffusiveFlux[i]-convectiveFlux[i];
			VecSetValues(rhs,1,&row,&value,ADD_VALUES);
			if (grid.face[f].bc==INTERNAL) { 
				row=(grid.myOffset+neighbor)*2+i;
				value*=-1.;
				VecSetValues(rhs,1,&row,&value,ADD_VALUES);
			}
		}
		
		// Calculate flux jacobians
		
		// Assumes k flux doesn't change with omega and vice versa 
		// This is true for convective flux (effect of mu_t in diffusive flux ignored)

		double leftRho=grid.cell[parent].rho;
		double rightRho=1.; // if right cell is not internal, jacobian is not inserted into the matrix, so this doesn't matter
		if (grid.face[f].bc==INTERNAL) rightRho=grid.cell[neighbor].rho;
		// dF_k/dk_left
		jacL[0]=weightL*grid.face[f].mdot*grid.face[f].area; // convective
		if (!extrapolated) jacL[0]+=(viscosity+mu_t*sigma_k)/(left2right.dot(grid.face[f].normal))*grid.face[f].area; // diffusive
		// dF_k/dk_right
		jacR[0]=weightR*grid.face[f].mdot*grid.face[f].area; // convective
		if (!extrapolated) jacR[0]-=(viscosity+mu_t*sigma_k)/(left2right.dot(grid.face[f].normal))*grid.face[f].area; // diffusive
		// dF_omega/dOmega_left
		jacL[1]=weightL*grid.face[f].mdot*grid.face[f].area; // convective
		if (!extrapolated) jacL[1]+=(viscosity+mu_t*sigma_omega)/(left2right.dot(grid.face[f].normal))*grid.face[f].area; // diffusive
		// dF_omega/dOmega_right
		jacR[1]=weightR*grid.face[f].mdot*grid.face[f].area; // convective
		if (!extrapolated) jacR[1]-=(viscosity+mu_t*sigma_omega)/(left2right.dot(grid.face[f].normal))*grid.face[f].area; // diffusive
		
		// Insert flux jacobians for the parent cell
		// left_k/left_k
		row=(grid.myOffset+parent)*2; col=row; value=jacL[0];
		MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
		// left_omega/left_omega
		row++; col++; value=jacL[1];
		MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
		if (grid.face[f].bc==INTERNAL) { 
			// left_k/right_k
			row=(grid.myOffset+parent)*2; col=(grid.myOffset+neighbor)*2; value=jacR[0];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			// left_omega/right_omega
			row++; col++; value=jacR[1];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			
			// Insert flux jacobians for the neighbor cell
			// right_k/right_k
			row=(grid.myOffset+neighbor)*2; col=row; value=-jacR[0];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			// right_omega/right_omega
			row++; col++; value=-jacR[1];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			// right_k/left_k
			row=(grid.myOffset+neighbor)*2; col=(grid.myOffset+parent)*2; value=-jacL[0];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			// right_omega/left_omega
			row++; col++; value=-jacL[1];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
		} else if (grid.face[f].bc==GHOST) { 
			// left_k/right_k
			row=(grid.myOffset+parent)*2; col=(grid.ghost[-1*neighbor-1].matrix_id)*2; value=jacR[0];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			// left_omega/right_omega
			row++; col++; value=jacR[1];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
		}

	} // for faces
	
	// Now to a cell loop to add the unsteady and source terms
	double divU,Prod_k;
	int count_prod_k_clip=0;
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
		MatSetValues(impOP,1,&row,1,&row,&value,ADD_VALUES);
		row++;
		MatSetValues(impOP,1,&row,1,&row,&value,ADD_VALUES);

		if (TURBULENCE_MODEL==BSL) {
			blending=get_blending(cell[c].k,cell[c].omega,grid.cell[c].rho,grid.cell[c].closest_wall_distance,cell[c].grad[0],cell[c].grad[1]);
		}
		
		sigma_omega=blending*komega.sigma_omega+(1.-blending)*kepsilon.sigma_omega;
		sigma_k=blending*komega.sigma_k+(1.-blending)*kepsilon.sigma_k;
		alpha=blending*komega.alpha+(1.-blending)*kepsilon.alpha;
		beta=blending*komega.beta+(1.-blending)*kepsilon.beta;
		beta_star=blending*komega.beta_star+(1.-blending)*kepsilon.beta_star;
		
		// Calculate the source terms
		mu_t=grid.cell[c].rho*cell[c].k/cell[c].omega;
		
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
		
		Prod_k=mu_t*(strainRate*strainRate-2./3.*divU*divU)-2./3.*grid.cell[c].rho*cell[c].k*divU;

		if (Prod_k>10./grid.cell[c].volume) {
			count_prod_k_clip++;
			Prod_k=10./grid.cell[c].volume;
		}
		
		cross_diffusion=0.;
		if (TURBULENCE_MODEL!=KOMEGA) cross_diffusion=cell[c].grad[0].dot(cell[c].grad[1]);
			
		source[0]=Prod_k; // k production
		source[0]-=beta_star*grid.cell[c].rho*cell[c].omega*cell[c].k; // k destruction
		
		source[1]=alpha/(mu_t/grid.cell[c].rho+EPS)*Prod_k; // omega production
		source[1]-=beta*grid.cell[c].rho*cell[c].omega*cell[c].omega; // omega destruction
				
		source[1]+=2.*(1.-blending)*grid.cell[c].rho*komega.sigma_omega*cross_diffusion/cell[c].omega;
		
		// Add source terms to rhs
		for (int i=0;i<2;++i) {
			row=(grid.myOffset+c)*2+i;
			value=source[i]*grid.cell[c].volume;
			VecSetValues(rhs,1,&row,&value,ADD_VALUES);
		}
		
		// Add source jacobians
		// Only include destruction terms
		
		// dS_k/dk
		row=(grid.myOffset+c)*2; col=row;
		value=beta_star*grid.cell[c].rho*cell[c].omega*grid.cell[c].volume; // approximate 
		MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
		
		// dS_k/dOmega
		col++; 
		value=beta_star*grid.cell[c].rho*cell[c].k*grid.cell[c].volume; 
		MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
		
		// dS_omega/dOmega
		row++; 
		value=2.*beta*grid.cell[c].rho*cell[c].omega*grid.cell[c].volume; // approximate
		// Add cross-diffusion term jacobian
		value+=2.*(1.-blending)*komega.sigma_omega*komega.sigma_omega*cross_diffusion/(cell[c].omega*cell[c].omega)*grid.cell[c].volume;
		MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
	
		
	} // end cell loop
	
	if (count_prod_k_clip>0) {
		cout << "[I Rank=" << Rank << "] k-production term is clipped for " <<  count_prod_k_clip << " cells" << endl;
	}
	
	return;
} // end RANS::terms

void get_kOmega() {
	
	// Get left k and omega
	double delta[2];

	if (order==SECOND) {
		for (unsigned int i=0;i<2;++i) {
			delta[i]=(grid.face[f].centroid-grid.cell[parent].centroid).dot(rans.cell[parent].grad[i]);
		}
	} else {
		for (unsigned int i=0;i<2;++i) delta[i]=0.;
	}
	
	leftK_center=rans.cell[parent].k;
	leftK=leftK_center+delta[0];
	leftOmega_center=rans.cell[parent].omega;
	leftOmega=leftOmega_center+delta[1];
	
	leftK=max(leftK,kLowLimit);
	leftK=min(leftK,kHighLimit);
	leftOmega=max(leftOmega,omegaLowLimit);

	// Find face averaged quantities
	faceRho=0.; faceGradK=0.; faceGradOmega=0.;
 	map<int,double>::iterator fit;
	for (fit=grid.face[f].average.begin();fit!=grid.face[f].average.end();fit++) {
		if ((*fit).first>=0) { // if contribution is coming from a real cell
			faceRho+=(*fit).second*grid.cell[(*fit).first].rho;
			faceGradK+=(*fit).second*rans.cell[(*fit).first].grad[0];
			faceGradOmega+=(*fit).second*rans.cell[(*fit).first].grad[1];
		} else { // if contribution is coming from a ghost cell
			faceRho+=(*fit).second*grid.ghost[-1*((*fit).first+1)].rho;
			faceGradK+=(*fit).second*rans.ghost[-1*((*fit).first+1)].grad[0];
			faceGradOmega+=(*fit).second*rans.ghost[-1*((*fit).first+1)].grad[1];
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
				delta[i]=(grid.face[f].centroid-grid.cell[neighbor].centroid).dot(rans.cell[neighbor].grad[i]);
			}
		} else {
			for (unsigned int i=0;i<2;++i) delta[i]=0.;
		}

		rightK_center=rans.cell[neighbor].k;
		rightK=rightK_center+delta[0];
		rightOmega_center=rans.cell[neighbor].omega;
		rightOmega=rightOmega_center+delta[1];
	
		rightK=max(rightK,kLowLimit);
		rightK=min(rightK,kHighLimit);
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
			extrapolated=true;
			weightL=1.;// This is true when rightK=leftK and rightOmega=leftOmega
		} else if (bc.region[grid.face[f].bc].type==SLIP) {
			rightK_center=leftK_center; 
			rightOmega_center=leftOmega_center;
			extrapolated=true;
			weightL=1.;
		} else if (bc.region[grid.face[f].bc].type==INLET) {
			rightK=bc.region[grid.face[f].bc].k;
			rightOmega=bc.region[grid.face[f].bc].omega;
			rightK_center=2.*rightK-leftK_center;
			rightOmega_center=2.*rightOmega-leftOmega_center;
			//weightL=0.; 
		} else if (bc.region[grid.face[f].bc].type==OUTLET) {
			rightK_center=leftK_center; 
			rightOmega_center=leftOmega_center;
			extrapolated=true;
			weightL=1.;
		}
		
		
		
	} else { // partition boundary

		int g=-1*grid.face[f].neighbor-1; // ghost cell index

		if (order==SECOND) {
			for (unsigned int i=0;i<2;++i) {
				delta[i]=(grid.face[f].centroid-grid.ghost[g].centroid).dot(rans.ghost[g].grad[i]);
			}
		}
		
		rightK_center=rans.ghost[g].k;
		rightK=rightK_center+delta[0];
		rightOmega_center=rans.ghost[g].omega;
		rightOmega=rightOmega_center+delta[1];
		
		rightK=max(rightK,kLowLimit);
		rightK=min(rightK,kHighLimit);
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
	double F,arg,arg1,arg2,arg3,CD,cross_diff;
	
	cross_diff=gradK[0]*gradOmega[0]+gradK[1]*gradOmega[1]+gradK[2]*gradOmega[2];
	CD=max(2.*rho*rans.komega.sigma_omega*cross_diff/omega,1.e-20);
	arg1=sqrt(k)/(0.09*omega*y);
	arg2=500.*viscosity/(rho*y*y*omega);
	arg3=4.*rho*rans.komega.sigma_omega*k/(CD*y*y);
	arg=min(max(arg1,arg2),arg3);
	F=tanh(pow(arg,4.));

	return F;
}

