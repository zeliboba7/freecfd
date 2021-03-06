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
#include "rans.h"

extern RANS_Model komega,kepsilon;

void get_kOmega(void);
double get_blending(double &k,double &omega,double &rho,double &y,double &visc,Vec3D &gradK,Vec3D &gradOmega);

int parent,neighbor,f;
double lam_visc,turb_visc;

double leftK,leftOmega, rightK, rightOmega, faceK, faceOmega, faceRho;
Vec3D faceGradK, faceGradOmega, left2right;

double mdot,weightL,weightR;
bool extrapolated;

void RANS::terms(void) {

	double blending;
	double sigma_k,sigma_omega,beta,beta_star,alpha;
	double value;
	double dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz;
	int row,col;
	double convectiveFlux[2],diffusiveFlux[2],source[2];
	double jacL[2],jacR[2];
	double cross_diffusion;
	double closest_wall_distance;

	MatZeroEntries(impOP); // Flush the implicit operator
	VecSet(rhs,0.); // Flush the right hand side

	// Loop through faces
	for (f=0;f<grid[gid].faceCount;++f) {
		extrapolated=false;
		if (model==KOMEGA) blending=1.;
		if (model==KEPSILON) blending=0.;
		
		parent=grid[gid].face[f].parent; neighbor=grid[gid].face[f].neighbor;
		
		lam_visc=ns[gid].material.viscosity(ns[gid].T.face(f));
		turb_visc=mu_t.face(f);
		
		// Fill in the yplus info for output
		int bcno=grid[gid].face[f].bc;
		if (bcno>=0) {
			Vec3D tau=ns[gid].tau.bc(bcno,f);
			double tau_w=fabs((tau-tau.dot(grid[gid].face[f].normal)*grid[gid].face[f].normal));
			double u_star=sqrt(tau_w/ns[gid].rho.face(f));
			double height=(grid[gid].face[f].centroid-grid[gid].cell[parent].centroid).dot(grid[gid].face[f].normal);
			yplus.bc(bcno,f)=ns[gid].rho.face(f)*u_star*height/lam_visc;
		}
		
		// Convective flux is based on mdot calculated through the Riemann solver right after
		// main flow was updated
		
		// The following weights are consistent with the splitting of the Riemann solver.
		mdot=ns[gid].mdot.face(f);
		weightL=ns[gid].weightL.face(f);
		weightR=1.-weightL;

		// Get left, right and face values of k and omega as well as the face normal gradients
		get_kOmega();

		convectiveFlux[0]=mdot*(weightL*leftK+weightR*rightK)*grid[gid].face[f].area;
		convectiveFlux[1]=mdot*(weightL*leftOmega+weightR*rightOmega)*grid[gid].face[f].area;

		// Diffusive k and omega fluxes	
		if (model==BSL || model==SST) {
			closest_wall_distance=grid[gid].cell[parent].closest_wall_distance;
			blending=get_blending(faceK,faceOmega,faceRho,closest_wall_distance,lam_visc,faceGradK,faceGradOmega);
		}
		
		sigma_omega=blending*komega.sigma_omega+(1.-blending)*kepsilon.sigma_omega;
		sigma_k=blending*komega.sigma_k+(1.-blending)*kepsilon.sigma_k;
		alpha=blending*komega.alpha+(1.-blending)*kepsilon.alpha;
		beta=blending*komega.beta+(1.-blending)*kepsilon.beta;
		beta_star=blending*komega.beta_star+(1.-blending)*kepsilon.beta_star;
	
		diffusiveFlux[0]=(lam_visc+turb_visc*sigma_k)*faceGradK.dot(grid[gid].face[f].normal)*grid[gid].face[f].area;
		diffusiveFlux[1]=(lam_visc+turb_visc*sigma_omega)*faceGradOmega.dot(grid[gid].face[f].normal)*grid[gid].face[f].area;

		// Fill in rhs vector for rans scalars
		for (int i=0;i<2;++i) {
			row=(grid[gid].myOffset+parent)*2+i;
			value=diffusiveFlux[i]-convectiveFlux[i];
			VecSetValues(rhs,1,&row,&value,ADD_VALUES);
			if (grid[gid].face[f].bc==INTERNAL_FACE) { 
				row=(grid[gid].myOffset+neighbor)*2+i;
				value*=-1.;
				VecSetValues(rhs,1,&row,&value,ADD_VALUES);
			}
		}
		
		// Calculate flux jacobians
		
		// Assumes k flux doesn't change with omega and vice versa 
		// This is true for convective flux (effect of mu_t in diffusive flux ignored)
		double l2rmag=fabs(left2right);
		double AoverH=grid[gid].face[f].area*grid[gid].face[f].normal.dot(left2right.norm())/(l2rmag);
		// dF_k/dk_left
		jacL[0]=weightL*mdot*grid[gid].face[f].area; // convective
		if (extrapolated) jacL[0]+=weightR*mdot*grid[gid].face[f].area; // convective
		if (!extrapolated) jacL[0]+=(lam_visc+turb_visc*sigma_k)*AoverH; // diffusive

		// dF_omega/dOmega_left
		jacL[1]=weightL*mdot*grid[gid].face[f].area; // convective
		if (extrapolated) jacL[1]+=weightR*mdot*grid[gid].face[f].area; // convective
		if (!extrapolated) jacL[1]+=(lam_visc+turb_visc*sigma_omega)*AoverH; // diffusive

		if (grid[gid].face[f].bc<0) {
			// dF_k/dk_right
			jacR[0]=weightR*mdot*grid[gid].face[f].area; // convective
			jacR[0]-=(lam_visc+turb_visc*sigma_k)*AoverH; // diffusive
			// dF_omega/dOmega_right
			jacR[1]=weightR*mdot*grid[gid].face[f].area; // convective
			jacR[1]-=(lam_visc+turb_visc*sigma_omega)*AoverH; // diffusive
		}
		
		// Insert flux jacobians for the parent cell
		// left_k/left_k
		row=(grid[gid].myOffset+parent)*2; col=row; value=jacL[0];
		MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
		// left_omega/left_omega
		row++; col++; value=jacL[1];
		MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
		if (grid[gid].face[f].bc==INTERNAL_FACE) { 
			// left_k/right_k
			row=(grid[gid].myOffset+parent)*2; col=(grid[gid].myOffset+neighbor)*2; value=jacR[0];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			// left_omega/right_omega
			row++; col++; value=jacR[1];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			
			// Insert flux jacobians for the neighbor cell
			// right_k/right_k
			row=(grid[gid].myOffset+neighbor)*2; col=row; value=-jacR[0];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			// right_omega/right_omega
			row++; col++; value=-jacR[1];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			// right_k/left_k
			row=(grid[gid].myOffset+neighbor)*2; col=(grid[gid].myOffset+parent)*2; value=-jacL[0];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			// right_omega/left_omega
			row++; col++; value=-jacL[1];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
		} else if (grid[gid].face[f].bc==PARTITION_FACE) { 
			// left_k/right_k
			row=(grid[gid].myOffset+parent)*2; col=(grid[gid].cell[neighbor].matrix_id)*2; value=jacR[0];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			// left_omega/right_omega
			row++; col++; value=jacR[1];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
		}

	} // for faces
	
	// Now do a cell loop to add the unsteady and source terms
	double divU,Prod_k,Dest_k;
	for (int c=0;c<grid[gid].cell.size();++c) {

		if (model==BSL || model==SST) {
			blending=get_blending(k.cell(c),omega.cell(c),ns[gid].rho.cell(c),grid[gid].cell[c].closest_wall_distance,lam_visc,gradk.cell(c),gradomega.cell(c));
		}

		sigma_omega=blending*komega.sigma_omega+(1.-blending)*kepsilon.sigma_omega;
		sigma_k=blending*komega.sigma_k+(1.-blending)*kepsilon.sigma_k;
		alpha=blending*komega.alpha+(1.-blending)*kepsilon.alpha;
		beta=blending*komega.beta+(1.-blending)*kepsilon.beta;
		beta_star=blending*komega.beta_star+(1.-blending)*kepsilon.beta_star;
		
		// Calculate the source terms
			
		dudx=ns[gid].gradu.cell(c)[0];
		dudy=ns[gid].gradu.cell(c)[1];
		dudz=ns[gid].gradu.cell(c)[2];
		
		dvdx=ns[gid].gradv.cell(c)[0];
		dvdy=ns[gid].gradv.cell(c)[1];
		dvdz=ns[gid].gradv.cell(c)[2];
		
		dwdx=ns[gid].gradw.cell(c)[0];
		dwdy=ns[gid].gradw.cell(c)[1];
		dwdz=ns[gid].gradw.cell(c)[2];
		
		divU=dudx+dvdy+dwdz;
		// sr or strainRate variable here is actually the term: sqrt(2*S_ij*S_ij)
		// Where S_ij=1/2*(du_i/dx_j+du_j/dx_i) is the proper strain rate tensor
		double sr=sqrt(2.*(dudx*dudx + dvdy*dvdy + dwdz*dwdz)
				+ (dudy+dvdx)*(dudy+dvdx) + (dudz+dwdx)*(dudz+dwdx)
				+ (dvdz+dwdy)*(dvdz+dwdy)
			       );
		
		// Store this for sst model
		strainRate.cell(c)=sr;
	
		if (c>=grid[gid].cellCount) continue;
	
		Prod_k=mu_t.cell(c)*(sr*sr-2./3.*divU*divU)-2./3.*ns[gid].rho.cell(c)*k.cell(c)*divU;
		Dest_k=beta_star*ns[gid].rho.cell(c)*omega.cell(c)*k.cell(c);
		
		// Limit production of k
		Prod_k=min(Prod_k,10.*Dest_k);
		
		cross_diffusion=0.;
		if (model!=KOMEGA) cross_diffusion=gradk.cell(c).dot(gradomega.cell(c));
			
		source[0]=Prod_k; // k production
		source[0]-=Dest_k; // k destruction
		
		source[1]=alpha*ns[gid].rho.cell(c)/mu_t.cell(c)*Prod_k; // omega production
		source[1]-=beta*ns[gid].rho.cell(c)*omega.cell(c)*omega.cell(c); // omega destruction
				
		source[1]+=2.*(1.-blending)*ns[gid].rho.cell(c)*komega.sigma_omega*cross_diffusion/omega.cell(c);
		
		// Add source terms to rhs
		for (int i=0;i<2;++i) {
			row=(grid[gid].myOffset+c)*2+i;
			value=source[i]*grid[gid].cell[c].volume;
			VecSetValues(rhs,1,&row,&value,ADD_VALUES);
		}
		
		// Add source jacobians
		// Only include destruction terms
		
		// dS_k/dk
		row=(grid[gid].myOffset+c)*2; col=row;
		value=beta_star*ns[gid].rho.cell(c)*omega.cell(c)*grid[gid].cell[c].volume; // approximate 
		MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
		
		// dS_k/dOmega
		col++; 
		value=beta_star*ns[gid].rho.cell(c)*k.cell(c)*grid[gid].cell[c].volume; 
		MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
		
		// dS_omega/dOmega
		row++; 
		value=2.*beta*ns[gid].rho.cell(c)*omega.cell(c)*grid[gid].cell[c].volume; // approximate
		// Add cross-diffusion term jacobian
		value+=2.*(1.-blending)*komega.sigma_omega*komega.sigma_omega*cross_diffusion/(omega.cell(c)*omega.cell(c))*grid[gid].cell[c].volume;
		MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);

		
	} // end cell loop
	
	MatAssemblyBegin(impOP,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(impOP,MAT_FINAL_ASSEMBLY);

	VecAssemblyBegin(rhs);
	VecAssemblyEnd(rhs);
	
	return;
} // end RANS::terms

void RANS::get_kOmega() {
	
	leftK=k.cell(parent);
	leftOmega=omega.cell(parent);
	rightK=k.cell(neighbor);
	rightOmega=omega.cell(neighbor);
	
	int bcno=grid[gid].face[f].bc;
	if (bcno>=0) {
		if (bc[gid][bcno].type==SYMMETRY || bc[gid][bcno].type==OUTLET) {
			extrapolated=true;
		} else if (bc[gid][bcno].type==WALL && bc[gid][bcno].kind==SLIP) {
			extrapolated=true;
		}
		faceGradK=gradk.cell(parent);
		faceGradOmega=gradomega.cell(parent);
	} else {
		faceGradK=gradk.face(f);
		faceGradOmega=gradomega.face(f);
	}

	faceK=k.face(f);
	faceOmega=omega.face(f);
	faceRho=ns[gid].rho.face(f);

	// Find distance between left and right centroids 
	Vec3D leftCentroid,rightCentroid;
	leftCentroid=grid[gid].cell[parent].centroid;
	rightCentroid=grid[gid].cell[neighbor].centroid;
	left2right=rightCentroid-leftCentroid;
	
	Vec3D l2rnormal=left2right;
	l2rnormal=l2rnormal.norm();
	double l2rmag=fabs(left2right);
	
	faceGradK-=faceGradK.dot(l2rnormal)*l2rnormal;
	faceGradK+=((rightK-leftK)/(l2rmag))*l2rnormal;
	
	faceGradOmega-=faceGradOmega.dot(l2rnormal)*l2rnormal;
	faceGradOmega+=((rightOmega-leftOmega)/(l2rmag))*l2rnormal;

	return;	
}

double get_blending(double &k,double &omega,double &rho,double &y,double &visc,Vec3D &gradK,Vec3D &gradOmega) {
	double arg,arg1,arg2,arg3,CD,cross_diff,F1;
	
	cross_diff=gradK[0]*gradOmega[0]+gradK[1]*gradOmega[1]+gradK[2]*gradOmega[2];
	CD=max(2.*rho*kepsilon.sigma_omega*cross_diff/omega,1.e-10);
	arg1=2.*sqrt(k)/(kepsilon.beta_star*omega*y);
	arg2=500.*visc/(rho*y*y*omega);
	arg3=4.*rho*kepsilon.sigma_omega*k/(CD*y*y);
	arg=min(max(arg1,arg2),arg3);
	F1=tanh(pow(arg,4.));
	return F1;
}

