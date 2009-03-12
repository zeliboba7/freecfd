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
#include "flamelet.h"
#include "rans.h"
extern Flamelet flamelet;
extern RANS rans;

void Flamelet::terms(void) {

	unsigned int parent,neighbor,f;
	double leftZ,leftZvar,rightZ,rightZvar;
	Vec3D faceGradZ, faceGradZvar, left2right;
	double mu_t=0.;
	double weightL,weightR;
	double weightL1,weightR1;
	double value;
	int row,col;
	double convectiveFlux[2],diffusiveFlux[2],source[2];
	double jacL[2],jacR[2];
	bool extrapolated;
	double mu,diff;
	map<int,double>::iterator fit;
	
	MatZeroEntries(impOP); // Flush the implicit operator
	VecSet(rhs,0.); // Flush the right hand side

	// Loop through faces
	Vec3D sum=0.;
	for (f=0;f<grid.faceCount;++f) {
		
		extrapolated=false;
		
		parent=grid.face[f].parent; neighbor=grid.face[f].neighbor;
		
		// Get left, right and face values of Z and Zvar as well as the face normal gradients
		get_Z_Zvar(parent,neighbor,f,
			   leftZ,leftZvar,rightZ,rightZvar,
      			   faceGradZ,faceGradZvar,left2right,weightL,extrapolated);
		
		// The following weights are consistent with rho splitting of the AUSM scheme.
		// The idea is that Z behaves a lot like rho
		// This is the only stable configuration out of many others tried.
		weightL1=double(grid.face[f].upstream);
		weightR1=1.-weightL1;
		
		// Convective flux is based on the mdot calculated through the Riemann solver right after
		// main flow was updated
		convectiveFlux[0]=grid.face[f].mdot*(weightL1*leftZ+weightR1*rightZ)*grid.face[f].area;
		convectiveFlux[1]=grid.face[f].mdot*(weightL1*leftZvar+weightR1*rightZvar)*grid.face[f].area;

		mu=flamelet.face[f].mu;
		mu_t=rans.face[f].mu_t;
		
		diff=0.;
		for (fit=grid.face[f].average.begin();fit!=grid.face[f].average.end();fit++) {
			if ((*fit).first>=0) { // if contribution is coming from a real cell
				diff+=(*fit).second*cell[(*fit).first].diffusivity;
			} else { // if contribution is coming from a ghost cell
				diff+=(*fit).second*ghost[-1*((*fit).first+1)].diffusivity;
			}
		}
		
		diffusiveFlux[0]=(diff+mu_t/constants.sigma_t)*faceGradZ.dot(grid.face[f].normal)*grid.face[f].area;
		diffusiveFlux[1]=(diff+mu_t/constants.sigma_t)*faceGradZvar.dot(grid.face[f].normal)*grid.face[f].area;
		
		// Fill in rhs vector
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
		// dF_Z/dZ_left
		jacL[0]=weightL1*grid.face[f].mdot*grid.face[f].area; // convective
		if (!extrapolated) jacL[0]+=(diff+mu_t/constants.sigma_t)/(left2right.dot(grid.face[f].normal))*grid.face[f].area; // diffusive
		// dF_Z/dZ_right
		jacR[0]=weightR1*grid.face[f].mdot*grid.face[f].area; // convective
		if (!extrapolated) jacR[0]-=(diff+mu_t/constants.sigma_t)/(left2right.dot(grid.face[f].normal))*grid.face[f].area; // diffusive
		// dF_Zvar/dZvar_left
		jacL[1]=weightL1*grid.face[f].mdot*grid.face[f].area; // convective
		if (!extrapolated) jacL[1]+=(diff+mu_t/constants.sigma_t)/(left2right.dot(grid.face[f].normal))*grid.face[f].area; // diffusive
		// dF_Zvar/dZvar_right
		jacR[1]=weightR1*grid.face[f].mdot*grid.face[f].area; // convective
		if (!extrapolated) jacR[1]-=(diff+mu_t/constants.sigma_t)/(left2right.dot(grid.face[f].normal))*grid.face[f].area; // diffusive
		
		// Insert flux jacobians for the parent cell
		// left_Z/left_Z
		row=(grid.myOffset+parent)*2; col=row; value=jacL[0];
		MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
		// left_Zvar/left_Zvar
		row++; col++; value=jacL[1];
		MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
		if (grid.face[f].bc==INTERNAL) { 
			// left_Z/right_Z
			row=(grid.myOffset+parent)*2; col=(grid.myOffset+neighbor)*2; value=jacR[0];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			// left_Zvar/right_Zvar
			row++; col++; value=jacR[1];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			
			// Insert flux jacobians for the neighbor cell
			// right_Z/right_Z
			row=(grid.myOffset+neighbor)*2; col=row; value=-jacR[0];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			// right_Zvar/right_Zvar
			row++; col++; value=-jacR[1];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			// right_Z/left_Z
			row=(grid.myOffset+neighbor)*2; col=(grid.myOffset+parent)*2; value=-jacL[0];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			// right_Zvar/left_Zvar
			row++; col++; value=-jacL[1];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
		} else if (grid.face[f].bc==GHOST) { 
			// left_Z/right_Z
			row=(grid.myOffset+parent)*2; col=(grid.ghost[-1*neighbor-1].matrix_id)*2; value=jacR[0];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			// left_Zvar/right_Zvar
			row++; col++; value=jacR[1];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
		}

	} // for faces
	
	// Now do a cell loop to add the unsteady and source terms
	double Prod_Zvar,Dest_Zvar;
	double drho_dZ,drho_dZvar;
	double Zplus,ZvarPlus,Chi,ChiPlus;
	for (unsigned int c=0;c<grid.cellCount;++c) {

		// Unsteady Z term should be:
		// (rho+Z*drho/dZ)dZ/dt
		
		Chi=2.*rans.cell[c].omega*cell[c].Zvar*rans.kepsilon.beta_star;
		if (cell[c].Z<(1.-sqrt_machine_error)){
			// forward differencing
			Zplus=cell[c].Z+sqrt_machine_error;
			drho_dZ=(table.get_rho(Zplus,cell[c].Zvar,Chi)-grid.cell[c].rho)/sqrt_machine_error;
				 
		} else {
			// backward differencing
			Zplus=cell[c].Z-sqrt_machine_error;
			drho_dZ=(grid.cell[c].rho-table.get_rho(Zplus,cell[c].Zvar,Chi))/sqrt_machine_error;
		}
		
		
		if (cell[c].Zvar<(0.25-sqrt_machine_error)){
			// forward differencing
			ZvarPlus=cell[c].Zvar+sqrt_machine_error;
			ChiPlus=2.*rans.cell[c].omega*ZvarPlus*rans.kepsilon.beta_star;
			drho_dZvar=(table.get_rho(cell[c].Z,ZvarPlus,ChiPlus)-grid.cell[c].rho)/sqrt_machine_error;
				 
		} else {
			// backward differencing
			ZvarPlus=cell[c].Zvar-sqrt_machine_error;
			ChiPlus=2.*rans.cell[c].omega*ZvarPlus*rans.kepsilon.beta_star;
			drho_dZvar=(grid.cell[c].rho-table.get_rho(cell[c].Z,ZvarPlus,ChiPlus))/sqrt_machine_error;
		}
		
		// Insert unsteady term
		row=(grid.myOffset+c)*2;
		value=grid.cell[c].volume*(grid.cell[c].rho+cell[c].Z*drho_dZ)/(grid.cell[c].dt*relaxation);
		MatSetValues(impOP,1,&row,1,&row,&value,ADD_VALUES);
		value=grid.cell[c].volume*(grid.cell[c].rho+cell[c].Zvar*drho_dZvar)/(grid.cell[c].dt*relaxation);
		row++;
		MatSetValues(impOP,1,&row,1,&row,&value,ADD_VALUES);
		
		// Calculate the source terms
		mu_t=rans.cell[c].mu_t;
		
		Prod_Zvar=constants.Cg*mu_t*cell[c].grad[0].dot(cell[c].grad[0]);
		Dest_Zvar=constants.Cd*grid.cell[c].rho*rans.kepsilon.beta_star*rans.cell[c].omega*cell[c].Zvar;

			
		source[0]=0.; // Z sources
		
		source[1]=Prod_Zvar-Dest_Zvar; // Zvar sources
		
		// Add source terms to rhs
		for (int i=0;i<2;++i) {
			row=(grid.myOffset+c)*2+i;
			value=source[i]*grid.cell[c].volume;
			VecSetValues(rhs,1,&row,&value,ADD_VALUES);
		}
		
		// Add source jacobians
		
		// Only includr the destruction term
		// dS_Zvar/dZvar
		row=(grid.myOffset+c)*2+1; col=row;
		value=(constants.Cd*grid.cell[c].rho*rans.kepsilon.beta_star*rans.cell[c].omega)*grid.cell[c].volume;
		MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
	
		
	} // end cell loop
	
	return;
} // end Flamelet::terms

void Flamelet::get_Z_Zvar(unsigned int &parent,unsigned int &neighbor,unsigned int &f,
			  double &leftZ,double &leftZvar,
			  double &rightZ,double &rightZvar,
    			  Vec3D &faceGradZ,Vec3D &faceGradZvar,Vec3D &left2right,
	 		  double &weightL,bool &extrapolated) {
	
	double faceZ,faceZvar;
	double leftZ_center,rightZ_center,leftZvar_center,rightZvar_center;
	
	// Get left Z and Zvar
	double delta[2];

	if (order==SECOND) {
		for (unsigned int i=0;i<2;++i) {
			delta[i]=(grid.face[f].centroid-grid.cell[parent].centroid).dot(cell[parent].grad[i]);
		}
	} else {
		for (unsigned int i=0;i<2;++i) delta[i]=0.;
	}
	
	leftZ_center=cell[parent].Z;
	leftZ=leftZ_center+delta[0];
	leftZvar_center=cell[parent].Zvar;
	leftZvar=leftZvar_center+delta[1];

	leftZ=max(0.,leftZ);
	leftZ=min(1.,leftZ);
	leftZvar=max(0.,leftZvar);
	leftZvar=min(0.25,leftZvar);
	
	// Find face averaged quantities
	faceGradZ=0.; faceGradZvar=0.;
 	map<int,double>::iterator fit;
	for (fit=grid.face[f].average.begin();fit!=grid.face[f].average.end();fit++) {
		if ((*fit).first>=0) { // if contribution is coming from a real cell
			faceGradZ+=(*fit).second*cell[(*fit).first].grad[0];
			faceGradZvar+=(*fit).second*cell[(*fit).first].grad[1];
		} else { // if contribution is coming from a ghost cell
			faceGradZ+=(*fit).second*ghost[-1*((*fit).first+1)].grad[0];
			faceGradZvar+=(*fit).second*ghost[-1*((*fit).first+1)].grad[1];
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
		
	// Get right Z and Zvar
	if (grid.face[f].bc==INTERNAL) {// internal face

		if (order==SECOND) {
			for (unsigned int i=0;i<2;++i) {
				delta[i]=(grid.face[f].centroid-grid.cell[neighbor].centroid).dot(cell[neighbor].grad[i]);
			}
		} else {
			for (unsigned int i=0;i<2;++i) delta[i]=0.;
		}

		rightZ_center=cell[neighbor].Z;
		rightZ=rightZ_center+delta[0];
		rightZvar_center=cell[neighbor].Zvar;
		rightZvar=rightZvar_center+delta[1];
		
		rightZ=max(0.,rightZ);
		rightZ=min(1.,rightZ);
		rightZvar=max(0.,rightZvar);
		rightZvar=min(0.25,rightZvar);
		
	} else if (grid.face[f].bc>=0) { // boundary face
		// Default is zero normal
		rightZ=leftZ; rightZvar=leftZvar;
		rightZ_center=leftZ_center;
		rightZvar_center=leftZvar_center;
		extrapolated=true;

		if (bc.region[grid.face[f].bc].type==INLET) {
			rightZ=bc.region[grid.face[f].bc].Z;
			rightZvar=bc.region[grid.face[f].bc].Zvar;
			rightZ_center=2.*rightZ-leftZ_center;
			rightZvar_center=2.*rightZvar-leftZvar_center;
		} 
		
	} else { // partition boundary

		int g=-1*grid.face[f].neighbor-1; // ghost cell index

		if (order==SECOND) {
			for (unsigned int i=0;i<2;++i) {
				delta[i]=(grid.face[f].centroid-grid.ghost[g].centroid).dot(ghost[g].grad[i]);
			}
		}
		
		rightZ_center=ghost[g].Z;
		rightZ=rightZ_center+delta[0];
		rightZvar_center=ghost[g].Zvar;
		rightZvar=rightZvar_center+delta[1];
		
		rightZ=max(0.,rightZ);
		rightZ=min(1.,rightZ);
		rightZvar=max(0.,rightZvar);
		rightZvar=min(0.25,rightZvar);

	}
	

	faceGradZ-=faceGradZ.dot(grid.face[f].normal)*grid.face[f].normal;
	faceGradZ+=(rightZ_center-leftZ_center)/(left2right.dot(grid.face[f].normal))*grid.face[f].normal;
	
	faceGradZvar-=faceGradZvar.dot(grid.face[f].normal)*grid.face[f].normal;
	faceGradZvar+=(rightZvar_center-leftZvar_center)/(left2right.dot(grid.face[f].normal))*grid.face[f].normal;

	faceZ=0.5*(leftZ+rightZ);
	faceZvar=0.5*(leftZvar+rightZvar);
	
	return;	
}

