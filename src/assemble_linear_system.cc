/************************************************************************
	
	Copyright 2007-2008 Emre Sozer & Patrick Clark Trizila

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
#include "state_cache.h"
#include "petsc_functions.h"

extern BC bc;
extern InputFile input;

void get_jacobians(const int var);
void convective_face_flux(Cell_State &left,Cell_State &right,Face_State &face,double flux[]);
void diffusive_face_flux(Cell_State &left,Cell_State &right,Face_State &face,double flux[]);
void left_state_update(Cell_State &left,Face_State &face);
void right_state_update(Cell_State &left,Cell_State &right,Face_State &face);
void face_geom_update(Face_State &face,unsigned int f);
void face_state_update(Cell_State &left,Cell_State &right,Face_State &face);
void state_perturb(Cell_State &state,Face_State &face,int var,double epsilon);
void face_state_adjust(Cell_State &left,Cell_State &right,Face_State &face,int var);
void sources(Cell_State &state,double source[],bool forJacobian=false);

namespace state {
	Cell_State left,right,leftPlus,rightPlus;
	Face_State face;
	Fluxes flux, fluxPlus;
	vector<double> jacobianLeft,jacobianRight,sourceJacLeft,sourceJacRight;
	vector<double> sourceLeft,sourceRight,sourceLeftPlus,sourceRightPlus;
	bool doLeftSourceJac,doRightSourceJac;
}

void assemble_linear_system(void) {

	using namespace state;
	using state::left;
	using state::right;
	
	unsigned int parent,neighbor,f;
	int row,col;
	
	vector<bool> cellVisited;
	for (unsigned int c=0;c<grid.cellCount;++c) cellVisited.push_back(false);
	PetscScalar value;

	flux.convective.resize(nSolVar);
	flux.diffusive.resize(nSolVar);
	fluxPlus.convective.resize(nSolVar);
	fluxPlus.diffusive.resize(nSolVar);
	jacobianLeft.resize(nSolVar);
	jacobianRight.resize(nSolVar);
	left.update.resize(nSolVar);
	right.update.resize(nSolVar);
	leftPlus.update.resize(nSolVar);
	rightPlus.update.resize(nSolVar);
	sourceLeft.resize(nSolVar);
	sourceRight.resize(nSolVar);
	sourceLeftPlus.resize(nSolVar);
	sourceRightPlus.resize(nSolVar);
	sourceJacLeft.resize(nSolVar);
	sourceJacRight.resize(nSolVar);
	
	for (int m=0;m<nSolVar;++m) flux.diffusive[m]=0.;
	
	bool implicit=true;
	if (TIME_INTEGRATOR==FORWARD_EULER) implicit=false;
	// Loop through faces
	for (f=0;f<grid.faceCount;++f) {

		doLeftSourceJac=false; doRightSourceJac=false;
		for (int m=0;m<nSolVar;++m) { 
			sourceLeft[m]=0.;
			sourceRight[m]=0.;
		}
		parent=grid.face[f].parent; neighbor=grid.face[f].neighbor;

		// Populate the state caches
		face_geom_update(face,f);
		left_state_update(left,face);
		right_state_update(left,right,face);
		if (EQUATIONS==NS) face_state_update(left,right,face);

		// Get unperturbed flux values
		convective_face_flux(left,right,face,&flux.convective[0]);
		if (EQUATIONS==NS) diffusive_face_flux(left,right,face,&flux.diffusive[0]);

		// Add Sources
		if (!cellVisited[parent]){
			sources(left,&sourceLeft[0]);
			cellVisited[parent]=true;
			doLeftSourceJac=true;
		}
		
		if(face.bc==INTERNAL && !cellVisited[neighbor]) {
			sources(right,&sourceRight[0]);
			cellVisited[neighbor]=true;
			doRightSourceJac=true;
		}
		
		// Integrate boundary fluxes
		if ((timeStep) % integrateBoundaryFreq == 0) {
			if (face.bc>=0) {
				bc.region[face.bc].mass+=flux.convective[0];
				for (int i=0;i<3;++i) bc.region[face.bc].momentum[i]+=(flux.convective[i+1]-flux.diffusive[i+1]);
				bc.region[face.bc].energy+=(flux.convective[4]-flux.diffusive[4]);
			}
		}

		// Fill in rhs vector
		for (int i=0;i<nSolVar;++i) {
			row=(grid.myOffset+parent)*nSolVar+i;
			value=flux.diffusive[i]-flux.convective[i]+sourceLeft[i];
			VecSetValues(rhs,1,&row,&value,ADD_VALUES);
			if (grid.face[f].bc==INTERNAL) { 
				row=(grid.myOffset+neighbor)*nSolVar+i;
				value=-1.*(flux.diffusive[i]-flux.convective[i])+sourceRight[i];
				VecSetValues(rhs,1,&row,&value,ADD_VALUES);
			}
		}

		if (implicit) {

			for (int i=0;i<nSolVar;++i) {
				for (int m=0;m<nSolVar;++m) { 
					sourceJacLeft[m]=0.;
					sourceJacRight[m]=0.;
				}
				get_jacobians(i);
		
				// Add change of flux (flux Jacobian) to implicit operator
				for (int j=0;j<nSolVar;++j) {
					col=(grid.myOffset+parent)*nSolVar+i; // Effect of parent ith var perturbation
					row=(grid.myOffset+parent)*nSolVar+j; // on parent jth flux
					value=-1.*jacobianLeft[j];
					if (doLeftSourceJac) value-=sourceJacLeft[j];
					MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
					if (face.bc==INTERNAL) { 
						row=(grid.myOffset+neighbor)*nSolVar+j; // on neighbor jth flux
						value=jacobianLeft[j]; 
						MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
					}
				} // for j
	
				if (face.bc==INTERNAL || face.bc==GHOST) { 
	
					// Add change of flux (flux Jacobian) to implicit operator
					for (int j=0;j<nSolVar;++j) {
						if (face.bc==INTERNAL) { 
							col=(grid.myOffset+neighbor)*nSolVar+i; // Effect of neighbor ith var perturbation
							row=(grid.myOffset+neighbor)*nSolVar+j; // on neighbor jth flux
							value=jacobianRight[j];
							if (doRightSourceJac) value-=sourceJacRight[j];
							MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
							row=(grid.myOffset+parent)*nSolVar+j; // on parent jth flux
							value=-1.*jacobianRight[j];
							MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
						} else  { // Ghost (only add effect on parent cell, effect on itself is taken care of in its own partition
							row=(grid.myOffset+parent)*nSolVar+j;
							col=(grid.ghost[-1*neighbor-1].matrix_id)*nSolVar+i;
							value=jacobianRight[j];
							MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
						}
					}
				} // if 
			} // for i
		} // if implicit
	} // for faces
	return;
} // end function

void get_jacobians(const int var) {

	using namespace state;
	using state::left;
	using state::right;
	double epsilon;
	double factor=0.01;

	leftPlus=left;
	rightPlus=right;
	
	for (int m=0;m<nSolVar;++m) {
		sourceLeftPlus[m]=sourceLeft[m];
		sourceRightPlus[m]=sourceRight[m];
		sourceJacLeft[m]=0.;
		sourceJacRight[m]=0.;
	}
	
	if (left.update[var]>0.) {epsilon=max(sqrt_machine_error,factor*left.update[var]);}
	else {epsilon=min(-1.*sqrt_machine_error,factor*left.update[var]); }
	
	if (var==5 || var==6 ) epsilon=sqrt_machine_error;
	// Perturb left state
	state_perturb(leftPlus,face,var,epsilon);
	// If right state is a boundary, correct the condition according to changes in left state
	if (face.bc>=0) right_state_update(leftPlus,rightPlus,face);
	convective_face_flux(leftPlus,rightPlus,face,&fluxPlus.convective[0]);
	
	if (EQUATIONS==NS) {
		face_state_adjust(leftPlus,rightPlus,face,var);
		diffusive_face_flux(leftPlus,rightPlus,face,&fluxPlus.diffusive[0]);
	} 
	
	
// 	if (doLeftSourceJac) {
// 		sources(left,&sourceLeft[0],true);
// 		sources(leftPlus,&sourceLeftPlus[0],true);
// 	}
	
	if (doLeftSourceJac && var==5) {
		sourceJacLeft[5]=-0.09*left.rho*left.omega*left.volume;
	}
	if (doLeftSourceJac && var==6) {
		sourceJacLeft[6]=-3./40.*2.*left.rho*left.omega*left.volume;
	}
	
	for (int j=0;j<nSolVar;++j){
		jacobianLeft[j]=(fluxPlus.diffusive[j]-flux.diffusive[j]-fluxPlus.convective[j]+flux.convective[j])/epsilon;
		//if (doLeftSourceJac) sourceJacLeft[j]=(sourceLeftPlus[j]-sourceLeft[j])/epsilon;
	}

	if (face.bc==INTERNAL || face.bc==GHOST) { 
		
		if (right.update[var]>0.) {
			epsilon=max(sqrt_machine_error,factor*right.update[var]); 
		} else {
			epsilon=min(-1.*sqrt_machine_error,factor*right.update[var]); 
		}

		if (var==5 || var==6 ) epsilon=sqrt_machine_error;
		
		state_perturb(rightPlus,face,var,epsilon);
			
		convective_face_flux(left,rightPlus,face,&fluxPlus.convective[0]);
		if (EQUATIONS==NS) {
			face_state_adjust(left,rightPlus,face,var);
			diffusive_face_flux(left,rightPlus,face,&fluxPlus.diffusive[0]);
		} 
		
// 		if (doRightSourceJac) {
// 			sources(right,&sourceRight[0],true);
// 			sources(rightPlus,&sourceRightPlus[0],true);
// 		}
		
		
		if (doRightSourceJac && var==5) {
			sourceJacRight[5]=-0.09*right.rho*right.omega*right.volume;
		}
		if (doRightSourceJac && var==6) {
			sourceJacRight[6]=-3./40.*2.*right.rho*right.omega*right.volume;
		}
		
		for (int j=0;j<nSolVar;++j){
			jacobianRight[j]=(fluxPlus.diffusive[j]-flux.diffusive[j]-fluxPlus.convective[j]+flux.convective[j])/epsilon;
			//if (doRightSourceJac) sourceJacRight[j]=(sourceRightPlus[j]-sourceRight[j])/epsilon;
		}
	}
	face_state_adjust(left,right,face,var);
}

void left_state_update(Cell_State &left,Face_State &face) {

	unsigned int parent;
	Vec3D deltaV;
	double delta[nSolVar];

	parent=grid.face[face.index].parent;

	if (order==SECOND) {
		for (unsigned int i=0;i<nSolVar;++i) {
			delta[i]=(grid.face[face.index].centroid-grid.cell[parent].centroid).dot(grid.cell[parent].limited_grad[i]);
		}
	} else {
		for (unsigned int i=0;i<nSolVar;++i) delta[i]=0.;
	}

	for (unsigned int i=0;i<nSolVar;++i) left.update[i]=grid.cell[parent].update[i];
	
	deltaV[0]=delta[1]; deltaV[1]=delta[2]; deltaV[2]=delta[3];
	
	delta[5]=0.; delta[6]=0.;
	
	// Set left primitive variables
	left.p=grid.cell[parent].p+delta[0];
	left.v_center=grid.cell[parent].v;
	left.v=left.v_center+deltaV;
	left.T_center=grid.cell[parent].T;
	left.T=left.T_center+delta[4];
	left.rho=eos.rho(left.p,left.T);
	left.a=sqrt(Gamma*(left.p+Pref)/left.rho);
	left.H=left.a*left.a/(Gamma-1.)+0.5*left.v.dot(left.v);
	left.mu=grid.cell[parent].mu;		
	left.vN[0]=left.v.dot(face.normal);
	left.vN[1]=left.v.dot(face.tangent1);
	left.vN[2]=left.v.dot(face.tangent2);
	left.volume=grid.cell[parent].volume;
	if (TURBULENCE_MODEL!=NONE) {
		left.k_center=grid.cell[parent].k;
		left.k=max(0.,left.k_center+delta[5]);
		left.omega_center=grid.cell[parent].omega;
		left.omega=max(omegaLowLimit,left.omega_center+delta[6]);
		left.gradU=grid.cell[parent].limited_grad[1];
		left.gradV=grid.cell[parent].limited_grad[2];
		left.gradW=grid.cell[parent].limited_grad[3];
		left.gradK=grid.cell[parent].limited_grad[5];
		left.gradOmega=grid.cell[parent].limited_grad[6];
	}
	
	return;
} // end left_state_update

void right_state_update(Cell_State &left,Cell_State &right,Face_State &face) {

	if (face.bc==INTERNAL) {// internal face

		Vec3D deltaV;
		double delta[nSolVar];
		unsigned int neighbor;

  		neighbor=grid.face[face.index].neighbor;

		if (order==SECOND) {
			for (unsigned int i=0;i<nSolVar;++i) {
				delta[i]=(grid.face[face.index].centroid-grid.cell[neighbor].centroid).dot(grid.cell[neighbor].limited_grad[i]);
			}
		} else {
			for (unsigned int i=0;i<nSolVar;++i) delta[i]=0.;
		}

		for (unsigned int i=0;i<nSolVar;++i) right.update[i]=grid.cell[neighbor].update[i];
		
		deltaV[0]=delta[1]; deltaV[1]=delta[2]; deltaV[2]=delta[3];
		delta[5]=0.; delta[6]=0.;

		// Set right primitive variables
		right.p=grid.cell[neighbor].p+delta[0];
		right.v_center=grid.cell[neighbor].v;
		right.v=right.v_center+deltaV;
		right.T_center=grid.cell[neighbor].T;
		right.T=right.T_center+delta[4];
		right.rho=eos.rho(right.p,right.T);
		right.mu=grid.cell[neighbor].mu;
		right.volume=grid.cell[neighbor].volume;
		if (TURBULENCE_MODEL!=NONE) {
			right.k_center=grid.cell[neighbor].k;
			right.k=right.k_center+delta[5];
			right.omega_center=grid.cell[neighbor].omega;
			right.omega=right.omega_center+delta[6];
			right.gradU=grid.cell[neighbor].limited_grad[1];
			right.gradV=grid.cell[neighbor].limited_grad[2];
			right.gradW=grid.cell[neighbor].limited_grad[3];
			right.gradK=grid.cell[neighbor].limited_grad[5];
			right.gradOmega=grid.cell[neighbor].limited_grad[6];
		}
		
	} else if (face.bc>=0) { // boundary face
		right.mu=left.mu;
		
		for (unsigned int i=0;i<nSolVar;++i) right.update[i]=0.;
		
		
		if (bc.region[face.bc].specified==BC_STATE) {
			right.p=bc.region[face.bc].p;
			right.T=bc.region[face.bc].T;
			right.T_center=right.T;
			right.rho=bc.region[face.bc].rho;
		} else if (bc.region[face.bc].specified==BC_P) {
			right.p=0.5*(left.p+bc.region[face.bc].p);
			right.T=left.T; // temperature is extrapolated
			right.T_center=left.T_center+2.*(right.T-left.T_center);
			right.rho=eos.rho(right.p,right.T); 
		} else if (bc.region[face.bc].specified==BC_T) {
			right.T=bc.region[face.bc].T;
			right.T_center=right.T;
			right.p=left.p; // pressure is extrapolated
			right.rho=eos.rho(right.p,right.T); 
		} else if (bc.region[face.bc].specified==BC_RHO) {
			right.rho=0.5*(left.rho+bc.region[face.bc].rho);
			right.p=left.p; // pressure is extrapolated
			right.T=eos.T(right.p,right.rho); 
			right.T_center=left.T_center+2.*(right.T-left.T_center);
		} else {
			// If nothing is specified, everything is extrapolated
			right.p=left.p;
			right.T=left.T; 
			right.T_center=left.T_center+2.*(right.T-left.T_center);
			right.rho=left.rho;
		}
		
		if (bc.region[face.bc].thermalType==ADIABATIC) right.T_center=left.T_center;
		
		if (bc.region[face.bc].type==OUTLET) {
			right.v=left.v;
			right.v_center=left.v_center+2.*(right.v-left.v_center); 
			right.k=left.k;
			right.k_center=left.k_center+2.*(right.k-left.k_center);
			right.omega=left.omega;
			right.omega_center=left.omega_center+2.*(right.omega-left.omega_center);
		} else if (bc.region[face.bc].type==SLIP) {
			right.v=left.v-2.*left.v.dot(face.normal)*face.normal;
			right.v_center=left.v_center-2.*left.v_center.dot(face.normal)*face.normal;
			right.k=0.;
			right.k_center=0.;
			right.omega=60.*viscosity/(right.rho*0.075*pow(0.5*fabs((face.left2right).dot(face.normal)),2.));
			//right.omega=left.omega;
			right.omega_center=left.omega_center+2.*(right.omega-left.omega_center);
		} else if (bc.region[face.bc].type==SYMMETRY) {
			right.v=left.v-2.*left.v.dot(face.normal)*face.normal;
			right.v_center=left.v_center-2.*left.v_center.dot(face.normal)*face.normal;
			right.k=left.k;
			right.k_center=left.k_center;
			right.omega=left.omega;
			right.omega_center=left.omega_center;
		} else if (bc.region[face.bc].type==NOSLIP) {
			right.v=-1.*left.v;
			right.v_center=-1.*left.v_center;
			right.k=0.;
			right.k_center=0.;
			right.omega=60.*viscosity/(right.rho*0.075*pow(0.5*fabs((face.left2right).dot(face.normal)),2.));
			right.omega_center=left.omega_center+2.*(right.omega-left.omega_center);
		} else if (bc.region[face.bc].type==INLET) {
			right.v=bc.region[face.bc].v;
			right.v_center=right.v;
			right.k=bc.region[face.bc].k;
			right.k_center=right.k;
			right.omega=bc.region[face.bc].omega;
			right.omega_center=right.omega;
		}
		
	} else { // partition boundary

		int g=-1*grid.face[face.index].neighbor-1; // ghost cell index
		Vec3D deltaV;
		double delta[nSolVar];
		if (order==SECOND) {
			for (unsigned int i=0;i<nSolVar;++i) {
				delta[i]=(grid.face[face.index].centroid-grid.ghost[g].centroid).dot(grid.ghost[g].limited_grad[i]);
			}
		}

		for (unsigned int i=0;i<nSolVar;++i) right.update[i]=grid.ghost[g].update[i];
		
		deltaV[0]=delta[1]; deltaV[1]=delta[2]; deltaV[2]=delta[3];
		delta[5]=0.; delta[6]=0.;
		right.p=grid.ghost[g].p+delta[0];	
		right.v_center=grid.ghost[g].v;
		right.v=right.v_center+deltaV;
		right.T_center=grid.ghost[g].T;
		right.T=right.T_center+delta[4];
		right.rho=eos.rho(right.p,right.T);
		right.mu=grid.ghost[g].mu;
		if (TURBULENCE_MODEL!=NONE) {
			right.k_center=grid.ghost[g].k;
			right.k=right.k_center+delta[5];
			right.omega_center=grid.ghost[g].omega;
			right.omega=max(omegaLowLimit,right.omega_center+delta[6]);
		}
	}

	right.a=sqrt(Gamma*(right.p+Pref)/right.rho);
	right.H=right.a*right.a/(Gamma-1.)+0.5*right.v.dot(right.v);
	right.vN[0]=right.v.dot(face.normal);
	right.vN[1]=right.v.dot(face.tangent1);
	right.vN[2]=right.v.dot(face.tangent2);
	
	right.k=max(0.,right.k);
	right.k_center=max(0.,right.k_center);
	right.omega=max(omegaLowLimit,right.omega);
	right.omega_center=max(omegaLowLimit,right.omega_center);
	
	return;
} // end right_state_update

void face_geom_update(Face_State &face,unsigned int f) {
	face.index=f;
	face.normal=grid.face[f].normal;
	face.tangent1=(grid.face[f].centroid-grid.face[f].node(0));
	face.tangent1/=fabs(face.tangent1);
	// Cross the tangent vector with the normal vector to get the second tangent
	face.tangent2=(face.normal).cross(face.tangent1);
	face.area=grid.face[f].area;
	face.bc=grid.face[f].bc;
	Vec3D leftCentroid,rightCentroid;
	leftCentroid=grid.cell[grid.face[f].parent].centroid;
	if (face.bc==INTERNAL) { rightCentroid=grid.cell[grid.face[f].neighbor].centroid;}
	else if (face.bc>=0) {
		rightCentroid=leftCentroid+2.*(grid.face[f].centroid-leftCentroid).dot(face.normal)*face.normal;
	} else { rightCentroid=grid.ghost[-1*grid.face[f].neighbor-1].centroid;}
	face.left2right=rightCentroid-leftCentroid;
	return;
} // end face_geom_update

void face_state_update(Cell_State &left,Cell_State &right,Face_State &face) {

	map<int,double>::iterator fit;
	// Find face averaged variables
	face.gradU=0.; face.gradV=0.; face.gradW=0.; face.gradT=0.; face.gradK=0.; face.gradOmega=0.;
	for (fit=grid.face[face.index].average.begin();fit!=grid.face[face.index].average.end();fit++) {
		if ((*fit).first>=0) { // if contribution is coming from a real cell
			face.gradU+=(*fit).second*grid.cell[(*fit).first].grad[1];
			face.gradV+=(*fit).second*grid.cell[(*fit).first].grad[2];
			face.gradW+=(*fit).second*grid.cell[(*fit).first].grad[3];
			face.gradT+=(*fit).second*grid.cell[(*fit).first].grad[4];
			face.gradK+=(*fit).second*grid.cell[(*fit).first].grad[5];
			face.gradOmega+=(*fit).second*grid.cell[(*fit).first].grad[6];
		} else { // if contribution is coming from a ghost cell
			face.gradU+=(*fit).second*grid.ghost[-1*((*fit).first+1)].grad[1];
			face.gradV+=(*fit).second*grid.ghost[-1*((*fit).first+1)].grad[2];
			face.gradW+=(*fit).second*grid.ghost[-1*((*fit).first+1)].grad[3];
			face.gradT+=(*fit).second*grid.ghost[-1*((*fit).first+1)].grad[4];
			face.gradK+=(*fit).second*grid.ghost[-1*((*fit).first+1)].grad[5];
			face.gradOmega+=(*fit).second*grid.ghost[-1*((*fit).first+1)].grad[6];
		}
	}
			
	face.gradU-=face.gradU.dot(face.normal)*face.normal;
	face.gradU+=((right.v_center[0]-left.v_center[0])/(face.left2right.dot(face.normal)))*face.normal;

	face.gradV-=face.gradV.dot(face.normal)*face.normal;
	face.gradV+=((right.v_center[1]-left.v_center[1])/(face.left2right.dot(face.normal)))*face.normal;

	face.gradW-=face.gradW.dot(face.normal)*face.normal;
	face.gradW+=((right.v_center[2]-left.v_center[2])/(face.left2right.dot(face.normal)))*face.normal;

	if (TURBULENCE_MODEL!=NONE) {
		face.gradK-=face.gradK.dot(face.normal)*face.normal;
		face.gradK+=((right.k_center-left.k_center)/(face.left2right.dot(face.normal)))*face.normal;
	
		face.gradOmega-=face.gradOmega.dot(face.normal)*face.normal;
		face.gradOmega+=((right.omega_center-left.omega_center)/(face.left2right.dot(face.normal)))*face.normal;
	}

	// Boundary conditions are already taken care of in right state update
	face.p=0.5*(left.p+right.p);
	face.v=0.5*(left.v+right.v);
	face.T=0.5*(left.T+right.T);
	face.k=0.5*(left.k+right.k);
	face.omega=0.5*(left.omega+right.omega);
	face.rho=0.5*(left.rho+right.rho);
	face.mu=0.5*(left.mu+right.mu);
	
	return;
} // end face_state_update

void state_perturb(Cell_State &state,Face_State &face,int var,double epsilon) {

	switch (var)
	{
		case 0 : // p
			state.p+=epsilon;
			state.rho=eos.rho(state.p,state.T);
			state.a=sqrt(Gamma*(state.p+Pref)/state.rho);
			state.H=state.a*state.a/(Gamma-1.)+0.5*state.v.dot(state.v);
			break;
		case 1 : // uleft,face,
			state.v[0]+=epsilon;
			state.v_center[0]+=epsilon;
			state.H=state.a*state.a/(Gamma-1.)+0.5*state.v.dot(state.v);
			state.vN[0]=state.v.dot(face.normal);
			state.vN[1]=state.v.dot(face.tangent1);
			state.vN[2]=state.v.dot(face.tangent2);
			break;
		case 2 : // v
			state.v[1]+=epsilon;
			state.v_center[1]+=epsilon;
			state.H=state.a*state.a/(Gamma-1.)+0.5*state.v.dot(state.v);
			state.vN[0]=state.v.dot(face.normal);
			state.vN[1]=state.v.dot(face.tangent1);
			state.vN[2]=state.v.dot(face.tangent2);
			break;
		case 3 : // w
			state.v[2]+=epsilon;
			state.v_center[2]+=epsilon;
			state.H=state.a*state.a/(Gamma-1.)+0.5*state.v.dot(state.v);
			state.vN[0]=state.v.dot(face.normal);
			state.vN[1]=state.v.dot(face.tangent1);
			state.vN[2]=state.v.dot(face.tangent2);
			break;
		case 4 : // T
			state.T+=epsilon;
			state.T_center+=epsilon;
			state.rho=eos.rho(state.p,state.T);
			state.a=sqrt(Gamma*(state.p+Pref)/state.rho);
			state.H=state.a*state.a/(Gamma-1.)+0.5*state.v.dot(state.v);
			break;
		case 5 : // k
			state.k+=epsilon;
			state.k_center+=epsilon;
			break;
		case 6 : // omega
			state.omega+=epsilon;
			state.omega_center+=epsilon;
			break;
	}
	
	return;
} // end state_perturb

void face_state_adjust(Cell_State &left,Cell_State &right,Face_State &face,int var) {

	switch (var)
	{
		case 0 : // p
			face.p=0.5*(left.p+right.p);
			face.rho=eos.rho(face.p,face.T);
			break;
		case 1 : // u
			face.v[0]=0.5*(left.v[0]+right.v[0]);
			face.gradU-=face.gradU.dot(face.normal)*face.normal;
			face.gradU+=((right.v_center[0]-left.v_center[0])/(face.left2right.dot(face.normal)))*face.normal;
			break;
		case 2 : // v
			face.v[1]=0.5*(left.v[1]+right.v[1]);
			face.gradV-=face.gradV.dot(face.normal)*face.normal;
			face.gradV+=((right.v_center[1]-left.v_center[1])/(face.left2right.dot(face.normal)))*face.normal;
			break;
		case 3 : // w
			face.v[2]=0.5*(left.v[2]+right.v[2]);
			face.gradW-=face.gradW.dot(face.normal)*face.normal;
			face.gradW+=((right.v_center[2]-left.v_center[2])/(face.left2right.dot(face.normal)))*face.normal;
			break;
		case 4 : // T
			face.T=0.5*(left.T+right.T);
			face.rho=eos.rho(face.p,face.T);
			face.gradT-=face.gradT.dot(face.normal)*face.normal;
			face.gradT+=((right.T_center-left.T_center)/(face.left2right.dot(face.normal)))*face.normal;
			break;
		case 5 : // k
			face.k=0.5*(left.k+right.k);
			face.gradK-=face.gradK.dot(face.normal)*face.normal;
			face.gradK+=((right.k_center-left.k_center)/(face.left2right.dot(face.normal)))*face.normal;
			break;
		case 6 : // omega
			face.omega=0.5*(left.omega+right.omega);
			face.gradOmega-=face.gradOmega.dot(face.normal)*face.normal;
			face.gradOmega+=((right.omega_center-left.omega_center)/(face.left2right.dot(face.normal)))*face.normal;
			break;
	}

	return;
} // end face_state_adjust
