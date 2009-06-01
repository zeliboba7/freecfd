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
#include "state_cache.h"
#include "petsc_functions.h"
#include "flamelet.h"
#include "rans.h"

extern BC bc;
extern InputFile input;
extern Flamelet flamelet;
extern RANS rans;

void get_jacobians(const int var);
void convective_face_flux(Cell_State &left,Cell_State &right,Face_State &face,double flux[]);
void apply_bcs(Cell_State &left,Cell_State &right,Face_State &face);
void diffusive_face_flux(Cell_State &left,Cell_State &right,Face_State &face,double flux[]);
void left_state_update(Cell_State &left,Face_State &face);
void right_state_update(Cell_State &left,Cell_State &right,Face_State &face);
void face_geom_update(Face_State &face,int f);
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
	
	int parent,neighbor,f;
	int row,col;
	
	vector<bool> cellVisited;
	for (int c=0;c<grid.cellCount;++c) cellVisited.push_back(false);
	PetscScalar value;

	flux.convective.resize(5);
	flux.diffusive.resize(5);
	fluxPlus.convective.resize(5);
	fluxPlus.diffusive.resize(5);
	jacobianLeft.resize(5);
	jacobianRight.resize(5);
	left.update.resize(5);
	right.update.resize(5);
	leftPlus.update.resize(5);
	rightPlus.update.resize(5);
	sourceLeft.resize(5);
	sourceRight.resize(5);
	sourceLeftPlus.resize(5);
	sourceRightPlus.resize(5);
	sourceJacLeft.resize(5);
	sourceJacRight.resize(5);
	
	for (int m=0;m<5;++m) flux.diffusive[m]=0.;
	
	bool implicit=true;
	if (TIME_INTEGRATOR==FORWARD_EULER) implicit=false;
	// Loop through faces
	for (f=0;f<grid.faceCount;++f) {
		
		doLeftSourceJac=false; doRightSourceJac=false;
		for (int m=0;m<5;++m) { 
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
		for (int i=0;i<5;++i) {
			row=(grid.myOffset+parent)*5+i;
			value=flux.diffusive[i]-flux.convective[i]+sourceLeft[i];
			VecSetValues(rhs,1,&row,&value,ADD_VALUES);
			if (grid.face[f].bc==INTERNAL) { 
				row=(grid.myOffset+neighbor)*5+i;
				value=-1.*(flux.diffusive[i]-flux.convective[i])+sourceRight[i];
				VecSetValues(rhs,1,&row,&value,ADD_VALUES);
			}
		}

		if (implicit && ps_timeStep==1) {

			for (int i=0;i<5;++i) { // perturb each variable
				
				for (int m=0;m<5;++m) { 
					sourceJacLeft[m]=0.;
					sourceJacRight[m]=0.;
				}
				
				get_jacobians(i);
				
				// Add change of flux (flux Jacobian) to implicit operator
				for (int j=0;j<5;++j) {
					col=(grid.myOffset+parent)*5+i; // Effect of parent ith var perturbation
					row=(grid.myOffset+parent)*5+j; // on parent jth flux
					value=-1.*jacobianLeft[j];
					//if (doLeftSourceJac) value-=sourceJacLeft[j];
					MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
					if (face.bc==INTERNAL) { 
						row=(grid.myOffset+neighbor)*5+j; // on neighbor jth flux
						value=jacobianLeft[j]; 
						MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
					}
				} // for j
				
				if (face.bc==INTERNAL || face.bc==GHOST) { 
	
					// Add change of flux (flux Jacobian) to implicit operator
					for (int j=0;j<5;++j) {
						if (face.bc==INTERNAL) { 
							col=(grid.myOffset+neighbor)*5+i; // Effect of neighbor ith var perturbation
							row=(grid.myOffset+neighbor)*5+j; // on neighbor jth flux
							value=jacobianRight[j];
							//if (doRightSourceJac) value-=sourceJacRight[j];
							MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
							row=(grid.myOffset+parent)*5+j; // on parent jth flux
							value=-1.*jacobianRight[j];
							MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
						} else  { // Ghost (only add effect on parent cell, effect on itself is taken care of in its own partition
							row=(grid.myOffset+parent)*5+j;
							col=(grid.ghost[-1*neighbor-1].matrix_id)*5+i;
							value=-1.*jacobianRight[j];
							MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
						}
					}
						
				} // if 
				
			} // for i (each perturbed variable)
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
	
	for (int m=0;m<5;++m) {
		sourceLeftPlus[m]=sourceLeft[m];
		sourceRightPlus[m]=sourceRight[m];
		sourceJacLeft[m]=0.;
		sourceJacRight[m]=0.;
	}
	
	if (left.update[var]>=0.) {epsilon=max(sqrt_machine_error,factor*left.update[var]);}
	else {epsilon=min(-1.*sqrt_machine_error,factor*left.update[var]); }

	// Perturb left state
	if (FLAMELET) {
		if (var==4) { // var=4 is Z for flamelet
			if ((left.Z+epsilon)>1.){
				epsilon=-1.*sqrt_machine_error;
			} else if ((left.Z+epsilon)<0.) {
				epsilon=sqrt_machine_error;
			}
		}
	}
	state_perturb(leftPlus,face,var,epsilon);
	// If right state is a boundary, correct the condition according to changes in left state
	if (face.bc>=0) right_state_update(leftPlus,rightPlus,face);

	convective_face_flux(leftPlus,rightPlus,face,&fluxPlus.convective[0]);
	
	if (EQUATIONS==NS) {
		face_state_adjust(leftPlus,rightPlus,face,var);
		diffusive_face_flux(leftPlus,rightPlus,face,&fluxPlus.diffusive[0]);
	} 
	
	if (doLeftSourceJac) {
		sources(left,&sourceLeft[0],true);
		sources(leftPlus,&sourceLeftPlus[0],true);
	}
	
	for (int j=0;j<5;++j){
		jacobianLeft[j]=(fluxPlus.diffusive[j]-flux.diffusive[j]-fluxPlus.convective[j]+flux.convective[j])/epsilon;
		if (doLeftSourceJac) sourceJacLeft[j]=(sourceLeftPlus[j]-sourceLeft[j])/epsilon;
	}
	
	if (face.bc==INTERNAL || face.bc==GHOST) { 
		
		if (right.update[var]>=0.) {
			epsilon=max(sqrt_machine_error,factor*right.update[var]); 
		} else {
			epsilon=min(-1.*sqrt_machine_error,factor*right.update[var]); 
		}
		
		if (FLAMELET) {
			if (var==4) {
				if ((right.Z+epsilon)>1.){
					epsilon=-1.*sqrt_machine_error;
				} else if ((right.Z+epsilon)<0.) {
					epsilon=sqrt_machine_error;
				}
			}
		}
		
		state_perturb(rightPlus,face,var,epsilon);
		
		convective_face_flux(left,rightPlus,face,&fluxPlus.convective[0]);
		
		if (EQUATIONS==NS) {
			face_state_adjust(left,rightPlus,face,var);
			diffusive_face_flux(left,rightPlus,face,&fluxPlus.diffusive[0]);
		} 
		
		if (doRightSourceJac) {
			sources(right,&sourceRight[0],true);
			sources(rightPlus,&sourceRightPlus[0],true);
		}
		
		for (int j=0;j<5;++j){
			jacobianRight[j]=(fluxPlus.diffusive[j]-flux.diffusive[j]-fluxPlus.convective[j]+flux.convective[j])/epsilon;
			if (doRightSourceJac) sourceJacRight[j]=(sourceRightPlus[j]-sourceRight[j])/epsilon;
		}
	
	}
	
	face_state_adjust(left,right,face,var);
}

void left_state_update(Cell_State &left,Face_State &face) {
	
	int parent;
	Vec3D deltaV;
	double delta[5];
	
	parent=grid.face[face.index].parent;

	if (order==SECOND) {
		for (int i=0;i<5;++i) {
			delta[i]=(grid.face[face.index].centroid-grid.cell[parent].centroid).dot(grid.cell[parent].grad[i]);
		}
	} else {
		for (int i=0;i<5;++i) delta[i]=0.;
	}

	for (int i=0;i<5;++i) left.update[i]=grid.cell[parent].update[i];
	
	deltaV[0]=delta[1]; deltaV[1]=delta[2]; deltaV[2]=delta[3];
	
	// Set left primitive variables
	left.p=grid.cell[parent].p+delta[0];
	left.v_center=grid.cell[parent].v;
	left.v=left.v_center+deltaV;
	left.T_center=grid.cell[parent].T;
	left.T=left.T_center+delta[4];
	left.rho=eos.rho(left.p,left.T);
	left.gamma=Gamma;
	
	if (FLAMELET) {
		delta[1]=(grid.face[face.index].centroid-grid.cell[parent].centroid).dot(flamelet.cell[parent].grad[0]);
		delta[2]=(grid.face[face.index].centroid-grid.cell[parent].centroid).dot(flamelet.cell[parent].grad[1]);
		left.Z_center=flamelet.cell[parent].Z;
		left.Z=left.Z_center+delta[1];
		left.Zvar=flamelet.cell[parent].Zvar+delta[2];
		left.Z=max(0.,left.Z);
		left.Z=min(1.,left.Z);
		left.Zvar=max(0.,left.Zvar);
		left.Zvar=min(0.25,left.Zvar);
		left.Chi=log10(2.0*rans.kepsilon.beta_star*rans.cell[parent].omega*left.Zvar);
		flamelet.table.get_rho_T_comp(left.p,left.Z,left.Zvar,left.Chi,left.rho,left.T);
		left.R=UNIV_GAS_CONST/flamelet.table.get_Mw(left.Z,left.Zvar,left.Chi);
		left.gamma=flamelet.table.get_gamma(left.Z,left.Zvar,left.Chi,false);
		left.update[4]=flamelet.cell[parent].update[0];
	}
	
	left.a=sqrt(left.gamma*(left.p+Pref)/left.rho);
	left.H=left.a*left.a/(left.gamma-1.)+0.5*left.v.dot(left.v);
	left.vN[0]=left.v.dot(face.normal);
	left.vN[1]=left.v.dot(face.tangent1);
	left.vN[2]=left.v.dot(face.tangent2);
	left.volume=grid.cell[parent].volume;
	
	return;
} // end left_state_update

void right_state_update(Cell_State &left,Cell_State &right,Face_State &face) {
	
	if (face.bc==INTERNAL) {// internal face

		Vec3D deltaV;
		double delta[5];
		int neighbor;

  		neighbor=grid.face[face.index].neighbor;

		if (order==SECOND) {
			for (int i=0;i<5;++i) {
				delta[i]=(grid.face[face.index].centroid-grid.cell[neighbor].centroid).dot(grid.cell[neighbor].grad[i]);
			}
		} else {
			for (int i=0;i<5;++i) delta[i]=0.;
		}

		for (int i=0;i<5;++i) right.update[i]=grid.cell[neighbor].update[i];
		
		deltaV[0]=delta[1]; deltaV[1]=delta[2]; deltaV[2]=delta[3];

		// Set right primitive variables
		right.p=grid.cell[neighbor].p+delta[0];
		right.v_center=grid.cell[neighbor].v;
		right.v=right.v_center+deltaV;
		right.T_center=grid.cell[neighbor].T;
		right.T=right.T_center+delta[4];
		right.rho=eos.rho(right.p,right.T);
		right.volume=grid.cell[neighbor].volume;
		right.gamma=Gamma;
		
		if (FLAMELET) {
			delta[1]=(grid.face[face.index].centroid-grid.cell[neighbor].centroid).dot(flamelet.cell[neighbor].grad[0]);
			delta[2]=(grid.face[face.index].centroid-grid.cell[neighbor].centroid).dot(flamelet.cell[neighbor].grad[1]);
			right.Z_center=flamelet.cell[neighbor].Z;
			right.Z=right.Z_center+delta[1];
			right.Zvar=flamelet.cell[neighbor].Zvar+delta[2];
			right.Z=max(0.,right.Z);
			right.Z=min(1.,right.Z);
			right.Zvar=max(0.,right.Zvar);
			right.Zvar=min(0.25,right.Zvar);
			right.Chi=log10(2.0*rans.kepsilon.beta_star*rans.cell[neighbor].omega*right.Zvar);
			flamelet.table.get_rho_T_comp(right.p,right.Z,right.Zvar,right.Chi,right.rho,right.T);
			right.R=UNIV_GAS_CONST/flamelet.table.get_Mw(right.Z,right.Zvar,right.Chi);
			right.gamma=flamelet.table.get_gamma(right.Z,right.Zvar,right.Chi,false);
			right.update[4]=flamelet.cell[neighbor].update[0];
		}
		
	} else if (face.bc>=0) { // boundary face
		
		right.gamma=Gamma;
		bool skip=false;
		if (FLAMELET) {
			// Default is zero normal
			right.Z=left.Z; right.Zvar=left.Zvar; right.Chi=left.Chi;
			right.Z_center=left.Z_center;
			right.R=UNIV_GAS_CONST/flamelet.table.get_Mw(right.Z,right.Zvar,right.Chi);
			right.gamma=flamelet.table.get_gamma(right.Z,right.Zvar,right.Chi,false);
			
			if (bc.region[face.bc].type==INLET) {
				right.Z=bc.region[face.bc].Z;
				right.Zvar=bc.region[face.bc].Zvar;
				right.Z_center=2.*right.Z-left.Z_center;
				double mdotNR=-1.*(bc.region[face.bc].mdot/bc.region[face.bc].area)*bc.region[face.bc].areaVec.norm().dot(face.normal);
				// extrapolate pressure
				right.p=left.p;
				flamelet.table.get_rho_T_comp(right.p,right.Z,right.Zvar,right.Chi,right.rho,right.T);
				right.gamma=flamelet.table.get_gamma(right.Z,right.Zvar,right.Chi);
				right.a=sqrt(right.gamma*right.R*(right.T+Tref));
				right.v=mdotNR/right.rho*face.normal;
				right.v_center=2.*right.v-left.v_center;
				skip=true;
			} 

		}
		
		for (int i=0;i<5;++i) right.update[i]=0.;
		if (!skip) apply_bcs(left,right,face);
		
	} else { // partition boundary

		int g=-1*grid.face[face.index].neighbor-1; // ghost cell index
		Vec3D deltaV;
		double delta[5];
		if (order==SECOND) {
			for (int i=0;i<5;++i) {
				delta[i]=(grid.face[face.index].centroid-grid.ghost[g].centroid).dot(grid.ghost[g].grad[i]);
			}
		} else {
			for (int i=0;i<5;++i) delta[i]=0.;
		}

		for (int i=0;i<5;++i) right.update[i]=grid.ghost[g].update[i];
		
		deltaV[0]=delta[1]; deltaV[1]=delta[2]; deltaV[2]=delta[3];

		right.p=grid.ghost[g].p+delta[0];	
		right.v_center=grid.ghost[g].v;
		right.v=right.v_center+deltaV;
		right.T_center=grid.ghost[g].T;
		right.T=right.T_center+delta[4];
		right.rho=eos.rho(right.p,right.T);
		right.gamma=Gamma;
		
		if (FLAMELET) {
			delta[1]=(grid.face[face.index].centroid-grid.ghost[g].centroid).dot(flamelet.ghost[g].grad[0]);
			delta[2]=(grid.face[face.index].centroid-grid.ghost[g].centroid).dot(flamelet.ghost[g].grad[1]);
			right.Z_center=flamelet.ghost[g].Z;
			right.Z=right.Z_center+delta[1];
			right.Zvar=flamelet.ghost[g].Zvar+delta[2];
			right.Z=max(0.,right.Z);
			right.Z=min(1.,right.Z);
			right.Zvar=max(0.,right.Zvar);
			right.Zvar=min(0.25,right.Zvar);
			right.Chi=log10(2.0*rans.kepsilon.beta_star*rans.ghost[g].omega*right.Zvar);
			flamelet.table.get_rho_T_comp(right.p,right.Z,right.Zvar,right.Chi,right.rho,right.T);
			right.R=UNIV_GAS_CONST/flamelet.table.get_Mw(right.Z,right.Zvar,right.Chi);
			right.gamma=flamelet.table.get_gamma(right.Z,right.Zvar,right.Chi,false);
			right.update[4]=sqrt_machine_error;
		}
	
	}
	
	right.a=sqrt(right.gamma*(right.p+Pref)/right.rho);
	right.H=right.a*right.a/(right.gamma-1.)+0.5*right.v.dot(right.v);
	right.vN[0]=right.v.dot(face.normal);
	right.vN[1]=right.v.dot(face.tangent1);
	right.vN[2]=right.v.dot(face.tangent2);
		
	return;
} // end right_state_update

void face_geom_update(Face_State &face,int f) {
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
	face.gradU=0.; face.gradV=0.; face.gradW=0.; face.gradT=0.; face.gradZ=0.; face.diff=0.;
	for (fit=grid.face[face.index].average.begin();fit!=grid.face[face.index].average.end();fit++) {
		if ((*fit).first>=0) { // if contribution is coming from a real cell
			face.gradU+=(*fit).second*grid.cell[(*fit).first].grad[1];
			face.gradV+=(*fit).second*grid.cell[(*fit).first].grad[2];
			face.gradW+=(*fit).second*grid.cell[(*fit).first].grad[3];
			face.gradT+=(*fit).second*grid.cell[(*fit).first].grad[4];
			if (FLAMELET) {
				face.gradZ+=(*fit).second*flamelet.cell[(*fit).first].grad[0];
				face.diff+=(*fit).second*flamelet.cell[(*fit).first].diffusivity;
			}
		} else { // if contribution is coming from a ghost cell
			face.gradU+=(*fit).second*grid.ghost[-1*(*fit).first-1].grad[1];
			face.gradV+=(*fit).second*grid.ghost[-1*(*fit).first-1].grad[2];
			face.gradW+=(*fit).second*grid.ghost[-1*(*fit).first-1].grad[3];
			face.gradT+=(*fit).second*grid.ghost[-1*(*fit).first-1].grad[4];
			if (FLAMELET) {
				face.gradZ+=(*fit).second*flamelet.ghost[-1*(*fit).first-1].grad[0];
				face.diff+=(*fit).second*flamelet.ghost[-1*(*fit).first-1].diffusivity;
			}
		}
	}
			
	face.gradU-=face.gradU.dot(face.normal)*face.normal;
	face.gradU+=((right.v_center[0]-left.v_center[0])/(face.left2right.dot(face.normal)))*face.normal;

	face.gradV-=face.gradV.dot(face.normal)*face.normal;
	face.gradV+=((right.v_center[1]-left.v_center[1])/(face.left2right.dot(face.normal)))*face.normal;

	face.gradW-=face.gradW.dot(face.normal)*face.normal;
	face.gradW+=((right.v_center[2]-left.v_center[2])/(face.left2right.dot(face.normal)))*face.normal;
	
	face.gradT-=face.gradT.dot(face.normal)*face.normal;
	face.gradT+=((right.T_center-left.T_center)/(face.left2right.dot(face.normal)))*face.normal;
	
	if (FLAMELET) {
		face.diff=max(0.,face.diff);
		face.gradZ-=face.gradZ.dot(face.normal)*face.normal;
		face.gradZ+=((right.Z_center-left.Z_center)/(face.left2right.dot(face.normal)))*face.normal;
	}

	// Boundary conditions are already taken care of in right state update
	face.p=0.5*(left.p+right.p);
	face.v=0.5*(left.v+right.v);
	face.T=0.5*(left.T+right.T);
	
	return;
} // end face_state_update

void state_perturb(Cell_State &state,Face_State &face,int var,double epsilon) {
	
	switch (var)
	{
		case 0 : // p
			state.p+=epsilon;
			state.rho+=state.rho/(state.p+Pref)*epsilon;
			state.a=sqrt(state.gamma*(state.p+Pref)/state.rho);
			break;
		case 1 : // u
			state.v[0]+=epsilon;
			state.v_center[0]+=epsilon;
			state.vN[0]=state.v.dot(face.normal);
			state.vN[1]=state.v.dot(face.tangent1);
			state.vN[2]=state.v.dot(face.tangent2);
			break;
		case 2 : // v
			state.v[1]+=epsilon;
			state.v_center[1]+=epsilon;
			state.vN[0]=state.v.dot(face.normal);
			state.vN[1]=state.v.dot(face.tangent1);
			state.vN[2]=state.v.dot(face.tangent2);
			break;
		case 3 : // w
			state.v[2]+=epsilon;
			state.v_center[2]+=epsilon;
			state.vN[0]=state.v.dot(face.normal);
			state.vN[1]=state.v.dot(face.tangent1);
			state.vN[2]=state.v.dot(face.tangent2);
			break;
		case 4 : // T
			if (FLAMELET) { // perturb Z
				state.Z+=epsilon;
				state.Z_center+=epsilon;
				flamelet.table.get_rho_T_comp(state.p,state.Z,state.Zvar,state.Chi,state.rho,state.T);
				state.R=UNIV_GAS_CONST/flamelet.table.get_Mw(state.Z,state.Zvar,state.Chi);
				state.gamma=flamelet.table.get_gamma(state.Z,state.Zvar,state.Chi,false);
			} else {
				state.T+=epsilon;
				state.T_center+=epsilon;
				state.rho-=state.rho/(state.T+Tref)*epsilon;
			}
			state.a=sqrt(state.gamma*(state.p+Pref)/state.rho);
			break;
 	}
	
	state.H=state.a*state.a/(state.gamma-1.)+0.5*state.v.dot(state.v);
	
	return;
} // end state_perturb

void face_state_adjust(Cell_State &left,Cell_State &right,Face_State &face,int var) {

	switch (var)
	{
		case 0 : // p
			face.p=0.5*(left.p+right.p);
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
			if (FLAMELET) {
				face.gradZ-=face.gradZ.dot(face.normal)*face.normal;
				face.gradZ+=((right.Z_center-left.Z_center)/(face.left2right.dot(face.normal)))*face.normal;
			} else {
				face.T=0.5*(left.T+right.T);
				face.gradT-=face.gradT.dot(face.normal)*face.normal;
				face.gradT+=((right.T_center-left.T_center)/(face.left2right.dot(face.normal)))*face.normal;
			}
			break;
	}

	return;
} // end face_state_adjust
