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
#include "hc.h"

namespace hc_state {
	HC_Cell_State left,right,leftPlus,rightPlus;
	HC_Face_State face;
	double flux, fluxPlus;
	double jacobianLeft,jacobianRight,sourceJacLeft,sourceJacRight;
	double sourceLeft,sourceRight,sourceLeftPlus,sourceRightPlus;
	bool doLeftSourceJac,doRightSourceJac;
}

void HeatConduction::assemble_linear_system(void) {

	using namespace hc_state;
	using hc_state::left;
	using hc_state::right;
	
	int parent,neighbor,f;
	int row,col;
	
	vector<bool> cellVisited;
	for (int c=0;c<grid[gid].cellCount;++c) cellVisited.push_back(false);
	
	sqrt_machine_error=sqrt(std::numeric_limits<double>::epsilon());
	
	PetscScalar value;
	
	flux=0.;
	
	bool implicit=true;

	// Loop through faces
	for (f=0;f<grid[gid].faceCount;++f) {

		doLeftSourceJac=false; doRightSourceJac=false;
		sourceLeft=0.; sourceRight=0.;

		parent=grid[gid].face[f].parent; neighbor=grid[gid].face[f].neighbor;

		// Populate the state caches
		face_geom_update(face,f);
		left_state_update(left,face);
		right_state_update(left,right,face);
		face_state_update(left,right,face);

		diffusive_face_flux(face,flux);

		// Add Sources
		if (!cellVisited[parent]){
			sources(left,sourceLeft);
			cellVisited[parent]=true;
			doLeftSourceJac=true;
		}

		if(face.bc==INTERNAL_FACE && !cellVisited[neighbor]) {
			sources(right,sourceRight);
			cellVisited[neighbor]=true;
			doRightSourceJac=true;
		}
	
		// Integrate boundary fluxes
		// TODO: Activate this
		/*
		if ((timeStep) % integrateBoundaryFreq == 0) {
			if (face.bc>=0) {
				bc.region[face.bc].energy-=flux;
			}
		}
		*/
			
		// Fill in rhs vector
		row=grid[gid].myOffset+parent;
		value=flux+sourceLeft;
		VecSetValues(rhs,1,&row,&value,ADD_VALUES);
		if (grid[gid].face[f].bc==INTERNAL_FACE) { 
			row=grid[gid].myOffset+neighbor;
			value=-1.*flux+sourceRight;
			VecSetValues(rhs,1,&row,&value,ADD_VALUES);
		}

		//if (implicit) { // TODO: Get this working

		sourceJacLeft=0.;
		sourceJacRight=0.;

		get_jacobians();
		
		col=grid[gid].myOffset+parent; // Effect of parent perturbation
		row=grid[gid].myOffset+parent; // on parent flux
		value=-1.*jacobianLeft;
		//if (doLeftSourceJac) value-=sourceJacLeft[j];
		MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
		if (face.bc==INTERNAL_FACE) { 
			row=grid[gid].myOffset+neighbor; // on neighbor flux
			value=jacobianLeft; 
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
		}
		
		if (face.bc==INTERNAL_FACE || face.bc==GHOST_FACE) { 
			
			// Add change of flux (flux Jacobian) to implicit operator
			if (face.bc==INTERNAL_FACE) { 
				col=grid[gid].myOffset+neighbor; // Effect of neighbor perturbation
				row=grid[gid].myOffset+neighbor; // on neighbor flux
				value=jacobianRight;
				//if (doRightSourceJac) value-=sourceJacRight[j];
				MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
				row=grid[gid].myOffset+parent; // on parent flux
				value=-1.*jacobianRight;
				MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			} else  { // Ghost (only add effect on parent cell, effect on itself is taken care of in its own partition
				row=grid[gid].myOffset+parent;
				col=grid[gid].ghost[-1*neighbor-1].matrix_id;
				value=-1.*jacobianRight;
				MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			}
			
		} // if 
		
		//} // if implicit
		
	} // for faces
	
	return;
} // end function

void HeatConduction::get_jacobians(void) {

	using namespace hc_state;
	using hc_state::left;
	using hc_state::right;
	double epsilon;
	double factor=0.01;

	leftPlus=left;
	rightPlus=right;
	
	sourceLeftPlus=sourceLeft;
	sourceRightPlus=sourceRight;
	sourceJacLeft=0.;
	sourceJacRight=0.;
	
	if (left.update>=0.) {epsilon=max(sqrt_machine_error,factor*left.update);}
	else {epsilon=min(-1.*sqrt_machine_error,factor*left.update); }

	// Perturb left state
	state_perturb(leftPlus,face,epsilon);
	// If right state is a boundary, correct the condition according to changes in left state
	if (face.bc>=0) right_state_update(leftPlus,rightPlus,face);
	face_state_adjust(leftPlus,rightPlus,face);
	diffusive_face_flux(face,fluxPlus);
	
	if (doLeftSourceJac) {
		sources(left,sourceLeft,true);
		sources(leftPlus,sourceLeftPlus,true);
	}

	jacobianLeft=(fluxPlus-flux)/epsilon;
	if (doLeftSourceJac) sourceJacLeft=(sourceLeftPlus-sourceLeft)/epsilon;
		
	if (face.bc==INTERNAL_FACE || face.bc==GHOST_FACE) { 
		
		if (right.update>=0.) epsilon=max(sqrt_machine_error,factor*right.update); 
		else epsilon=min(-1.*sqrt_machine_error,factor*right.update); 
				
		state_perturb(rightPlus,face,epsilon);
		face_state_adjust(left,rightPlus,face);
		diffusive_face_flux(face,fluxPlus);
		
		if (doRightSourceJac) {
			sources(right,sourceRight,true);
			sources(rightPlus,sourceRightPlus,true);
		}
		
		jacobianRight=(fluxPlus-flux)/epsilon;
		if (doRightSourceJac) sourceJacRight=(sourceRightPlus-sourceRight)/epsilon;		

	}
	
	face_state_adjust(left,right,face);
}

void HeatConduction::left_state_update(HC_Cell_State &left,HC_Face_State &face) {
	
	int parent;
	
	parent=grid[gid].face[face.index].parent;

	Vec3D cell2face=grid[gid].face[face.index].centroid-grid[gid].cell[parent].centroid;

	left.T_center=T.cell(parent);
	left.T=left.T_center+cell2face.dot(gradT.cell(parent));
	left.update=update.cell(parent);
	left.volume=grid[gid].cell[parent].volume;
	
	return;
}

void HeatConduction::right_state_update(HC_Cell_State &left,HC_Cell_State &right,HC_Face_State &face) {

	if (face.bc==INTERNAL_FACE) {// internal face

		int neighbor=grid[gid].face[face.index].neighbor;
		
		Vec3D cell2face=grid[gid].face[face.index].centroid-grid[gid].cell[neighbor].centroid;

		right.T_center=T.cell(neighbor);
		right.T=right.T_center+cell2face.dot(gradT.cell(neighbor));
		right.update=update.cell(neighbor);
		right.volume=grid[gid].cell[neighbor].volume;
		
	} else if (face.bc>=0) { // boundary face

		right.update=0.;
		if (bc[gid][face.bc].thermalType==FIXED_T) {
			right.T=T.bc(face.bc,face.index);
			right.T_center=2.*right.T-left.T_center;
		}
		// qdot BC's are taken care of in the diffusive_face_flux fuction
		
	} else { // partition boundary

		int g=-1*grid[gid].face[face.index].neighbor-1; // ghost cell index
		
		Vec3D cell2face=grid[gid].face[face.index].centroid-grid[gid].ghost[g].centroid;

		right.T_center=T.ghost(g);
		right.T=right.T_center+cell2face.dot(gradT.ghost(g));
		right.update=update.ghost(g);
		right.volume=grid[gid].ghost[g].volume;
	}
	
	return;
}

void HeatConduction::face_geom_update(HC_Face_State &face,int f) {
	face.index=f;
	face.normal=grid[gid].face[f].normal;
	face.tangent1=(grid[gid].face[f].centroid-grid[gid].faceNode(f,0));
	face.tangent1/=fabs(face.tangent1);
	// Cross the tangent vector with the normal vector to get the second tangent
	face.tangent2=(face.normal).cross(face.tangent1);
	face.area=grid[gid].face[f].area;
	face.bc=grid[gid].face[f].bc;
	Vec3D leftCentroid,rightCentroid;
	leftCentroid=grid[gid].cell[grid[gid].face[f].parent].centroid;
	if (face.bc==INTERNAL_FACE) { rightCentroid=grid[gid].cell[grid[gid].face[f].neighbor].centroid;}
	else if (face.bc>=0) {
		rightCentroid=leftCentroid+2.*(grid[gid].face[f].centroid-leftCentroid).dot(face.normal)*face.normal;
	} else { rightCentroid=grid[gid].ghost[-1*grid[gid].face[f].neighbor-1].centroid;}
	face.left2right=rightCentroid-leftCentroid;
	return;
} // end face_geom_update

void HeatConduction::face_state_update(HC_Cell_State &left,HC_Cell_State &right,HC_Face_State &face) {

	face.gradT=gradT.face(face.index);
	
	face.gradT-=face.gradT.dot(face.normal)*face.normal;
	face.gradT+=((right.T_center-left.T_center)/(face.left2right.dot(face.normal)))*face.normal;

	// Boundary conditions are already taken care of in right state update
	// TODO: In first order, this won't be a good averaging
	face.T=0.5*(left.T+right.T);
	face.lambda=material.therm_cond(face.T);
	
	return;
} // end face_state_update

void HeatConduction::state_perturb(HC_Cell_State &state,HC_Face_State &face,double epsilon) {
	
	state.T+=epsilon;
	state.T_center+=epsilon;
	
	return;
} // end state_perturb

void HeatConduction::face_state_adjust(HC_Cell_State &left,HC_Cell_State &right,HC_Face_State &face) {

	face.T=0.5*(left.T+right.T);
	face.gradT-=face.gradT.dot(face.normal)*face.normal;
	face.gradT+=((right.T_center-left.T_center)/(face.left2right.dot(face.normal)))*face.normal;
	face.lambda=material.therm_cond(face.T);
	
	return;
} // end face_state_adjust
