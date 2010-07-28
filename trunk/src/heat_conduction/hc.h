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
#ifndef HC_H
#define HC_H

#include "petscksp.h"

#include "vec3d.h"
#include "grid.h"
#include "inputs.h"
#include "variable.h"
#include "utilities.h"
#include "bc.h"
#include "hc_state_cache.h"
#include "commons.h"
#include "bc_interface.h"
#include "material.h"

extern InputFile input;
extern vector<Grid> grid;
extern vector<vector<BCregion> > bc;
extern vector<Variable<double> > dt;
extern vector<vector<BC_Interface> > interface; // for each grid

// Class for Navier-Stokes equations
class HeatConduction {
public:
	int gid; // Grid id
	int nVars;
	int Rank,np; // Current processors index and total number of processors
	
	// Inputs
	double rtol,abstol;
	int maxits;
	double sqrt_machine_error;
	
	// Total residuals
	double resT;
	
	// Scalar variables
	Variable<double> T,update,qdot;
	// Vector variables
	Variable<Vec3D> gradT;
	
	MATERIAL material;
	
	// PETSC variables
	KSP ksp; // linear solver context
	Vec deltaU,rhs; // solution, residual vectors
	Mat impOP; // implicit operator matrix
	
	HeatConduction (void); // Empty constructor

	void initialize(void);
	void create_vars(void);
	void apply_initial_conditions(void);
	void mpi_init(void);
	void mpi_update_ghost_primitives(void);
	void mpi_update_ghost_gradients(void);
	void calc_cell_grads (void);
	void set_bcs(void);
	
	void petsc_init(void);
	void petsc_solve(int &nIter,double &rNorm);
	void petsc_destroy(void);
	
	void solve(int timeStep);
	
	void initialize_linear_system();
	void assemble_linear_system(void);

	void get_jacobians(void);
	void diffusive_face_flux(HC_Face_State &face,double &flux);
	void left_state_update(HC_Cell_State &left,HC_Face_State &face);
	void right_state_update(HC_Cell_State &left,HC_Cell_State &right,HC_Face_State &face);
	void face_geom_update(HC_Face_State &face,int f);
	void face_state_update(HC_Cell_State &left,HC_Cell_State &right,HC_Face_State &face);
	void state_perturb(HC_Cell_State &state,HC_Face_State &face,double epsilon);
	void face_state_adjust(HC_Cell_State &left,HC_Cell_State &right,HC_Face_State &face);
	void sources(HC_Cell_State &state,double &source,bool forJacobian=false);

	void apply_bcs(HC_Cell_State &left,HC_Cell_State &right,HC_Face_State &face);
	void velocity_inlet(HC_Cell_State &left,HC_Cell_State &right,HC_Face_State &face);
	void mdot_inlet(HC_Cell_State &left,HC_Cell_State &right,HC_Face_State &face);
	void outlet(HC_Cell_State &left,HC_Cell_State &right,HC_Face_State &face);
	void noslip(HC_Cell_State &left,HC_Cell_State &right,HC_Face_State &face);
	void slip(HC_Cell_State &left,HC_Cell_State &right,HC_Face_State &face);
	void symmetry(HC_Cell_State &left,HC_Cell_State &right,HC_Face_State &face);
	
	void update_variables(void);
	void write_restart(int timeStep);
	void read_restart(int restart_step,vector<vector<int> > &partitionMap);
};

#endif
