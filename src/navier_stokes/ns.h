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
#ifndef NS_H
#define NS_H

#include "petscksp.h"

#include "vec3d.h"
#include "grid.h"
#include "inputs.h"
#include "variable.h"
#include "utilities.h"
#include "material.h"
#include "bc.h"
#include "ns_state_cache.h"
#include "commons.h"
#include "bc_interface.h"
#include "loads.h"

#define NONE -1
// Options for limiter
#define VK 1
#define BJ 2
#define MINMOD 3
// Options for order of accuracy
#define FIRST 1
#define SECOND 2
// Options for convective_flux_function
#define ROE 1
#define VAN_LEER 2
#define AUSM_PLUS_UP 3
#define SD_SLAU 4
#define SW 5
// Options for preconditioner
#define WS95 1

extern InputFile input;
extern vector<Grid> grid;
extern vector<vector<BCregion> > bc;
extern vector<Variable<double> > dt;
extern vector<Variable<double> > dtau;
extern vector<vector<BC_Interface> > interface; // for each grid
extern vector<bool> turbulent;
extern vector<Loads> loads;

// Class for Navier-Stokes equations
class NavierStokes {
public:
	int gid; // Grid id
	int nVars;
	int Rank,np; // Current processors index and total number of processors
	int timeStep;
	int ps_step_max;
	int ps_step;
	int nIter;
	double rNorm,res,ps_res;
	double qmax[5],qmin[5];
      
	// MPI exchange buffer
	vector<double> sendBuffer;
	vector<double> recvBuffer;
	vector<int> mpi_send_offset,mpi_recv_offset;
	int send_req_count,recv_req_count;
	vector<MPI_Request> send_request, recv_request;
	
	// Inputs
	double rtol,abstol;
	int maxits;
	int order,jac_order;
	int limiter_function;
	int convective_flux_function;
	double limiter_threshold;
	double Minf;
	int preconditioner;
	double wdiss,bl_height;
	
	double small_number;
	double order_factor;
	
	// Total residuals
	vector<double> first_residuals,first_ps_residuals;
	
	// Scalar variables
	Variable<double> rho,p,T,qdot,mdot,weightL,p_total,T_total;
	// Vector variables
	Variable<Vec3D> V,gradu,gradv,gradw,gradp,gradT,tau;
	vector<Variable<double> > update,limiter;

	MATERIAL material;
	
	// PETSC variables
	KSP ksp; // linear solver context
	PC pc; // preconditioner context
	Vec deltaU,rhs; // solution, residual vectors
	Mat impOP; // implicit operator matrix
	Vec soln_n; // Solution at n level
	Vec pseudo_delta; // u^k-u^n
	Vec pseudo_right;
	
	NavierStokes (void); // Empty constructor
	// TODO: sort the following list of functions in the proper order of application
	void initialize(int ps_step_max);
	void create_vars(void);
	void apply_initial_conditions(void);
	void mpi_init(void);
	void mpi_update_ghost_primitives(void);
	void mpi_update_ghost_gradients(void);
	void calc_cell_grads (void);
	void set_bcs(void);
	void set_interfaces(void);

	void petsc_init(void);
	void petsc_solve(void);
	void petsc_destroy(void);
	
	void calc_limiter(void);
	void venkatakrishnan_limiter(void); 
	void barth_jespersen_limiter(void);
	void minmod_limiter(void);
		
	void solve(int timeStep,int ps_step);
	
	void cons2prim(int cid,vector<vector<double> > &P);
	void preconditioner_ws95(int c,vector<vector<double> > &P);
	void assemble_linear_system(void);
	void time_terms(void);
	void left_state_update(NS_Cell_State &left,NS_Face_State &face);
	void right_state_update(NS_Cell_State &left,NS_Cell_State &right,NS_Face_State &face);
	void face_geom_update(NS_Face_State &face,int f);
	void face_state_update(NS_Cell_State &left,NS_Cell_State &right,NS_Face_State &face);
	void face_state_adjust(NS_Cell_State &left,NS_Cell_State &right,NS_Face_State &face,int var);
	void state_perturb(NS_Cell_State &state,NS_Face_State &face,int var,double epsilon);
	void convective_face_flux(NS_Cell_State &left,NS_Cell_State &right,NS_Face_State &face,double flux[]);
	void diffusive_face_flux(NS_Cell_State &left,NS_Cell_State &right,NS_Face_State &face,double flux[]);
	void sources(NS_Cell_State &state,double source[],bool forJacobian=false);
	void get_jacobians(const int var);
	void apply_bcs(NS_Cell_State &left,NS_Cell_State &right,NS_Face_State &face);
	void velocity_inlet(NS_Cell_State &left,NS_Cell_State &right,NS_Face_State &face);
	void mdot_inlet(NS_Cell_State &left,NS_Cell_State &right,NS_Face_State &face);
	void stagnation_inlet(NS_Cell_State &left,NS_Cell_State &right,NS_Face_State &face);
	void outlet(NS_Cell_State &left,NS_Cell_State &right,NS_Face_State &face);
	void wall(NS_Cell_State &left,NS_Cell_State &right,NS_Face_State &face,bool slip=false);
	void symmetry(NS_Cell_State &left,NS_Cell_State &right,NS_Face_State &face);
	void update_variables(void);
	void update_boundaries(void);
	void write_restart(int timeStep);
	void read_restart(int restart_step,vector<vector<int> > &partitionMap);
	void find_min_max (void);

};

#endif
