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
#ifndef RANS_H
#define RANS_H

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
#include "ns.h"

extern InputFile input;
extern vector<Grid> grid;
extern vector<vector<BCregion> > bc;
extern vector<Variable<double> > dt;
extern vector<vector<BC_Interface> > interface; 
extern vector<NavierStokes> ns;

// Options for order of accuracy
#define FIRST 1
#define SECOND 2
// Options for turbulence model
#define KEPSILON 1
#define KOMEGA 2
#define BSL 3
#define SST 4

class RANS_Model {
public:
	double sigma_k,sigma_omega,beta,beta_star,kappa,alpha;
};

class RANS {
public:
	int gid; // Grid id
	int order,model;
	int Rank,np;
	int nVars;
	Variable<double> k,omega,mu_t,strainRate,yplus;
	Variable<Vec3D> gradk,gradomega;
	vector<Variable<double> > update,limiter;
	double kLowLimit,kHighLimit,omegaLowLimit,viscosityRatioLimit;
	double Pr_t;
	double rtol,abstol;
	int maxits;
	int timeStep;
	int ps_step_max;
	int ps_step;	
	int nIter;
	double rNorm,res,ps_res;
	// Total residuals
	vector<double> first_residuals,first_ps_residuals;

    // MPI exchange buffer
    vector<double> sendBuffer;
    vector<double> recvBuffer;
    vector<int> mpi_send_offset,mpi_recv_offset;
	int send_req_count,recv_req_count;
	vector<MPI_Request> send_request, recv_request;
	
	// PETSC variables
	KSP ksp; // linear solver context
	PC pc; // preconditioner context
	Vec deltaU,rhs; // solution, residual vectors
	Mat impOP; // implicit operator matrix
	Vec soln_n; // Solution at n level
	Vec pseudo_delta; // u^k-u^n
	Vec pseudo_right;
	
	MATERIAL material;
	
	RANS (void); 
	void initialize(int ps_step_max);
	void mpi_init(void);
	void create_vars (void);
	void apply_initial_conditions (void);
	void set_bcs (void);
	void mpi_update_ghost_primitives(void);
	void calc_cell_grads(void);
	void mpi_update_ghost_gradients(void);
	void petsc_init(void);
	void petsc_solve(void);
	void petsc_destroy(void);
	void calc_limiter(void);
	void barth_jespersen_limiter(void);
	void solve (int timeStep,int ps_step);
	void initialize_linear_system(void);
	void get_kOmega(void);
	void terms(void);
	void time_terms(void);
	void update_eddy_viscosity(void);
	void update_variables(void);
	void write_restart(int timeStep);
	void read_restart(int restart_step,vector<vector<int> > &partitionMap);
	
};

#endif
