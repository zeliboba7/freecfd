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
#ifndef RANS_H
#define RANS_H

#include "petscksp.h"
#include <mpi.h>

#include "commons.h"
#include "bc.h"

extern IndexMaps maps;
extern BC bc;

class RANS_Model {
	public:
	double sigma_k,sigma_omega,beta,beta_star,kappa,alpha;
};

class RANS_Face {
	public:
	double mu_t;
};

class RANS_Cell {
	public:
	double k,omega,mu_t,strainRate,filterFunction;
	double update[2];
	Vec3D grad[2]; // k and omega gradients
};

class RANS_Ghost {
	public:
	double k,omega,mu_t;
	Vec3D grad[2]; // k and omega gradients
};

class RANS {
	public:
	std::vector<RANS_Face> face;
	std::vector<RANS_Cell> cell;
	std::vector<RANS_Ghost> ghost;
	RANS_Model kepsilon,komega;
	KSP ksp; // linear solver context
	Vec deltaU,rhs; // solution, residual vectors
	Mat impOP; // implicit operator matrix
	MPI_Datatype MPI_GHOST;
	MPI_Datatype MPI_GRAD;
	RANS(void); // constructor
	void allocate(void); // allocate memory for cell and face variables
	void petsc_init(double rtol,double abstol,int maxits);
	void petsc_solve(int &nIter, double &rNorm);
	void petsc_destroy(void);
	void mpi_init(void);
	void mpi_update_ghost(void);
	void mpi_update_ghost_gradients(void);
	void update_cell_eddy_viscosity(void); // Updates eddy viscosity stored at the cell centers
	void update_face_eddy_viscosity(void); // Updates eddy viscosity stored at the faces
	void gradients(void);
	void limit_gradients(void);
	void terms(void);
	void update(double &resK, double &resOmega);
};

#endif
