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
#ifndef FLAMELET_H
#define FLAMELET_H

#include "petscksp.h"
#include <mpi.h>

#include "commons.h"
#include "bc.h"

extern IndexMaps maps;
extern BC bc;

class Flamelet_Constants {
	public:
	double sigma_t,Cg,Cd;
};

class Flamelet_Cell {
	public:
	double Z,Zvar;
	double update[2];
	Vec3D grad[2]; // k and omega gradients
};

class Flamelet_Ghost {
	public:
	double Z,Zvar;
	Vec3D grad[2]; // k and omega gradients
};

class Flamelet {
	public:
	std::vector<Flamelet_Cell> cell;
	std::vector<Flamelet_Ghost> ghost;
	Flamelet_Constants constants;
	KSP ksp; // linear solver context
	Vec deltaU,rhs; // solution, residual vectors
	Mat impOP; // implicit operator matrix
	MPI_Datatype MPI_GHOST;
	MPI_Datatype MPI_GRAD;
	Flamelet(void); // constructor
	void allocate(void); // allocate memory for cell and face variables
	void petsc_init(double rtol,double abstol,int maxits);
	void petsc_solve(int &nIter, double &rNorm);
	void petsc_destroy(void);
	void mpi_init(void);
	void mpi_update_ghost(void);
	void mpi_update_ghost_gradients(void);
	void gradients(void);
	void limit_gradients(void);
	void terms(void);
	void get_Z_Zvar(unsigned int &parent,unsigned int &neighbor,unsigned int &f,
				double &leftZ,double &leftZvar,
      				double &rightZ,double &rightZvar,
      				Vec3D &faceGradZ,Vec3D &faceGradZvar,Vec3D &left2right);
	void update(double &resZ, double &resZvar);
};

#endif
