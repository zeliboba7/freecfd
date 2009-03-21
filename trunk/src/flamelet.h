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

class Species {
	public:
	string name;
	double Mw;
	Species(string name_in);
};

class Flamelet_Table {
	public:
	int nZ,nZvar,nChi;
	int i1,i2,i3;  // indeicies for Z,Zvar and Chi
	vector<double> weights;
	std::vector<double> Z;
	std::vector<double> Zvar;
	std::vector<double> Chi;
	std::vector<vector<vector<double> > > rho;
	std::vector<vector<vector<double> > > T;
	std::vector<vector<vector<double> > > mu;
	std::vector<vector<vector<double> > > diffusivity;
	std::vector<vector<vector<vector<double> > > > Y; // mass fractions
	std::vector<Species> species;
	void read(string fileName);
	void get_weights(double &Z_in, double &Zvar_in, double &Chi_in);
	double get_rho(double &Z_in, double &Zvar_in, double &Chi_in,bool refreshWeights=true);
	double get_T(double &Z_in, double &Zvar_in, double &Chi_in,bool refreshWeights=true);
	void get_rho_T_comp(double &p_in, double &Z_in, double &Zvar_in, double &Chi_in,double &rho_out, double &T_out);
	double get_mu(double &Z_in, double &Zvar_in, double &Chi_in,bool refreshWeights=true);
	double get_diffusivity(double &Z_in, double &Zvar_in, double &Chi_in,bool refreshWeights=true);
	double get_Mw(double &Z_in, double &Zvar_in, double &Chi_in,bool refreshWeights=true);
	double get_drho_dZ(double &Z_in, double &Zvar_in, double &Chi_in);
	double get_drho_dZvar(double &Z_in, double &Zvar_in, double &Chi_in);
};

class Flamelet_Constants {
	public:
	double sigma_t,Cg,Cd;
};

class Flamelet_Face {
	public:
	double mu,RL,RR; // R is the gas constant
};

class Flamelet_Cell {
	public:
	double Z,Zvar,Chi,mu,diffusivity,R; // R is the gas constant
	double update[2];
	Vec3D grad[2]; // Z and Zvar gradients
};

class Flamelet_Ghost {
	public:
	double Z,Zvar,Chi,mu,diffusivity,R; // R is the gas constant
	Vec3D grad[2]; // Z and Zvar gradients
};

class Flamelet {
	public:
	std::vector<Flamelet_Cell> cell;
	std::vector<Flamelet_Face> face;
	std::vector<Flamelet_Ghost> ghost;
	Flamelet_Constants constants;
	Flamelet_Table table;
	double relaxation;
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
      				Vec3D &faceGradZ,Vec3D &faceGradZvar,Vec3D &left2right,
	  			bool &extrapolated);
	void update(double &resZ, double &resZvar);
	void lookup(void);
	void update_face_properties(void);
};

#endif
