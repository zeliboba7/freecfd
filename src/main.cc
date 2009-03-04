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

#include <iostream>
#include <fstream>
#include <map>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <limits>

#include "commons.h"
#include "rans.h"
#include "flamelet.h"
#include "inputs.h"
#include "bc.h"
#include "mpi_functions.h"
#include "petsc_functions.h"
#include "probe.h"
using namespace std;

// Function prototypes
void read_inputs(InputFile &input);
void initialize(InputFile &input);
void write_output(double time, InputFile input);
void write_restart(double time);
void read_restart(double &time);
void setBCs(InputFile &input, BC &bc);
bool withinBox(Vec3D point,Vec3D corner_1,Vec3D corner_2);
bool withinCylinder(Vec3D point,Vec3D center,double radius,Vec3D axisDirection,double height);
bool withinSphere(Vec3D point,Vec3D center,double radius);
double minmod(double a, double b);
double maxmod(double a, double b);
double doubleMinmod(double a, double b);
double harmonic(double a, double b);
double superbee(double a, double b);
void update(void);
string int2str(int number) ;
void assemble_linear_system(void);
void initialize_linear_system(void);
void update_face_mdot(void);
void linear_system_turb(void);
void get_dt(void);
void get_CFLmax(void);
void set_probes(void);
void set_integrateBoundary(void);
inline double signof(double a) { return (a == 0.) ? 0. : (a<0. ? -1. : 1.); }
double temp;

// Allocate container for boundary conditions
BC bc;
InputFile input;

// Turbulence model needs to be declared regardless of being used or not
// But memory won't be allocated for cell,face,etc.. variables for laminar cases
// So it is OK
RANS rans;
// Same goes for the flamelet model
Flamelet flamelet;

vector<Probe> probes;
vector<BoundaryFlux> boundaryFluxes;

double resP,resV,resT,resK,resOmega,resZ,resZvar;

int main(int argc, char *argv[]) {
	
	// Set global parameters
	sqrt_machine_error=sqrt(std::numeric_limits<double>::epsilon());
	
	// Initialize mpi
	mpi_init(argc,argv);

	string inputFileName, gridFileName;
	inputFileName.assign(argv[1]);
	gridFileName.assign(argv[1]);
	inputFileName+=".in";
	gridFileName+=".cgns";
	//gridFileName+=".dat";

	restart=0;
	if (argc>2) restart=atoi(argv[2]);

	input.setFile(inputFileName);
	read_inputs(input);
	// Set equation of state based on the inputs
	eos.set();

	// Physical time
	double time=0.;

	if (restart!=0 && ramp) {
		ramp=false;
	}
	int timeStepMax=input.section("timeMarching").get_int("numberOfSteps");
	
	grid.read(gridFileName);
	// Hand shake with other processors
	// (agree on which cells' data to be shared)
	mpi_handshake();
	mpi_get_ghost_centroids();
	
	if (TURBULENCE_MODEL!=NONE) {
		rans.allocate();
		rans.mpi_init();
		if (Rank==0) cout << "[I] Initialized turbulence model" << endl;
	}
	
	if (FLAMELET) {
		flamelet.allocate();
		flamelet.mpi_init();
		if (Rank==0) cout << "[I] Initialized flamelet model" << endl;
		flamelet.table.read(input.section("flamelet").get_string("tableFile"));
		exit(1);
		if (Rank==0) cout << "[I] Read flamelet pdf table" << endl;
	}
	
	initialize(input);
	if (Rank==0) cout << "[I] Applied initial conditions" << endl;

	setBCs(input,bc);
	if (Rank==0) cout << "[I] Set boundary conditions" << endl;

	set_probes();
	set_integrateBoundary();

	grid.lengthScales();
	
	if (Rank==0) cout << "[I] Calculating node averaging metrics, this might take a while..." << endl;
	grid.nodeAverages(); // Linear triangular (tetrahedral) + idw blended mode
	//grid.nodeAverages_idw(); // Inverse distance based mode
	grid.faceAverages();
	grid.gradMaps();
			
	if (restart!=0) {
		for (int p=0;p<np;++p) {
			if (Rank==p) read_restart(time);
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}

	// Initialize petsc
	petsc_init(argc,argv,
			input.section("linearSolver").get_double("relTolerance"),
			input.section("linearSolver").get_double("absTolerance"),
			input.section("linearSolver").get_int("maxIterations"));
	
	if (TURBULENCE_MODEL!=NONE) {
		// Initialize petsc for turbulence model
		rans.petsc_init(input.section("linearSolver").get_double("relTolerance"),
					input.section("linearSolver").get_double("absTolerance"),
					input.section("linearSolver").get_int("maxIterations"));
	}
	
	if (FLAMELET) {
		// Initialize petsc for turbulence model
		flamelet.petsc_init(input.section("linearSolver").get_double("relTolerance"),
				input.section("linearSolver").get_double("absTolerance"),
				input.section("linearSolver").get_int("maxIterations"));
	}
	
	// Initialize time step or CFL size if ramping
	if (ramp) {
		if (TIME_STEP_TYPE==FIXED) dt=ramp_initial;
		if (TIME_STEP_TYPE==CFL_MAX) CFLmax=ramp_initial;
		if (TIME_STEP_TYPE==CFL_LOCAL) CFLlocal=ramp_initial;
	}

	if (Rank==0) cout << "[I] Beginning time loop" << endl;
	// Open residual file 
	ofstream iterFile;
	iterFile.open ("linear_solver.log", ios::out);
	iterFile.precision(3); iterFile << left << scientific;
	if (Rank==0) cout << "[I] Linear solver details are logged in linear_solver.log" << endl;
	// Begin timing
	MPI_Barrier(MPI_COMM_WORLD);
	double timeRef, timeEnd;
	timeRef=MPI_Wtime();

	/*****************************************************************************************/
	// Begin time loop
	/*****************************************************************************************/
	bool firstTimeStep=true;
	for (timeStep=restart+1;timeStep<=timeStepMax+restart;++timeStep) {

		int nIter,nIterTurb,nIterFlame; // Number of linear solver iterations
		double rNorm,rNormTurb,rNormFlame; // Residual norm of the linear solver

		// Flush boundary forces
		if ((timeStep) % integrateBoundaryFreq == 0) {
			for (int i=0;i<bcCount;++i) {
				bc.region[i].mass=0.;
				bc.region[i].momentum=0.;
				bc.region[i].energy=0.;
			}
		}
		
		if (firstTimeStep) {
			mpi_update_ghost_primitives();
			// Gradient are calculated for NS and/or second order schemes
			if (EQUATIONS==NS | order==SECOND ) {
			// Calculate all the cell gradients for each variable
				grid.gradients();
			// Update gradients of the ghost cells
				mpi_update_ghost_gradients();
			}
		
			if (LIMITER!=NONE) {
			// Limit gradients
				grid.limit_gradients(); 
			// Update limited gradients of the ghost cells
				mpi_update_ghost_gradients();
			}
		}


		// Get time step (will handle fixed, CFL and CFLramp)
		get_dt();
		if (TIME_STEP_TYPE==FIXED) { get_CFLmax(); } 
		else if (TIME_STEP_TYPE==CFL_LOCAL) { CFLmax=CFLlocal;}
		
		if (FLAMELET) {
			if (firstTimeStep) {
				flamelet.mpi_update_ghost();
				flamelet.gradients();
				flamelet.mpi_update_ghost_gradients();
				flamelet.limit_gradients();
				flamelet.mpi_update_ghost_gradients(); // Called again to update the ghost gradients
				update_face_mdot();
			} else {
				update_face_mdot();
				flamelet.terms();
				flamelet.petsc_solve(nIterFlame,rNormFlame);
				flamelet.update(resZ,resZvar);
				flamelet.update_rho(); // TODO Update the ghost rho's in flamelet mpi structure
				mpi_update_ghost_primitives();
				flamelet.mpi_update_ghost();
				flamelet.gradients();
				flamelet.mpi_update_ghost_gradients();
				flamelet.limit_gradients();
				flamelet.mpi_update_ghost_gradients(); // Called again to update the ghost gradients
			}
		}
		
		
		if (TURBULENCE_MODEL!=NONE) {
			if (firstTimeStep) {
				rans.mpi_update_ghost();
				rans.update_eddy_viscosity();
			} else {
				rans.gradients();
				rans.mpi_update_ghost_gradients();
				rans.limit_gradients();
				rans.mpi_update_ghost_gradients(); // Called again to update the ghost gradients
				if (!FLAMELET) update_face_mdot(); // Already done if flamelet model is being used
				rans.terms();
				rans.petsc_solve(nIterTurb,rNormTurb);
				rans.update(resK,resOmega);
				rans.mpi_update_ghost();
				rans.update_eddy_viscosity();
			}
		}
				
		initialize_linear_system();
		assemble_linear_system();
		petsc_solve(nIter,rNorm);
		update();
		mpi_update_ghost_primitives();
		// Gradient are calculated for NS and/or second order schemes
		if (EQUATIONS==NS | order==SECOND ) {
			// Calculate all the cell gradients for each variable
			grid.gradients();
			// Update gradients of the ghost cells
			mpi_update_ghost_gradients();
		}
		
		if (LIMITER!=NONE) {
			// Limit gradients
			grid.limit_gradients(); 
			// Update limited gradients of the ghost cells
			mpi_update_ghost_gradients();
		}

		// Advance physical time
		time += dt;
		if (Rank==0) {
			// Output residual stream labels
			cout.precision(3);
			cout << left << scientific;
			
 			if (timeStep==(restart+1))  {
				// Write screen labels
				cout << setw(8) << "step";
				cout << setw(12) << "time";
				cout << setw(12) << "dt";
				cout << setw(12) << "CFLmax";
				cout << setw(12) << "resP";
				cout << setw(12) << "resV";
				cout << setw(12) << "resT";
				if (TURBULENCE_MODEL!=NONE) cout << setw(12) << "resK" << setw(12) << "resOmega";
				if (FLAMELET) cout << setw(12) << "resZ" << setw(12) << "resZvar";
				cout << endl;	
 			}

			cout << setw(8) << timeStep;
			cout << setw(12) << time;
			cout << setw(12) << dt;
			cout << setw(12) << CFLmax;
			cout << setw(12) << resP;
			cout << setw(12) << resV;
			cout << setw(12) << resT;
			iterFile << timeStep << "\t" << nIter << "\t" << rNorm;
			
			if (!firstTimeStep && TURBULENCE_MODEL!=NONE) {
				cout << setw(12) << resK << setw(12) << resOmega;
				iterFile << "\t" << nIterTurb << "\t" << rNormTurb;
			}
			if (!firstTimeStep && FLAMELET) {
				cout << setw(12) << resZ << setw(12) << resZvar;
				iterFile << "\t" << nIterFlame << "\t" << rNormFlame;
			}
			cout << endl; iterFile << endl;
			
			
			
		}
		
		// Ramp-up if needed
		if (ramp) {
			if (TIME_STEP_TYPE==FIXED) dt=min(dt*ramp_growth,dtTarget);
			if (TIME_STEP_TYPE==CFL_MAX) CFLmax=min(CFLmax*ramp_growth,CFLmaxTarget);
			if (TIME_STEP_TYPE==CFL_LOCAL) CFLlocal=min(CFLlocal*ramp_growth,CFLlocalTarget);
		}
			
		if ((timeStep) % restartFreq == 0) {
			for (int p=0;p<np;p++) {
				if (p==Rank) write_restart(time);
				// Syncronize the processors
				MPI_Barrier(MPI_COMM_WORLD);
			}
		}

		if ((timeStep) % outFreq == 0) write_output(time,input);
		
		if ((timeStep) % probeFreq == 0) {

			for (int p=0;p<probes.size();++p) {
				Cell &c=grid.cell[probes[p].nearestCell];
				ofstream file;
				file.open((probes[p].fileName).c_str(),ios::app);
				file << timeStep << "\t" << setw(16) << setprecision(8)  << time << "\t" << c.p << "\t" << c.v[0] << "\t" << c.v[1] << "\t" << c.v[2] << "\t" << c.T << "\t" << c.rho << "\t" ;
				if (TURBULENCE_MODEL!=NONE) {
					unsigned int cc=probes[p].nearestCell;
					file << rans.cell[cc].k << "\t" << rans.cell[cc].omega << "\t" << c.rho*rans.cell[cc].k/rans.cell[cc].omega ;
				}
				if (FLAMELET) {
					unsigned int cc=probes[p].nearestCell;
					file << flamelet.cell[cc].Z << "\t" << flamelet.cell[cc].Zvar ;
				}
				file << endl;
				file.close();
			}
		}

		if ((timeStep) % integrateBoundaryFreq == 0) {
			for (int n=0;n<boundaryFluxes.size();++n) {
				// Sum integrated boundary fluxes accross partitions
				double myFlux[5],integratedFlux[5];
				myFlux[0]=bc.region[boundaryFluxes[n].bc-1].mass;
				for (int i=0;i<3;++i) myFlux[i+1]=bc.region[boundaryFluxes[n].bc-1].momentum[i];
				myFlux[4]=bc.region[boundaryFluxes[n].bc-1].energy;				
				if (np!=1) {
					MPI_Reduce(&myFlux,&integratedFlux,5, MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
					MPI_Barrier(MPI_COMM_WORLD);
				} else {for (int i=0;i<5;++i) integratedFlux[i]=myFlux[i]; }
				// TODO Instead of data transfer between my and integrated 
				// Try MPI_IN_PLACE, also is the Barrier need here?
				if (Rank==0) {
					ofstream file;
					file.open((boundaryFluxes[n].fileName).c_str(),ios::app);
					file << timeStep << setw(16) << setprecision(8)  << "\t" << time;
					for (int i=0;i<5;++i) file << "\t" << integratedFlux[i];
					file << endl;
					file.close();
				}
			}
		}
		firstTimeStep=false;
	
	}

	iterFile.close();
	/*****************************************************************************************/
	// End time loop
	/*****************************************************************************************/
	// Syncronize the processors
	MPI_Barrier(MPI_COMM_WORLD);

	// Report the wall time
	if (Rank==0) {
		timeEnd=MPI_Wtime();
		cout << "* Wall time: " << timeEnd-timeRef << " seconds" << endl;
	}
	
	if (TURBULENCE_MODEL!=NONE) rans.petsc_destroy();
	if (FLAMELET) flamelet.petsc_destroy();
	petsc_finalize();

	return 0;
}

double minmod(double a, double b) {
	if ((a*b)<0.) {
		return 0.;
	} else if (fabs(a) < fabs(b)) {
		return a;
	} else {
		return b;
	}
}

double maxmod(double a, double b) {
	if ((a*b)<0.) {
		return 0.;
	} else if (fabs(a) < fabs(b)) {
		return b;
	} else {
		return a;
	}
}

double doubleMinmod(double a, double b) {
	if ((a*b)<0.) {
		return 0.;
	} else {
		return signof(a+b)*min(fabs(2.*a),min(fabs(2.*b),fabs(0.5*(a+b))));
	}
	
}

double harmonic(double a, double b) {
	if ((a*b)<0.) {
		return 0.;
	} else {
		if (fabs(a+b)<1.e-12) return 0.5*(a+b);
		return signof(a+b)*min(fabs(0.5*(a+b)),fabs(2.*a*b/(a+b)));
	}
}

double superbee(double a, double b) {
	return minmod(maxmod(a,b),minmod(2.*a,2.*b));
}


void update(void) {

	resP=0.; resV=0.; resT=0.;
	for (unsigned int c = 0;c < grid.cellCount;++c) {
		grid.cell[c].p +=grid.cell[c].update[0];
		grid.cell[c].v[0] +=grid.cell[c].update[1];
		grid.cell[c].v[1] +=grid.cell[c].update[2];
		grid.cell[c].v[2] +=grid.cell[c].update[3];
		if (!FLAMELET) {
			grid.cell[c].T += grid.cell[c].update[4];
			grid.cell[c].rho=eos.rho(grid.cell[c].p,grid.cell[c].T);
		}
		
		resP+=grid.cell[c].update[0]*grid.cell[c].update[0];
		resV+=grid.cell[c].update[1]*grid.cell[c].update[1]+grid.cell[c].update[2]*grid.cell[c].update[2]+grid.cell[c].update[3]*grid.cell[c].update[3];
		resT+=grid.cell[c].update[4]*grid.cell[c].update[4];
	} // cell loop
	double residuals[3],totalResiduals[3];
	residuals[0]=resP; residuals[1]=resV; residuals[2]=resT;

        if (np!=1) {
        	MPI_Reduce(&residuals,&totalResiduals,3, MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		resP=totalResiduals[0]; resV=totalResiduals[1]; resT=totalResiduals[2];
        }
	resP=sqrt(resP); resV=sqrt(resV); resT=sqrt(resT);
	return;
}

string int2str(int number) {
	std::stringstream ss;
	ss << number;
	return ss.str();
}

void get_dt(void) {
	
	if (TIME_STEP_TYPE==CFL_MAX) {
		// Determine time step with CFL condition
		dt=1.E20;
		for (cit=grid.cell.begin();cit<grid.cell.end();cit++) {
			double a=sqrt(Gamma*((*cit).p+Pref)/(*cit).rho);
			for (int i=0;i<3;++i) dt=min(dt,CFLmax*(*cit).lengthScale/(fabs((*cit).v[i])+a));
		}
		MPI_Allreduce(&dt,&dt,1, MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
	} 

	return;
} // end get_dt

void get_CFLmax(void) {

	CFLmax=0.;
	for (cit=grid.cell.begin();cit<grid.cell.end();cit++) {
		double a=sqrt(Gamma*((*cit).p+Pref)/(*cit).rho);
		//cout << Gamma <<  "\t" << (*cit).p << "\t" << Pref << "\t" << (*cit).rho << endl;
		for (int i=0;i<3;++i) CFLmax=max(CFLmax,(fabs((*cit).v[0])+a)*dt/(*cit).lengthScale);
	}
	MPI_Allreduce(&CFLmax,&CFLmax,1, MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

	return;
}

bool withinBox(Vec3D point,Vec3D corner_1,Vec3D corner_2) {
	for (int i=0;i<3;++i) {
		if (corner_1[i]>corner_2[i]) {
			if (point[i]>corner_1[i]) return false;
			if (point[i]<corner_2[i]) return false;
		} else {
			if (point[i]<corner_1[i]) return false;
			if (point[i]>corner_2[i]) return false;
		}
	}
	return true;
}

bool withinCylinder(Vec3D point,Vec3D center,double radius,Vec3D axisDirection,double height) {
	
	Vec3D onAxis=(point-center).dot(axisDirection)*axisDirection;
	if (fabs(onAxis)>0.5*height) return false;
	
	Vec3D offAxis=(point-center)-onAxis;
	if (fabs(offAxis)>radius) return false;
	
	return true;
}

bool withinSphere(Vec3D point,Vec3D center,double radius) {
	return (fabs(point-center)<=radius);
}
