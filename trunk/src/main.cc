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

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <mpi.h>
#include <vector>
#include "inputs.h"
using namespace std;
#include "grid.h"
#include "ns.h"
#include "hc.h"
#include "rans.h"
#include "bc.h"
#include "commons.h"
#include "bc_interface.h"
#include "loads.h"

// Function prototypes
void read_inputs(void);
void set_bcs(int gid);
void write_volume_output(int gid, int step);
void write_surface_output(int gid, int step);
void write_restart(int gid,int timeStep,double time);
void write_loads(int gid,int timeStep,double time);
void read_restart(int gid,int restart_step,double &time);
void set_lengthScales(int gid);

void set_time_step_options(void);
void update_time_step_options(void);
void update_time_step(int timeStep,double &time,double &max_cfl,int gid);

void set_pseudo_time_step_options(void);
void update_pseudo_time_step_options(void);
void update_pseudo_time_step(int ps_tep,double &max_cfl,int gid);

void bc_interface_sync(void);
void face_interpolation_weights(int gid);
void node_interpolation_weights(int gid);
void gradient_maps(int gid);
	
// Global declerations
InputFile input;
vector<InputFile> material_input;
vector<Grid> grid;
vector<vector<BCregion> > bc;
vector<NavierStokes> ns;
vector<HeatConduction> hc;
vector<RANS> rans;
vector<bool> turbulent;
// Time step for each grid
vector<Variable<double> > dt;
vector<Variable<double> > dtau;
vector<int> equations;
vector<vector<BC_Interface> > interface; // for each grid
vector<Loads> loads;

int Rank,np;
int gradient_test;
double min_x,max_x;
bool pseudo_time_active;

static char help[] = "Free CFD\n - A free general purpose computational fluid dynamics code";

// This is to suppress PETSc's error output
int PetscPrintError(const char error[],...){
	if (Rank==0) cerr << "PETSc Error ... exiting" << endl;
	exit(1);
	return 0;
}

int main(int argc, char *argv[]) {

	// Initialize mpi
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	
	// Initialize petsc
	PetscInitialize(&argc,&argv,(char *)0,help);
	PetscErrorPrintf = PetscPrintError;
	
	string inputFileName;
	inputFileName.assign(argv[1]);
	inputFileName+=".in";

	int restart_step=0;	
	bool PREP=false;
	bool OUTPUT_ONLY=false;
	if (argc>2) {
		string str_arg;
		str_arg=argv[2];
		if (str_arg=="restart") {
			restart_step=atoi(argv[3]);
		} else if (str_arg=="output") {
			OUTPUT_ONLY=true;
			if (argc>3) restart_step=atoi(argv[3]);
			else restart_step=0;
		} else if (str_arg=="prep") {
			PREP=true;
			remove("grid.raw");
		}
	}
	// Read the input file
	input.setFile(inputFileName);
	read_inputs();
	
	equations.resize(input.section("grid",0).count);
	turbulent.resize(input.section("grid",0).count);
	for (int gid=0;gid<input.section("grid",0).count;++gid) {
		if (input.section("grid",gid).get_string("equations")=="navierstokes") {
			equations[gid]=NS;
		} else if (input.section("grid",gid).get_string("equations")=="heatconduction") {
			equations[gid]=HEAT;
		} 
		if (input.section("grid",gid).subsection("turbulence").is_found) turbulent[gid]=true;
	}

	// Allocate space for each grid and equations and BC's on each grid
	grid.resize(input.section("grid",0).count);
	ns.resize(grid.size());
	rans.resize(grid.size());
	hc.resize(grid.size());
	bc.resize(grid.size());
	interface.resize(grid.size());
	
	// Read the grid and initialize
	for (int gid=0;gid<grid.size();++gid) {
		grid[gid].dimension=input.section("grid",gid).get_int("dimension");
		grid[gid].gid=gid;
		// Read the grid raw data from file
		grid[gid].read(input.section("grid",gid).get_string("file"),input.section("grid",gid).get_string("format"));
		if (PREP && Rank==0) {grid[gid].write_raw(); continue;}
		// Do the transformations
		int tcount=input.section("grid",gid).subsection("transform",0).count;
		
		for (int t=0;t<tcount;++t) {
			if (input.section("grid",gid).subsection("transform",t).get_string("function")=="translate") {
				Vec3D begin=input.section("grid",gid).subsection("transform",t).get_Vec3D("anchor");
				Vec3D end=input.section("grid",gid).subsection("transform",t).get_Vec3D("end");
				grid[gid].translate(begin,end);
			} else if (input.section("grid",gid).subsection("transform",t).get_string("function")=="scale") {
				Vec3D anchor=input.section("grid",gid).subsection("transform",t).get_Vec3D("anchor");
				Vec3D factor=input.section("grid",gid).subsection("transform",t).get_Vec3D("factor");
				grid[gid].scale(anchor,factor);
			} else if (input.section("grid",gid).subsection("transform",t).get_string("function")=="rotate") {
				Vec3D anchor=input.section("grid",gid).subsection("transform",t).get_Vec3D("anchor");
				Vec3D axis=input.section("grid",gid).subsection("transform",t).get_Vec3D("axis");
				double angle=input.section("grid",gid).subsection("transform",t).get_double("angle");
				grid[gid].rotate(anchor,axis,angle);
			} 
		}

		// Establish connectivity, area, volume etc... all the needed information
		grid[gid].setup();
		set_bcs(gid);
		
		set_lengthScales(gid);
		if (Rank==0) cout << "[I grid=" << gid+1 << " ] Calculating face averaging metrics" << endl;

		face_interpolation_weights(gid);
		node_interpolation_weights(gid);
		gradient_maps(gid);
	}

	if (PREP) return 0;
	
	for (int gid=0;gid<grid.size();++gid) {
		for (int i=0;i<interface[gid].size();++i) {
			interface[gid][i].setup();
		}
	}

	vector<double> time (grid.size(),0.);
	vector<double> max_cfl (grid.size(),0.);
	vector<double> ps_max_cfl (grid.size(),0.);
	
	int ps_step_max=input.section("pseudotime").get_int("numberofsteps");
	int ps_step;
	
	for (int gid=0;gid<grid.size();++gid) {
		if (equations[gid]==NS) {
			if (Rank==0) cout << "[I grid=" << gid+1 << " ] Initializing Navier Stokes solver" << endl; 
			ns[gid].gid=gid;
			ns[gid].initialize(ps_step_max);
			if (turbulent[gid]) {
				if (Rank==0) cout << "[I grid=" << gid+1 << " ] Initializing RANS solver" << endl; 
				rans[gid].gid=gid;
				rans[gid].initialize(ps_step_max);
			}
		}
		if (equations[gid]==HEAT) {
			if (Rank==0) cout << "[I grid=" << gid+1 << " ] Initializing Heat Conduction solver" << endl;
			hc[gid].gid=gid;
			hc[gid].initialize();
		}
		if (restart_step>0) read_restart(gid,restart_step,time[gid]);
	}

	set_time_step_options();
	set_pseudo_time_step_options();
	int timeStepMax=input.section("timemarching").get_int("numberofsteps");
	int time_step_update_freq=input.section("timemarching").get_int("updatefrequency");
	int pseudo_time_step_update_freq=input.section("pseudotime").get_int("updatefrequency");
	double ps_tolerance=input.section("pseudotime").get_double("residualtolerance");
	
	// Get the output frequency for each grid
	vector<int> volume_plot_freq,surface_plot_freq,restart_freq;
	loads.resize(grid.size());
	for (int gid=0;gid<grid.size();++gid) {
		volume_plot_freq.push_back(input.section("grid",gid).subsection("writeoutput").get_int("volumeplotfrequency"));
		surface_plot_freq.push_back(input.section("grid",gid).subsection("writeoutput").get_int("surfaceplotfrequency"));
		restart_freq.push_back(input.section("grid",gid).subsection("writeoutput").get_int("restartfrequency"));
		loads[gid].moment_center=input.section("grid",gid).subsection("writeoutput").get_Vec3D("momentcenter");
		loads[gid].frequency=input.section("grid",gid).subsection("writeoutput").get_int("loadfrequency");
		vector<string> temp;
		temp=input.section("grid",gid).subsection("writeoutput").get_stringList("includebcs");
		loads[gid].force.resize(temp.size());
		loads[gid].moment.resize(temp.size());
		for (int i=0;i<temp.size();++i) {
			loads[gid].include_bcs.push_back(atoi(temp[i].c_str())-1);
			loads[gid].force[i]=0.;
			loads[gid].moment[i]=0.;
		}
	}
	
	if (equations[0]==NS && gradient_test!=NONE) {
		// Dump gradient errors and exit
		if (Rank==0) cout << "\n[I] Writing gradient error output" << endl;
		input.section("grid",0).subsection("writeoutput").stringLists["volumevariables"].value.clear();
		input.section("grid",0).subsection("writeoutput").stringLists["volumevariables"].value.push_back("rank");
		input.section("grid",0).subsection("writeoutput").stringLists["volumevariables"].value.push_back("grad");
		input.section("grid",0).subsection("writeoutput").stringLists["volumevariables"].value.push_back("percent_grad_error");
		input.section("grid",0).subsection("writeoutput").stringLists["volumevariables"].value.push_back("volume");
		write_volume_output(0,0);
		exit(1);
	}
	
	// If turbulence model is not on, clear the turbulence variables from the output option list
	for (int gid=0;gid<grid.size();++gid) {
		if (!turbulent[gid]) {
			vector<string> *list;
			list= &input.section("grid",0).subsection("writeoutput").stringLists["volumevariables"].value;
			for (int i=0;i<(*list).size();++i) {
				if ((*list)[i]=="k" || (*list)[i]=="omega" || (*list)[i]=="mu_t" || (*list)[i]=="gradk" || (*list)[i]=="gradomega" || (*list)[i]=="yplus") {
						(*list)[i]="null";
				}
			}
			list= &input.section("grid",0).subsection("writeoutput").stringLists["surfacevariables"].value;
			for (int i=0;i<(*list).size();++i) {
				if ((*list)[i]=="k" || (*list)[i]=="omega" || (*list)[i]=="mu_t" || (*list)[i]=="gradk" || (*list)[i]=="gradomega" || (*list)[i]=="yplus") {
					(*list)[i]="null";
				}
			}
		}
	}
	
	if (OUTPUT_ONLY==true) {
		for (int gid=0;gid<grid.size();++gid) {
			if (Rank==0) cout << "[I] Writing surface output for grid=" << gid+1 << endl;
			write_surface_output(gid,restart_step);
			if (Rank==0) cout << "[I] Writing volume output for grid=" << gid+1 << endl;
			write_volume_output(gid, restart_step);
		}
		exit(0);
	}

	if (Rank==0) cout << "[I] Beginning time loop\n" << endl; 
	// Write out label for residuals
	if (Rank==0) {
		cout << "=============================================================================================================" << endl;
		cout << "step -- for each grid -> [grid-no  time  -- for each equation -> {cfl-max linear-iterations total-residual} ]" << endl;
		cout << "=============================================================================================================" << endl;
	}
	cout << setprecision(3) << scientific;
	fstream convergence;
	if (fexists("convergence.dat") && restart_step>0) convergence.open("convergence.dat",fstream::out | fstream::app);
	else convergence.open("convergence.dat",fstream::out);
	convergence << setprecision(3) << scientific;

	/*****************************************************************************************/
	// Begin time loop
	/*****************************************************************************************/
	MPI_Barrier(MPI_COMM_WORLD);
      	double timeRef,timeEnd;
	timeRef=MPI_Wtime();


	bool lastTimeStep=false;
	
	for (int timeStep=restart_step+1;timeStep<=timeStepMax+restart_step;++timeStep) {
		if (timeStep==(timeStepMax+restart_step)) lastTimeStep=true;
		for (int gid=0;gid<grid.size();++gid) {
			update_time_step(timeStep,time[gid],max_cfl[gid],gid);
			if (equations[gid]==NS) {
				for (ps_step=1;ps_step<=ps_step_max;++ps_step) {
					for (int b=0;b<loads[gid].include_bcs.size();++b) {
						loads[gid].force[b]=0.;
						loads[gid].moment[b]=0.;
					}
					if (ps_step_max>1) update_pseudo_time_step(ps_step,ps_max_cfl[gid],gid);
					 ns[gid].solve(timeStep,ps_step);
					// Write screen output for pseudo time iteration
					if (Rank==0 && ps_step_max>1) {
						cout        << "\t" << ps_step << "\t" << ps_max_cfl[gid] << "\t" << ns[gid].nIter << "\t" << ns[gid].ps_res;
						convergence << "\t" << ps_step << "\t" << ps_max_cfl[gid] << "\t" << ns[gid].nIter << "\t" << ns[gid].ps_res;
						if (turbulent[gid]) {
							cout        << "\t" << rans[gid].nIter << "\t" << rans[gid].ps_res;
							convergence << "\t" << rans[gid].nIter << "\t" << rans[gid].ps_res;

						}
						cout        << endl;
						convergence << endl;
					}
					if (ps_step_max>1) {
						if (ns[gid].ps_res < ps_tolerance) {
							if (!turbulent[gid]) break;
							else if (rans[gid].ps_res < ps_tolerance) break;
						}
					}
				}
			}
			if (equations[gid]==HEAT) hc[gid].solve(timeStep);
			bc_interface_sync();
			// Screen output
			if (Rank==0) {
				cout        << timeStep << "\t" << gid+1 << "\t" << time[gid];
				convergence << timeStep << "\t" << gid+1 << "\t" << time[gid];
				if (equations[gid]==NS) {
					cout        << "\t" << max_cfl[gid] << "\t" << ns[gid].nIter << "\t" << ns[gid].res;
					convergence << "\t" << max_cfl[gid] << "\t" << ns[gid].nIter << "\t" << ns[gid].res;
					if (turbulent[gid]) {
						cout        << "\t" << rans[gid].nIter << "\t" << rans[gid].res;
						convergence << "\t" << rans[gid].nIter << "\t" << rans[gid].res;
					}
				}
				if (equations[gid]==HEAT) {
					cout        << "\t" << hc[gid].nIter << "\t" << hc[gid].res;
					convergence << "\t" << hc[gid].nIter << "\t" << hc[gid].res;
				}
				cout        << endl;
				convergence << endl;
			}
			if (timeStep%volume_plot_freq[gid]==0 || lastTimeStep) {
				if (Rank==0) cout << "[I] Writing volume output for grid=" << gid+1 << endl;
				write_volume_output(gid,timeStep);
			} // end if
			if (timeStep%surface_plot_freq[gid]==0 || lastTimeStep) {
				if (Rank==0) cout << "[I] Writing surface output for grid=" << gid+1 << endl;
				write_surface_output(gid,timeStep);
			} // end if
			if (timeStep%restart_freq[gid]==0 || lastTimeStep) {
				if (Rank==0) cout << "[I] Writing restart for grid=" << gid+1 << endl;
				write_restart(gid,timeStep,time[gid]);
			} // end if
			if (timeStep%loads[gid].frequency==0) {
				write_loads(gid,timeStep,time[gid]);
			} // end if
			if (timeStep%time_step_update_freq==0) {
				update_time_step_options();
				time_step_update_freq=input.section("timemarching").get_int("updatefrequency");
				timeStepMax=input.section("timemarching").get_int("numberofsteps");
			} // end if
			if (timeStep%pseudo_time_step_update_freq==0) {
				update_pseudo_time_step_options();
				ps_tolerance=input.section("pseudotime").get_double("residualtolerance");
				ps_step_max=input.section("pseudotime").get_int("numberofsteps");
				pseudo_time_step_update_freq=input.section("pseudotime").get_int("updatefrequency");
				if (equations[gid]==NS) {
					if (input.section("pseudotime").get_string("preconditioner")=="none") {
						ns[gid].preconditioner=NONE;
					} else if (input.section("pseudotime").get_string("preconditioner")=="ws95") {
						ns[gid].preconditioner=WS95;
					}
				}
			} // end if	
			
		} // end grid loop
		
		if (fexists("dump_restart")) {
			if (Rank==0) cout << "[I] Writing restart for all grids" << endl;
			for (int gid=0;gid<grid.size();++gid) write_restart(gid,timeStep,time[gid]);
			MPI_Barrier(MPI_COMM_WORLD);
			if (Rank==0) remove("dump_restart");
		} 
		if (fexists("dump_volume")) {
			if (Rank==0) cout << "[I] Writing volume output for all grids" << endl;
			for (int gid=0;gid<grid.size();++gid) write_volume_output(gid,timeStep);
			MPI_Barrier(MPI_COMM_WORLD);
			if (Rank==0) remove("dump_volume");
		} 
		if (fexists("dump_surface")) {
			if (Rank==0) cout << "[I] Writing surface output for all grids" << endl;
			for (int gid=0;gid<grid.size();++gid) write_surface_output(gid,timeStep);
			MPI_Barrier(MPI_COMM_WORLD);
			if (Rank==0) remove("dump_surface");
		} 
		if (fexists("dump_all")) {
			if (Rank==0) {
				cout << "[I] Writing restart for all grids" << endl;
				cout << "[I] Writing volume output for all grids" << endl;
				cout << "[I] Writing surface output for all grids" << endl;
			}
			for (int gid=0;gid<grid.size();++gid) {
				write_restart(gid,timeStep,time[gid]);
				write_volume_output(gid,timeStep);
				write_surface_output(gid,timeStep);
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (Rank==0) remove("dump_all");
		} 	

	}
	/*****************************************************************************************/
	// End time loop
	/*****************************************************************************************/	
	convergence.close();	
	MPI_Barrier(MPI_COMM_WORLD);

      	if (Rank==0) {
      	      timeEnd=MPI_Wtime();
      	      cout << "[I] Wall time = " << timeEnd-timeRef << " second" << endl;
      	}

	PetscFinalize();
	MPI_Finalize();
	
	return 0;
}                	
                 	
void set_lengthScales(int gid) {
	// Loop through the cells and calculate length scales for each cell
	for (int c=0;c<grid[gid].cellCount;++c) {
		grid[gid].cell[c].lengthScale=1.e20;
		double height;
		int f;
		for (int cf=0;cf<grid[gid].cell[c].faceCount;++cf) {
			f=grid[gid].cell[c].faces[cf];
			height=fabs(grid[gid].face[f].normal.dot(grid[gid].face[f].centroid-grid[gid].cell[c].centroid));
			bool skipScale=false;
			if (grid[gid].face[f].bc>=0) {
				if (bc[gid][grid[gid].face[f].bc].type==SYMMETRY) {
					skipScale=true;
				}
			}
			if (!skipScale) grid[gid].cell[c].lengthScale=min(grid[gid].cell[c].lengthScale,height);
		}
	}
	
	// Find out the global, grid cell length scale
	if (grid[gid].dimension==3) {
		grid[gid].lengthScale=pow(grid[gid].globalTotalVolume,1./3.);
		if (Rank==0) cout << "[I] Grid length scale is set to " << grid[gid].lengthScale << endl;
	} else {
		
		// If the problem is 2D, finding the grid length scale is a bit more challenging
		// Find the symmetry BC region with the largest area, and take the square root
	        grid[gid].lengthScale=0.;	
		for (int b=0;b<bc[gid].size();++b) {
			if (bc[gid][b].type==SYMMETRY) grid[gid].lengthScale=max(grid[gid].lengthScale,sqrt(0.5*bc[gid][b].total_area));
		}
		if (Rank==0) cout << "[I] Grid length scale is set to " << grid[gid].lengthScale << endl;
	}
	
	return;
} 
