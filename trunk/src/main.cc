#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include <iomanip>
using namespace std;

#include "grid.h"
#include "inputs.h"
#include "sparse.h"
#include "bc.h"

// Function prototypes
void read_inputs(InputFile &input);
void initialize(Grid &grid, InputFile &input);
void write_tec(string fileName, double time);
void read_tec(string fileName, double &time);
void writeTecplotMacro(int restart, int timeStepMax, int outFreq);
void set_bcs(Grid& grid, InputFile &input, BC &bc);
bool within_box(Vec3D centroid, Vec3D box_1, Vec3D box_2);
double minmod(double a, double b);
double maxmod(double a, double b);
double superbee(double a, double b);
void update(double dt, double gamma);
string get_filename(string begin, int number, string ext) ;

void fou(double gamma);
void hancock_predictor(double gamma, double dt, string limiter);
void hancock_corrector(double gamma, string limiter);
void diff_flux(double mu);

// Initiate grid class
Grid grid;
// Allocate container for boundary conditions
BC bc;

int np, rank;

int main (int argc, char *argv[]) {

	// Initialize mpi
	MPI_Init (&argc,&argv);
	// Find the number of processors
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	// Find current processor rank
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); 

	string inputFileName, gridFileName;
	inputFileName.assign(argv[1]);
	gridFileName.assign(argv[1]);
	inputFileName+=".in";
	gridFileName+=".cgns";
	
	int restart=0;
	if (argc>2) restart=atoi(argv[2]);
	
	InputFile input(inputFileName);
	read_inputs(input);
	grid.read(gridFileName);
/*
	initialize(grid,input);

	cout << "* Applied initial conditions" << endl;
	
	set_bcs(grid,input,bc);

	cout << "* Set boundary conditions" << endl;	
	
	double gamma = input.section["fluidProperties"].doubles["gamma"];
	double time = 0.;

	double dt = input.section["timeMarching"].doubles["step"];
	double CFL= input.section["timeMarching"].doubles["CFL"];
	int timeStepMax = input.section["timeMarching"].ints["numberOfSteps"];
	int outFreq = input.section["timeMarching"].ints["outFreq"];
	double mu;
	if (input.section["fluidProperties"].subsections["viscosity"].strings["type"]=="fixed") {
		mu=input.section["fluidProperties"].subsections["viscosity"].doubles["value"];
	}
	string limiter=input.section["numericalOptions"].strings["limiter"];

	// Calculate and store face average contributions from each cell
	unsigned int n,c;
	double weight,weightSum;
	Vec3D node2cell;
	VecSparse nodeAverage;
	
	for (unsigned int f=0;f<grid.faceCount;++f) {
		grid.face[f].cellContributions.flush();
		for (unsigned int fn=0;fn<grid.face[f].nodeCount;++fn) {
			n=grid.face[f].nodes[fn];
			weightSum=0.;
			nodeAverage.flush();
			for (unsigned int nc=0;nc<grid.node[n].cells.size();++nc) {
				c=grid.node[n].cells[nc];
				node2cell=grid.node[n]-grid.cell[c].centroid;
				weight=1./(node2cell.dot(node2cell));
				weightSum+=weight;
				nodeAverage.put(c,weight);
			}
			grid.face[f].cellContributions+=(nodeAverage/weightSum)/grid.face[f].nodeCount;
		}
	}

	cout << "* Calculated cell gradient contributions" << endl;
	
	if (restart!=0) {
		string fileName=get_filename("out", restart, "dat");
		fstream file;
		file.open(fileName.c_str());
		if (file.is_open()) {
			cout << "* Restarting from " << fileName  << endl;
			file.close();
			read_tec(fileName,time);
		} else {
			cerr << "[!!] Restart "<< fileName << " could not be found." << endl;
			exit(0);
		}
	}


	cout << "* Beginning time loop" << endl;
	
	// Begin time loop
	for (int timeStep = restart;timeStep < input.section["timeMarching"].ints["numberOfSteps"]+restart;++timeStep) {
		
		if (input.section["timeMarching"].strings["type"]=="CFL") {
			// Determine time step with CFL condition
			double cellScaleX, cellScaleY, cellScaleZ;
			dt=1.E20;
			for (unsigned int c=0;c<grid.cellCount;++c) {
				double a=sqrt(gamma*grid.cell[c].p/grid.cell[c].rho);
				cellScaleX=cellScaleY=cellScaleZ=0.;
				for (unsigned int cn=0;cn<grid.cell[c].nodeCount;++cn) {
					for (unsigned int cn2=0;cn2<grid.cell[c].nodeCount;++cn2) {
						cellScaleX=max(cellScaleX,fabs(grid.cell[c].node(cn).comp[0]-grid.cell[c].node(cn2).comp[0]));
						cellScaleY=max(cellScaleY,fabs(grid.cell[c].node(cn).comp[1]-grid.cell[c].node(cn2).comp[1]));
						cellScaleZ=max(cellScaleZ,fabs(grid.cell[c].node(cn).comp[2]-grid.cell[c].node(cn2).comp[2]));					
					}
				}
				dt=min(dt,CFL*cellScaleX/(fabs(grid.cell[c].v.comp[0])+a));
				dt=min(dt,CFL*cellScaleY/(fabs(grid.cell[c].v.comp[1])+a));
				dt=min(dt,CFL*cellScaleZ/(fabs(grid.cell[c].v.comp[2])+a));			
			}
		}
				
		if (input.section["numericalOptions"].strings["order"]=="first") {	
			fou(gamma);
			if (input.section["equations"].strings["set"]=="NS") grid.gradients();
		} else {
			
			// Calculate all the cell gradients for each variable
			grid.gradients();
			
			double rhoOld[grid.cellCount] ,pOld[grid.cellCount] ; 
			Vec3D vOld[grid.cellCount];		
			// Backup variables
			for (unsigned int c = 0;c < grid.cellCount;++c) {
				rhoOld[c]=grid.cell[c].rho;
				vOld[c]=grid.cell[c].v;
				pOld[c]=grid.cell[c].p;
			}
			//hancock_predictor(gamma,dt,limiter);
			hancock_corrector(gamma,limiter);
			// Restore variables
			for (unsigned int c = 0;c < grid.cellCount;++c) {
				grid.cell[c].rho=rhoOld[c];
				grid.cell[c].v=vOld[c];
				grid.cell[c].p=pOld[c];
			}
		}
		
		if (input.section["equations"].strings["set"]=="NS") diff_flux(mu);
				
		update(dt,gamma);
		
		cout << timeStep << "\t" << setprecision(4) << scientific << time << "\t" << dt << endl;
		time += dt;
		if ((timeStep + 1) % outFreq == 0) {
			string fileName;
			fileName=get_filename("out",timeStep+1, "dat") ;
			// Write tecplot output file
			write_tec(fileName,time);
		}
	}
	writeTecplotMacro(restart,timeStepMax, outFreq);
*/
	MPI_Finalize();
	return 0;
}

double minmod(double a, double b) {
	if (a*b < 0 ) {
		return 0.;
	} else if ( fabs(a)<fabs(b)) {
		return a;
	} else {
		return b;
	}
}

double maxmod(double a, double b) {
	if (a*b < 0 ) {
		return 0.;
	} else if ( fabs(a)<fabs(b)) {
		return b;
	} else {
		return a;
	}
}

double superbee(double a, double b) {
	return minmod(maxmod(a,b),minmod(a,b));
}

bool within_box(Vec3D centroid, Vec3D box_1, Vec3D box_2) {
	double tocorners_x=fabs(centroid.comp[0]-box_1.comp[0])+fabs(centroid.comp[0]-box_2.comp[0]);
	double tocorners_y=fabs(centroid.comp[1]-box_1.comp[1])+fabs(centroid.comp[1]-box_2.comp[1]);
	double tocorners_z=fabs(centroid.comp[2]-box_1.comp[2])+fabs(centroid.comp[2]-box_2.comp[2]);				
	if (tocorners_x<=fabs(box_1.comp[0]-box_2.comp[0]) ) {
		if (tocorners_y<=fabs(box_1.comp[1]-box_2.comp[1]) ) {
			if (tocorners_z<=fabs(box_1.comp[2]-box_2.comp[2]) ) {
				// The cell centroid is inside the box region
				return true;						
			}
		}
	}
	return false;
}

void update(double dt, double gamma) {
	
	double conservative[5];	
	for (unsigned int c = 0;c < grid.cellCount;++c) {
		conservative[0] = grid.cell[c].rho;
		conservative[1] = grid.cell[c].rho * grid.cell[c].v.comp[0];
		conservative[2] = grid.cell[c].rho * grid.cell[c].v.comp[1];
		conservative[3] = grid.cell[c].rho * grid.cell[c].v.comp[2];
		conservative[4] = 0.5 * grid.cell[c].rho * grid.cell[c].v.dot(grid.cell[c].v) + grid.cell[c].p / (gamma - 1.);
		for (int i = 0;i < 5;++i) {
			conservative[i] -= dt / grid.cell[c].volume * grid.cell[c].flux[i];
			grid.cell[c].flux[i] = 0.;
		}
		grid.cell[c].rho = conservative[0];
		grid.cell[c].v.comp[0] = conservative[1] / conservative[0];
		grid.cell[c].v.comp[1] = conservative[2] / conservative[0];
		grid.cell[c].v.comp[2] = conservative[3] / conservative[0];
		grid.cell[c].p = (conservative[4] - 0.5 * conservative[0] * grid.cell[c].v.dot(grid.cell[c].v)) * (gamma - 1.);
	} // cell loop
	return;
}

string get_filename(string begin, int number, string ext) {
	char dummy[12];				
	// Print timeStep integer to character
	sprintf(dummy, "%12d", number);
	// Convert character to string and erase leading whitespaces
	string fileName = dummy;
	fileName.erase(0, fileName.rfind(" ", fileName.length()) + 1);
	fileName = begin + fileName + "."  + ext;
	return fileName;
}