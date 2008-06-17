#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include <iomanip>
#include <sys/stat.h>
#include <sys/types.h>
using namespace std;

#include "grid.h"
#include "inputs.h"
#include "sparse.h"
#include "bc.h"

#include "petscksp.h"


// Function prototypes
void read_inputs(InputFile &input);
void initialize(Grid &grid, InputFile &input);
void write_tec(int timeStep, double time);
void write_vtk(int timeStep);
void write_vtk_parallel(int timeStep);
void read_tec(int restart, int global2local[], double &time);
void writeTecplotMacro(int restart, int timeStepMax, int outFreq);
void set_bcs(Grid& grid, InputFile &input, BC &bc);
bool within_box(Vec3D centroid, Vec3D box_1, Vec3D box_2);
double minmod(double a, double b);
double maxmod(double a, double b);
double superbee(double a, double b);
void update(double dt);
void updateImp(double dt);
string int2str(int number) ;
void fou(void);
void hancock_predictor(double dt, string limiter);
void hancock_corrector(string limiter);
void diff_flux(double mu);
void assembleConvJac(Mat impOP,unsigned int f, unsigned int p);
void jac(Mat impOP);

// PETSC 
static char help[] = "Free CFD\n";
//extern PetscErrorCode petscRHS(SNES,Vec,Vec,void*);
//extern PetscErrorCode petscMAT(SNES,Vec,Mat*,Mat*,MatStructure*,void*);

// Initiate grid class
Grid grid;
// Allocate container for boundary conditions
BC bc;

int np, rank;
double Gamma,dt;

int main(int argc, char *argv[]) {

	// Initialize mpi
	MPI_Init(&argc,&argv);
	// Find the number of processors
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	// Find current processor rank
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Initialize PETSC
	KSP ksp; // linear solver context
	PC pc; // preconditioner context
	Vec deltaU,rhs; // solution, residual vectors
	Mat impOP; // implicit operator matrix
	PetscErrorCode ierr;

	PetscInitialize(&argc,&argv,(char *)0,help);

	//Create nonlinear solver context
	KSPCreate(PETSC_COMM_WORLD,&ksp);

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
	
	vector<unsigned int> sendCells[np];
	unsigned int recvCount[np];
	{
		int maxGhost=grid.globalCellCount/np*2;		
		unsigned int ghosts2receive[np][maxGhost],ghosts2send[np][maxGhost];
		for (unsigned int p=0;p<np;++p) {
			ghosts2receive[p][0]=0;
		}
		for (unsigned int g=0; g<grid.ghostCount; ++g) {
			unsigned int p=grid.ghost[g].partition;
			ghosts2receive[p][ghosts2receive[p][0]+1]=grid.ghost[g].globalId;
			++ghosts2receive[p][0];
		}
		//if (rank==0) for (int p=0;p<np;++p) for (int g=0;g<ghosts2receive[p][0];++g) cout << p << "\t" << ghosts2receive[p][0] << "\t" << ghosts2receive[p][g+1] << endl;
		for (unsigned int p=0;p<np;++p) {
			MPI_Alltoall(ghosts2receive,maxGhost,MPI_UNSIGNED,ghosts2send,maxGhost,MPI_UNSIGNED,MPI_COMM_WORLD);
		}
		//if (rank==0) for (int p=0;p<np;++p) for (int g=0;g<ghosts2send[p][0];++g) cout << p << "\t" << ghosts2send[p][0] << "\t" << ghosts2send[p][g+1] << endl;
		// Transfer data to more efficient containers
		for (unsigned int p=0;p<np;++p) {
			for (unsigned int i=1;i<=ghosts2send[p][0];++i) sendCells[p].push_back(ghosts2send[p][i]);
			recvCount[p]=ghosts2receive[p][0];
		}
	} // end scope
	
	initialize(grid,input);
	if (rank==0) cout << "[I] Applied initial conditions" << endl;

	set_bcs(grid,input,bc);
	if (rank==0) cout << "[I] Set boundary conditions" << endl;
	
	grid.nodeAverages();
	grid.faceAverages();
	grid.gradMaps();
	
	double time = 0.;
	Gamma = input.section["fluidProperties"].doubles["gamma"];
	dt = input.section["timeMarching"].doubles["step"];
	double CFL= input.section["timeMarching"].doubles["CFL"];
	int timeStepMax = input.section["timeMarching"].ints["numberOfSteps"];
	int outFreq = input.section["timeMarching"].ints["outFreq"];
	double mu;
	if (input.section["fluidProperties"].subsections["viscosity"].strings["type"]=="fixed") {
		mu=input.section["fluidProperties"].subsections["viscosity"].doubles["value"];
	}
	string limiter=input.section["numericalOptions"].strings["limiter"];
	double sharpeningFactor=input.section["numericalOptions"].doubles["sharpeningFactor"];
	
	// Need a conversion map from globalId to local index
	int global2local[grid.globalCellCount];
	for (unsigned int c=0;c<grid.globalCellCount;++c) global2local[c]=-1;
	for (unsigned int c=0;c<grid.cellCount;++c) {
		global2local[grid.cell[c].globalId]=c;
	}
	
	if (restart!=0) {
		for (int p=0;p<np;++p) {
			if (rank==p) read_tec(restart,global2local,time);
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
	
	struct mpiGhost {
		unsigned int partition,globalId;
		double vars[5];
	};
	
	int array_of_block_lengths[2]={2,5};
	MPI_Aint extent;
	MPI_Type_extent(MPI_UNSIGNED,&extent);
	MPI_Aint array_of_displacements[2]={0,2*extent};
	MPI_Datatype array_of_types[2]={MPI_UNSIGNED,MPI_DOUBLE};
	MPI_Datatype MPI_GHOST;
	MPI_Type_struct(2,array_of_block_lengths,array_of_displacements,array_of_types,&MPI_GHOST);
	MPI_Type_commit(&MPI_GHOST);
	
	struct mpiGrad {
		unsigned int globalId;
		double grads[15];
	};
	array_of_block_lengths[0]=1;array_of_block_lengths[1]=15;
	array_of_displacements[1]=extent;
	MPI_Datatype MPI_GRAD;
	MPI_Type_struct(2,array_of_block_lengths,array_of_displacements,array_of_types,&MPI_GRAD);
	MPI_Type_commit(&MPI_GRAD);
	
	for (unsigned int g=0;g<grid.ghostCount;++g) {
		global2local[grid.ghost[g].globalId]=g;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	double timeRef, timeEnd;
	timeRef=MPI_Wtime();

	PetscScalar *dU,*ff,value;
	VecCreateMPI(PETSC_COMM_WORLD,grid.cellCount*5,grid.globalCellCount*5,&rhs);
	VecSetFromOptions(rhs);
	VecDuplicate(rhs,&deltaU);
	VecSet(deltaU,0.);
	//VecView(deltaU,PETSC_VIEWER_STDOUT_WORLD);
	MatCreate(PETSC_COMM_WORLD,&impOP);
	MatSetSizes(impOP,grid.cellCount*5,grid.cellCount*5,grid.globalCellCount*5,grid.globalCellCount*5);
	MatSetFromOptions(impOP);
	
	KSPSetOperators(ksp,impOP,impOP,SAME_NONZERO_PATTERN);
	KSPSetFromOptions(ksp);
	
	if (rank==0) cout << "[I] Beginning time loop" << endl;
	// Begin time loop
	for (int timeStep=restart+1;timeStep<=input.section["timeMarching"].ints["numberOfSteps"]+restart;++timeStep) {
		int nIter;
		// TODO partition variable need not be send and received
		for (unsigned int p=0;p<np;++p) {
			if (rank!=p) {
				mpiGhost sendBuffer[sendCells[p].size()];
				mpiGhost recvBuffer[recvCount[p]];
				int id;
				for (unsigned int g=0;g<sendCells[p].size();++g)	{
					id=global2local[sendCells[p][g]];
					sendBuffer[g].partition=rank;
					sendBuffer[g].globalId=grid.cell[id].globalId;
					sendBuffer[g].vars[0]=grid.cell[id].rho;
					sendBuffer[g].vars[1]=grid.cell[id].v.comp[0];
					sendBuffer[g].vars[2]=grid.cell[id].v.comp[1];
					sendBuffer[g].vars[3]=grid.cell[id].v.comp[2];
					sendBuffer[g].vars[4]=grid.cell[id].p;
				}

				int tag=rank; // tag is set to source
				MPI_Sendrecv(sendBuffer,sendCells[p].size(),MPI_GHOST,p,0,recvBuffer,recvCount[p],MPI_GHOST,p,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		
				for (unsigned int g=0;g<recvCount[p];++g) {
					id=global2local[recvBuffer[g].globalId];
					grid.ghost[id].partition=recvBuffer[g].partition;
					// TODO is it necessary to overwrite gloabalId?
					grid.ghost[id].globalId=recvBuffer[g].globalId;
					grid.ghost[id].rho=recvBuffer[g].vars[0];
					grid.ghost[id].v.comp[0]=recvBuffer[g].vars[1];
					grid.ghost[id].v.comp[1]=recvBuffer[g].vars[2];
					grid.ghost[id].v.comp[2]=recvBuffer[g].vars[3];
					grid.ghost[id].p=recvBuffer[g].vars[4];
				}
				
			}
		}

// 		cerr << "[DEBUG rank=" << rank << " ] updated ghosts" << endl;
		
		if (input.section["timeMarching"].strings["type"]=="CFL") {
			// Determine time step with CFL condition
			//double cellScaleX, cellScaleY, cellScaleZ;
			double lengthScale;
			dt=1.E20;
			for (unsigned int c=0;c<grid.cellCount;++c) {
				double a=sqrt(Gamma*grid.cell[c].p/grid.cell[c].rho);
				lengthScale=grid.cell[c].lengthScale;
				dt=min(dt,CFL*lengthScale/(fabs(grid.cell[c].v.comp[0])+a));
				dt=min(dt,CFL*lengthScale/(fabs(grid.cell[c].v.comp[1])+a));
				dt=min(dt,CFL*lengthScale/(fabs(grid.cell[c].v.comp[2])+a));
			}
			MPI_Allreduce(&dt,&dt,1, MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
		}


		if (input.section["numericalOptions"].strings["order"]=="first") {

			int index;
			fou();

			if (input.section["timeMarching"].strings["integrator"]=="backwardEuler") {
				MatZeroEntries(impOP);

				unsigned int counter=0;
				for (unsigned int c=0;c<grid.cellCount;++c) {
					for (int i=0;i<5;++i) {
						value=-1.*grid.cell[c].flux[i]; 
				//if (grid.cell[c].centroid.comp[0]>0.23 && grid.cell[c].centroid.comp[0]<0.26) cout << c << "\t" << ff << endl;
						index=grid.cell[c].globalId*5+i;
						VecSetValues(rhs,1,&index,&value,INSERT_VALUES);

						value=grid.cell[c].volume/dt;
	// 					if (i==4) {
	// 						double Mach=fabs(grid.cell[c].v)/sqrt(Gamma*grid.cell[c].p/grid.cell[c].rho);
	// 						double MachTom2=1./(Mach*Mach);
	// // 						value+=0.5*grid.cell[c].v.dot(grid.cell[c].v)*(MachTom2-1.)*grid.cell[c].rho;
	// // 						value+=grid.cell[c].v.comp[0]*(1.-MachTom2)*grid.cell[c].rho*grid.cell[c].v.comp[0];
	// // 						value+=grid.cell[c].v.comp[1]*(1.-MachTom2)*grid.cell[c].rho*grid.cell[c].v.comp[1];
	// // 						value+=grid.cell[c].v.comp[2]*(1.-MachTom2)*grid.cell[c].rho*grid.cell[c].v.comp[2];
	// // 						value+=MachTom2*0.5*grid.cell[c].rho*(grid.cell[c].v.dot(grid.cell[c].v))+grid.cell[c].p/(Gamma-1.);
	// 						value+=0.5*grid.cell[c].v.dot(grid.cell[c].v)*(MachTom2-1.)*grid.cell[c].rho;
	// 						value+=grid.cell[c].v.comp[0]*(1.-MachTom2)*grid.cell[c].rho*grid.cell[c].v.comp[0];
	// 						value+=grid.cell[c].v.comp[1]*(1.-MachTom2)*grid.cell[c].rho*grid.cell[c].v.comp[1];
	// 						value+=grid.cell[c].v.comp[2]*(1.-MachTom2)*grid.cell[c].rho*grid.cell[c].v.comp[2];
	// 						value+=MachTom2*0.5*grid.cell[c].rho*(grid.cell[c].v.dot(grid.cell[c].v))+grid.cell[c].p/(Gamma-1.);
	// 
	// 					}
						MatSetValues(impOP,1,&index,1,&index,&value,INSERT_VALUES);
						counter++;
					}
				}

				MatAssemblyBegin(impOP,MAT_FLUSH_ASSEMBLY);
				MatAssemblyEnd(impOP,MAT_FLUSH_ASSEMBLY);

				//cout << "before" << endl;
				jac(impOP);
				//cout << "after" << endl;
	// 			for (unsigned int f=0;f<grid.faceCount;++f) {
	// 				if (grid.face[f].parent!=grid.cellCount-1) assembleConvJac(impOP,f,grid.face[f].parent);
	// 				if (grid.face[f].neighbor!=grid.cellCount-1) assembleConvJac(impOP,f,grid.face[f].neighbor);
	// 			}

				//VecRestoreArray(f,&ff);
				
				MatAssemblyBegin(impOP,MAT_FINAL_ASSEMBLY);
				MatAssemblyEnd(impOP,MAT_FINAL_ASSEMBLY);
				
				//MatView(impOP,PETSC_VIEWER_STDOUT_WORLD);
				//MatView(impOP,PETSC_VIEWER_DRAW_WORLD);
				
				VecAssemblyBegin(rhs);
				VecAssemblyEnd(rhs);

				KSPSolve(ksp,rhs,deltaU);

				KSPGetIterationNumber(ksp,&nIter);
		
				
				//VecView(deltaU,PETSC_VIEWER_STDOUT_WORLD);
				VecGetArray(deltaU,&dU);//VecGetArray(rhs,&ff);
				counter=0;
				for (unsigned int c=0;c<grid.cellCount;++c) {
					for (int i=0;i<5;++i) {
						//if (grid.cell[c].centroid.comp[0]>0.23 && grid.cell[c].centroid.comp[0]<0.26) cout << c << "\t" << dU[counter] << "\t" << ff[counter] << "\t" << grid.cell[c].volume/dt << endl;
						grid.cell[c].flux[i]=dU[counter];
						//cout << counter << "\t" << dU[counter] << endl;
						counter++;
					}
				}
				VecRestoreArray(deltaU,&dU);//VecRestoreArray(rhs,&ff);
				updateImp(dt);

			} // if backwardEuler

			if (input.section["timeMarching"].strings["integrator"]=="forwardEuler") update(dt);
			
			//if (input.section["equations"].strings["set"]=="NS") grid.gradients();
		} else { // if not first order

			// Calculate all the cell gradients for each variable
			grid.gradients();

			// Update ghost gradients
			for (unsigned int p=0;p<np;++p) {
				mpiGrad sendBuffer[sendCells[p].size()];
				mpiGrad recvBuffer[recvCount[p]];
				int id;
				for (unsigned int g=0;g<sendCells[p].size();++g) {
					id=global2local[sendCells[p][g]];
					sendBuffer[g].globalId=grid.cell[id].globalId;
					int count=0;
					for (unsigned int var=0;var<5;++var) {
						for (unsigned int comp=0;comp<3;++comp) {
							sendBuffer[g].grads[count]=grid.cell[id].grad[var].comp[comp];
							count++;
						}
					}
				}

				int tag=rank; // tag is set to source
				MPI_Sendrecv(sendBuffer,sendCells[p].size(),MPI_GRAD,p,0,recvBuffer,recvCount[p],MPI_GRAD,p,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		
				for (unsigned int g=0;g<recvCount[p];++g) {
					id=global2local[recvBuffer[g].globalId];
					int count=0;
					for (unsigned int var=0;var<5;++var) {
						for (unsigned int comp=0;comp<3;++comp) {
							grid.ghost[id].grad[var].comp[comp]=recvBuffer[g].grads[count];
							count++;
						}
					}
				}
			}

			// Limit gradients
			grid.limit_gradients(limiter,sharpeningFactor);
			
			// Send limited gradients
			for (unsigned int p=0;p<np;++p) {
				mpiGrad sendBuffer[sendCells[p].size()];
				mpiGrad recvBuffer[recvCount[p]];
				int id;
				for (unsigned int g=0;g<sendCells[p].size();++g) {
					id=global2local[sendCells[p][g]];
					sendBuffer[g].globalId=grid.cell[id].globalId;
					int count=0;
					for (unsigned int var=0;var<5;++var) {
						for (unsigned int comp=0;comp<3;++comp) {
							sendBuffer[g].grads[count]=grid.cell[id].limited_grad[var].comp[comp];
							count++;
						}
					}
				}

				int tag=rank; // tag is set to source
				MPI_Sendrecv(sendBuffer,sendCells[p].size(),MPI_GRAD,p,0,recvBuffer,recvCount[p],MPI_GRAD,p,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		
				for (unsigned int g=0;g<recvCount[p];++g) {
					id=global2local[recvBuffer[g].globalId];
					int count=0;
					for (unsigned int var=0;var<5;++var) {
						for (unsigned int comp=0;comp<3;++comp) {
							grid.ghost[id].limited_grad[var].comp[comp]=recvBuffer[g].grads[count];
							count++;
						}
					}
				}
			}

			hancock_corrector(limiter);
			//fou();
			update(dt);
		}

		//if (input.section["equations"].strings["set"]=="NS") diff_flux(mu);
			


		if (rank==0) cout << timeStep << "\t" << setprecision(4) << scientific << time << "\t" << dt << "\t" << nIter << endl;
		time += dt;
		if ((timeStep) % outFreq == 0) {
			mkdir("./output",S_IRWXU);
			if (input.section["solutionFormat"].strings["output"]=="vtk") {
				// Write vtk output file
				if (rank==0) write_vtk_parallel(timeStep);
				write_vtk(timeStep);
			} else if(input.section["solutionFormat"].strings["output"]=="tecplot") {
				// Write tecplot output file
				for (int p=0;p<np;++p) {
					if(rank==p) write_tec(timeStep,time);
					MPI_Barrier(MPI_COMM_WORLD);
				}
			}
			ofstream file;
			string fileName;
			for (int p=0;p<np;++p) {
				if (rank==p) {
					fileName="./output/partitionMap"+int2str(timeStep)+".dat";
					if (rank==0) { file.open(fileName.c_str(), ios::out); file << np << endl;}
					else { file.open(fileName.c_str(), ios::app); }
					file << grid.nodeCount << "\t" << grid.cellCount << endl;
					for (unsigned int c=0;c<grid.cellCount;++c) file << grid.cell[c].globalId << endl;
					file.close();
				}
				MPI_Barrier(MPI_COMM_WORLD);
			}
		}
	}
	//writeTecplotMacro(restart,timeStepMax, outFreq);
	// Syncronize the processors
	MPI_Barrier(MPI_COMM_WORLD);

	// Report the wall time
	if (rank==0) {
		timeEnd=MPI_Wtime();
		cout << "* Wall time: " << timeEnd-timeRef << " seconds" << endl;
	}
	//MPI_Finalize();

	KSPDestroy(ksp);
	MatDestroy(impOP);
	VecDestroy(rhs);
	VecDestroy(deltaU);
	PetscFinalize();
		
	return 0;
}

double minmod(double a, double b) {
	if ((a*b)<0) {
		return 0.;
	} else if (fabs(a) < fabs(b)) {
		return a;
	} else {
		return b;
	}
}

double maxmod(double a, double b) {
	if ((a*b)<0) {
		return 0.;
	} else if (fabs(a) < fabs(b)) {
		return b;
	} else {
		return a;
	}
}

double superbee(double a, double b) {
	return minmod(maxmod(a,b),minmod(a,b));
}

bool within_box(Vec3D centroid, Vec3D box_1, Vec3D box_2) {
	double tocorners_x=fabs(centroid.comp[0]-box_1.comp[0]) +fabs(centroid.comp[0]-box_2.comp[0]);
	double tocorners_y=fabs(centroid.comp[1]-box_1.comp[1]) +fabs(centroid.comp[1]-box_2.comp[1]);
	double tocorners_z=fabs(centroid.comp[2]-box_1.comp[2]) +fabs(centroid.comp[2]-box_2.comp[2]);
	if (tocorners_x<=fabs(box_1.comp[0]-box_2.comp[0])) {
		if (tocorners_y<=fabs(box_1.comp[1]-box_2.comp[1])) {
			if (tocorners_z<=fabs(box_1.comp[2]-box_2.comp[2])) {
				// The cell centroid is inside the box region
				return true;
			}
		}
	}
	return false;
}

void update(double dt) {

	double conservative[5];
	for (unsigned int c = 0;c < grid.cellCount;++c) {
		conservative[0] = grid.cell[c].rho;
		conservative[1] = grid.cell[c].rho * grid.cell[c].v.comp[0];
		conservative[2] = grid.cell[c].rho * grid.cell[c].v.comp[1];
		conservative[3] = grid.cell[c].rho * grid.cell[c].v.comp[2];
		conservative[4] = 0.5 * grid.cell[c].rho * grid.cell[c].v.dot(grid.cell[c].v) + grid.cell[c].p / (Gamma - 1.);
		for (int i = 0;i < 5;++i) {
			conservative[i] -= dt / grid.cell[c].volume * grid.cell[c].flux[i];
			grid.cell[c].flux[i] = 0.;
		}
		grid.cell[c].rho = conservative[0];
		grid.cell[c].v.comp[0] = conservative[1] / conservative[0];
		grid.cell[c].v.comp[1] = conservative[2] / conservative[0];
		grid.cell[c].v.comp[2] = conservative[3] / conservative[0];
		grid.cell[c].p = (conservative[4] - 0.5 * conservative[0] * grid.cell[c].v.dot(grid.cell[c].v)) * (Gamma - 1.);
	} // cell loop
	return;
}


void updateImp(double dt) {

	double conservative[5];
	for (unsigned int c = 0;c < grid.cellCount;++c) {
		conservative[0] = grid.cell[c].rho;
		conservative[1] = grid.cell[c].rho * grid.cell[c].v.comp[0];
		conservative[2] = grid.cell[c].rho * grid.cell[c].v.comp[1];
		conservative[3] = grid.cell[c].rho * grid.cell[c].v.comp[2];
		conservative[4] = 0.5 * grid.cell[c].rho * grid.cell[c].v.dot(grid.cell[c].v) + grid.cell[c].p / (Gamma - 1.);
		for (int i = 0;i < 5;++i) {
			conservative[i] += grid.cell[c].flux[i];
			grid.cell[c].flux[i] = 0.;
		}
		grid.cell[c].rho = conservative[0];
		grid.cell[c].v.comp[0] = conservative[1] / conservative[0];
		grid.cell[c].v.comp[1] = conservative[2] / conservative[0];
		grid.cell[c].v.comp[2] = conservative[3] / conservative[0];
		grid.cell[c].p = (conservative[4] - 0.5 * conservative[0] * grid.cell[c].v.dot(grid.cell[c].v)) * (Gamma - 1.);
	} // cell loop
	return;
}

string int2str(int number) {
	char dummy[12];
	// Print integer to character
	sprintf(dummy, "%12d", number);
	// Convert character to string and erase leading whitespaces
	string name = dummy;
	name.erase(0, name.rfind(" ", name.length()) + 1);
	return name;
}


void assembleConvJac(Mat impOP, unsigned int f, unsigned int c) {
			
	double a2=Gamma-1.;
	double a3=Gamma-2.;
	Vec3D normal;
	normal=grid.face[f].normal;
	if ((grid.face[f].centroid-grid.cell[c].centroid).dot(normal)<0) normal*=-1;
	
	double faceVel=grid.cell[c].v.dot(normal);
	double phi=0.5*(Gamma-1.)*grid.cell[c].v.dot(grid.cell[c].v);
	double a1=Gamma/(Gamma-1.)*(grid.cell[c].p/grid.cell[c].rho+phi)-phi;
	int diag=grid.cell[c].globalId*5;
	int row=diag;
	int col;
	double value;

	double Ac[5][5];
	
	Ac[0][0]=0.;
	Ac[0][1]=normal.comp[0];
	Ac[0][2]=normal.comp[1];
	Ac[0][3]=normal.comp[2];
	Ac[0][4]=0.;
	
	Ac[1][0]=normal.comp[0]*phi-grid.cell[c].v.comp[0]*faceVel;
	Ac[1][1]=faceVel-a3*normal.comp[0]*grid.cell[c].v.comp[0];
	Ac[1][2]=normal.comp[1]*grid.cell[c].v.comp[0]-a2*normal.comp[0]*grid.cell[c].v.comp[1];
	Ac[1][3]=normal.comp[2]*grid.cell[c].v.comp[0]-a2*normal.comp[0]*grid.cell[c].v.comp[2];
	Ac[1][4]=a2*normal.comp[0];

	Ac[2][0]=normal.comp[1]*phi-grid.cell[c].v.comp[1]*faceVel;
	Ac[2][1]=normal.comp[0]*grid.cell[c].v.comp[1]-a2*normal.comp[1]*grid.cell[c].v.comp[0];
	Ac[2][2]=faceVel-a3*normal.comp[1]*grid.cell[c].v.comp[1];
	Ac[2][3]=normal.comp[2]*grid.cell[c].v.comp[1]-a2*normal.comp[1]*grid.cell[c].v.comp[2];
	Ac[2][4]=a2*normal.comp[1];

	Ac[3][0]=normal.comp[2]*phi-grid.cell[c].v.comp[2]*faceVel;
	Ac[3][1]=normal.comp[0]*grid.cell[c].v.comp[2]-a2*normal.comp[2]*grid.cell[c].v.comp[0];
	Ac[3][2]=normal.comp[1]*grid.cell[c].v.comp[2]-a2*normal.comp[2]*grid.cell[c].v.comp[1];
	Ac[3][3]=faceVel-a3*normal.comp[2]*grid.cell[c].v.comp[2];
	Ac[3][4]=a2*normal.comp[2];

	Ac[4][0]=faceVel*(phi-a1);
	Ac[4][1]=normal.comp[0]*a1-a2*grid.cell[c].v.comp[0]*faceVel;
	Ac[4][2]=normal.comp[1]*a1-a2*grid.cell[c].v.comp[1]*faceVel;
	Ac[4][3]=normal.comp[2]*a1-a2*grid.cell[c].v.comp[2]*faceVel;
	Ac[4][4]=Gamma*faceVel;
	
	for (int i=0;i<5;++i) {
		row=diag+i;
		for (int j=0;j<5;++j) {
			col=row+j;
			value=0.5*grid.face[f].area*grid.cell[c].volume*Ac[i][j];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
		}
	}

	return;
}
