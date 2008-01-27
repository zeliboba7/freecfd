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

int main(int argc, char *argv[]) {

	// Initialize mpi
	MPI_Init(&argc,&argv);
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
		
	// We want to move this to grid.cc cause we already do this search
	//TODO There may be no need for communication
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
	
	initialize(grid,input);
	//cout << "* Applied initial conditions" << endl;

	set_bcs(grid,input,bc);
	//cout << "* Set boundary conditions" << endl;

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

	/*
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
	//cout << "* Calculated cell gradient contributions" << endl;


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
	*/

	//cout << "* Beginning time loop" << endl;

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
	
	// Need a conversion map from globalId to local index
	unsigned int global2local[grid.globalCellCount];
	for (unsigned int c=0;c<grid.cellCount;++c) {
		global2local[grid.cell[c].globalId]=c;
	}
	for (unsigned int g=0;g<grid.ghostCount;++g) {
		global2local[grid.ghost[g].globalId]=g;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	double timeRef, timeEnd;
	timeRef=MPI_Wtime();
	
	// Begin time loop
	for (int timeStep = restart;timeStep < input.section["timeMarching"].ints["numberOfSteps"]+restart;++timeStep) {

		for (unsigned int p=0;p<np;++p) {
			if (rank!=p) {
				mpiGhost sendBuffer[ghosts2send[p][0]];
				mpiGhost recvBuffer[ghosts2receive[p][0]];
				int id;
				for (unsigned int g=0;g<ghosts2send[p][0];++g)	{
					id=global2local[ghosts2send[p][g+1]];
					sendBuffer[g].partition=rank;
					sendBuffer[g].globalId=grid.cell[id].globalId;
					sendBuffer[g].vars[0]=grid.cell[id].rho;
					sendBuffer[g].vars[1]=grid.cell[id].v.comp[0];
					sendBuffer[g].vars[2]=grid.cell[id].v.comp[1];
					sendBuffer[g].vars[3]=grid.cell[id].v.comp[2];
					sendBuffer[g].vars[4]=grid.cell[id].p;
				}
				//if (rank==1) cout << sendBuffer[2].globalId << "\t" << sendBuffer[2].vars[0] << "\t" << sendBuffer[2].vars[4] << endl;
				int tag=rank; // tag is set to source
				//sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status
				MPI_Sendrecv(sendBuffer,ghosts2send[p][0],MPI_GHOST,p,0,recvBuffer,ghosts2receive[p][0],MPI_GHOST,p,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		
				for (unsigned int g=0;g<ghosts2receive[p][0];++g) {
					id=global2local[recvBuffer[g].globalId];
					grid.ghost[id].partition=recvBuffer[g].partition;
					grid.ghost[id].globalId=recvBuffer[g].globalId;
					grid.ghost[id].rho=recvBuffer[g].vars[0];
					grid.ghost[id].v.comp[0]=recvBuffer[g].vars[1];
					grid.ghost[id].v.comp[1]=recvBuffer[g].vars[2];
					grid.ghost[id].v.comp[2]=recvBuffer[g].vars[3];
					grid.ghost[id].p=recvBuffer[g].vars[4];
				}
				
				//MPI_Send(sendBuffer,ghosts2send[p][0],MPI_GHOST,p,tag,MPI_COMM_WORLD);
			}
		}

		/*
		for (unsigned int p=0;p<np;++p) {
			if (rank!=p) {
				mpiGhost recvBuffer[ghosts2receive[p][0]];
				int tag=p; // again tag is set to source
				MPI_Recv(recvBuffer,ghosts2receive[p][0],MPI_GHOST,p,tag,MPI_COMM_WORLD,&status);
				int id;
				for (unsigned int g=0;g<ghosts2receive[p][0];++g) {
					id=global2local[recvBuffer[g].globalId];
					grid.ghost[id].partition=recvBuffer[g].partition;
					grid.ghost[id].globalId=recvBuffer[g].globalId;
					grid.ghost[id].rho=recvBuffer[g].vars[0];
					grid.ghost[id].v.comp[0]=recvBuffer[g].vars[1];
					grid.ghost[id].v.comp[1]=recvBuffer[g].vars[2];
					grid.ghost[id].v.comp[2]=recvBuffer[g].vars[3];
					grid.ghost[id].p=recvBuffer[g].vars[4];
				}
			}
		}
		*/

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
			MPI_Allreduce(&dt,&dt,1, MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
		}
		
		//if (input.section["numericalOptions"].strings["order"]=="first") {
		fou(gamma);
		
			//if (input.section["equations"].strings["set"]=="NS") grid.gradients();
			/*
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
		*/
			
		update(dt,gamma);

		if (rank==0) cout << timeStep << "\t" << setprecision(4) << scientific << time << "\t" << dt << endl;
		time += dt;
		if ((timeStep + 1) % outFreq == 0) {
			string fileName;
			fileName=get_filename("out",timeStep+1,"dat") ;
			// Write tecplot output file
			for (int p=0;p<np;++p) {
				if(rank==p) write_tec(fileName,time);
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
	MPI_Finalize();
	return 0;
}

double minmod(double a, double b) {
	if ((a*b)<0) {
		return 0.;
	} else if (fabs(a) <fabs(b)) {
		return a;
	} else {
		return b;
	}
}

double maxmod(double a, double b) {
	if ((a*b)<0) {
		return 0.;
	} else if (fabs(a) <fabs(b)) {
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
	sprintf(dummy, "%12d", rank);
	fileName = begin + fileName + "."  + ext;
	return fileName;
}
