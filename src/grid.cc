/************************************************************************
	
	Copyright 2007-2008 Emre Sozer & Patrick Clark Trizila

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
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <iomanip>
#include<cmath>
using namespace std;
#include <cgnslib.h>
#include <parmetis.h>


#include "grid.h"
#include "bc.h"

#define EPS 1e-10

				 
extern Grid grid;
extern BC bc;
extern int np, Rank;
extern double Gamma;
extern double Pref;

double superbee(double a, double b);
double minmod(double a, double b);
int gelimd(double **a,double *b,double *x, int n);

Grid::Grid() {
	;
}

int Grid::read(string fname) {
	fstream file;
	fileName=fname;
	file.open(fileName.c_str());
	if (file.is_open()) {
		if (Rank==0) cout << "[I] Found grid file " << fileName  << endl;
		file.close();
		ReadCGNS();
		return 1;
	} else {
		if (Rank==0) cerr << "[E] Grid file "<< fileName << " could not be found." << endl;
		return 0;
	}
}

int Grid::ReadCGNS() {

	int fileIndex,baseIndex,zoneIndex,sectionIndex,nBases,nZones,nSections,nBocos;
	char zoneName[20],sectionName[20]; //baseName[20]
	//  int nBases,cellDim,physDim;
	int size[3];
	globalNodeCount=0;
	globalCellCount=0;
	globalFaceCount=0;
	
	// Open the grid file for reading
	cg_open(fileName.c_str(),MODE_READ,&fileIndex);

	// Read number of bases
	cg_nbases(fileIndex,&nBases);

	// Number of bases is typically 1 
	if (Rank==0) cout << "[I] Number of Bases= " << nBases << endl;

	//for (int baseIndex=1;baseIndex<=nBases;++baseIndex) {
	// TODO For now assuming there is only one base as I don't know in which cases there would be more
	baseIndex=1;
	// Read number of zones (number of blocks in the grid)
	cg_nzones(fileIndex,baseIndex,&nZones);
	int zoneNodeCount[nZones],zoneCellCount[nZones];
	if (Rank==0) cout << "[I] Number of Zones= " << nZones << endl;
	std::vector<int> elemConnIndex,elemConnectivity;
	std::vector<ElementType_t> elemTypes;
	std::vector<double> coordX[nZones],coordY[nZones],coordZ[nZones];
	std::vector<int> zoneCoordMap[nZones];

	// Add to the total number of boundary condition regions in the whole domain
	int totalnBocos=0;
	for (int zoneIndex=1;zoneIndex<=nZones;++zoneIndex) {
		cg_nbocos(fileIndex,baseIndex,zoneIndex,&nBocos);
		totalnBocos+=nBocos;
	}
	std::vector<int> bocos[totalnBocos],bocoConnectivity[totalnBocos];
	int bocoCount=0;
	for (int zoneIndex=1;zoneIndex<=nZones;++zoneIndex) {
		// Read the zone
		cg_zone_read(fileIndex,baseIndex,zoneIndex,zoneName,size);
		// These are the number of cells and nodes in that zone
		zoneNodeCount[zoneIndex-1]=size[0];
		zoneCellCount[zoneIndex-1]=size[1];
		// Read number of sections
		cg_nsections(fileIndex,baseIndex,zoneIndex,&nSections);
		cg_nbocos(fileIndex,baseIndex,zoneIndex,&nBocos);
		if (Rank==0) cout << "[I] In Zone " << zoneName << endl;
		if (Rank==0) cout << "[I] ...Number of Nodes= " << size[0] << endl;
		if (Rank==0) cout << "[I] ...Number of Cells= " << size[1] << endl;
		if (Rank==0) cout << "[I] ...Number of Sections= " << nSections << endl;
		if (Rank==0) cout << "[I] ...Number of Boundary Conditions= " << nBocos << endl;
			
		// Read the node coordinates
		int nodeStart[3],nodeEnd[3];
		nodeStart[0]=nodeStart[1]=nodeStart[2]=1;
		nodeEnd[0]=nodeEnd[1]=nodeEnd[2]=size[0];
				
		double x[size[0]],y[size[0]],z[size[0]];
		cg_coord_read(fileIndex,baseIndex,zoneIndex,"CoordinateX",RealDouble,nodeStart,nodeEnd,&x);
		cg_coord_read(fileIndex,baseIndex,zoneIndex,"CoordinateY",RealDouble,nodeStart,nodeEnd,&y);
		cg_coord_read(fileIndex,baseIndex,zoneIndex,"CoordinateZ",RealDouble,nodeStart,nodeEnd,&z);
		for (int i=0;i<size[0];++i) {
			coordX[zoneIndex-1].push_back(x[i]);
			coordY[zoneIndex-1].push_back(y[i]);
			coordZ[zoneIndex-1].push_back(z[i]);
			// This map takes in zone index and node number within that zone
			// Returns the global index of that node
			// It is initialized to -1 for now
			zoneCoordMap[zoneIndex-1].push_back(-1);
		}

		// In case there are multiple connected zones, collapse the repeated nodes and fix the node numbering
		if (zoneIndex==1) { // If the first zone
			for (int c=0;c<coordX[0].size();++c) {
				// Global node count is incremented as new nodes are found.
				// When in the first zone, every node is new.
				zoneCoordMap[0][c]=globalNodeCount;
				globalNodeCount++;
			}
		}
		// Scan the coordinates of all the other zones before this one for duplicates
		for (int z=0;z<zoneIndex-1;++z) {
			for (int c=0;c<coordX[zoneIndex-1].size();++c) {
				bool foundFlag=false;
				for (int c2=0;c2<coordX[z].size();++c2) {
					if (fabs(coordX[zoneIndex-1][c]-coordX[z][c2])<1.e-7 && fabs(coordY[zoneIndex-1][c]-coordY[z][c2])<1.e-7 && fabs(coordZ[zoneIndex-1][c]-coordZ[z][c2])<1.e-7) {
						zoneCoordMap[zoneIndex-1][c]=zoneCoordMap[z][c2];
						foundFlag=true;
					}
				}
				if (!foundFlag) {
					zoneCoordMap[zoneIndex-1][c]=globalNodeCount;
					globalNodeCount++;
				}
			}
		}

		// Read eement ranges for all the boundary condition regions within the current zone
		int bc_range[nBocos][2];
		for (int bocoIndex=1;bocoIndex<=nBocos;++bocoIndex) {
			int dummy;
			cg_boco_read(fileIndex,baseIndex,zoneIndex,bocoIndex,bc_range[bocoIndex-1],&dummy);
		} // for boco
			
		// Loop sections within the zone
		// These include connectivities of cells and bonudary faces
		for (int sectionIndex=1;sectionIndex<=nSections;++sectionIndex) {
			ElementType_t elemType;
			int elemNodeCount,elemStart,elemEnd,nBndCells,parentFlag;
			// Read the section
			cg_section_read(fileIndex,baseIndex,zoneIndex,sectionIndex,sectionName,&elemType,&elemStart,&elemEnd,&nBndCells,&parentFlag);

			switch (elemType) {
				case TRI_3:
					elemNodeCount=3; break;
				case QUAD_4:
					elemNodeCount=4; break;
				case TETRA_4:
					elemNodeCount=4; break;
				case PENTA_6:
					elemNodeCount=6; break;
				case HEXA_8:
					elemNodeCount=8; break;
			} //switch
			int elemNodes[elemEnd-elemStart+1][elemNodeCount];

			// Read element node connectivities
			cg_elements_read(fileIndex,baseIndex,zoneIndex,sectionIndex,*elemNodes,0);

			// Only pick the volume elements
			if (elemType==TETRA_4 | elemType==PENTA_6 | elemType==HEXA_8 ) {
				if (Rank==0) cout << "[I]    ...Found Volume Section " << sectionName << endl;
				// elements array serves as a start index for connectivity list elemConnectivity
				for (int elem=0;elem<=(elemEnd-elemStart);++elem) {
					elemConnIndex.push_back(elemConnectivity.size());
					elemTypes.push_back(elemType);
					for (int n=0;n<elemNodeCount;++n) elemConnectivity.push_back(zoneCoordMap[zoneIndex-1][elemNodes[elem][n]-1]);
				}
				globalCellCount+=(elemEnd-elemStart+1);
			} else { // If not a volume element
				// Check if a boundary condition section
				bool bcFlag=false;
				for (int nbc=0;nbc<nBocos;++nbc) {
					if (elemStart==bc_range[nbc][0] && elemEnd==bc_range[nbc][1]) {
						bcFlag=true;
						break;
					}
				}
				if (bcFlag) {
					if (Rank==0) cout << "[I]    ...Found BC Section " << sectionName << endl;
					for (int elem=0;elem<=(elemEnd-elemStart);++elem) {
						bocos[bocoCount].push_back(bocoConnectivity[bocoCount].size());
						for (int n=0;n<elemNodeCount;++n) bocoConnectivity[bocoCount].push_back(zoneCoordMap[zoneIndex-1][elemNodes[elem][n]-1]);
					}
					bocoCount+=1;
				}
			}// if
		} // for section
	} // for zone
	//} // for base

	if (Rank==0) cout << "[I] Total Node Count= " << globalNodeCount << endl;
	// Merge coordinates of the zones
	double x[globalNodeCount],y[globalNodeCount],z[globalNodeCount];
	int counter=0;
	// for zone 0
	for (int n=0;n<coordX[0].size();++n) {
		x[counter]=coordX[0][n];
		y[counter]=coordY[0][n];
		z[counter]=coordZ[0][n];
		counter++;
	}
	for (int zone=1;zone<nZones;++zone) {
		for (int n=0;n<coordX[zone].size();++n) {
			if (zoneCoordMap[zone][n]>zoneCoordMap[zone-1][zoneCoordMap[zone-1].size()-1]) {
				x[counter]=coordX[zone][n];
				y[counter]=coordY[zone][n];
				z[counter]=coordZ[zone][n];
				counter++;
			}
		}
	}
	if (counter!=globalNodeCount) cerr << "[E] counter is different from globalNodeCount" << endl;
	if (Rank==0) cout << "[I] Total Cell Count= " << globalCellCount << endl;
	// Store element node counts
	int elemNodeCount[globalCellCount];
	for (unsigned int c=0;c<globalCellCount-1;++c) {
		elemNodeCount[c]=elemConnIndex[c+1]-elemConnIndex[c];
	}
	elemNodeCount[globalCellCount-1]=elemConnectivity.size()-elemConnIndex[globalCellCount-1];
	
	// Initialize the partition sizes
	cellCount=floor(globalCellCount/np);
	int baseCellCount=cellCount;
	unsigned int offset=Rank*cellCount;
	if (Rank==np-1) cellCount=cellCount+globalCellCount-np*cellCount;

	//Implementing Parmetis
	/* ParMETIS_V3_PartMeshKway(idxtype *elmdist, idxtype *eptr, idxtype *eind, idxtype *elmwgt, int *wgtflag, int *numflag, int *ncon, int * ncommonnodes, int *nparts, float *tpwgts, float *ubvec, int *options, int *edgecut, idxtype *part, MPI_Comm) */

	/*  Definining variables
	elmdist- look into making the 5 arrays short int (for performance
		on 64 bit arch)
	eptr- like xadf
	eind- like adjncy
	elmwgt- null (element weights)
	wgtflag- 0 (no weights, can take value of 0,1,2,3 see documentation)
	numflag- 0 C-style numbers, 1 Fortran-style numbers
	ncon- 1  ( # of constraints)
	ncommonnodes- 4 ( we can probably put this to 3)
	nparts- # of processors (Note: BE CAREFUL if != to # of proc)
	tpwgts- 
	ubvec-  (balancing constraints,if needed 1.05 is a good value)
	options- [0 1 15] for default
	edgecut- output, # of edges cut (measure of communication)
	part- output, where our elements should be
	comm- most likely MPI_COMM_WORLD
	*/

	idxtype elmdist[np+1];
	idxtype *eptr;
	eptr = new idxtype[cellCount+1];
	idxtype *eind;
	int eindSize=0;
	if ((offset+cellCount)==globalCellCount) {
		eindSize=elemConnectivity.size()-elemConnIndex[offset];
	}
	else {
		eindSize=elemConnIndex[offset+cellCount]-elemConnIndex[offset]+1;
	}
	eind = new idxtype[eindSize];
	idxtype* elmwgt = NULL;
	int wgtflag=0; // no weights associated with elem or edges
	int numflag=0; // C-style numbering
	int ncon=1; // # of weights or constraints
	int ncommonnodes; ncommonnodes=3; // set to 3 for tetrahedra or mixed type

	float tpwgts[np];
	for (unsigned int p=0; p<np; ++p) tpwgts[p]=1./float(np);
	float ubvec=1.02;
	int options[3]; // default values for timing info set 0 -> 1

	options[0]=0; options[1]=1; options[2]=15;
	int edgecut ; // output
	idxtype* part = new idxtype[cellCount];

	for (unsigned int p=0;p<np;++p) elmdist[p]=p*floor(globalCellCount/np);
	elmdist[np]=globalCellCount;// Note this is because #elements mod(np) are all on last proc

	for (unsigned int c=0; c<cellCount;++c) {
		eptr[c]=elemConnIndex[offset+c]-elemConnIndex[offset];
	}
	if ((offset+cellCount)==globalCellCount) {
		 eptr[cellCount]=elemConnectivity.size()-elemConnIndex[offset];
	} else {
		eptr[cellCount]=elemConnIndex[offset+cellCount]-elemConnIndex[offset];
	}

	for (unsigned int i=0; i<eindSize; ++i) {
			eind[i]=elemConnectivity[elemConnIndex[offset]+i];
	}

	ompi_communicator_t* commWorld=MPI_COMM_WORLD;

	ParMETIS_V3_PartMeshKway(elmdist,eptr,eind, elmwgt,
	                         &wgtflag, &numflag, &ncon, &ncommonnodes,
	                         &np, tpwgts, &ubvec, options, &edgecut,
	                         part,&commWorld) ;
	
	delete[] eptr;
	delete[] eind;


	// Distribute the part list to each proc
	// Each proc has an array of length globalCellCount which says the processor number that cell belongs to [cellMap]
	int recvCounts[np];
	int displs[np];
	for (int p=0;p<np;++p) {
		recvCounts[p]=baseCellCount;
		displs[p]=p*baseCellCount;
	}
	recvCounts[np-1]=baseCellCount+globalCellCount-np*baseCellCount;
	int cellMap[globalCellCount];
	//cellMap of a cell returns which processor it is assigned to
	MPI_Allgatherv(part,cellCount,MPI_INT,cellMap,recvCounts,displs,MPI_INT,MPI_COMM_WORLD);

	// Find new local cellCount after ParMetis distribution
	cellCount=0.;
	int otherCellCounts[np]; 
	for (unsigned int p=0;p<np;p++) otherCellCounts[p]=0; 
	
	for (unsigned int c=0;c<globalCellCount;++c) {
		otherCellCounts[cellMap[c]]+=1;
		if (cellMap[c]==Rank) ++cellCount;
	}
	cout << "[I Rank=" << Rank << "] Number of Cells= " << cellCount << endl;
		
	//node.reserve(nodeCount/np);
	face.reserve(cellCount);
	cell.reserve(cellCount);

	// Create the nodes and cells for each partition
	// Run through the list of cells and check if it belongs to current partition
	// Loop through the cell's nodes
	// Mark the visited nodes so that no duplicates are created (nodeFound array).
	bool nodeFound[globalNodeCount];
	unsigned int nodeMap[globalNodeCount];

	for (unsigned int n=0;n<globalNodeCount;++n) nodeFound[n]=false;
	nodeCount=0;
	
	for (unsigned int c=0;c<globalCellCount;++c) {
		if (cellMap[c]==Rank) {
			unsigned int cellNodes[elemNodeCount[c]];
			for (unsigned int n=0;n<elemNodeCount[c];++n) {

				if (!nodeFound[elemConnectivity[elemConnIndex[c]+n]]) {
					Node temp;
					temp.id=nodeCount;
					temp.globalId=elemConnectivity[elemConnIndex[c]+n];
					temp.comp[0]=x[temp.globalId];
					temp.comp[1]=y[temp.globalId];
					temp.comp[2]=z[temp.globalId];
					node.push_back(temp);
					nodeFound[temp.globalId]=true;
					nodeMap[temp.globalId]=temp.id;
					++nodeCount;
				}
				cellNodes[n]=nodeMap[elemConnectivity[elemConnIndex[c]+n]];
			}

			Cell temp;
			temp.Construct(elemTypes[c],cellNodes);
			temp.globalId=c;
			cell.push_back(temp);

		}
	}

	cout << "[I Rank=" << Rank << "] created cells and nodes" << endl;

	// Create a mapping from node globalId's to node local id's
	int nodeGlobal2local [grid.globalNodeCount];
	// Initialize the array to -1 (meaning not in current partition)
	for (unsigned int n=0;n<grid.globalNodeCount;++n) nodeGlobal2local[n]=-1;
	// If the queried globalId is not on this partition, returns -1 as localId
	for (unsigned int n=0;n<grid.nodeCount;++n) {
		nodeGlobal2local[grid.node[n].globalId]=n;
	}

	// This stores the bc region numbers that 'some' nodes touch to
	// Those 'some' nodes happen to be the first nodes in each bc face connectivity list
	std::map<unsigned int, vector<int> > bNodeRegions;
	// For those 'some' nodes, this stores the conn list of the corresponsing bc faces
	std::map<unsigned int, vector<set<unsigned int> > > bNodeRegionNodeSets;
	// Mark nodes that touch boundaries and store which boundary(s) in node.bc set
	for (int breg=0;breg<totalnBocos;++breg) { // for each boundary condition region defined in grid file
		int bfnodeCount;
		for (int bf=0;bf<bocos[breg].size();++bf) { // for each boundary face in current region
			if (bocoConnectivity[breg][bocos[breg][bf]]!=-1) {
				// Get number of nodes of the boundary face
				if (bf==bocos[breg].size()-1) {
					bfnodeCount=bocoConnectivity[breg].size()-bocos[breg][bf];
				} else {
					bfnodeCount=bocos[breg][bf+1]-bocos[breg][bf];
				}
				unsigned int faceNode;
				// See if the face is on the current partition
				bool onCurrent=true;
				for (unsigned int i=0;i<bfnodeCount;++i) {
					faceNode=nodeGlobal2local[ bocoConnectivity[breg][bocos[breg][bf]+i] ];
					if (faceNode<0) {
						onCurrent=false;
						break;
					}
				}
				if (onCurrent) {
					// Grab the first node of the the boundary face
					faceNode=nodeGlobal2local[ bocoConnectivity[breg][bocos[breg][bf]] ];
					bNodeRegions[faceNode].push_back(breg);
					std::set<unsigned int> nodeSet;
					nodeSet.clear();
					for (unsigned int i=0;i<bfnodeCount;++i) {
						nodeSet.insert(nodeGlobal2local[ bocoConnectivity[breg][bocos[breg][bf]+i] ]);
						
					}
					bNodeRegionNodeSets[faceNode].push_back(nodeSet);
				}
			}
		}
	}

	//Create the Mesh2Dual inputs
	//idxtype elmdist [np+1] (stays the same size)
	eptr = new idxtype[cellCount+1];
	eindSize=0;
	for (unsigned int c=0;c<cellCount;++c) {
		eindSize+=cell[c].nodeCount;
	}
	eind = new idxtype[eindSize];
	// numflag and ncommonnodes previously defined
	ncommonnodes=1;
	idxtype* xadj;
	idxtype* adjncy;

	elmdist[0]=0;
	for (unsigned int p=1;p<=np;p++) elmdist[p]=otherCellCounts[p-1]+elmdist[p-1];
	eptr[0]=0;
	for (unsigned int c=1; c<=cellCount;++c) eptr[c]=eptr[c-1]+cell[c-1].nodeCount;
	int eindIndex=0;
	for (unsigned int c=0; c<cellCount;c++){
		for (unsigned int cn=0; cn<cell[c].nodeCount; ++cn) {
			eind[eindIndex]=cell[c].node(cn).globalId;
			++eindIndex;
		}	
	}

	ParMETIS_V3_Mesh2Dual(elmdist, eptr, eind, &numflag, &ncommonnodes, &xadj, &adjncy, &commWorld);
	
	// Construct the list of cells for each node
	bool flag;
	for (unsigned int c=0;c<cellCount;++c) {
		unsigned int n;
		for (unsigned int cn=0;cn<cell[c].nodeCount;++cn) {
			n=cell[c].nodes[cn];
			flag=false;
			for (unsigned int i=0;i<node[n].cells.size();++i) {
				if (node[n].cells[i]==c) {
					flag=true;
					break;
				}
			}
			if (!flag) {
				node[n].cells.push_back(c);
			}
		}
	}

	cout << "[I Rank=" << Rank << "] computed node-cell connectivity" << endl;
	
	// Construct the list of neighboring cells for each cell
	int c2;
	for (unsigned int c=0;c<cellCount;++c) {
		unsigned int n;
		for (unsigned int cn=0;cn<cell[c].nodeCount;++cn) {
			n=cell[c].nodes[cn];
			for (unsigned int nc=0;nc<node[n].cells.size();++nc) {
				c2=node[n].cells[nc];
				flag=false;
				for (unsigned int cc=0;cc<cell[c].neighborCells.size();++cc) {
					if(cell[c].neighborCells[cc]==c2) {
						flag=true;
						break;
					}
				}
				if (!flag) cell[c].neighborCells.push_back(c2);
			} // end node cell loop
		} // end cell node loop
		cell[c].neighborCellCount=cell[c].neighborCells.size();
	} // end cell loop
	
	// Set face connectivity lists
	int hexaFaces[6][4]= {
		{0,3,2,1},
		{4,5,6,7},
		{1,2,6,5},
		{0,4,7,3},
		{1,5,4,0},
		{2,3,7,6}
	};

	int prismFaces[5][4]= {
		{0,2,1,0},
		{3,4,5,0},
		{0,3,5,2},
		{1,2,5,4},
		{0,1,4,3},
	};

	int tetraFaces[4][3]= {
		{0,2,1},
		{1,2,3},
		{0,3,2},
		{0,1,3}
	};


	// Search and construct faces
	faceCount=0;

	int boundaryFaceCount[totalnBocos];
	for (int boco=0;boco<totalnBocos;++boco) boundaryFaceCount[boco]=0;
	
	
	double timeRef, timeEnd;
	if (Rank==0) timeRef=MPI_Wtime();

	// Loop through all the cells
	for (unsigned int c=0;c<cellCount;++c) {
		// Loop through the faces of the current cell
		for (unsigned int cf=0;cf<cell[c].faceCount;++cf) {
			Face tempFace;
			unsigned int *tempNodes;
			switch (elemTypes[c]) {
				case TETRA_4:
					tempFace.nodeCount=3;
					tempNodes= new unsigned int[3];
					break;
				case PENTA_6:
					if (cf<2) {
						tempFace.nodeCount=3;
						tempNodes= new unsigned int[3];
					} else {
						tempFace.nodeCount=4;
						tempNodes= new unsigned int[4];
					}
					break;
				case HEXA_8:
					tempFace.nodeCount=4;
					tempNodes= new unsigned int[4];
					break;
			}
			tempFace.id=faceCount;
			// Assign current cell as the parent cell
			tempFace.parent=c;
			// Assign boundary type as internal by default, will be overwritten later
			tempFace.bc=-1;
			// Store the nodes of the current face
			//if (faceCount==face.capacity()) face.reserve(int (face.size() *0.10) +face.size()) ; //TODO check how the size grows by default

			for (unsigned int fn=0;fn<tempFace.nodeCount;++fn) {
				switch (elemTypes[c]) {
					case TETRA_4: tempNodes[fn]=cell[c].node(tetraFaces[cf][fn]).id; break;
					case PENTA_6: tempNodes[fn]=cell[c].node(prismFaces[cf][fn]).id; break;
					case HEXA_8: tempNodes[fn]=cell[c].node(hexaFaces[cf][fn]).id; break;
				}
			}
			// Find the neighbor cell
			bool internal=false;
			bool unique=true;
			// Loop cells neighboring the first node of the current face
			for (unsigned int nc=0;nc<node[tempNodes[0]].cells.size();++nc) {
				// i is the neighbor cell index
				unsigned int i=node[tempNodes[0]].cells[nc];
				// if neighbor cell is not the current cell itself, and it has the same nodes as the face
				if (i!=c && cell[i].HaveNodes(tempFace.nodeCount,tempNodes)) {
					// if the neighbor cell index is smaller then the current cell index,
					// it has already been processed
					if (i>c) {
						tempFace.neighbor=i;
						internal=true;
					} else {
						unique=false;
					}
				}
			}
			if (unique) {
				for (unsigned int fn=0;fn<tempFace.nodeCount;++fn) tempFace.nodes.push_back(tempNodes[fn]);
				if (!internal) { // the face is either at inter-partition or boundary
					// Unassigned boundary face, as the bc faces and inter-partition faces are found,
					// this will be overwritten. At the end, there should be no face with this value.
					tempFace.bc=-2; 
					// Loop the current face nodes
					for (unsigned int i=0;i<tempFace.nodeCount;++i) {
						unsigned int nn=tempNodes[i];
						bool match;
						// Find if this node has the boundary region data
						// If the face is at a boundary, one of the nodes should have this
						if(bNodeRegions.find(nn)!=bNodeRegions.end()) {
							// Loop through the different boundary regions assigned for that node					
							for (int j=0;j<bNodeRegions[nn].size();++j) {
								// Loop the current face nodes
								match=true;
								for (int k=0;k<tempFace.nodeCount;++k) {
									// if the current face node is not in the current boundary face's node list
									if (bNodeRegionNodeSets[nn][j].find(tempNodes[k])==bNodeRegionNodeSets[nn][j].end()) {
 										match=false;
 										break;
									}
								}
								if (match) {
									tempFace.bc=bNodeRegions[nn][j];
									boundaryFaceCount[bNodeRegions[nn][j]]++;
									break;
								}
							} // For each boundary regions that node belongs to
							if (match) break;
						} // if node has the boundary region data
					} // for tempFace nodes
				} // if not internal
				
				face.push_back(tempFace);
				cell[c].faces.push_back(tempFace.id);
				if (internal) cell[tempFace.neighbor].faces.push_back(tempFace.id);
				++faceCount;
			}
			delete [] tempNodes;
		} //for face cf
	} // for cells c

	cout << "[I Rank=" << Rank << "] Number of Faces=" << faceCount << endl;
	for (int boco=0;boco<totalnBocos;++boco) cout << "[I Rank=" << Rank << "] Number of Faces on BC_" << boco+1 << "=" << boundaryFaceCount[boco] << endl;
	
	if (Rank==0) {
		timeEnd=MPI_Wtime();
		cout << "[I Rank=" << Rank << "] Time spent on finding faces= " << timeEnd-timeRef << " sec" << endl;
	}
	
	// Determine and mark faces adjacent to other partitions
	// Create ghost elemets to hold the data from other partitions
	ghostCount=0;

	if (np==1) {
		for (unsigned int c=0;c<cellCount;++c) {
			cell[c].ghostCount=cell[c].ghosts.size();
		} // end cell loop
	} else {

		int counter=0;
		int cellCountOffset[np];
		
		for (unsigned int p=0;p<np;++p) {
			cellCountOffset[p]=counter;
			counter+=otherCellCounts[p];
		}
		
		// Now find the metis2global mapping
		int metis2global[globalCellCount];
		int counter2[np];
		for (unsigned int p=0;p<np;++p) counter2[p]=0;
		for (unsigned int c=0;c<globalCellCount;++c) {
			metis2global[cellCountOffset[cellMap[c]]+counter2[cellMap[c]]]=c;
			counter2[cellMap[c]]++;
		}

		int foundFlag[globalCellCount];
		for (unsigned int c=0; c<globalCellCount; ++c) foundFlag[c]=0;

		unsigned int parent, metisIndex, gg, matchCount;
		map<unsigned int,unsigned int> ghostGlobal2local;
		
		map<int,set<int> > nodeCellSet;
		Vec3D nodeVec;
		
		for (unsigned int f=0; f<faceCount; ++f) {
			if (face[f].bc==-2) { // if an assigned boundary face
				parent=face[f].parent;
				for (unsigned int adjCount=0;adjCount<(xadj[parent+1]-xadj[parent]);++adjCount)  {
					metisIndex=adjncy[xadj[parent]+adjCount];
					gg=metis2global[metisIndex];

					if (metisIndex<cellCountOffset[Rank] || metisIndex>=(cellCount+cellCountOffset[Rank])) {
						matchCount=0;
						for (unsigned int fn=0;fn<face[f].nodeCount;++fn) {
							set<int> tempSet;
							nodeCellSet.insert(pair<unsigned int,set<int> >(face[f].nodes[fn],tempSet) );
							for (unsigned int gn=0;gn<elemNodeCount[gg];++gn) {
								if (elemConnectivity[elemConnIndex[gg]+gn]==face[f].node(fn).globalId) {
									nodeCellSet[face[f].nodes[fn]].insert(gg);
									++matchCount;
								}
							}
						}
						if (matchCount>0 && foundFlag[gg]==0) {
							if (matchCount>=3) foundFlag[gg]=3;
							if (matchCount<3) foundFlag[gg]=matchCount;
							Ghost temp;
							temp.globalId=gg;
							temp.partition=cellMap[gg];
							// Calculate the centroid
							temp.centroid=0.;
							for (unsigned int gn=0;gn<elemNodeCount[gg];++gn) {
								nodeVec.comp[0]=x[elemConnectivity[elemConnIndex[gg]+gn]];
								nodeVec.comp[1]=y[elemConnectivity[elemConnIndex[gg]+gn]];
								nodeVec.comp[2]=z[elemConnectivity[elemConnIndex[gg]+gn]];
								temp.centroid+=nodeVec;
							}
							temp.centroid/=double(elemNodeCount[gg]);
							ghost.push_back(temp);
							ghostGlobal2local.insert(pair<unsigned int,unsigned int>(gg,ghostCount));
							++ghostCount;
						}
						if (matchCount>=3) {
							foundFlag[gg]=3;
							face[f].bc=-1*ghostGlobal2local[gg]-3;
						} 
					}
				}	
			}
		}

		map<int,set<int> >::iterator mit;
		set<int>::iterator sit;
		for ( mit=nodeCellSet.begin() ; mit != nodeCellSet.end(); mit++ ) {
			for ( sit=(*mit).second.begin() ; sit != (*mit).second.end(); sit++ ) {
				node[(*mit).first].ghosts.push_back(ghostGlobal2local[*sit]);
				//cout << (*mit).first << "\t" << *sit << endl; // DEBUG
			}	
		}
		
		// Construct the list of neighboring ghosts for each cell
		int g;
		bool flag;
		for (unsigned int c=0;c<cellCount;++c) {
			unsigned int n;
			for (unsigned int cn=0;cn<cell[c].nodeCount;++cn) {
				n=cell[c].nodes[cn];
				for (unsigned int ng=0;ng<node[n].ghosts.size();++ng) {
					g=node[n].ghosts[ng];
					flag=false;
					for (unsigned int cg=0;cg<cell[c].ghosts.size();++cg) {
						if(cell[c].ghosts[cg]==g) {
							flag=true;
							break;
						}
					} // end cell ghost loop
					if (flag==false) {
						cell[c].ghosts.push_back(g);
						ghost[g].cells.push_back(c);
					}
				} // end node ghost loop
			} // end cell node loop
			cell[c].ghostCount=cell[c].ghosts.size();
		} // end cell loop
	} // if (np!=1) 

	cout << "[I Rank=" << Rank << "] Number of Inter-Partition Ghost Cells= " << ghostCount << endl;
	
	// Now loop through faces and calculate centroids and areas
	for (unsigned int f=0;f<faceCount;++f) {
		Vec3D centroid=0.;
		Vec3D areaVec=0.;
		for (unsigned int n=0;n<face[f].nodeCount;++n) {
			centroid+=face[f].node(n);
		}
		centroid/=double(face[f].nodeCount);
		face[f].centroid=centroid;
		for (unsigned int n=0;n<face[f].nodeCount-1;++n) {
			areaVec+=0.5* (face[f].node(n)-centroid).cross(face[f].node(n+1)-centroid);
		}
		areaVec+=0.5*(face[f].node(face[f].nodeCount-1)-centroid).cross(face[f].node(0)-centroid);
		if (areaVec.dot(centroid-cell[face[f].parent].centroid) <0.) {
			// [TBM] Need to swap the face and reflect the area vector
			cout << "face " << f << " should be swapped" << endl;
		}
		face[f].area=fabs(areaVec);
		face[f].normal=areaVec/face[f].area;
	}
			
	// Loop through the cells and calculate the volumes and length scales
	double totalVolume=0.;
	for (unsigned int c=0;c<cellCount;++c) {
		cell[c].lengthScale=1.e20;
		double volume=0.,height;
		unsigned int f;
		for (unsigned int cf=0;cf<cell[c].faceCount;++cf) {
			f=cell[c].faces[cf];
			height=fabs(face[f].normal.dot(face[f].centroid-cell[c].centroid));
			volume+=1./3.*face[f].area*height;
			cell[c].lengthScale=min(cell[c].lengthScale,height);
		}
		cell[c].volume=volume;
		totalVolume+=volume;
	}
	
	cout << "[I Rank=" << Rank << "] Total Volume= " << setw(16) << setprecision(8) << scientific << totalVolume << endl;
	return 0;

}


Node::Node(double x, double y, double z) {
	comp[0]=x;
	comp[1]=y;
	comp[2]=z;
}

Cell::Cell(void) {
	;
}

int Cell::Construct(const ElementType_t elemType, unsigned int nodeList[]) {

	switch (elemType) {
		case TETRA_4:
			faceCount=4;
			nodeCount=4;
			break;
		case PENTA_6:
			faceCount=5;
			nodeCount=6;
			break;
		case HEXA_8:
			faceCount=6;
			nodeCount=8;
			break;
	}
	type=elemType;
	nodes.reserve(nodeCount);
	if (nodeCount==0) {
		cerr << "[!! proc " << Rank << " ] Number of nodes of the cell must be specified before allocation" << endl;
		return -1;
	} else {
		centroid=0.;
		for (unsigned int i=0;i<nodeCount;++i) {
			nodes.push_back(nodeList[i]);
			centroid+=node(i);
		}
		centroid/=nodeCount;
		return 0;
	}
}

bool Cell::HaveNodes(unsigned int &nodelistsize, unsigned int nodelist []) {
	unsigned int matchCount=0,nodeId;
	for (unsigned int j=0;j<nodeCount;++j) {
		nodeId=node(j).id;
		for (unsigned int i=0;i<nodelistsize;++i) {
			if (nodelist[i]==nodeId) ++matchCount;
			if (matchCount==3) return true;
		}
	}
	return false;
}

int Grid::face_exists(int &parentCell) {
	//int s=face.size();
	for (int f=0;f<faceCount;++f) {
		if (face[f].parent==parentCell) return 1;
	}
	return 0;
}

Node& Cell::node(int n) {
	return grid.node[nodes[n]];
};

Face& Cell::face(int f) {
	return grid.face[faces[f]];
};

Node& Face::node(int n) {
	return grid.node[nodes[n]];
};

void Grid::nodeAverages() {
	
	double **a;
	double *b;
	double *weights;
	
	for (unsigned int n=0;n<nodeCount;++n) {
		set<int> stencil; //interpolation stencil
		set<int>::iterator sit,sit1,sit2;
		int c1,c2,c3,c4;
		// Initialize stencil to nearest neighbor cells
		for (int nc=0;nc<node[n].cells.size();++nc) stencil.insert(node[n].cells[nc]);
		string method;
		Vec3D planeNormal;
		// if stencil doesn't have at least 4 points, expand it to include 2nd nearest neighbor cells
		if (stencil.size()<4) {
			for (int nc=0;nc<node[n].cells.size();++nc) {
				int ncell=node[n].cells[nc];
				for (int cc=0;cc<cell[ncell].neighborCells.size();++cc) {
					stencil.insert(cell[ncell].neighborCells[cc]);
				}
			}
		}
		if (stencil.size()>=4) { // Candidate for tetra interpolation
			// Pick first 3 points in the stencil
			sit=stencil.begin();
			c1=*sit; sit++;
			c2=*sit; sit++;
			c3=*sit;
			// Calculate the normal vector of the plane they form
			planeNormal=((cell[c2].centroid-cell[c1].centroid).cross(cell[c3].centroid-cell[c2].centroid)).norm();
			// Check the remaining stencil cell centroids to see if they are on the same plane too
			method="tri";
			for (sit=sit;sit!=stencil.end();sit++) {
				// If the centroid lies on the plane fomed by first three
				if (fabs(planeNormal.dot((cell[*sit].centroid-cell[*(stencil.begin())].centroid).norm()))>1.e-7) {
					method="tetra";
					// Finding on off-plane point is enough
					// That means tetra method can be used
					break;
				}
			}
		}
		else if (stencil.size()==3) { method="tri";}
		else if (stencil.size()==2) { method="line";}
		else if (stencil.size()==1) { method="point";}

		// TODO implement tetra method
		if (method=="tri") {
			// Find out best combination neigboring cells for triangulation
			double edgeRatio,bestEdgeRatio=0.; // smallest edge length / largest edge length
			double edge1,edge2,edge3,maxEdge,minEdge;
			for (sit=stencil.begin();sit!=stencil.end();sit++) {
				sit1=sit; sit1++;
				for (sit1=sit1;sit1!=stencil.end();sit1++) {
					sit2=sit; sit2++;
					for (sit2=sit2;sit2!=stencil.end();sit2++) {
						edge1=fabs(cell[*sit].centroid-cell[*sit1].centroid);
						edge2=fabs(cell[*sit].centroid-cell[*sit2].centroid);
						edge3=fabs(cell[*sit1].centroid-cell[*sit2].centroid);
						maxEdge=max(edge1,edge2); maxEdge=max(maxEdge,edge3);
						minEdge=min(edge1,edge2); minEdge=min(minEdge,edge3);
						edgeRatio=minEdge/maxEdge; // always smaller than 1, ideal is 1
						// TODO this criteria is not good, even a zero area triangle might return 0.5 (3 points on the same line)
						// Switch to area/perimeter^2
						if (edgeRatio>bestEdgeRatio) {
							bestEdgeRatio=edgeRatio;
							c1=*sit;
							c2=*sit1;
							c3=*sit2;
						}
					}
				}
			}
			// 
			// Work on plane coordinates
			// pp is the node point projected onto the plane
			Vec3D pp,origin,basis1,basis2;
			origin=cell[c1].centroid;
			basis1=(cell[c2].centroid-origin).norm();
			planeNormal=(basis1.cross((cell[c3].centroid-origin).norm()));
			// If the plane normal magnitude is zero, then
			// the 3 points do not form a plane but a line
			if (fabs(planeNormal)<1.e-7) {
				method="line";
			}
			if (method!="line") { // If method is not switched to line
				// Normalize the plane normal vector
				planeNormal/=fabs(planeNormal);
				basis2=-1.*(basis1.cross(planeNormal)).norm();
				// Project the node point to the plane
				pp=node[n]-origin;
				pp-=pp.dot(planeNormal)*planeNormal;
				pp+=origin;
				// Form the linear system
				a = new double* [3];
				for (int i=0;i<3;++i) a[i]=new double[3];
				b = new double[3];
				weights= new double [3];
				a[0][0]=(cell[c1].centroid).dot(basis1);
				a[0][1]=(cell[c2].centroid).dot(basis1);
				a[0][2]=(cell[c3].centroid).dot(basis1);
				a[1][0]=(cell[c1].centroid).dot(basis2);
				a[1][1]=(cell[c2].centroid).dot(basis2);
				a[1][2]=(cell[c3].centroid).dot(basis2);
				a[2][0]=1.;
				a[2][1]=1.;
				a[2][2]=1.;
				b[0]=pp.dot(basis1);
				b[1]=pp.dot(basis2);
				b[2]=1.;
				// Solve the 3x3 linear system by Gaussion Elimination
				gelimd(a,b,weights,3);
				// Insert the weights
				node[n].average.clear();
				node[n].average.insert(pair<int,double>(c1,weights[0]));
				node[n].average.insert(pair<int,double>(c2,weights[1]));
				node[n].average.insert(pair<int,double>(c3,weights[2]));
			} // end if method is not switched to line
		} // end if method is tri
		if (method=="line") {
			// Chose the closest 2 points in stencil
			double distanceMin=1.e20;
			double distance;
			for (sit=stencil.begin();sit!=stencil.end();sit++) {
				distance=fabs(node[n]-cell[*sit].centroid);
				if (distance<distanceMin) {
					distanceMin=distance;
					sit2=sit1;
					sit1=sit;
				}
			}
			Vec3D pp,p1,p2,direction;
			double a,weight1,weight2;
			p1=cell[*sit1].centroid;
			p2=cell[*sit2].centroid;
			direction=(p2-p1).norm();
			a=(node[n]-p1).dot(direction)/fabs(p2-p1);
			weight1=a;
			weight2=1.-a;
			node[n].average.clear();
			node[n].average.insert(pair<int,double>(*sit1,weight1));
			node[n].average.insert(pair<int,double>(*sit2,weight2));
		} // end if method is line
		if (method=="point") {
			node[n].average.clear();
			node[n].average.insert(pair<int,double>(*(stencil.begin()),1.));
		} // end if method is point
	} // node loop
	
	return;
} // end Grid::nodeAverages()

void Grid::faceAverages() {

	unsigned int n;
	map<int,double>::iterator it;
	map<int,double>::iterator fit;
	bool simple;
	set<int>::iterator bcit;
	double cell2face;
	double weightSum;
	double weight,weightSumCheck;

 	for (unsigned int f=0;f<faceCount;++f) {
		face[f].average.clear();
		simple=false;
		// if one of the face nodes is touching a boundary other than slip,
		// switch to simpler averaging
		for (unsigned int fn=0;fn<face[f].nodeCount;++fn) {
			n=face[f].nodes[fn];
			for (bcit=node[n].bcs.begin();bcit!=node[n].bcs.end();bcit++) {
				if (bc.region[*bcit].type!="symmetry") {
					simple=true;
					break;
				}
			}
		}
		simple=false;
		// This method averages the parent and neighbor cell values
		if (simple) {
			unsigned int c;
			weightSum=0.;
			for (int i=0;i<2;++i) {
				if (i==0) { c=face[f].parent; } else { c=face[f].neighbor;}
				cell2face=(cell[c].centroid-face[f].centroid).dot(face[f].normal);
				//cell2face=fabs(cell[c].centroid-face[f].centroid);
				//weight=1./(cell2face*cell2face);
				//weight=cell[c].volume;
				weightSum+=weight;
				face[f].average.insert(pair<int,double>(c,weight));
			}
			for (fit=face[f].average.begin();fit!=face[f].average.end();fit++) {
				face[f].average[(*fit).first]/=weightSum;
			}
			
		} else { // This method uses node based averaging
			weightSum=0.;
			for (unsigned int fn=0;fn<face[f].nodeCount;++fn) {
				n=face[f].nodes[fn];
				weight=fabs(node[n]-face[f].centroid);
				weight*=weight;
				weightSum+=weight;
				for ( it=node[n].average.begin() ; it != node[n].average.end(); it++ ) {
					if (face[f].average.find((*it).first)!=face[f].average.end()) { // if the cell contributing to the node average is already contained in the face average map
						face[f].average[(*it).first]+=weight*(*it).second;
					} else {
						face[f].average.insert(pair<int,double>((*it).first,weight*(*it).second));
					}
				}
			}
			weightSumCheck=0.;
			for ( fit=face[f].average.begin() ; fit != face[f].average.end(); fit++ ) {
				(*fit).second/=weightSum;
				weightSumCheck+=(*fit).second;
			}
			//cout << setw(16) << setprecision(8) << scientific << weightSumCheck << endl;
		} // end face node loop

				
	} // end face loop

} // end Grid::faceAverages()

void Grid::gradMaps() {
	
	unsigned int f;
	map<int,double>::iterator it;
	Vec3D areaVec;
	for (unsigned int c=0;c<cellCount;++c) {
		cell[c].gradMap.clear();
		for (unsigned int cf=0;cf<cell[c].faceCount;++cf) {
			f=cell[c].faces[cf];
			if (face[f].bc<0 ) { // if internal or interpartition face
				areaVec=face[f].normal*face[f].area/cell[c].volume;
				if (face[f].parent!=c) areaVec*=-1.;
				for (it=face[f].average.begin();it!=face[f].average.end();it++) {
					if (cell[c].gradMap.find((*it).first)!=cell[c].gradMap.end()) { // if the cell contributing to the face average is already contained in the cell gradient map
						cell[c].gradMap[(*it).first]+=(*it).second*areaVec;
					} else {
						cell[c].gradMap.insert(pair<int,Vec3D>((*it).first,(*it).second*areaVec));
					}
				}
			} // end if internal face
		} // end cell face loop
	} // end cell loop
	
} // end Grid::gradMaps()
	

void Grid::gradients(void) {

	// Calculate cell gradients

	map<int,Vec3D>::iterator it;
	map<int,double>::iterator fit;
	unsigned int f;
	Vec3D faceVel,areaVec;
	double faceRho,faceP,faceK,faceOmega;
	double Mach;
	
	for (unsigned int c=0;c<cellCount;++c) {
		// Initialize all gradients to zero
		for (unsigned int i=0;i<7;++i) cell[c].grad[i]=0.;
		// Add internal and interpartition face contributions
		for (it=cell[c].gradMap.begin();it!=cell[c].gradMap.end(); it++ ) {
			if ((*it).first>=0) { // if contribution is coming from a real cell
				cell[c].grad[0]+=(*it).second*cell[(*it).first].rho;
				for (unsigned int i=1;i<4;++i) cell[c].grad[i]+=(*it).second*cell[(*it).first].v.comp[i-1];
				cell[c].grad[4]+=(*it).second*cell[(*it).first].p;
				cell[c].grad[5]+=(*it).second*cell[(*it).first].k;
				cell[c].grad[6]+=(*it).second*cell[(*it).first].omega;
			} else { // if contribution is coming from a ghost cell
				cell[c].grad[0]+=(*it).second*ghost[-1*((*it).first+1)].rho;
				for (unsigned int i=1;i<4;++i) cell[c].grad[i]+=(*it).second*ghost[-1*((*it).first+1)].v.comp[i-1];
				cell[c].grad[4]+=(*it).second*ghost[-1*((*it).first+1)].p;
				cell[c].grad[5]+=(*it).second*ghost[-1*((*it).first+1)].k;
				cell[c].grad[6]+=(*it).second*ghost[-1*((*it).first+1)].omega;
			}
		} // end gradMap loop

		// Add boundary face contributions
		for (unsigned int cf=0;cf<cell[c].faceCount;++cf) {
			f=cell[c].faces[cf];
			if (face[f].bc>=0) { // if a boundary face
				areaVec=face[f].area*face[f].normal/cell[c].volume;
				if (face[f].parent!=c) areaVec*=-1.;

				faceVel=0.;faceRho=0.;faceP=0.;faceK=0.;faceOmega=0.;
				for (fit=face[f].average.begin();fit!=face[f].average.end();fit++) {
					if ((*fit).first>=0) { // if contribution is coming from a real cell
						faceRho+=(*fit).second*cell[(*fit).first].rho;
						faceVel+=(*fit).second*cell[(*fit).first].v;
						faceP+=(*fit).second*cell[(*fit).first].p;
						faceK+=(*fit).second*cell[(*fit).first].k;
						faceOmega+=(*fit).second*cell[(*fit).first].omega;
					} else { // if contribution is coming from a ghost cell
						faceRho+=(*fit).second*ghost[-1*((*fit).first+1)].rho;
						faceVel+=(*fit).second*ghost[-1*((*fit).first+1)].v;
						faceP+=(*fit).second*ghost[-1*((*fit).first+1)].p;
						faceK+=(*fit).second*ghost[-1*((*fit).first+1)].k;
						faceOmega+=(*fit).second*ghost[-1*((*fit).first+1)].omega;
					}
				}
				
				if (bc.region[face[f].bc].type=="inlet") {
					faceRho=bc.region[face[f].bc].rho;
					faceVel=bc.region[face[f].bc].v;
					faceK=bc.region[face[f].bc].k;
					faceOmega=bc.region[face[f].bc].omega;
					//faceP=cell[face[f].parent].p;
					//faceP=bc.region[face[f].bc].p;
				}
				
				if (bc.region[grid.face[f].bc].type=="outlet" &&
					bc.region[grid.face[f].bc].kind=="fixedPressure") {
					// find Mach number
// 					Mach=(cell[c].v.dot(face[f].normal))/sqrt(Gamma*(cell[c].p+Pref)/cell[c].rho);
// 					if (Mach<1.) faceP=bc.region[face[f].bc].p;
				}
				// Kill the wall normal component for slip or symmetry, pressure and rho is extrapolated
				if (bc.region[face[f].bc].type=="slip" | bc.region[face[f].bc].type=="symmetry") faceVel-=faceVel.dot(face[f].normal)*face[f].normal;
				// Kill the velocity for no-slip, pressure and rho is extrapolated
				if (bc.region[face[f].bc].type=="noslip") faceVel=0.;

				cell[c].grad[0]+=faceRho*areaVec;
				cell[c].grad[1]+=faceVel.comp[0]*areaVec;
				cell[c].grad[2]+=faceVel.comp[1]*areaVec;
				cell[c].grad[3]+=faceVel.comp[2]*areaVec;
				cell[c].grad[4]+=faceP*areaVec;

			} // end if a boundary face
		} // end cell face loop
		
	} // end cell loop
 	
} // end Grid::gradients(void)

void Grid::limit_gradients(string limiter, double sharpeningFactor) {
	
	unsigned int neighbor,g;
	Vec3D maxGrad[7],minGrad[7];
	
	for (unsigned int c=0;c<cellCount;++c) {

		// Initialize min and max to current cells values
		for (unsigned int i=0;i<7;++i) maxGrad[i]=minGrad[i]=cell[c].grad[i];
		// Find extremes in neighboring real cells
		for (unsigned int cc=0;cc<cell[c].neighborCellCount;++cc) {
			neighbor=cell[c].neighborCells[cc];
			for (unsigned int var=0;var<7;++var) {
				for (unsigned int comp=0;comp<3;++comp) {
					maxGrad[var].comp[comp]=max(maxGrad[var].comp[comp],(1.-sharpeningFactor)*cell[neighbor].grad[var].comp[comp]+sharpeningFactor*cell[c].grad[var].comp[comp]);
					minGrad[var].comp[comp]=min(minGrad[var].comp[comp],(1.-sharpeningFactor)*cell[neighbor].grad[var].comp[comp]+sharpeningFactor*cell[c].grad[var].comp[comp]);
				}
			}
		}
		// Find extremes in neighboring ghost cells
		for (unsigned int cg=0;cg<cell[c].ghostCount;++cg) {
			g=cell[c].ghosts[cg];
			for (unsigned int var=0;var<7;++var) {
				for (unsigned int comp=0;comp<3;++comp) {
					maxGrad[var].comp[comp]=max(maxGrad[var].comp[comp],(1.-sharpeningFactor)*ghost[g].grad[var].comp[comp]+sharpeningFactor*cell[c].grad[var].comp[comp]);
					minGrad[var].comp[comp]=min(minGrad[var].comp[comp],(1.-sharpeningFactor)*ghost[g].grad[var].comp[comp]+sharpeningFactor*cell[c].grad[var].comp[comp]);
				}
			}
		}
		if(limiter=="superbee") for (unsigned int var=0;var<7;++var) for (unsigned int comp=0;comp<3;++comp) cell[c].limited_grad[var].comp[comp]=superbee(maxGrad[var].comp[comp],minGrad[var].comp[comp]);
		if(limiter=="minmod") for (unsigned int var=0;var<7;++var) for (unsigned int comp=0;comp<3;++comp) cell[c].limited_grad[var].comp[comp]=minmod(maxGrad[var].comp[comp],minGrad[var].comp[comp]);

	}

} // end Grid::limit_gradients(string limiter, double sharpeningFactor)

int gelimd(double **a,double *b,double *x, int n)
{
	double tmp,pvt,*t;
	int i,j,k,itmp;

	for (i=0;i<n;i++) {             // outer loop on rows
		pvt = a[i][i];              // get pivot value
		if (fabs(pvt) < EPS) {
			for (j=i+1;j<n;j++) {
				if(fabs(pvt = a[j][i]) >= EPS) break;
			}
			if (fabs(pvt) < EPS) return 1;     // nowhere to run!
			t=a[j];                 // swap matrix rows...
			a[j]=a[i];
			a[i]=t;
			tmp=b[j];               // ...and result vector
			b[j]=b[i];
			b[i]=tmp;
		}
// (virtual) Gaussian elimination of column
		for (k=i+1;k<n;k++) {       // alt: for (k=n-1;k>i;k--)
			tmp = a[k][i]/pvt;
			for (j=i+1;j<n;j++) {   // alt: for (j=n-1;j>i;j--)
				a[k][j] -= tmp*a[i][j];
			}
			b[k] -= tmp*b[i];
		}
	}
// Do back substitution
	for (i=n-1;i>=0;i--) {
		x[i]=b[i];
		for (j=n-1;j>i;j--) {
			x[i] -= a[i][j]*x[j];
		}
		x[i] /= a[i][i];
	}
	return 0;
}
