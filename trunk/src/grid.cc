#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <iomanip>
using namespace std;
#include<cmath>
#include <cgnslib.h>
#include <parmetis.h>


#include "grid.h"
#include "bc.h"

extern Grid grid;
extern BC bc;
extern int np, rank;

Grid::Grid() {
	;
}

int Grid::read(string fname) {
	fstream file;
	fileName=fname;
	file.open(fileName.c_str());
	if (file.is_open()) {
		if (rank==0) cout << "* Found grid file " << fileName  << endl;
		file.close();
		ReadCGNS();
		return 1;
	} else {
		if (rank==0) cerr << "[!!] Grid file "<< fileName << " could not be found." << endl;
		return 0;
	}
}

int Grid::ReadCGNS() {

	int fileIndex,baseIndex,zoneIndex,sectionIndex;
	char zoneName[20],sectionName[20]; //baseName[20]
	//  int nBases,cellDim,physDim;
	int size[3];

	cg_open(fileName.c_str(),MODE_READ,&fileIndex);
	baseIndex=1;
	zoneIndex=1;

	cg_zone_read(fileIndex,baseIndex,zoneIndex,zoneName,size);

	globalNodeCount=size[0];
	globalCellCount=size[1];
	globalFaceCount=0;

	if (rank==0) {
		cout << "* Total number of nodes: " << globalNodeCount << endl;
		cout << "* Total number of cells: " << globalCellCount << endl;
	}

	// Initialize the partition sizes
	cellCount=floor(globalCellCount/np);
	int baseCellCount=cellCount;
	unsigned int offset=rank*cellCount;
	if (rank==np-1) cellCount=cellCount+globalCellCount-np*cellCount;

	int nodeStart[3],nodeEnd[3];
	nodeStart[0]=nodeStart[1]=nodeStart[2]=1;
	nodeEnd[0]=nodeEnd[1]=nodeEnd[2]=globalNodeCount;

	double x[globalNodeCount],y[globalNodeCount],z[globalNodeCount];

	cg_coord_read(fileIndex,baseIndex,zoneIndex,"CoordinateX",RealDouble,nodeStart,nodeEnd,&x);
	cg_coord_read(fileIndex,baseIndex,zoneIndex,"CoordinateY",RealDouble,nodeStart,nodeEnd,&y);
	cg_coord_read(fileIndex,baseIndex,zoneIndex,"CoordinateZ",RealDouble,nodeStart,nodeEnd,&z);

	// Determine the number of sections in the zone
	int sectionCount;
	cg_nsections(fileIndex,baseIndex,zoneIndex, &sectionCount);

	if (rank==0) cout << "* Number of sections found in zone " << zoneName << ": " << sectionCount << endl;

	ElementType_t elemType;
	int elemNodeCount;
	int elemStart,elemEnd,nBndCells,parentFlag;

	unsigned int section=0;
	cg_section_read(fileIndex,baseIndex,zoneIndex,section+1,sectionName,&elemType,&elemStart,&elemEnd,&nBndCells,&parentFlag);
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
	}
	int elemNodes[elemEnd-elemStart+1][elemNodeCount];
	cg_elements_read(fileIndex,baseIndex,zoneIndex,section+1,*elemNodes,0);

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

	idxtype elmdist[np + 1];
	idxtype *eptr;
	eptr = new idxtype[cellCount+1];
	idxtype *eind;
	eind = new idxtype[cellCount*elemNodeCount];

	idxtype* elmwgt = NULL;
	int wgtflag=0; // no weights associated with elem or edges
	int numflag=0; // C-style numbering
	int ncon=1; // # of weights or constraints
	int ncommonnodes; ncommonnodes=3; // set to 3 for tetrahedra or mixed type

	float tpwgts[np];
	for (unsigned int p=0; p<np; ++p) tpwgts[p]=1./float(np);
	float ubvec=1.05;
	int options[3]; // default values for timing info set 0 -> 1

	options[0]=0; options[1]=1; options[2]=15;
	int edgecut ; // output
	idxtype* part = new idxtype[cellCount];

	for (unsigned int p=0;p<np;++p) elmdist[p]=p*floor(globalCellCount/np);
	elmdist[np]=globalCellCount;// Note this is because #elements mod(np) are all on last proc

	eptr[0]=0;
	for (unsigned int c=1; c<=cellCount;++c) eptr[c]=eptr[c-1]+elemNodeCount;

	int index=0;
	for (unsigned int c=0; c<cellCount; ++c) {
		for (unsigned int nc=0; nc<elemNodeCount; ++nc) {
			eind[index]=elemNodes[c+offset][nc]-1;
			++index;
		}
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
		if (cellMap[c]==rank) ++cellCount;
	}
	cout << "* Processor " << rank << " has " << cellCount << " cells" << endl;
		
	node.reserve(nodeCount/np);
	face.reserve(cellCount);
	cell.reserve(cellCount);

	// Create the nodes and cells for each partition
	// Run through the list of cells and check if it belongs to current partition
	// Loop through the cell's nodes
	// Mark the visited nodes so that no duplicates are created (nodeFound array).
	bool nodeFound[globalNodeCount];
	unsigned int nodeMap[globalNodeCount];
	unsigned int cellNodes[elemNodeCount];
	for (unsigned int n=1;n<=globalNodeCount;++n) nodeFound[n]=false;
	nodeCount=0;
	for (unsigned int c=0;c<globalCellCount;++c) {
		if (cellMap[c]==rank) {
			for (unsigned int n=0;n<elemNodeCount;++n) {
				if (!nodeFound[elemNodes[c][n]]) {
					Node temp;
					temp.id=nodeCount;
					temp.globalId=elemNodes[c][n]-1;
					temp.comp[0]=x[elemNodes[c][n]-1];
					temp.comp[1]=y[elemNodes[c][n]-1];
					temp.comp[2]=z[elemNodes[c][n]-1];
					node.push_back(temp);
					nodeFound[elemNodes[c][n]]=true;
					nodeMap[elemNodes[c][n]]=temp.id;
					++nodeCount;
				}
				cellNodes[n]=nodeMap[elemNodes[c][n]];
			}
			Cell temp;
			temp.Construct(elemType,cellNodes);
			temp.globalId=c;
			cell.push_back(temp);
		}
	}

	cout << "* Processor " << rank << " has created its cells and nodes" << endl;

	//Create the Mesh2Dual inputs
	//idxtype elmdist [np+1] (stays the same size)
	eptr = new idxtype[cellCount+1];
	eind = new idxtype[cellCount*elemNodeCount];
	// numflag and ncommonnodes previously defined
	idxtype* xadj;
	idxtype* adjncy;
	
	elmdist[0]=0;
	for (unsigned int p=1;p<=np;p++) elmdist[p]=otherCellCounts[p-1]+elmdist[p-1];
	eptr[0]=0;
	for (unsigned int c=1; c<=cellCount;++c) eptr[c]=eptr[c-1]+elemNodeCount;
	index=0;
	for (unsigned int c=0; c<cellCount;c++){
		for (unsigned int cn=0; cn<cell[c].nodeCount; ++cn) {
			eind[index]=cell[c].node(cn).globalId;
			++index;
		}	
	}
	
	ParMETIS_V3_Mesh2Dual(elmdist, eptr, eind, &numflag, &ncommonnodes, &xadj, &adjncy, &commWorld);

	// Construct the list of cells for each node
	int flag;
	for (unsigned int c=0;c<cellCount;++c) {
		int n;
		for (unsigned int cn=0;cn<cell[c].nodeCount;++cn) {
			n=cell[c].nodes[cn];
			flag=0;
			for (unsigned int i=0;i<node[n].cells.size();++i) {
				if (node[n].cells[i]==c) {
					flag=1;
					break;
				}
			}
			if (!flag) {
				node[n].cells.push_back(c);
			}
		}
	}

	cout << "* Processor " << rank << " has computed its node-cell connectivity" << endl;
	
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

	int boundaryFaceCount=0;
	// Loop through all the cells

	double timeRef, timeEnd;
	if (rank==0) timeRef=MPI_Wtime();
	
	for (unsigned int c=0;c<cellCount;++c) {
		// Loop through the faces of the current cell
		for (unsigned int cf=0;cf<cell[c].faceCount;++cf) {
			Face tempFace;
			unsigned int *tempNodes;
			switch (elemType) {
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
				switch (elemType) {
					case TETRA_4: tempNodes[fn]=cell[c].node(tetraFaces[cf][fn]).id; break;
					case PENTA_6: tempNodes[fn]=cell[c].node(prismFaces[cf][fn]).id; break;
					case HEXA_8: tempNodes[fn]=cell[c].node(hexaFaces[cf][fn]).id; break;
				}
				//tempFace.nodes.push_back(tempNodes[fn]); // FIXME I think this is unnecessary
			}
			// Find the neighbor cell
			bool internal=false;
			bool unique=true;
			for (unsigned int nc=0;nc<node[tempNodes[0]].cells.size();++nc) {
				unsigned int i=node[tempNodes[0]].cells[nc];
				if (i!=c && cell[i].HaveNodes(tempFace.nodeCount,tempNodes)) {
					if (i>c) {
						tempFace.neighbor=i;
						internal=true;
					} else {
						unique=false;
					}
				}
			}
			if (unique) {
				if (!internal) {
					tempFace.bc=0;
					++boundaryFaceCount;
				}
				for (unsigned int fn=0;fn<tempFace.nodeCount;++fn) tempFace.nodes.push_back(tempNodes[fn]);
				face.push_back(tempFace);
				cell[c].faces.push_back(tempFace.id);
				if (internal) cell[tempFace.neighbor].faces.push_back(tempFace.id);
				++faceCount;
			}
			delete [] tempNodes;
		} //for face cf
	} // for cells c

	cout << "* Processor " << rank << " has " << faceCount << " faces, " << boundaryFaceCount << " are on boundaries" << endl;
	
	if (rank==0) {
		timeEnd=MPI_Wtime();
		cout << "* Processor 0 took " << timeEnd-timeRef << " sec to find its faces" << endl;
	}
	
	// Determine and mark faces adjacent to other partitions
	// Create ghost elemets to hold the data from other partitions
	ghostCount=0;
	
	if (np!=1) {
		int counter=0;
		int cellCountOffset[np];
		
		for (unsigned int p=0;p<np;++p) {
			cellCountOffset[p]=counter;
			counter+=otherCellCounts[p];
		}
		
		// Now find the metis2global mapping
		index=0;
		int metis2global[globalCellCount];
		int counter2[np];
		for (unsigned int p=0;p<np;++p) counter2[p]=0;
		for (unsigned int c=0;c<globalCellCount;++c) {
			metis2global[cellCountOffset[cellMap[c]]+counter2[cellMap[c]]]=c;
			counter2[cellMap[c]]++;
		}

		bool foundFlag[globalCellCount];
		for (unsigned int c=0; c<globalCellCount; ++c) foundFlag[c]=false;

		unsigned int parent, metisIndex, gg, matchCount;
		map<unsigned int,unsigned int> ghostGlobal2local;
		for (unsigned int f=0; f<faceCount; ++f) {
			if (face[f].bc>=0) { // if not an internal face or if not found before as partition boundary
				parent=face[f].parent;
				for (unsigned int adjCount=0;adjCount<(xadj[parent+1]-xadj[parent]);++adjCount)  {
					metisIndex=adjncy[xadj[parent]+adjCount];
					gg=metis2global[metisIndex];

						if (metisIndex<cellCountOffset[rank] || metisIndex>=(cellCount+cellCountOffset[rank])) {
							matchCount=0;
							for (unsigned int fn=0;fn<face[f].nodeCount;++fn) {
								for (unsigned int gn=0;gn<elemNodeCount;++gn) {
									if ((elemNodes[gg][gn]-1)==face[f].node(fn).globalId) ++matchCount;
								}
							}
							if (matchCount==face[f].nodeCount) {
								if (!foundFlag[gg]) {
									foundFlag[gg]=true;
									Ghost temp;
									temp.globalId=gg;
									temp.partition=cellMap[gg];
									ghost.push_back(temp);
									face[f].bc=-1*ghostCount-2;
									ghostGlobal2local.insert(pair<unsigned int,unsigned int>(gg,ghostCount));
									++ghostCount;
								}
								face[f].bc=-1*ghostGlobal2local[gg]-2;
								break;
							}
						}
					
				}
			}
		}
				
//  				for (unsigned int c=0; c<globalCellCount; ++c) { //TODO mesh2dual could narrow this search significantly
//  					if (cellMap[c]!=rank && !foundFlag[c]) {
//  						unsigned int matchCount=0;
//  						for (unsigned int fn=0;fn<face[f].nodeCount;++fn) {
//  							for (unsigned int cn=0;cn<elemNodeCount;++cn) {
//  								if ((elemNodes[c][cn]-1) ==face[f].node(fn).globalId) ++matchCount;
//  							}
//  						}
//  						if (matchCount==face[f].nodeCount) {
//  							foundFlag[c]=1;
//  							Ghost temp;
//  							temp.partition=cellMap[c];
//  							temp.globalId=c;
//  							ghost.push_back(temp);
//  							face[f].bc=-1*ghostCount-2;
//  							++ghostCount;
//  							break;
//  						}
//  					}
//  				}

	}
	
	cout << "* Processor " << rank << " has " << ghostCount << " ghost cells at partition boundaries" << endl;
	
	/*
		for (section=1; section<sectionCount; ++section) {
			cg_section_read(fileIndex,baseIndex,zoneIndex,section+1,sectionName,&elemType,&elemStart,&elemEnd,&nBndCells,&parentFlag);
			switch (elemType)  {
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
			}
			int elemNodes[elemEnd-elemStart+1][elemNodeCount];
			cg_elements_read(fileIndex,baseIndex,zoneIndex,section+1,*elemNodes,0);

			//cout << "* Began applying " << sectionName << " boundary conditions" << endl;
			// Now all the cells and faces are constructed, read the boundary conditions
			int elemCount=elemEnd-elemStart+1;
			int mark[elemCount];
			for (unsigned int e=0;e<elemCount;++e) {
				for (unsigned int n=0;n<elemNodeCount;++n) --elemNodes[e][n];
				mark[e]=0;
			}
			int flag, count=0.;
			unsigned int f,e,n,n2;
			for (f=0; f<faceCount; ++f) {
				if (face[f].bc==-2 && elemNodeCount==face[f].nodeCount) {// meaning a boundary face of unknown type
					for (e=0;e<elemCount;++e) {
						if (mark[e]==0) {
							for (n=0;n<elemNodeCount;++n) {
								flag=0;
								for (n2=0;n2<elemNodeCount;++n2) {
									if (elemNodes[e][n]==face[f].nodes[n2]) {
										flag=1; break;
									}
								}
								if (flag==0) break;
							}
							if (flag==1) {
								mark[e]=1;
								face[f].bc=section-1;
								++count;
							}
						}
					}
				}
			}
			cout << "* Boundary condition section found and applied: " << sectionName << "\t" << count << endl;
			if (count!=elemCount) cout << "!!! Something is terribly wrong here !!!" << endl;
		} // end loop over sections
	*/

// Now loop through faces and calculate centroids and areas
	for (unsigned int f=0;f<faceCount;++f) {
		Vec3D centroid=0.;
		Vec3D areaVec=0.;
		for (unsigned int n=0;n<face[f].nodeCount;++n) {
			centroid+=face[f].node(n);
		}
		centroid/=face[f].nodeCount;
		face[f].centroid=centroid;
		for (unsigned int n=0;n<face[f].nodeCount-1;++n) {
			areaVec+=0.5* (face[f].node(n)-centroid).cross(face[f].node(n+1)-centroid);
		}
		areaVec+=0.5* (face[f].node(face[f].nodeCount-1)-centroid).cross(face[f].node(0)-centroid);
		if (areaVec.dot(centroid-cell[face[f].parent].centroid) <0.) {
			// [TBM] Need to swap the face and reflect the area vector
			cout << "face " << f << " should be swapped" << endl;
		}
		face[f].area=fabs(areaVec);
		face[f].normal=areaVec/face[f].area;
	}

// Loop through the cells and calculate the volumes
	double totalVolume=0.;
	for (unsigned int c=0;c<cellCount;++c) {
		double volume=0.;
		for (unsigned int f=0;f<cell[c].faceCount;++f) {
			// [TBM] Area of a prism. Good for HEXA cells only
			volume+=1./3.*cell[c].face(f).area*fabs(cell[c].face(f).normal.dot(cell[c].face(f).centroid-cell[c].centroid));
		}
		cell[c].volume=volume;
		totalVolume+=volume;
	}
	cout << "* Processor " << rank << " Total Volume: " << totalVolume << endl;
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
		cerr << "[!! proc " << rank << " ] Number of nodes of the cell must be specified before allocation" << endl;
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

void Grid::gradients(void) {

	// Calculate cell gradients

	// Initialize all gradients to zero
	for (unsigned int c = 0;c<grid.cellCount;++c) {
		for (unsigned int i=0;i<5;++i) cell[c].grad[i]=0.;
	}

	int parent,neighbor;
	double average[5],factor;

	// Loop faces
	Vec3D areaVec, averageVel;
	for (unsigned int f=0;f<faceCount;++f) {
		parent=face[f].parent; neighbor=face[f].neighbor;
		areaVec=face[f].normal*face[f].area;
		for (unsigned int i=0;i<5;++i) average[i]=0.;
		for (unsigned int i=0;i<face[f].cellContributions.indices.size();++i) {
			factor=face[f].cellContributions.data[i];
			average[0]+=cell[face[f].cellContributions.indices[i]].rho*factor;
			average[1]+=cell[face[f].cellContributions.indices[i]].v.comp[0]*factor;
			average[2]+=cell[face[f].cellContributions.indices[i]].v.comp[1]*factor;
			average[3]+=cell[face[f].cellContributions.indices[i]].v.comp[2]*factor;
			average[4]+=cell[face[f].cellContributions.indices[i]].p*factor;
		}
		if (face[f].bc==-1) { // Internal face
			for (unsigned int i=0;i<5;++i) cell[neighbor].grad[i]-=average[i]*areaVec/cell[neighbor].volume;
		} else {	// Boundary face
			if (bc.region[face[f].bc].type=="slip") {
				for (unsigned int i=0;i<3;++i) averageVel.comp[i]=average[i+1];
				averageVel-=averageVel.dot(face[f].normal) *face[f].normal;
				for (unsigned int i=0;i<3;++i) average[i+1]=averageVel.comp[i];
				//average[0]=grid.cell[parent].rho;
				//average[4]=grid.cell[parent].p;
			} else if (bc.region[face[f].bc].type=="noslip") {
				for (unsigned int i=1;i<4;++i) average[i]=0.;
			} else if (bc.region[face[f].bc].type=="inlet") {
				average[0]=bc.region[face[f].bc].rho;
				average[1]=bc.region[face[f].bc].v.comp[0];
				average[2]=bc.region[face[f].bc].v.comp[1];
				average[3]=bc.region[face[f].bc].v.comp[2];
				average[4]=bc.region[face[f].bc].p;
			} else if (bc.region[face[f].bc].type=="outlet") {
				/*average[0]=grid.cell[parent].rho;
				average[1]=grid.cell[parent].v.comp[0];
				average[2]=grid.cell[parent].v.comp[1];
				average[3]=grid.cell[parent].v.comp[2];
				average[4]=grid.cell[parent].p; */
			}
		}

		cell[parent].grad[0]+=average[0]*areaVec/cell[parent].volume;
		cell[parent].grad[4]+=average[4]*areaVec/cell[parent].volume;
		for (unsigned int i=1;i<4;++i) cell[parent].grad[i]+=average[i]*areaVec/cell[parent].volume;
	}
}
