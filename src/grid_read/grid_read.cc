#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include<cmath>
using namespace std;
#include <cgnslib.h>


int main() {

	int fileIndex,baseIndex,zoneIndex,sectionIndex,nBases,nZones,nSections,nBocos;
	char zoneName[20],sectionName[20]; //baseName[20]
	//  int nBases,cellDim,physDim;
	int size[3];
	int globalCellCount=0;

	string fileName;
	fileName="input.cgns";
	// Open the grid file for reading
	cg_open(fileName.c_str(),MODE_READ,&fileIndex);

	// Read number of bases
	cg_nbases(fileIndex,&nBases);

	cout << "Number of Bases= " << nBases << endl;

	//for (int baseIndex=1;baseIndex<=nBases;++baseIndex) {
	// For now assuming there is only one base as I don't know in which cases there would be more
	baseIndex=1;
		// Read number of zones
		cg_nzones(fileIndex,baseIndex,&nZones);
		int zoneNodeCount[nZones],zoneCellCount[nZones];
		cout << "Number of Zones= " << nZones << endl;
		std::vector<int> cells,cellConnectivity;
		std::vector<double> coordX[nZones],coordY[nZones],coordZ[nZones];
		std::vector<int> zoneCoordMap[nZones];
		int uniqNodeCount=0;

		// Get total number of boundary conditions
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
			cout << "In Zone " << zoneName << endl;
			cout << "...Number of Nodes= " << size[0] << endl;
			cout << "...Number of Cells= " << size[1] << endl;
			cout << "...Number of Sections= " << nSections << endl;
			cout << "...Number of Boundary Conditions= " << nBocos << endl;
			
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
				zoneCoordMap[zoneIndex-1].push_back(-1);
			}

			// In case there are multiple connected zones, collapse the repeated nodes and fix the node numbering
			if (zoneIndex==1) {
				for (int c=0;c<coordX[0].size();++c) {
					uniqNodeCount++;
					zoneCoordMap[0][c]=uniqNodeCount;
				}
			}
			for (int z=0;z<=zoneIndex-2;++z) {
				for (int c=0;c<coordX[zoneIndex-1].size();++c) {
					bool foundFlag=false;
					for (int c2=0;c2<coordX[z].size();++c2) {
						if (fabs(coordX[zoneIndex-1][c]-coordX[z][c2])<1.e-7 && fabs(coordY[zoneIndex-1][c]-coordY[z][c2])<1.e-7 && fabs(coordZ[zoneIndex-1][c]-coordZ[z][c2])<1.e-7) {
							zoneCoordMap[zoneIndex-1][c]=zoneCoordMap[z][c2];
							foundFlag=true;
						}
					}
					if (!foundFlag) {
						uniqNodeCount++;
						zoneCoordMap[zoneIndex-1][c]=uniqNodeCount;
					}
				}
			}

			int bc_range[nBocos][2];
 			for (int bocoIndex=1;bocoIndex<=nBocos;++bocoIndex) {
 				int dummy;
 				cg_boco_read(fileIndex,baseIndex,zoneIndex,bocoIndex,bc_range[bocoIndex-1],&dummy);
  			} // for boco
			
			// Loop sections within the zone
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
					cout << "   ...Found Volume Section " << sectionName << endl;
					// elements array serves as an start index for connectivity list elemConnectivity
					for (int elem=0;elem<=(elemEnd-elemStart);++elem) {
						cells.push_back(cellConnectivity.size());
						for (int n=0;n<elemNodeCount;++n) cellConnectivity.push_back(zoneCoordMap[zoneIndex-1][elemNodes[elem][n]-1]);
					}
					globalCellCount+=(elemEnd-elemStart+1);
				} else {
					bool bcFlag=false;
					for (int nbc=0;nbc<nBocos;++nbc) {
						if (elemStart==bc_range[nbc][0] && elemEnd==bc_range[nbc][1]) {
							bcFlag=true;
							break;
						}
					}
					if (bcFlag) {
						cout << "   ...Found BC Section " << sectionName << endl;
						for (int elem=0;elem<=(elemEnd-elemStart);++elem) {
							bocos[bocoCount].push_back(bocoConnectivity[bocoCount].size());
							for (int n=0;n<elemNodeCount;++n) bocoConnectivity[bocoCount].push_back(zoneCoordMap[zoneIndex-1][elemNodes[elem][n]-1]);
						}
						bocoCount+=1;
					}
				}// if
			} // for section
			cout << "Total Node Count= " << uniqNodeCount << endl;
		} // for zone
	//} // for base


		fstream file;

		file.open("tecplot.dat", ios::out);

		file << "VARIABLES = \"x\", \"y\", \"z\" " << endl;
		file << "ZONE T=\"VOLUME MESH\"  N=" << uniqNodeCount << "  E=" << globalCellCount << endl;
		file << "DATAPACKING=BLOCK, ZONETYPE=FEBRICK"<< endl;

		for (int zone=0;zone<nZones;++zone) {
			for (int n=0;n<coordX[zone].size();++n) {
				if (zone==0 | zoneCoordMap[zone][n]>1200) {
					file << setw(16) << setprecision(8) << scientific << coordX[zone][n] << endl;
				}
			}
		}
		for (int zone=0;zone<nZones;++zone) {
			for (int n=0;n<coordX[zone].size();++n) {
				if (zone==0 | zoneCoordMap[zone][n]>1200) {
					file << setw(16) << setprecision(8) << scientific << coordY[zone][n] << endl;
				}
			}
		}
		for (int zone=0;zone<nZones;++zone) {
			for (int n=0;n<coordX[zone].size();++n) {
				if (zone==0 | zoneCoordMap[zone][n]>1200) {
					file << setw(16) << setprecision(8) << scientific << coordZ[zone][n] << endl;
				}
			}
		}
		int ni,nn;
		for (int c=0;c<globalCellCount;++c) {
			ni=cells[c];
			if (c==globalCellCount-1) {
				nn=cellConnectivity.size();
			} else {
				nn=cells[c+1];
			}
			if ((nn-ni)==8) { // hexa
				for (int n=0;n<8;++n) {
					ni=cells[c]+n;
					file << cellConnectivity[ni] << "  ";
				}
				file << endl;
			}
			else if ((nn-ni)==6) { //prism
				for (int n=0;n<=2;++n) {
					ni=cells[c]+n;
					file << cellConnectivity[ni] << "  ";
				}
				ni=cells[c]+2;
				file << cellConnectivity[ni] << "  ";
				for (int n=3;n<=5;++n) {
					ni=cells[c]+n;
					file << cellConnectivity[ni] << "  ";
				}
				ni=cells[c]+5;
				file << cellConnectivity[ni] << "  ";
				file << endl;
			}
		}

		// Write Boundary conditions
		for (int bc=0;bc<totalnBocos;++bc) {
			file << "ZONE T=\"BC_"<< bc+1 << "\", N=" << uniqNodeCount <<   " , E=" << bocos[bc].size() <<  ", ZONETYPE=";
			if ((bocos[bc][1]-bocos[bc][0])==4) {
				file << "FEQUADRILATERAL, VARSHARELIST = ([1,2,3]=1)"<< endl;
				for (int c=0;c<bocos[bc].size();++c) {
					ni=bocos[bc][c];
					if (c==bocos[bc].size()-1) {
						nn=bocoConnectivity[bc].size();
					} else {
						nn=bocos[bc][c+1];
					}
					for (int n=0;n<4;++n) {
						ni=bocos[bc][c]+n;
						file << bocoConnectivity[bc][ni] << "  ";
					}
				file << endl;
				}
			}
			if ((bocos[bc][1]-bocos[bc][0])==3) {
				file << "FETRIANGLE, VARSHARELIST = ([1,2,3]=1)"<< endl;
				for (int c=0;c<bocos[bc].size();++c) {
					ni=bocos[bc][c];
					if (c==bocos[bc].size()-1) {
						nn=bocoConnectivity[bc].size();
					} else {
						nn=bocos[bc][c+1];
					}
					for (int n=0;n<3;++n) {
						ni=bocos[bc][c]+n;
						file << bocoConnectivity[bc][ni] << "  ";
					}
					file << endl;
				}
			}
		}
		file.close();




// 
// 	Initialize the partition sizes
// 	cellCount=floor(globalCellCount/np);
// 	int baseCellCount=cellCount;
// 	unsigned int offset=rank*cellCount;
// 	if (rank==np-1) cellCount=cellCount+globalCellCount-np*cellCount;
// 
// 	int nodeStart[3],nodeEnd[3];
// 	nodeStart[0]=nodeStart[1]=nodeStart[2]=1;
// 	nodeEnd[0]=nodeEnd[1]=nodeEnd[2]=globalNodeCount;
// 
// 	double x[globalNodeCount],y[globalNodeCount],z[globalNodeCount];
// 
// 	cg_coord_read(fileIndex,baseIndex,zoneIndex,"CoordinateX",RealDouble,nodeStart,nodeEnd,&x);
// 	cg_coord_read(fileIndex,baseIndex,zoneIndex,"CoordinateY",RealDouble,nodeStart,nodeEnd,&y);
// 	cg_coord_read(fileIndex,baseIndex,zoneIndex,"CoordinateZ",RealDouble,nodeStart,nodeEnd,&z);
// 
// 	Determine the number of sections in the zone
// 	int sectionCount;
// 	cg_nsections(fileIndex,baseIndex,zoneIndex, &sectionCount);
// 
// 	if (rank==0) cout << "* Number of sections found in zone " << zoneName << ": " << sectionCount << endl;
// 
// 	ElementType_t elemType;
// 	int elemNodeCount;
// 	int elemStart,elemEnd,nBndCells,parentFlag;
// 
// 	unsigned int section=0;
// 	cg_section_read(fileIndex,baseIndex,zoneIndex,section+1,sectionName,&elemType,&elemStart,&elemEnd,&nBndCells,&parentFlag);
// 	switch (elemType) {
// 		case TRI_3:
// 			elemNodeCount=3; break;
// 		case QUAD_4:
// 			elemNodeCount=4; break;
// 		case TETRA_4:
// 			elemNodeCount=4; break;
// 		case PENTA_6:
// 			elemNodeCount=6; break;
// 		case HEXA_8:
// 			elemNodeCount=8; break;
// 	}
// 	int elemNodes[elemEnd-elemStart+1][elemNodeCount];
// 	cg_elements_read(fileIndex,baseIndex,zoneIndex,section+1,*elemNodes,0);
// 
// 	Implementing Parmetis
// 	/* ParMETIS_V3_PartMeshKway(idxtype *elmdist, idxtype *eptr, idxtype *eind, idxtype *elmwgt, int *wgtflag, int *numflag, int *ncon, int * ncommonnodes, int *nparts, float *tpwgts, float *ubvec, int *options, int *edgecut, idxtype *part, MPI_Comm) */
// 
// 	/*  Definining variables
// 	elmdist- look into making the 5 arrays short int (for performance
// 		on 64 bit arch)
// 	eptr- like xadf
// 	eind- like adjncy
// 	elmwgt- null (element weights)
// 	wgtflag- 0 (no weights, can take value of 0,1,2,3 see documentation)
// 	numflag- 0 C-style numbers, 1 Fortran-style numbers
// 	ncon- 1  ( # of constraints)
// 	ncommonnodes- 4 ( we can probably put this to 3)
// 	nparts- # of processors (Note: BE CAREFUL if != to # of proc)
// 	tpwgts- 
// 	ubvec-  (balancing constraints,if needed 1.05 is a good value)
// 	options- [0 1 15] for default
// 	edgecut- output, # of edges cut (measure of communication)
// 	part- output, where our elements should be
// 	comm- most likely MPI_COMM_WORLD
// 	*/
// 
// 	idxtype elmdist[np + 1];
// 	idxtype *eptr;
// 	eptr = new idxtype[cellCount+1];
// 	idxtype *eind;
// 	eind = new idxtype[cellCount*elemNodeCount];
// 
// 	idxtype* elmwgt = NULL;
// 	int wgtflag=0; // no weights associated with elem or edges
// 	int numflag=0; // C-style numbering
// 	int ncon=1; // # of weights or constraints
// 	int ncommonnodes; ncommonnodes=3; // set to 3 for tetrahedra or mixed type
// 
// 	float tpwgts[np];
// 	for (unsigned int p=0; p<np; ++p) tpwgts[p]=1./float(np);
// 	float ubvec=1.02;
// 	int options[3]; // default values for timing info set 0 -> 1
// 
// 	options[0]=0; options[1]=1; options[2]=15;
// 	int edgecut ; // output
// 	idxtype* part = new idxtype[cellCount];
// 
// 	for (unsigned int p=0;p<np;++p) elmdist[p]=p*floor(globalCellCount/np);
// 	elmdist[np]=globalCellCount;// Note this is because #elements mod(np) are all on last proc
// 
// 	eptr[0]=0;
// 	for (unsigned int c=1; c<=cellCount;++c) eptr[c]=eptr[c-1]+elemNodeCount;
// 
// 	int index=0;
// 	for (unsigned int c=0; c<cellCount; ++c) {
// 		for (unsigned int nc=0; nc<elemNodeCount; ++nc) {
// 			eind[index]=elemNodes[c+offset][nc]-1;
// 			++index;
// 		}
// 	}
// 
// 	ompi_communicator_t* commWorld=MPI_COMM_WORLD;
// 
// 	ParMETIS_V3_PartMeshKway(elmdist,eptr,eind, elmwgt,
// 	                         &wgtflag, &numflag, &ncon, &ncommonnodes,
// 	                         &np, tpwgts, &ubvec, options, &edgecut,
// 	                         part,&commWorld) ;
// 
// 	delete[] eptr;
// 	delete[] eind;
// 
// 	Distribute the part list to each proc
// 	Each proc has an array of length globalCellCount which says the processor number that cell belongs to [cellMap]
// 	int recvCounts[np];
// 	int displs[np];
// 	for (int p=0;p<np;++p) {
// 		recvCounts[p]=baseCellCount;
// 		displs[p]=p*baseCellCount;
// 	}
// 	recvCounts[np-1]=baseCellCount+globalCellCount-np*baseCellCount;
// 	int cellMap[globalCellCount];
// 	cellMap of a cell returns which processor it is assigned to
// 	MPI_Allgatherv(part,cellCount,MPI_INT,cellMap,recvCounts,displs,MPI_INT,MPI_COMM_WORLD);
// 
// 	Find new local cellCount after ParMetis distribution
// 	cellCount=0.;
// 	int otherCellCounts[np]; 
// 	for (unsigned int p=0;p<np;p++) otherCellCounts[p]=0; 
// 	
// 	for (unsigned int c=0;c<globalCellCount;++c) {
// 		otherCellCounts[cellMap[c]]+=1;
// 		if (cellMap[c]==rank) ++cellCount;
// 	}
// 	cout << "* Processor " << rank << " has " << cellCount << " cells" << endl;
// 		
// 	node.reserve(nodeCount/np);
// 	face.reserve(cellCount);
// 	cell.reserve(cellCount);
// 
// 	Create the nodes and cells for each partition
// 	Run through the list of cells and check if it belongs to current partition
// 	Loop through the cell's nodes
// 	Mark the visited nodes so that no duplicates are created (nodeFound array).
// 	bool nodeFound[globalNodeCount];
// 	unsigned int nodeMap[globalNodeCount];
// 	unsigned int cellNodes[elemNodeCount];
// 	for (unsigned int n=1;n<=globalNodeCount;++n) nodeFound[n]=false;
// 	nodeCount=0;
// 	fstream file;
// 	if (rank==0) file.open("connectivity.dat", ios::out);
// 	for (unsigned int c=0;c<globalCellCount;++c) {
// 		if (rank==0) {
// 			if (elemType==PENTA_6) {
// 				file << elemNodes[c][0] << "\t" ; 
// 				file << elemNodes[c][1] << "\t" ; 
// 				file << elemNodes[c][2] << "\t" ; 
// 				file << elemNodes[c][2] << "\t" ; 
// 				file << elemNodes[c][3] << "\t" ; 
// 				file << elemNodes[c][4] << "\t" ; 
// 				file << elemNodes[c][5] << "\t" ; 	
// 				file << elemNodes[c][5] << "\t" ; 
// 			} else {
// 				for (unsigned int i=0;i<elemNodeCount;++i) {
// 					file << elemNodes[c][i] << "\t";
// 				}
// 			}
// 		}
// 		if (cellMap[c]==rank) {
// 			for (unsigned int n=0;n<elemNodeCount;++n) {
// 				if (!nodeFound[elemNodes[c][n]]) {
// 					Node temp;
// 					temp.id=nodeCount;
// 					temp.globalId=elemNodes[c][n]-1;
// 					temp.comp[0]=x[elemNodes[c][n]-1];
// 					temp.comp[1]=y[elemNodes[c][n]-1];
// 					temp.comp[2]=z[elemNodes[c][n]-1];
// 					node.push_back(temp);
// 					nodeFound[elemNodes[c][n]]=true;
// 					nodeMap[elemNodes[c][n]]=temp.id;
// 					++nodeCount;
// 				}
// 				cellNodes[n]=nodeMap[elemNodes[c][n]];
// 			}
// 			Cell temp;
// 			temp.Construct(elemType,cellNodes);
// 			temp.globalId=c;
// 			cell.push_back(temp);
// 		}
// 		if (rank==0) file << "\n" ;
// 	}
// 	if (rank==0) file.close();
// 	
// 	cout << "* Processor " << rank << " has created its cells and nodes" << endl;
// 
// 	Create the Mesh2Dual inputs
// 	idxtype elmdist [np+1] (stays the same size)
// 	eptr = new idxtype[cellCount+1];
// 	eind = new idxtype[cellCount*elemNodeCount];
// 	numflag and ncommonnodes previously defined
// 	ncommonnodes=1;
// 	idxtype* xadj;
// 	idxtype* adjncy;
// 	
// 	elmdist[0]=0;
// 	for (unsigned int p=1;p<=np;p++) elmdist[p]=otherCellCounts[p-1]+elmdist[p-1];
// 	eptr[0]=0;
// 	for (unsigned int c=1; c<=cellCount;++c) eptr[c]=eptr[c-1]+elemNodeCount;
// 	index=0;
// 	for (unsigned int c=0; c<cellCount;c++){
// 		for (unsigned int cn=0; cn<cell[c].nodeCount; ++cn) {
// 			eind[index]=cell[c].node(cn).globalId;
// 			++index;
// 		}	
// 	}
// 	
// 	ParMETIS_V3_Mesh2Dual(elmdist, eptr, eind, &numflag, &ncommonnodes, &xadj, &adjncy, &commWorld);
// 	
// 	Construct the list of cells for each node
// 	bool flag;
// 	for (unsigned int c=0;c<cellCount;++c) {
// 		unsigned int n;
// 		for (unsigned int cn=0;cn<cell[c].nodeCount;++cn) {
// 			n=cell[c].nodes[cn];
// 			flag=false;
// 			for (unsigned int i=0;i<node[n].cells.size();++i) {
// 				if (node[n].cells[i]==c) {
// 					flag=true;
// 					break;
// 				}
// 			}
// 			if (!flag) {
// 				node[n].cells.push_back(c);
// 			}
// 		}
// 	}
// 
// 	cout << "* Processor " << rank << " has computed its node-cell connectivity" << endl;
// 	
// 	Construct the list of neighboring cells for each cell
// 	int c2;
// 	for (unsigned int c=0;c<cellCount;++c) {
// 		unsigned int n;
// 		for (unsigned int cn=0;cn<cell[c].nodeCount;++cn) {
// 			n=cell[c].nodes[cn];
// 			for (unsigned int nc=0;nc<node[n].cells.size();++nc) {
// 				c2=node[n].cells[nc];
// 				flag=false;
// 				for (unsigned int cc=0;cc<cell[c].neighborCells.size();++cc) {
// 					if(cell[c].neighborCells[cc]==c2) {
// 						flag=true;
// 						break;
// 					}
// 				}
// 				if (!flag) cell[c].neighborCells.push_back(c2);
// 			} // end node cell loop
// 		} // end cell node loop
// 		cell[c].neighborCellCount=cell[c].neighborCells.size();
// 	} // end cell loop
// 	
// 	Set face connectivity lists
// 	int hexaFaces[6][4]= {
// 		{0,3,2,1},
// 		{4,5,6,7},
// 		{1,2,6,5},
// 		{0,4,7,3},
// 		{1,5,4,0},
// 		{2,3,7,6}
// 	};
// 
// 	int prismFaces[5][4]= {
// 		{0,2,1,0},
// 		{3,4,5,0},
// 		{0,3,5,2},
// 		{1,2,5,4},
// 		{0,1,4,3},
// 	};
// 
// 	int tetraFaces[4][3]= {
// 		{0,2,1},
// 		{1,2,3},
// 		{0,3,2},
// 		{0,1,3}
// 	};
// 
// 
// 	Search and construct faces
// 	faceCount=0;
// 
// 	int boundaryFaceCount=0;
// 	Loop through all the cells
// 
// 	double timeRef, timeEnd;
// 	if (rank==0) timeRef=MPI_Wtime();
// 	
// 	for (unsigned int c=0;c<cellCount;++c) {
// 		Loop through the faces of the current cell
// 		for (unsigned int cf=0;cf<cell[c].faceCount;++cf) {
// 			Face tempFace;
// 			unsigned int *tempNodes;
// 			switch (elemType) {
// 				case TETRA_4:
// 					tempFace.nodeCount=3;
// 					tempNodes= new unsigned int[3];
// 					break;
// 				case PENTA_6:
// 					if (cf<2) {
// 						tempFace.nodeCount=3;
// 						tempNodes= new unsigned int[3];
// 					} else {
// 						tempFace.nodeCount=4;
// 						tempNodes= new unsigned int[4];
// 					}
// 					break;
// 				case HEXA_8:
// 					tempFace.nodeCount=4;
// 					tempNodes= new unsigned int[4];
// 					break;
// 			}
// 			tempFace.id=faceCount;
// 			Assign current cell as the parent cell
// 			tempFace.parent=c;
// 			Assign boundary type as internal by default, will be overwritten later
// 			tempFace.bc=-1;
// 			Store the nodes of the current face
// 			if (faceCount==face.capacity()) face.reserve(int (face.size() *0.10) +face.size()) ; //TODO check how the size grows by default
// 
// 			for (unsigned int fn=0;fn<tempFace.nodeCount;++fn) {
// 				switch (elemType) {
// 					case TETRA_4: tempNodes[fn]=cell[c].node(tetraFaces[cf][fn]).id; break;
// 					case PENTA_6: tempNodes[fn]=cell[c].node(prismFaces[cf][fn]).id; break;
// 					case HEXA_8: tempNodes[fn]=cell[c].node(hexaFaces[cf][fn]).id; break;
// 				}
// 				tempFace.nodes.push_back(tempNodes[fn]); // FIXME I think this is unnecessary
// 			}
// 			Find the neighbor cell
// 			bool internal=false;
// 			bool unique=true;
// 			for (unsigned int nc=0;nc<node[tempNodes[0]].cells.size();++nc) {
// 				unsigned int i=node[tempNodes[0]].cells[nc];
// 				if (i!=c && cell[i].HaveNodes(tempFace.nodeCount,tempNodes)) {
// 					if (i>c) {
// 						tempFace.neighbor=i;
// 						internal=true;
// 					} else {
// 						unique=false;
// 					}
// 				}
// 			}
// 			if (unique) {
// 				if (!internal) {
// 					tempFace.bc=0;
// 					++boundaryFaceCount;
// 				}
// 				for (unsigned int fn=0;fn<tempFace.nodeCount;++fn) tempFace.nodes.push_back(tempNodes[fn]);
// 				face.push_back(tempFace);
// 				cell[c].faces.push_back(tempFace.id);
// 				if (internal) cell[tempFace.neighbor].faces.push_back(tempFace.id);
// 				++faceCount;
// 			}
// 			delete [] tempNodes;
// 		} //for face cf
// 	} // for cells c
// 
// 	cout << "* Processor " << rank << " has " << faceCount << " faces, " << boundaryFaceCount << " are on boundaries" << endl;
// 	
// 	if (rank==0) {
// 		timeEnd=MPI_Wtime();
// 		cout << "* Processor 0 took " << timeEnd-timeRef << " sec to find its faces" << endl;
// 	}
// 	
// 	Determine and mark faces adjacent to other partitions
// 	Create ghost elemets to hold the data from other partitions
// 	ghostCount=0;
// 	
// 	if (np!=1) {
// 		int counter=0;
// 		int cellCountOffset[np];
// 		
// 		for (unsigned int p=0;p<np;++p) {
// 			cellCountOffset[p]=counter;
// 			counter+=otherCellCounts[p];
// 		}
// 		
// 		Now find the metis2global mapping
// 		index=0;
// 		int metis2global[globalCellCount];
// 		int counter2[np];
// 		for (unsigned int p=0;p<np;++p) counter2[p]=0;
// 		for (unsigned int c=0;c<globalCellCount;++c) {
// 			metis2global[cellCountOffset[cellMap[c]]+counter2[cellMap[c]]]=c;
// 			counter2[cellMap[c]]++;
// 		}
// 
// 		int foundFlag[globalCellCount];
// 		for (unsigned int c=0; c<globalCellCount; ++c) foundFlag[c]=0;
// 
// 		unsigned int parent, metisIndex, gg, matchCount;
// 		map<unsigned int,unsigned int> ghostGlobal2local;
// 		
// 		map<int,set<int> > nodeCellSet;
// 		Vec3D nodeVec;
// 		
// 		for (unsigned int f=0; f<faceCount; ++f) {
// 			if (face[f].bc>=0) { // if not an internal face or if not found before as partition boundary
// 				parent=face[f].parent;
// 				for (unsigned int adjCount=0;adjCount<(xadj[parent+1]-xadj[parent]);++adjCount)  {
// 					metisIndex=adjncy[xadj[parent]+adjCount];
// 					gg=metis2global[metisIndex];
// 
// 					if (metisIndex<cellCountOffset[rank] || metisIndex>=(cellCount+cellCountOffset[rank])) {
// 						matchCount=0;
// 						for (unsigned int fn=0;fn<face[f].nodeCount;++fn) {
// 							set<int> tempSet;
// 							nodeCellSet.insert(pair<unsigned int,set<int> >(face[f].nodes[fn],tempSet) );
// 							for (unsigned int gn=0;gn<elemNodeCount;++gn) {
// 								if ((elemNodes[gg][gn]-1)==face[f].node(fn).globalId) {
// 									nodeCellSet[face[f].nodes[fn]].insert(gg);
// 									++matchCount;
// 								}
// 							}
// 						}
// 						if (matchCount>0 && foundFlag[gg]==0) {
// 							if (matchCount>=3) foundFlag[gg]=3;
// 							if (matchCount<3) foundFlag[gg]=matchCount;
// 							Ghost temp;
// 							temp.globalId=gg;
// 							temp.partition=cellMap[gg];
// 							Calculate the centroid
// 							temp.centroid=0.;
// 							for (unsigned int gn=0;gn<elemNodeCount;++gn) {
// 								nodeVec.comp[0]=x[elemNodes[gg][gn]-1];
// 								nodeVec.comp[1]=y[elemNodes[gg][gn]-1];
// 								nodeVec.comp[2]=z[elemNodes[gg][gn]-1];
// 								temp.centroid+=nodeVec;
// 							}
// 							temp.centroid/=double(elemNodeCount);
// 							ghost.push_back(temp);
// 							ghostGlobal2local.insert(pair<unsigned int,unsigned int>(gg,ghostCount));
// 							++ghostCount;
// 						}
// 						if (matchCount>=3) {
// 							foundFlag[gg]=3;
// 							face[f].bc=-1*ghostGlobal2local[gg]-2;
// 						} 
// 					}
// 				}
// 					
// 			}
// 		}
// 
// 		map<int,set<int> >::iterator mit;
// 		set<int>::iterator sit;
// 		for ( mit=nodeCellSet.begin() ; mit != nodeCellSet.end(); mit++ ) {
// 			for ( sit=(*mit).second.begin() ; sit != (*mit).second.end(); sit++ ) {
// 				node[(*mit).first].ghosts.push_back(ghostGlobal2local[*sit]);
// 				cout << (*mit).first << "\t" << *sit << endl; // DEBUG
// 			}	
// 		}
// 		
// 		Construct the list of neighboring ghosts for each cell
// 		int g;
// 		bool flag;
// 		for (unsigned int c=0;c<cellCount;++c) {
// 			unsigned int n;
// 			for (unsigned int cn=0;cn<cell[c].nodeCount;++cn) {
// 				n=cell[c].nodes[cn];
// 				for (unsigned int ng=0;ng<node[n].ghosts.size();++ng) {
// 					g=node[n].ghosts[ng];
// 					flag=false;
// 					for (unsigned int cg=0;cg<cell[c].ghosts.size();++cg) {
// 						if(cell[c].ghosts[cg]==g) {
// 							flag=true;
// 							break;
// 						}
// 					} // end cell ghost loop
// 					if (flag==false) {
// 						cell[c].ghosts.push_back(g);
// 						ghost[g].cells.push_back(c);
// 					}
// 				} // end node ghost loop
// 			} // end cell node loop
// 			cell[c].ghostCount=cell[c].ghosts.size();
// 		} // end cell loop
// 		
// 		
// 		
// 	} // if (np!=1) 
// 	
// 	if (rank==2) { // DEBUG
// 		for (int n=0;n<nodeCount;++n) {
// 			for (int nc=0;nc<node[n].cells.size();++nc) cout << n << "\t" << "cells\t" << node[n].cells[nc] << endl;
// 			for (int ng=0;ng<node[n].ghosts.size();++ng) cout << n << "\t" << "ghosts\t" << node[n].ghosts[ng] << "\t" << ghost[node[n].ghosts[ng]].globalId <<  endl;
// 		}
// 	}
// 	
// 	
// 	cout << "* Processor " << rank << " has " << ghostCount << " ghost cells at partition boundaries" << endl;
// 	
// 	/*
// 		for (section=1; section<sectionCount; ++section) {
// 			cg_section_read(fileIndex,baseIndex,zoneIndex,section+1,sectionName,&elemType,&elemStart,&elemEnd,&nBndCells,&parentFlag);
// 			switch (elemType)  {
// 				case TRI_3:
// 					elemNodeCount=3; break;
// 				case QUAD_4:
// 					elemNodeCount=4; break;
// 				case TETRA_4:
// 					elemNodeCount=4; break;
// 				case PENTA_6:
// 					elemNodeCount=6; break;
// 				case HEXA_8:
// 					elemNodeCount=8; break;
// 			}
// 			int elemNodes[elemEnd-elemStart+1][elemNodeCount];
// 			cg_elements_read(fileIndex,baseIndex,zoneIndex,section+1,*elemNodes,0);
// 
// 			cout << "* Began applying " << sectionName << " boundary conditions" << endl;
// 			Now all the cells and faces are constructed, read the boundary conditions
// 			int elemCount=elemEnd-elemStart+1;
// 			int mark[elemCount];
// 			for (unsigned int e=0;e<elemCount;++e) {
// 				for (unsigned int n=0;n<elemNodeCount;++n) --elemNodes[e][n];
// 				mark[e]=0;
// 			}
// 			int flag, count=0.;
// 			unsigned int f,e,n,n2;
// 			for (f=0; f<faceCount; ++f) {
// 				if (face[f].bc==-2 && elemNodeCount==face[f].nodeCount) {// meaning a boundary face of unknown type
// 					for (e=0;e<elemCount;++e) {
// 						if (mark[e]==0) {
// 							for (n=0;n<elemNodeCount;++n) {
// 								flag=0;
// 								for (n2=0;n2<elemNodeCount;++n2) {
// 									if (elemNodes[e][n]==face[f].nodes[n2]) {
// 										flag=1; break;
// 									}
// 								}
// 								if (flag==0) break;
// 							}
// 							if (flag==1) {
// 								mark[e]=1;
// 								face[f].bc=section-1;
// 								++count;
// 							}
// 						}
// 					}
// 				}
// 			}
// 			cout << "* Boundary condition section found and applied: " << sectionName << "\t" << count << endl;
// 			if (count!=elemCount) cout << "!!! Something is terribly wrong here !!!" << endl;
// 		} // end loop over sections
// 	*/
// 
// Now loop through faces and calculate centroids and areas
// 	for (unsigned int f=0;f<faceCount;++f) {
// 		Vec3D centroid=0.;
// 		Vec3D areaVec=0.;
// 		for (unsigned int n=0;n<face[f].nodeCount;++n) {
// 			centroid+=face[f].node(n);
// 		}
// 		centroid/=face[f].nodeCount;
// 		face[f].centroid=centroid;
// 		for (unsigned int n=0;n<face[f].nodeCount-1;++n) {
// 			areaVec+=0.5* (face[f].node(n)-centroid).cross(face[f].node(n+1)-centroid);
// 		}
// 		areaVec+=0.5* (face[f].node(face[f].nodeCount-1)-centroid).cross(face[f].node(0)-centroid);
// 		if (areaVec.dot(centroid-cell[face[f].parent].centroid) <0.) {
// 			[TBM] Need to swap the face and reflect the area vector
// 			cout << "face " << f << " should be swapped" << endl;
// 		}
// 		face[f].area=fabs(areaVec);
// 		face[f].normal=areaVec/face[f].area;
// 	}
// 
// Loop through the cells and calculate the volumes and length scales
// 	double totalVolume=0.;
// 	for (unsigned int c=0;c<cellCount;++c) {
// 		cell[c].lengthScale=1.e20;
// 		double volume=0.,height;
// 		unsigned int f;
// 		for (unsigned int cf=0;cf<cell[c].faceCount;++cf) {
// 			FIXME Is this a generic volume formula?
// 			f=cell[c].faces[cf];
// 			height=fabs(face[f].normal.dot(face[f].centroid-cell[c].centroid));
// 			volume+=1./3.*face[f].area*height;
// 			cell[c].lengthScale=min(cell[c].lengthScale,height);
// 		}
// 		cell[c].volume=volume;
// 		totalVolume+=volume;
// 	}
// 	cout << "* Processor " << rank << " Total Volume: " << totalVolume << endl;*/
	return 0;

}

