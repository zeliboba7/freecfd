#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
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

Grid::Grid () {
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
    unsigned int offset=rank*cellCount;
    if (rank==np-1) cellCount=cellCount+globalCellCount-np*cellCount;

    int nodeStart[3],nodeEnd[3];
    nodeStart[0]=nodeStart[1]=nodeStart[2]=1;
    nodeEnd[0]=nodeEnd[1]=nodeEnd[2]=nodeCount;

    double x[nodeCount],y[nodeCount],z[nodeCount];

    cg_coord_read(fileIndex,baseIndex,zoneIndex,"CoordinateX",RealDouble,nodeStart,nodeEnd,&x);
    cg_coord_read(fileIndex,baseIndex,zoneIndex,"CoordinateY",RealDouble,nodeStart,nodeEnd,&y);
    cg_coord_read(fileIndex,baseIndex,zoneIndex,"CoordinateZ",RealDouble,nodeStart,nodeEnd,&z);

    //cout << "* Read coordinates"  << endl;
    //node.reserve(nodeCount);
    //face.reserve(cellCount);
    cell.reserve(cellCount);

    //cout << "* Reserved memory for the nodes"  << endl;

    /*
    for (unsigned int i=0;i<nodeCount;++i) {
        Node temp;
        temp.id=i;
        temp.comp[0]=x[i];
        temp.comp[1]=y[i];
        temp.comp[2]=z[i];
        node.push_back(temp);
    }

    cout << "* Nodes created"  << endl;
    */    

    // Determine the number of sections in the zone
    int sectionCount;
    cg_nsections(fileIndex,baseIndex,zoneIndex, &sectionCount);
    
    if (rank==0) cout << "* Number of sections found in zone " << zoneName << ": " << sectionCount << endl;

    ElementType_t elemType;
    int elemNodeCount;
    int elemStart,elemEnd,nBndCells,parentFlag;
    
    // Set face connectivity lists
    int hexaFaces[6][4]=
        {
            {0,3,2,1},
            {4,5,6,7},
            {1,2,6,5},
            {0,4,7,3},
            {1,5,4,0},
            {2,3,7,6}
        };
	
    int prismFaces[5][4]=
        {
            {0,2,1,0},
            {3,4,5,0},
            {0,3,5,2},
            {1,2,5,4},
            {0,1,4,3},
        };	

    int tetraFaces[4][3]=
        {
            {0,2,1},
            {1,2,3},
            {0,3,2},
            {0,1,3}
        };


	unsigned int section=0;
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

	// Create the nodes for each partition
	bool nodeFound[globalNodeCount];
	unsigned int nodeMap[globalNodeCount];
	for (unsigned int n=1;n<=globalNodeCount;++n) nodeFound[n]=false;
	nodeCount=0;
	for (unsigned int c=0;c<cellCount;++c) {
		for (unsigned int n=0;n<elemNodeCount;++n) {
			if (!nodeFound[elemNodes[c+offset][n]]) {
				++nodeCount;
        			Node temp;
        			temp.id=nodeCount-1;
				temp.globalId=elemNodes[c+offset][n]-1;
        			temp.comp[0]=x[elemNodes[c+offset][n]];
        			temp.comp[1]=y[elemNodes[c+offset][n]];
	        		temp.comp[2]=z[elemNodes[c+offset][n]];
        			node.push_back(temp);
				nodeFound[elemNodes[c+offset][n]]=true;
				nodeMap[elemNodes[c+offset][n]]=temp.id;
			}
		}
	}

	// Create the cells
	for (unsigned int c=0;c<cellCount;++c) {
		Cell temp;
		for (unsigned int n=0;n<elemNodeCount;++n) elemNodes[c+offset][n]=nodeMap[elemNodes[c+offset][n]];
		temp.Construct(elemType,elemNodes[c+offset]);
		temp.globalId=c+offset;
		cell.push_back(temp);
	}

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
	ncon- 0 ( # of constraints)
	ncommonnodes- 4 ( we can probably put this to 3)
	nparts- # of processors (Note: BE CAREFUL if != to # of proc)
	tpwgts- null
	ubvec- null (balancing constraints,if needed 1.05 is a good value)
	options- [0 1 15] for default
	edgecut- output, # of edges cut (measure of communication)
	part- output, where our elements should be
	comm- most likely MPI_COMM_WORLD
	*/

	idxtype elmdist[np + 1]; 
	idxtype eptr[cellCount+1];

	int count=0; for (unsigned int c=0; c<cellCount; ++c) count+=cell[c].nodeCount;
	idxtype  eind[count];

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
	for (unsigned int c=1; c<=cellCount;++c) eptr[c]=eptr[c-1]+cell[c-1].nodeCount;

	int index=0;
	for (unsigned int c=0; c<cellCount; ++c) {
		for (unsigned int nc=0; nc<cell[c].nodeCount; ++nc) {
		eind[index]=cell[c].node(nc).globalId;
		++index;
		}
	}

	ompi_communicator_t* commWorld=MPI_COMM_WORLD; 

	 ParMETIS_V3_PartMeshKway(elmdist,eptr,eind, elmwgt,
				 &wgtflag, &numflag, &ncon, &ncommonnodes,
				 &np, tpwgts, &ubvec, options, &edgecut,
				 part,&commWorld) ;

	//XXX Delete the stuff we created that is no longer needed by ParMetis
	// delete[] somestuffwedontneed;

/*
	// Construct the list of cells for each node
	int flag;
	for (unsigned int c=0;c<cellCount;++c) {
		for (unsigned int n=0;n<cell[c].nodeCount;++n) {
			flag=0;
			for (unsigned int i=0;i<cell[c].node(n).cells.size();++i) {
				if (cell[c].node(n).cells[i]==int(c)) {
					flag=1;
					break;
				}
			}
			if (!flag) {
				node[cell[c].nodes[n]].cells.push_back(c);
			}
		}
	}
	// Search and construct faces
	faceCount=0;
	unsigned int *tempNodes;
	int boundaryFaceCount=0;
	// Loop through all the cells

	for (unsigned int i=0;i<cellCount;++i) {
		// Loop through the faces of the current cell
		for (unsigned int f=0;f<cell[i].faceCount;++f) {
			Face tempFace;
			switch (elemType)  {						
				case TETRA_4:
					tempFace.nodeCount=3;
					tempNodes= new unsigned int[3];
					break;
				case PENTA_6:
					if (f<2) { 
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
			tempFace.parent=i;
			// Assign boundary type as internal by default, will be overwritten later
			tempFace.bc=-1;
			// Store the nodes of the current face
			if (faceCount==face.capacity()) face.reserve(int(face.size()*0.10)+face.size()) ;
			for (unsigned int n=0;n<tempFace.nodeCount;++n) {
				switch (elemType)  {
					case TETRA_4: tempNodes[n]=cell[i].node(tetraFaces[f][n]).id; break;
					case PENTA_6: tempNodes[n]=cell[i].node(prismFaces[f][n]).id; break;
					case HEXA_8: tempNodes[n]=cell[i].node(hexaFaces[f][n]).id; break;
				}
				tempFace.nodes.push_back(tempNodes[n]);
			}
			// Find the neighbor cell
			int flagInternal=0;
			int flagExists=1;
			for (unsigned int j=0;j<node[tempNodes[0]].cells.size();++j) {
				int c=node[tempNodes[0]].cells[j];
				if (int(i)!=c && cell[c].HaveNodes(tempFace.nodeCount,tempNodes)) {
					tempFace.neighbor=c;
					flagInternal=1;
					if (int(i)<c) flagExists=0;
					break;
				}
			}
			if (!flagInternal) {
				++boundaryFaceCount;
			} else if (!flagExists) {
				tempFace.bc=-1;
				face.push_back(tempFace);
				for (unsigned int fn=0;fn<tempFace.nodeCount;++fn) face[tempFace.id].nodes.push_back(tempNodes[fn]);
				cell[i].faces.push_back(tempFace.id);
				cell[tempFace.neighbor].faces.push_back(tempFace.id);
				++faceCount;
			}
		} //for face
	} // for cells
	//cout << "* Number of Faces: " << faceCount << endl;
	//cout << "* Number of Faces at Boundaries: " << boundaryFaceCount << endl;


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
            areaVec+=0.5*(face[f].node(n)-centroid).cross(face[f].node(n+1)-centroid);
        }
	areaVec+=0.5*(face[f].node(face[f].nodeCount-1)-centroid).cross(face[f].node(0)-centroid);
        if (areaVec.dot(centroid-cell[face[f].parent].centroid)<0.) {
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
    cout << "* Total Volume: " << totalVolume << endl;
            volume+=1./3.*cell[c].face(f).area*fabs(cell[c].face(f).normal.dot(cell[c].face(f).centroid-cell[c].centroid));
        }
        cell[c].volume=volume;
        totalVolume+=volume;
    }
    cout << "* Total Volume: " << totalVolume << endl;
    return 0;
*/
}


Node::Node (double x, double y, double z) {
    comp[0]=x;
    comp[1]=y;
    comp[2]=z;
}

Cell::Cell(void) {
    ;
}

int Cell::Construct(const ElementType_t elemType, const int nodeList[]) {

	switch (elemType)  {
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

int Cell::HaveNodes(unsigned int &nodelistsize, unsigned int nodelist []) {

    unsigned int matchCount=0;
    for (unsigned int i=0;i<nodelistsize;++i) {
        for (unsigned int j=0;j<nodeCount;++j) {
            if (nodelist[i]==node(j).id) ++matchCount;
	    if (matchCount==nodelistsize) return 1;		    
        }
    }

    return 0;
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
				averageVel-=averageVel.dot(face[f].normal)*face[f].normal;
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
