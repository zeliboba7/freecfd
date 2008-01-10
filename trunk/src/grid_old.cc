#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
#include<cmath>
#include <cgnslib.h>

#include "grid.h"
#include "bc.h"

extern Grid grid;
extern BC bc;

Grid::Grid() {
	;
}

int Grid::read(string fname) {
	fstream file;
	fileName=fname;
	file.open(fileName.c_str());
	if (file.is_open()) {
		cout << "* Found grid file " << fileName  << endl;
		file.close();
		ReadCGNS();
		return 1;
	} else {
		cerr << "[!!] Grid file "<< fileName << " could not be found." << endl;
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
	sectionIndex=1; // Assumed [TBM]

	cg_zone_read(fileIndex,baseIndex,zoneIndex,zoneName,size);

	nodeCount=size[0];
	cellCount=size[1];
	faceCount=0;

	cout << "* Number of nodes: " << nodeCount << endl;
	cout << "* Number of cells: " << cellCount << endl;

	int nodeStart[3],nodeEnd[3];
	nodeStart[0]=nodeStart[1]=nodeStart[2]=1;
	nodeEnd[0]=nodeEnd[1]=nodeEnd[2]=nodeCount;

	double x[nodeCount],y[nodeCount],z[nodeCount];

	cg_coord_read(fileIndex,baseIndex,zoneIndex,"CoordinateX",RealDouble,nodeStart,nodeEnd,&x);
	cg_coord_read(fileIndex,baseIndex,zoneIndex,"CoordinateY",RealDouble,nodeStart,nodeEnd,&y);
	cg_coord_read(fileIndex,baseIndex,zoneIndex,"CoordinateZ",RealDouble,nodeStart,nodeEnd,&z);

	cout << "* Read coordinates"  << endl;

	node.reserve(nodeCount);
	face.reserve(cellCount);
	cell.reserve(cellCount);

	cout << "* Reserved memory for the nodes"  << endl;

	for (unsigned int i=0;i<nodeCount;++i) {
		Node temp;
		temp.id=i;
		temp.comp[0]=x[i];
		temp.comp[1]=y[i];
		temp.comp[2]=z[i];
		node.push_back(temp);
	}

	cout << "* Nodes created"  << endl;

	ElementType_t elemType;
	int elemStart,elemEnd,nBndCells,parentFlag;

	cg_section_read(fileIndex,baseIndex,zoneIndex,sectionIndex,sectionName,&elemType,&elemStart,&elemEnd,&nBndCells,&parentFlag);

	int elemNodeCount;
	if (elemType==17) elemNodeCount=8;   // HEXA
	if (elemType==14) elemNodeCount=6;   // WEDGE
	if (elemType==10) elemNodeCount=4;   // TETRA


	cout << "* Cell node count="<< elemNodeCount  << endl;

	int elemNodes[elemEnd-elemStart+1][elemNodeCount];
	cg_elements_read(fileIndex,baseIndex,zoneIndex,sectionIndex,*elemNodes,0);

	cout << "* Read cell nodes"  << endl;

	if (elemType==17) {
		for (unsigned int c=0;c<cellCount;++c) {
			Cell temp;
			for (unsigned int n=0;n<elemNodeCount;++n) --elemNodes[c][n];
			temp.Construct(HEXA,elemNodeCount,elemNodes[c]);
			cell.push_back(temp);
		}
		cout << "* Constructed HEXA cells"  << endl;
	}

	if (elemType==14) {
		for (unsigned int c=0;c<cellCount;++c) {
			Cell temp;
			for (int n=0;n<elemNodeCount;++n) --elemNodes[c][n];
			temp.Construct(WEDGE,elemNodeCount,elemNodes[c]);
			cell.push_back(temp);
		}
		cout << "* Constructed WEDGE cells"  << endl;
	}

	if (elemType==10) {
		for (unsigned int c=0;c<cellCount;++c) {
			Cell temp;
			for (int n=0;n<elemNodeCount;++n) --elemNodes[c][n];
			temp.Construct(TETRA,elemNodeCount,elemNodes[c]);
			cell.push_back(temp);
		}
		cout << "* Constructed TETRA cells"  << endl;
	}



	for (unsigned int n=0;n<nodeCount;++n) {
		node[n].cells.reserve(8);
	}

	for (unsigned int c=0;c<cellCount;++c) {
		for (unsigned int n=0;n<cell[c].nodeCount;++n) {
			int flag=0;
			for (unsigned int i=0;i<cell[c].node(n).cells.size();++i) {
				if (cell[c].node(n).cells[i]==int (c)) {
					flag=1;
					break;
				}
			}
			if (!flag) {
				node[cell[c].nodes[n]].cells.push_back(c);
			}
		}
	}

	// Find out how many faces are there in the domain

	// Set hexa cell face connectivity list
	int hexaFaces[6][4]= {
		{0,3,2,1},
		{4,5,6,7},
		{1,2,6,5},
		{0,4,7,3},
		{1,5,4,0},
		{2,3,7,6}
	};

	int wedgeFaces[5][4]= {
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

	faceCount=0;
	unsigned int tempNodesSize;
	unsigned int *tempNodes;
	int boundaryFaceCount=0;
	// Loop through all the cells

	for (unsigned int i=0;i<cellCount;++i) {
		if (elemType==17) { // HEXA
			tempNodesSize=4;
			tempNodes= new unsigned int[4];
		}

		if (elemType==10) { // TETRA
			tempNodesSize=3;
			tempNodes= new unsigned int[3];
		}
		// Loop through the faces of the current cell
		for (unsigned int f=0;f<cell[i].faceCount;++f) {
			Face tempFace;
			if (elemType==17) tempFace.nodeCount=4;
			if (elemType==14) {
				if (f<2) {
					tempFace.nodeCount=3;
					tempNodesSize=3;
					tempNodes= new unsigned int[3];
				} else {
					tempFace.nodeCount=4;
					tempNodesSize=4;
					tempNodes= new unsigned int[4];
				}
			}
			if (elemType==10) tempFace.nodeCount=3;
			tempFace.id=faceCount;
			// Assign current cell as the parent cell
			tempFace.parent=i;
			// Assign boundary type as internal by default, will be overwritten when BCs applied
			tempFace.bc=-1;
			// Store the nodes of the current face
			if (faceCount==face.capacity()) face.reserve(int (face.size() *0.10) +face.size()) ;
			for (unsigned int n=0;n<tempFace.nodeCount;++n) {
				if (elemType==17) tempNodes[n]=cell[i].node(hexaFaces[f][n]).id;
				if (elemType==14) tempNodes[n]=cell[i].node(wedgeFaces[f][n]).id;
				if (elemType==10) tempNodes[n]=cell[i].node(tetraFaces[f][n]).id;
				tempFace.nodes.push_back(tempNodes[n]);
			}
			// Find the neighbor cell
			int flagInternal=0;
			int flagExists=1;
			for (unsigned int j=0;j<node[tempNodes[0]].cells.size();++j) {
				int c=node[tempNodes[0]].cells[j];
				if (int (i) !=c && cell[c].HaveNodes(tempNodesSize,tempNodes)) {
					tempFace.neighbor=c;
					flagInternal=1;
					// Check if any other face was found before with parent=c
					//flagExists=face_exists(c);
					if (int (i) <c) flagExists=0;
					break;
				}
			}
			if (!flagInternal) {
				++boundaryFaceCount;
				tempFace.neighbor=-1*boundaryFaceCount;
				face.push_back(tempFace);
				for (unsigned int fn=0;fn<tempNodesSize;++fn) face[tempFace.id].nodes.push_back(tempNodes[fn]);
				cell[i].faces.push_back(tempFace.id);
				++faceCount;
			} else if (!flagExists) {
				face.push_back(tempFace);
				for (unsigned int fn=0;fn<tempNodesSize;++fn) face[tempFace.id].nodes.push_back(tempNodes[fn]);
				cell[i].faces.push_back(tempFace.id);
				cell[tempFace.neighbor].faces.push_back(tempFace.id);
				++faceCount;
			}
		} //for face
	} // for cells

	cout << "* Number of Faces: " << faceCount << endl;
	cout << "* Number of Faces at Boundaries: " << boundaryFaceCount << endl;

// Now loop through faces and calculate centroids and areas
	for (unsigned int f=0;f<faceCount;++f) {
		Vec3D centroid=0.;
		Vec3D areaVec=0.;
		for (unsigned int n=0;n<face[f].nodeCount;++n) {
			centroid+=face[f].node(n);
		}
		centroid/=face[f].nodeCount;
		face[f].centroid=centroid;
		for (unsigned int n=0;n<face[f].nodeCount;n++) {
			areaVec+=0.5* (face[f].node(n)-centroid).cross(face[f].node(n+1)-centroid);
		}
		if (areaVec.dot(centroid-cell[face[f].parent].centroid) <0.) {
			// [TBM] Need to swap the face and reflect the area vector
			cout << "face " << f << " should be swapped" << endl;
		}
		face[f].area=fabs(areaVec);
		//cout << face[f].area << endl;
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

int Cell::Construct(const int celltype, const int nc, const int nodeList[]) {

	type=celltype;
	if (type==HEXA) faceCount=6;
	if (type==WEDGE) faceCount=5;
	if (type==TETRA) faceCount=4;
	nodeCount=nc;
	nodes.reserve(nodeCount);
	if (nodeCount==0) {
		cerr << "[!!] Number of nodes of the cell must be specified before allocation" << endl;
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
		areaVec=face[f].normal*face[f].area/cell[parent].volume;
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
			for (unsigned int i=0;i<5;++i) cell[neighbor].grad[i]-=average[i]*areaVec;
		} else {	// Boundary face
			if (bc.region[face[f].bc].type=="wall") {
				for (unsigned int i=0;i<3;++i) averageVel.comp[i]=average[i+1];
				averageVel-=averageVel.dot(face[f].normal) *face[f].normal;
				for (unsigned int i=0;i<3;++i) average[i+1]=averageVel.comp[i];
				//wallCellMark.push_back(parent);
				//wallCellFaceNormal.push_back(grid.face[f].normal);
			} else if (bc.region[face[f].bc].type=="inlet") {
				average[0]=bc.region[grid.face[f].bc].rho;
				average[1]=bc.region[grid.face[f].bc].v.comp[0];
				average[2]=bc.region[grid.face[f].bc].v.comp[1];
				average[3]=bc.region[grid.face[f].bc].v.comp[2];
				average[4]=bc.region[grid.face[f].bc].p;
			} else if (bc.region[face[f].bc].type=="outlet") {
				average[0]=grid.cell[parent].rho;
				average[1]=grid.cell[parent].v.comp[0];
				average[2]=grid.cell[parent].v.comp[1];
				average[3]=grid.cell[parent].v.comp[2];
				average[4]=grid.cell[parent].p;
			}
		}
		cell[parent].grad[0]+=average[0]*areaVec;
		cell[parent].grad[4]+=average[4]*areaVec;
		for (unsigned int i=1;i<4;++i) cell[parent].grad[i]+=average[i]*areaVec;
	}
}

