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
#include <cmath>
using namespace std;

#include "grid.h"

extern Grid grid;
extern int np, Rank;

extern GridRawData raw;
extern IndexMaps maps;

int Grid::create_nodes_cells() {

	// Reserve the right amount of memory beforehand
	node.reserve(nodeCount);
	cell.reserve(cellCount);
	
	// This stores the total node count in the current partition
	nodeCount=0;
	
	for (unsigned int c=0;c<globalCellCount;++c) {
		if (maps.cellOwner[c]==Rank) { // If the cell belongs to current proc
			int cellNodeCount; // Find the number of nodes of the cell from raw grid data
			if (c<globalCellCount-1) {
				cellNodeCount=raw.cellConnIndex[c+1]-raw.cellConnIndex[c];
			} else {
				cellNodeCount=raw.cellConnectivity.size()-raw.cellConnIndex[globalCellCount-1];
			}
			unsigned int cellNodes[cellNodeCount];
			for (unsigned int n=0;n<cellNodeCount;++n) { // Loop the cell  nodes
				if (maps.nodeGlobal2Local.find(raw.cellConnectivity[raw.cellConnIndex[c]+n])== maps.nodeGlobal2Local.end() ) { // If the node is not already found
					// Create the node
					Node temp;
					temp.id=nodeCount;
					temp.globalId=raw.cellConnectivity[raw.cellConnIndex[c]+n];
					temp.comp[0]=raw.x[temp.globalId];
					temp.comp[1]=raw.y[temp.globalId];
					temp.comp[2]=raw.z[temp.globalId];
					node.push_back(temp);
					maps.nodeGlobal2Local[temp.globalId]=temp.id;
					++nodeCount;
				}
				// Fill in cell nodes temp array with local node id's
				cellNodes[n]=maps.nodeGlobal2Local[raw.cellConnectivity[raw.cellConnIndex[c]+n]];
			} // end for each cell node
			// Create the cell
			Cell temp;
			temp.nodeCount=cellNodeCount;
			switch (cellNodeCount) {
				case 4: // Tetra
					temp.faceCount=4;
					break;
				case 5: // Pyramid
					temp.faceCount=5;
					break;
				case 6: // Prism
					temp.faceCount=5;
					break;
				case 8: // Hexa
					temp.faceCount=6;
					break;
			}
			temp.nodes.reserve(cellNodeCount);
			// Fill in the node list and calculate the centroid
			temp.centroid=0.;
			for (int n=0;n<temp.nodeCount;++n) {
				temp.nodes.push_back(cellNodes[n]);
				temp.centroid+=temp.node(n);
			}
			temp.centroid/=double(temp.nodeCount);
			temp.globalId=c;
			maps.cellGlobal2Local[temp.globalId]=cell.size();
			cell.push_back(temp);
		} // end if cell is in current proc
	} // end loop global cell count

	cout << "[I Rank=" << Rank << "] Created cells and nodes" << endl;

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

	cout << "[I Rank=" << Rank << "] Computed node-cell connectivity" << endl;
	
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

	cout << "[I Rank=" << Rank << "] Computed cell-cell connectivity" << endl;
	
	// Mark nodes that touch boundaries and store which boundary(s) in node.bc set
	for (int breg=0;breg<raw.bocoNameMap.size();++breg) { // for each boundary condition region defined in grid file
		int bfnodeCount; // Boundary face node count
		for (int bf=0;bf<raw.bocoConnIndex[breg].size();++bf) { // for each boundary face in current region
			if (raw.bocoConnectivity[breg][raw.bocoConnIndex[breg][bf]]!=-1) {
				// Get number of nodes of the boundary face
				if (bf==raw.bocoConnIndex[breg].size()-1) {
					bfnodeCount=raw.bocoConnectivity[breg].size()-raw.bocoConnIndex[breg][bf];
				} else {
					bfnodeCount=raw.bocoConnIndex[breg][bf+1]-raw.bocoConnIndex[breg][bf];
				}
				unsigned int faceNode;
				// See if the face is on the current partition
				bool onCurrent=true;
				for (unsigned int i=0;i<bfnodeCount;++i) {
					if (maps.nodeGlobal2Local.find(raw.bocoConnectivity[breg][raw.bocoConnIndex[breg][bf]+i]) == maps.nodeGlobal2Local.end()) {
						onCurrent=false;
						break;
					}
				}
				if (onCurrent) {
					// Grab the first node of the the boundary face
					faceNode=maps.nodeGlobal2Local[raw.bocoConnectivity[breg][raw.bocoConnIndex[breg][bf]] ];
					maps.nodeBCregions[faceNode].push_back(breg);
					std::set<unsigned int> nodeSet;
					nodeSet.clear();
					for (unsigned int i=0;i<bfnodeCount;++i) {
						nodeSet.insert(maps.nodeGlobal2Local[raw.bocoConnectivity[breg][raw.bocoConnIndex[breg][bf]+i] ]);
						
					}
					maps.nodeBCregionFaceConn[faceNode].push_back(nodeSet);
				}
			}
		}
	}
	
	
} //end Grid::create_nodes_cells

int Grid::create_faces() {

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
	int pyraFaces[5][4]= {
		{0,1,2,3},
		{1,0,4,0},
  		{0,3,4,0},
  		{3,2,4,0},
  		{2,1,4,0}
	};
	int tetraFaces[4][3]= {
		{0,2,1},
		{1,2,3},
  		{0,3,2},
  		{0,1,3}
	};

	// Search and construct faces
	faceCount=0;

	// Keep track of how many faces are on each boundary
	int boundaryFaceCount[raw.bocoNameMap.size()];
	for (int boco=0;boco<raw.bocoNameMap.size();++boco) boundaryFaceCount[boco]=0;
	
	// Time the face search
	double timeRef, timeEnd;
	if (Rank==0) timeRef=MPI_Wtime();

	// Loop through all the cells
	for (unsigned int c=0;c<cellCount;++c) {
		// Loop through the faces of the current cell
		for (unsigned int cf=0;cf<cell[c].faceCount;++cf) {
			Face tempFace;
			unsigned int *tempNodes;
			switch (cell[c].nodeCount) {
				case 4: // Tetrahedra
					tempFace.nodeCount=3;
					tempNodes= new unsigned int[3];
					break;
				case 5: // Pyramid
						if (cf<1) {
							tempFace.nodeCount=4;
							tempNodes= new unsigned int[4];
						} else {
							tempFace.nodeCount=3;
							tempNodes= new unsigned int[3];
						}
						break;
				case 6: // Prism
					if (cf<2) {
						tempFace.nodeCount=3;
						tempNodes= new unsigned int[3];
					} else {
						tempFace.nodeCount=4;
						tempNodes= new unsigned int[4];
					}
					break;
				case 8: // Brick 
					tempFace.nodeCount=4;
					tempNodes= new unsigned int[4];
					break;
			}
			// Face count is incremented everytime a new face is found
			tempFace.id=faceCount;
			// Assign current cell as the parent cell
			tempFace.parent=c;
			// Assign boundary type as internal by default, will be overwritten later
			tempFace.bc=-1;
			// Store the node local ids of the current face	
			for (unsigned int fn=0;fn<tempFace.nodeCount;++fn) {
				switch (cell[c].nodeCount) {
					case 4: tempNodes[fn]=cell[c].node(tetraFaces[cf][fn]).id; break;
					case 5: tempNodes[fn]=cell[c].node(pyraFaces[cf][fn]).id; break;
					case 6: tempNodes[fn]=cell[c].node(prismFaces[cf][fn]).id; break;
					case 8: tempNodes[fn]=cell[c].node(hexaFaces[cf][fn]).id; break;
				}
			}
			// Find the neighbor cell
			bool internal=false;
			bool unique=true;
			// Loop cells neighboring the first node of the current face
			for (unsigned int nc=0;nc<node[tempNodes[0]].cells.size();++nc) {
				// i is the neighbor cell index
				unsigned int i=node[tempNodes[0]].cells[nc];
				// If neighbor cell is not the current cell itself, and it has the same nodes as the face
				if (i!=c && cell[i].HaveNodes(tempFace.nodeCount,tempNodes)) {
					// If the neighbor cell index is smaller then the current cell index,
					// it has already been processed so skip it
					if (i>c) {
						tempFace.neighbor=i;
						internal=true;
					} else {
						unique=false;
					}
				}
			}
			if (unique) { // If a new face
				// Insert the node list
				for (unsigned int fn=0;fn<tempFace.nodeCount;++fn) tempFace.nodes.push_back(tempNodes[fn]);
				if (!internal) { // If the face is either at inter-partition or boundary
					tempFace.bc=-2; // Unassigned boundary type
					// As the bc faces and inter-partition faces are found, this will be overwritten.
					// At the end, there should be no face with this value.
					// Loop the current face nodes
					for (unsigned int i=0;i<tempFace.nodeCount;++i) {
						unsigned int nn=tempNodes[i];
						bool match;
						// Find 
						if(maps.nodeBCregions.find(nn)!=maps.nodeBCregions.end()) { // If this node touches a boundary
							// Loop through the different boundary regions assigned for that node
							for (int j=0;j<maps.nodeBCregions[nn].size();++j) {
								// For each boundary region type, loop the current face nodes to see it matches the boundary face 
								// Loop the current face nodes
								match=true;
								for (int k=0;k<tempFace.nodeCount;++k) {
									// if the current face node is not in the current boundary face's node list
									if (maps.nodeBCregionFaceConn[nn][j].find(tempNodes[k])==maps.nodeBCregionFaceConn[nn][j].end()) {
										match=false;
										break;
									}
								}
								if (match) {
									tempFace.bc=maps.nodeBCregions[nn][j];
									boundaryFaceCount[maps.nodeBCregions[nn][j]]++;
									break;
								}
							} // For each boundary region that the node belongs to
							// If a match is found, exit the loop for other bc regions
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
	for (int boco=0;boco<raw.bocoNameMap.size();++boco) cout << "[I Rank=" << Rank << "] Number of Faces on BC_" << boco+1 << "=" << boundaryFaceCount[boco] << endl;
	
	if (Rank==0) {
		timeEnd=MPI_Wtime();
		cout << "[I Rank=" << Rank << "] Time spent on finding faces= " << timeEnd-timeRef << " sec" << endl;
	}

	
} // end Grid::create_faces

int Grid::create_ghosts() {

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

		// Find out other partition's cell counts
		int otherCellCounts[np];
		for (int i=0;i<np;++i) otherCellCounts[i]=0;
		for (unsigned int c=0;c<globalCellCount;++c) otherCellCounts[maps.cellOwner[c]]+=1;
		
		for (unsigned int p=0;p<np;++p) {
			cellCountOffset[p]=counter;
			counter+=otherCellCounts[p];
		}
		
		// Now find metis2global index mapping
		int metis2global[globalCellCount];
		int counter2[np];
		for (unsigned int p=0;p<np;++p) counter2[p]=0;
		for (unsigned int c=0;c<globalCellCount;++c) {
			metis2global[cellCountOffset[maps.cellOwner[c]]+counter2[maps.cellOwner[c]]]=c;
			counter2[maps.cellOwner[c]]++;
		}

		int foundFlag[globalCellCount];
		for (unsigned int c=0; c<globalCellCount;++c) foundFlag[c]=0;

		unsigned int parent, metisIndex, gg, matchCount;
		
		map<int,set<int> > nodeCellSet;
		Vec3D nodeVec;
		
		for (unsigned int f=0; f<faceCount; ++f) {
			if (face[f].bc==-2) { // if an assigned boundary face
				parent=face[f].parent;
				for (unsigned int adjCount=0;adjCount<(maps.adjIndex[parent+1]-maps.adjIndex[parent]);++adjCount)  {
					metisIndex=maps.adjacency[maps.adjIndex[parent]+adjCount];
					gg=metis2global[metisIndex];
					int cellNodeCount; // Find the number of nodes of the cell from raw grid data
					if (gg<globalCellCount-1) {
						cellNodeCount=raw.cellConnIndex[gg+1]-raw.cellConnIndex[gg];
					} else {
						cellNodeCount=raw.cellConnectivity.size()-raw.cellConnIndex[globalCellCount-1];
					}
					if (metisIndex<cellCountOffset[Rank] || metisIndex>=(cellCount+cellCountOffset[Rank])) {
						matchCount=0;
						for (unsigned int fn=0;fn<face[f].nodeCount;++fn) {
							set<int> tempSet;
							nodeCellSet.insert(pair<unsigned int,set<int> >(face[f].nodes[fn],tempSet) );
							for (unsigned int gn=0;gn<cellNodeCount;++gn) {
								if (raw.cellConnectivity[raw.cellConnIndex[gg]+gn]==face[f].node(fn).globalId) {
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
							temp.partition=maps.cellOwner[gg];
							// Calculate the centroid
							temp.centroid=0.;

							for (unsigned int gn=0;gn<cellNodeCount;++gn) {
								nodeVec.comp[0]=raw.x[raw.cellConnectivity[raw.cellConnIndex[gg]+gn]];
								nodeVec.comp[1]=raw.y[raw.cellConnectivity[raw.cellConnIndex[gg]+gn]];
								nodeVec.comp[2]=raw.z[raw.cellConnectivity[raw.cellConnIndex[gg]+gn]];
								temp.centroid+=nodeVec;
							}
							temp.centroid/=double(cellNodeCount);
							maps.ghostGlobal2Local[temp.globalId]=ghost.size();
							ghost.push_back(temp);			
							++ghostCount;
						}
						if (matchCount>=3) {
							foundFlag[gg]=3;
							face[f].bc=-1*maps.ghostGlobal2Local[gg]-3;
						}
					}
				}
			}
		}

		map<int,set<int> >::iterator mit;
		set<int>::iterator sit;
		for ( mit=nodeCellSet.begin() ; mit != nodeCellSet.end(); mit++ ) {
			for ( sit=(*mit).second.begin() ; sit != (*mit).second.end(); sit++ ) {
				node[(*mit).first].ghosts.push_back(maps.ghostGlobal2Local[*sit]);
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
	
} // end int Grid::create_ghosts

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