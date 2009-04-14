/************************************************************************
	
	Copyright 2007-2009 Emre Sozer & Patrick Clark Trizila

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
#include "commons.h"
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <iomanip>
#include <cmath>
using namespace std;

extern GridRawData raw;
extern IndexMaps maps;

int Grid::create_nodes_cells() {

	// Reserve the right amount of memory beforehand
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
				unsigned int gid=raw.cellConnectivity[raw.cellConnIndex[c]+n]; // node globalId
				if (maps.nodeGlobal2Local.find(gid)==maps.nodeGlobal2Local.end() ) { // If the node is not already found
					// Create the node
					Node temp;
					temp.id=nodeCount;
					temp.globalId=gid;
					temp.comp[0]=raw.x[temp.globalId];
					temp.comp[1]=raw.y[temp.globalId];
					temp.comp[2]=raw.z[temp.globalId];
					maps.nodeGlobal2Local[temp.globalId]=temp.id;
					node.push_back(temp);
					++nodeCount;
				}
				// Fill in cell nodes temp array with local node id's
				cellNodes[n]=maps.nodeGlobal2Local[gid];
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
			
			// Fill in the node list
			for (int n=0;n<temp.nodeCount;++n) {
				temp.nodes.push_back(cellNodes[n]);
			}
			
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
	
	// Construct the list of neighboring cells (node neighbors) for each cell
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
	
	for (int nbc=0;nbc<raw.bocoNodes.size();++nbc) {
		set<int> temp;
		set<int>::iterator sit;
		for (sit=raw.bocoNodes[nbc].begin();sit!=raw.bocoNodes[nbc].end();sit++) {
			if (maps.nodeGlobal2Local.find(*sit)!=maps.nodeGlobal2Local.end()) {
				temp.insert(maps.nodeGlobal2Local[*sit]);
			}
		}
		raw.bocoNodes[nbc].swap(temp);
		temp.clear();
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
		{0,1,4,3}
	};
	int pyraFaces[5][4]= {
		{0,3,2,1},
		{0,1,4,0},
  		{1,2,4,0},
  		{3,4,2,0},
  		{0,4,3,0}
	};
	int tetraFaces[4][3]= {
		{0,2,1},
		{1,2,3},
  		{0,3,2},
  		{0,1,3}
	};

	// Search and construct faces
	faceCount=0;
	
	// Time the face search
	double timeRef, timeEnd;
	if (Rank==0) timeRef=MPI_Wtime();
	vector<unsigned int> unique_nodes;
	set<unsigned int> repeated_node_cells;
	// Loop through all the cells
	for (unsigned int c=0;c<cellCount;++c) {
		int degenerate_face_count=0;
		// Loop through the faces of the current cell
		for (unsigned int cf=0;cf<cell[c].faceCount;++cf) {
			bool degenerate=false;
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
			tempFace.bc=INTERNAL;
			// Store the node local ids of the current face	
			for (unsigned int fn=0;fn<tempFace.nodeCount;++fn) {
				switch (cell[c].nodeCount) {
					case 4: tempNodes[fn]=cell[c].node(tetraFaces[cf][fn]).id; break;
					case 5: tempNodes[fn]=cell[c].node(pyraFaces[cf][fn]).id; break;
					case 6: tempNodes[fn]=cell[c].node(prismFaces[cf][fn]).id; break;
					case 8: tempNodes[fn]=cell[c].node(hexaFaces[cf][fn]).id; break;
				}
			}
			// Check if there is a repeated node
			unique_nodes.clear();
			bool skip;
			for (unsigned int fn=0;fn<tempFace.nodeCount;++fn) {
				skip=false;
				for (unsigned int i=0;i<fn;++i) {
					if (tempNodes[fn]==tempNodes[i]) {
						skip=true;
						break;
					}
				}
				if (!skip) unique_nodes.push_back(tempNodes[fn]);
			}
			if (unique_nodes.size()!=tempFace.nodeCount) {
				repeated_node_cells.insert(c); // mark the owner cell (it has repeated nodes)
				if (unique_nodes.size()==2) { // If a face only has two unique nodes, mark as degenerate
					degenerate=true;
					degenerate_face_count++;
				}
				tempFace.nodeCount=unique_nodes.size();
				for (unsigned int fn=0;fn<tempFace.nodeCount;++fn) tempNodes[fn]=unique_nodes[fn];
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

			if (unique && !degenerate) { // If a new face
				// Insert the node list
				for (unsigned int fn=0;fn<tempFace.nodeCount;++fn) tempFace.nodes.push_back(tempNodes[fn]);
				if (!internal) { // If the face is either at inter-partition or boundary
					tempFace.bc=UNASSIGNED; // yet
					vector<int> face_matched_bcs;
					int cell_matched_bc=-1;
					bool match;
					for (int nbc=0;nbc<raw.bocoNameMap.size();++nbc) { // For each boundary condition region
						match=true;
						for (unsigned int i=0;i<tempFace.nodeCount;++i) { // For each node of the current face
							if (raw.bocoNodes[nbc].find(tempNodes[i])==raw.bocoNodes[nbc].end()) {
								match=false;
								break;
							}
						}
						if (match) { // This means that all the face nodes are on the current bc node list
							face_matched_bcs.push_back(nbc);
						}
						// There can be situations like back and front symmetry BC's in which
						// face nodes will match more than one boundary condition
						// Check if the owner cell has all its nodes on one of those bc's
						// and eliminate those
						if (cell_matched_bc==-1) {
							match=true;
							for (unsigned int i=0;i<cell[c].nodeCount;++i) { 
								if (raw.bocoNodes[nbc].find(cell[c].nodes[i])==raw.bocoNodes[nbc].end()) {
									match=false;
									break;
								}
							}
							if (match) { // This means that all the cell nodes are on the current bc node list
								cell_matched_bc=nbc;
							}
						}
						
						
					}
					if (face_matched_bcs.size()>1) {
						for (int fbc=0;fbc<face_matched_bcs.size();++fbc) {
							if(face_matched_bcs[fbc]!=cell_matched_bc) {
								tempFace.bc=face_matched_bcs[fbc];
								break;
							}
						}
					} else if (face_matched_bcs.size()==1) {
						tempFace.bc=face_matched_bcs[0];
					}
					// Some of these bc values will be overwritten later if the face is at a partition interface

				} // if not internal
				tempFace.parentIndex=cf;
				face.push_back(tempFace);
				cell[c].faces.push_back(tempFace.id);
				if (internal) cell[tempFace.neighbor].faces.push_back(tempFace.id);
				++faceCount;
			}
			delete [] tempNodes;
		} //for face cf
		cell[c].faceCount-=degenerate_face_count;
	} // for cells c

	// Loop cells that has repeated nodes and fix the node list
	set<unsigned int>::iterator sit,sit2;
	set<unsigned int> unique_cell_nodes;
	for (sit=repeated_node_cells.begin();sit!=repeated_node_cells.end();sit++) {
		for (int cn=0;cn<cell[(*sit)].nodeCount;++cn) unique_cell_nodes.insert(cell[(*sit)].nodes[cn]);
		cell[(*sit)].nodes.clear();
		for (sit2=unique_cell_nodes.begin();sit2!=unique_cell_nodes.end();sit2++) cell[(*sit)].nodes.push_back((*sit2));
		cell[(*sit)].nodeCount=unique_cell_nodes.size();
		unique_cell_nodes.clear();
	}
	repeated_node_cells.clear();
	
	for (unsigned int c=0; c<cellCount; ++c) {
		if (Rank==0) if (cell[c].faceCount != cell[c].faces.size() ) cout << "no match" << "\t" << c << "\t" << cell[c].faceCount << "\t" << cell[c].faces.size() << "\t" << cell[c].nodeCount << endl;
	}
	
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
			cell[c].ghostCount=0;
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
		
		map<int,set<int> > nodeGhostSet; // Stores global id's of ghost cells near each node
		Vec3D nodeVec;
		
		// Loop faces
		for (unsigned int f=0; f<faceCount; ++f) {
			if (face[f].bc==UNASSIGNED || face[f].bc>=0) { // if an unassigned boundary face
				parent=face[f].parent;
				// Loop through the cells that are adjacent to the current face's parent
				for (unsigned int adjCount=0;adjCount<(maps.adjIndex[parent+1]-maps.adjIndex[parent]);++adjCount)  {
					metisIndex=maps.adjacency[maps.adjIndex[parent]+adjCount];
					// Get global id of the adjacent cell
					gg=metis2global[metisIndex];
					int cellNodeCount; // Find the number of nodes of the cell from raw grid data
					if (gg<globalCellCount-1) {
						cellNodeCount=raw.cellConnIndex[gg+1]-raw.cellConnIndex[gg];
					} else {
						cellNodeCount=raw.cellConnectivity.size()-raw.cellConnIndex[globalCellCount-1];
					}
					// If that cell is not on the current partition
					if (metisIndex<cellCountOffset[Rank] || metisIndex>=(cellCount+cellCountOffset[Rank])) {
						// Count number of matches in node lists of the current face and the adjacent cell
						matchCount=0;
						for (unsigned int fn=0;fn<face[f].nodeCount;++fn) {
							set<int> tempSet;
							nodeGhostSet.insert(pair<unsigned int,set<int> >(face[f].nodes[fn],tempSet) );
							for (unsigned int gn=0;gn<cellNodeCount;++gn) {
								if (raw.cellConnectivity[raw.cellConnIndex[gg]+gn]==face[f].node(fn).globalId) {
									nodeGhostSet[face[f].nodes[fn]].insert(gg);
									++matchCount;
								}
							}
						}
						// foundFlag is 0 by default
						// 0 means that particular adjacent cell wasn't discovered as a ghost before
						// If so, create a new ghost
						if (matchCount>0 && foundFlag[gg]==0) {
							foundFlag[gg]=matchCount;
							Ghost temp;
							temp.globalId=gg;
							temp.partition=maps.cellOwner[gg];
							maps.ghostGlobal2Local[temp.globalId]=ghost.size();
							ghost.push_back(temp);	
							++ghostCount;
						}
						// If that ghost was found before, now we discovered another face also neighbors the same ghost
						if (matchCount==face[f].nodeCount) {
							foundFlag[gg]=matchCount;
							face[f].bc=GHOST;
							face[f].neighbor=-1*maps.ghostGlobal2Local[gg]-1;
						}
					}
				}
			}
		}
		
		// Store the local id's of ghosts touching each node
		map<int,set<int> >::iterator mit;
		set<int>::iterator sit;
		for ( mit=nodeGhostSet.begin() ; mit != nodeGhostSet.end(); mit++ ) {
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
	
	maps.nodeGlobal2Output.resize(globalNodeCount);
	for (unsigned int ng=0;ng<globalNodeCount;++ng) maps.nodeGlobal2Output[ng]=-1;
	nodeCountOffset=0;
	int nodeCountOffsetPlus=0;
	// Receive my nodeCountOffset from processor Rank-1

	// Set tag to destination
	if (Rank!=0) MPI_Recv(&nodeCountOffset,1,MPI_INT,Rank-1,Rank,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	if (Rank!=0) MPI_Recv(&maps.nodeGlobal2Output[0],globalNodeCount,MPI_INT,Rank-1,Rank,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	for (unsigned int n=0;n<nodeCount;++n) {
		if (maps.nodeGlobal2Output[node[n].globalId]==-1) {
			maps.nodeGlobal2Output[node[n].globalId]=nodeCountOffset+nodeCountOffsetPlus;
			nodeCountOffsetPlus++;
		}
	}
	nodeCountOffsetPlus+=nodeCountOffset;

 	if (Rank!=np-1) MPI_Send(&nodeCountOffsetPlus,1,MPI_INT,Rank+1,Rank+1,MPI_COMM_WORLD);
	if (Rank!=np-1) MPI_Send(&maps.nodeGlobal2Output[0],globalNodeCount,MPI_INT,Rank+1,Rank+1,MPI_COMM_WORLD);
	
} // end int Grid::create_ghosts

bool Cell::HaveNodes(unsigned int &nodelistsize, unsigned int nodelist []) {	
	bool match;
	for (unsigned int i=0;i<nodelistsize;++i) {
		match=false;
		for (unsigned int j=0;j<nodeCount;++j) if (nodelist[i]==nodes[j]) {match=true; break;}
		if (!match) return false;
	}
	return true;
}
