/************************************************************************
	
	Copyright 2007-2010 Emre Sozer

	Contact: emresozer@freecfd.com

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
#include "grid.h"

using namespace std;

int Grid::create_nodes_cells() {

	// Reserve the right amount of memory beforehand
	cell.reserve(cellCount);
	
	// This stores the total node count in the current partition
	nodeCount=0;
	
	for (int c=0;c<globalCellCount;++c) {
		if (maps.cellOwner[c]==Rank) { // If the cell belongs to current proc
			int cellNodeCount; // Find the number of nodes of the cell from raw grid data
			if (c<globalCellCount-1) {
				cellNodeCount=raw.cellConnIndex[c+1]-raw.cellConnIndex[c];
			} else {
				cellNodeCount=raw.cellConnectivity.size()-raw.cellConnIndex[globalCellCount-1];
			}
			int cellNodes[cellNodeCount];
			for (int n=0;n<cellNodeCount;++n) { // Loop the cell  nodes
				int ngid=raw.cellConnectivity[raw.cellConnIndex[c]+n]; // node globalId
				if (maps.nodeGlobal2Local.find(ngid)==maps.nodeGlobal2Local.end() ) { // If the node is not already found
					// Create the node
					Node temp;
					temp.id=nodeCount;
					temp.globalId=ngid;
					temp.comp[0]=raw.node[temp.globalId][0];
					temp.comp[1]=raw.node[temp.globalId][1];
					temp.comp[2]=raw.node[temp.globalId][2];
					maps.nodeGlobal2Local[temp.globalId]=temp.id;
					node.push_back(temp);
					++nodeCount;
				}
				// Fill in cell nodes temp array with local node id's
				cellNodes[n]=maps.nodeGlobal2Local[ngid];
			} // end for each cell node
			// Create the cell
			Cell temp;
			temp.nodeCount=cellNodeCount;
			if (raw.type==CELL) {
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
			} else {
				temp.faceCount=0;
				// Skip the face count for now.
				// It is done just after the cell loop for efficiency
				// Loop the left and right data and count
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

	cellCount=cell.size();
	
	if (raw.type==FACE) {
		int c;
		for (int i=0;i<raw.left.size();++i) {
			c=raw.left[i];
			if (c>=0 && maps.cellOwner[c]==Rank) { // If the cell belongs to current proc
				cell[maps.cellGlobal2Local[c]].faceCount++;
			}
		}
		for (int i=0;i<raw.right.size();++i) {
			c=raw.right[i];
			if (c>=0 && maps.cellOwner[c]==Rank) { // If the cell belongs to current proc
				cell[maps.cellGlobal2Local[c]].faceCount++;
			}
		}
	}
		
	if (Rank==0) cout << "[I] Created cells and nodes" << endl;

	// Construct the list of cells for each node
	bool flag;
	for (int c=0;c<cellCount;++c) {
		int n;
		for (int cn=0;cn<cell[c].nodeCount;++cn) {
			n=cell[c].nodes[cn];
			flag=false;
			for (int i=0;i<node[n].cells.size();++i) {
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

	if (Rank==0) cout << "[I] Computed node-cell connectivity" << endl;
	
	// Construct the list of neighboring cells (node neighbors) for each cell
	int c2;
	for (int c=0;c<cellCount;++c) {
		int n;
		for (int cn=0;cn<cell[c].nodeCount;++cn) { // Loop nodes of the cell
			n=cell[c].nodes[cn];
			for (int nc=0;nc<node[n].cells.size();++nc) { // Loop neighboring cells of the node
				c2=node[n].cells[nc];
				flag=false;
				for (int cc=0;cc<cell[c].neighborCells.size();++cc) { // Check if the cell was found before
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

	if (Rank==0) cout << "[I] Computed cell-cell connectivity" << endl;
	
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
	
	return 0;
	
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
	vector<int> unique_nodes;
	set<int> repeated_node_cells;
	// Loop through all the cells
	for (int c=0;c<cellCount;++c) {
		int degenerate_face_count=0;
		// Loop through the faces of the current cell
		for (int cf=0;cf<cell[c].faceCount;++cf) {
			bool degenerate=false;
			Face tempFace;
			int *tempNodes;
			switch (cell[c].nodeCount) {
				case 4: // Tetrahedra
					tempFace.nodeCount=3;
					tempNodes= new int[3];
					break;
				case 5: // Pyramid
					if (cf<1) {
						tempFace.nodeCount=4;
						tempNodes= new int[4];
					} else {
						tempFace.nodeCount=3;
						tempNodes= new int[3];
					}
					break;
				case 6: // Prism
					if (cf<2) {
						tempFace.nodeCount=3;
						tempNodes= new int[3];
					} else {
						tempFace.nodeCount=4;
						tempNodes= new int[4];
					}
					break;
				case 8: // Brick 
					tempFace.nodeCount=4;
					tempNodes= new int[4];
					break;
			}
			// Face count is incremented everytime a new face is found
			tempFace.id=faceCount;
			// Assign current cell as the parent cell
			tempFace.parent=c;
			// Assign boundary type as internal by default, will be overwritten later
			tempFace.bc=INTERNAL_FACE;
			// Store the node local ids of the current face	
			for (int fn=0;fn<tempFace.nodeCount;++fn) {
				switch (cell[c].nodeCount) {
					case 4: tempNodes[fn]=cellNode(c,tetraFaces[cf][fn]).id; break;
					case 5: tempNodes[fn]=cellNode(c,pyraFaces[cf][fn]).id; break;
					case 6: tempNodes[fn]=cellNode(c,prismFaces[cf][fn]).id; break;
					case 8: tempNodes[fn]=cellNode(c,hexaFaces[cf][fn]).id; break;
				}
			}
			// Check if there is a repeated node
			unique_nodes.clear();
			bool skip;
			for (int fn=0;fn<tempFace.nodeCount;++fn) {
				skip=false;
				for (int i=0;i<fn;++i) {
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
				for (int fn=0;fn<tempFace.nodeCount;++fn) tempNodes[fn]=unique_nodes[fn];
			}
			// Find the neighbor cell
			bool internal=false;
			bool unique=true;
			// Loop cells neighboring the first node of the current face
			for (int nc=0;nc<node[tempNodes[0]].cells.size();++nc) {
				// i is the neighbor cell index
				int i=node[tempNodes[0]].cells[nc];
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
				for (int fn=0;fn<tempFace.nodeCount;++fn) tempFace.nodes.push_back(tempNodes[fn]);
				if (!internal) { // If the face is either at inter-partition or boundary
					tempFace.bc=UNASSIGNED_FACE; // yet
					vector<int> face_matched_bcs;
					int cell_matched_bc=-1;
					bool match;
					for (int nbc=0;nbc<raw.bocoNameMap.size();++nbc) { // For each boundary condition region
						match=true;
						for (int i=0;i<tempFace.nodeCount;++i) { // For each node of the current face
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
							for (int i=0;i<cell[c].nodeCount;++i) { 
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
	set<int>::iterator sit;
	vector<int> repeated_nodes;
	for (sit=repeated_node_cells.begin();sit!=repeated_node_cells.end();sit++) {
		// Find repeated nodes
		repeated_nodes.clear();
		for (int cn=0;cn<cell[(*sit)].nodeCount;++cn) {
			for (int cn2=0;cn2<cn;++cn2) {
				if (cell[(*sit)].nodes[cn]==cell[(*sit)].nodes[cn2]) repeated_nodes.push_back(cell[(*sit)].nodes[cn]);
			}
		}
		if (cell[(*sit)].nodeCount==8 && repeated_nodes.size()==2) { // TODO Only Hexa to Penta mapping is handled for now
			cell[(*sit)].nodes.clear();
			// Loop triangular cell faces
			int rindex=-1;
			for (int cf=0;cf<cell[(*sit)].faceCount;++cf) {
				if (cellFace((*sit),cf).nodeCount==3) {
					// Loop the face nodes and see if the repeated node apears
					int fn;
					for (fn=0;fn<3;++fn) {
						if (cellFace((*sit),cf).nodes[fn]==repeated_nodes[0]) { rindex=0; break; }
						if (cellFace((*sit),cf).nodes[fn]==repeated_nodes[1]) { rindex=1; break; }
					}
					// Start from fn and fill the new cell node list
					if (fn==0) {
						cell[(*sit)].nodes.push_back(cellFace((*sit),cf).nodes[0]);
						cell[(*sit)].nodes.push_back(cellFace((*sit),cf).nodes[1]);
						cell[(*sit)].nodes.push_back(cellFace((*sit),cf).nodes[2]);
					} else if (fn==1) {
						cell[(*sit)].nodes.push_back(cellFace((*sit),cf).nodes[1]);
						cell[(*sit)].nodes.push_back(cellFace((*sit),cf).nodes[2]);
						cell[(*sit)].nodes.push_back(cellFace((*sit),cf).nodes[0]);
					} else if (fn==2) {
						cell[(*sit)].nodes.push_back(cellFace((*sit),cf).nodes[2]);
						cell[(*sit)].nodes.push_back(cellFace((*sit),cf).nodes[0]);
						cell[(*sit)].nodes.push_back(cellFace((*sit),cf).nodes[1]);
					}
					
				}
				
			}
			cell[(*sit)].nodeCount=6;
		}
	}		
 	repeated_node_cells.clear();
	
	for (int c=0; c<cellCount; ++c) {
		if (Rank==0) if (cell[c].faceCount != cell[c].faces.size() ) cout << "no match" << "\t" << c << "\t" << cell[c].faceCount << "\t" << cell[c].faces.size() << "\t" << cell[c].nodeCount << endl;
	}
	
	if (Rank==0) {
		timeEnd=MPI_Wtime();
		cout << "[I] Time spent on finding faces= " << timeEnd-timeRef << " sec" << endl;
	}

	for (int f=0;f<faceCount;++f) {
		for (int n=0;n<face[f].nodes.size();++n) faceNode(f,n).faces.push_back(f);	
		face[f].symmetry=false; // by default
	}
	
	return 0;

} // end Grid::create_faces

int Grid::create_faces2() {

	// Search and construct faces
	faceCount=0;
	// Loop the global faces
	int parent,neighbor;
	bool owner,inter_partition,swap;
	int ghostFaceCount=0;
	for (int f=0;f<globalFaceCount;++f) {
		// If the left or right cell is owned by the current processor, create the face
		parent=raw.left[f];
		neighbor=raw.right[f];
		owner=false; 
		inter_partition=true;
		swap=true;
		if (maps.cellOwner[parent]==Rank) {
			owner=true;
			swap=false;
			parent=maps.cellGlobal2Local[parent];
		}
		if (neighbor>=0 && maps.cellOwner[neighbor]==Rank) {
			if (owner) inter_partition=false;
			owner=true;
			neighbor=maps.cellGlobal2Local[neighbor];
		}
		
		if (owner) { // This face needs to be created in the current processor
			Face tempFace;
			tempFace.id=faceCount;
			tempFace.nodeCount=raw.faceNodeCount[f];
			if (neighbor<0) tempFace.bc=UNASSIGNED_FACE;
			else if (inter_partition) {tempFace.bc=GHOST_FACE; ghostFaceCount++;}
			else tempFace.bc=INTERNAL_FACE;

			tempFace.parent=parent;
			tempFace.neighbor=neighbor;	
			tempFace.nodes.resize(tempFace.nodeCount);
			for (int fn=0;fn<tempFace.nodeCount;++fn) {
				tempFace.nodes[fn]=maps.nodeGlobal2Local[raw.faceConnectivity[raw.faceConnIndex[f]+fn]];
			}
			// If ghost face, parent may not be owned by the current processor. Check for that
			if (swap) {
				// Swap parent and neighbor
				tempFace.parent=neighbor;
				tempFace.neighbor=parent;
				// Swap tempFace node ordering
				vector<int> tempv=tempFace.nodes;
				vector<int>::reverse_iterator vit;
				int temp=0;
				for (vit=tempFace.nodes.rbegin();vit<tempFace.nodes.rend();++vit) {
					tempv[temp]=*vit;
					temp++;
				}
				tempFace.nodes.swap(tempv);
			}

			if (tempFace.bc==UNASSIGNED_FACE) { // If the face is at a boundary
				vector<int> face_matched_bcs;
				int cell_matched_bc=-1;
				bool match;
				for (int nbc=0;nbc<raw.bocoNameMap.size();++nbc) { // For each boundary condition region
					match=true;
					for (int i=0;i<tempFace.nodeCount;++i) { // For each node of the current face
						if (raw.bocoNodes[nbc].find(tempFace.nodes[i])==raw.bocoNodes[nbc].end()) {
							match=false;
							break;
						}
					}
					if (match) { // This means that all the face nodes are on the current bc node list
						face_matched_bcs.push_back(nbc);
					}
					// There can be situations like back and front symmetry BC's in which
					// face nodes will match more than one boundary condition
					// Check if the parent cell has all its nodes on one of those bc's
					// and eliminate those
					
					if (cell_matched_bc==-1) {
						match=true;
						for (int i=0;i<cell[tempFace.parent].nodeCount;++i) { 
							if (raw.bocoNodes[nbc].find(cell[tempFace.parent].nodes[i])==raw.bocoNodes[nbc].end()) {
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
				} else {
					cout << "[E] Couldn't find BC for face " << f << endl;
					exit(1);
				}

			} // if face is on boundary
			face.push_back(tempFace);
			cell[tempFace.parent].faces.push_back(tempFace.id);
			if (tempFace.bc==INTERNAL_FACE) cell[tempFace.neighbor].faces.push_back(tempFace.id);
			++faceCount;
		} // end if owner
	} // end global face loop

	// Set faceCount for each cell
	for (int c=0;c<cellCount;++c) cell[c].faceCount=cell[c].faces.size();

	for (int f=0;f<faceCount;++f) {
		for (int n=0;n<face[f].nodes.size();++n) faceNode(f,n).faces.push_back(f);	
		face[f].symmetry=false; // by default
	}

	//cout << "[I Rank=" << Rank << "] Number of inter-partition faces = " << ghostFaceCount << endl;

	return 0;
}

int Grid::create_ghosts() {
	// Determine and mark faces adjacent to other partitions
	// Create ghost elemets to hold the data from other partitions
	ghostCount=0;

	if (np==1) {
		for (int c=0;c<cellCount;++c) {
			cell[c].ghostCount=0;
		} // end cell loop
	} else {

		int counter=0;
		int cellCountOffset[np];

		// Find out other partition's cell counts
		int otherCellCounts[np];
		for (int i=0;i<np;++i) otherCellCounts[i]=0;
		for (int c=0;c<globalCellCount;++c) otherCellCounts[maps.cellOwner[c]]++;
		
		for (int p=0;p<np;++p) {
			cellCountOffset[p]=counter;
			counter+=otherCellCounts[p];
		}
		// Now find metis2global index mapping
		int metis2global[globalCellCount];
		int counter2[np];
		for (int p=0;p<np;++p) counter2[p]=0;
		for (int c=0;c<globalCellCount;++c) {
			metis2global[cellCountOffset[maps.cellOwner[c]]+counter2[maps.cellOwner[c]]]=c;
			counter2[maps.cellOwner[c]]++;
		}

		int foundFlag[globalCellCount];
		for (int c=0; c<globalCellCount;++c) foundFlag[c]=0;

		int parent, metisIndex, gg, matchCount;
		
		map<int,set<int> > nodeGhostSet; // Stores global id's of ghost cells near each node
		Vec3D nodeVec;
		
		// Loop faces
		for (int f=0; f<faceCount; ++f) {
			if ((face[f].bc==UNASSIGNED_FACE || face[f].bc==GHOST_FACE) || face[f].bc>=0) { // if an unassigned boundary face
				parent=face[f].parent;
				// Loop through the cells that are adjacent to the current face's parent
				for (int adjCount=0;adjCount<(maps.adjIndex[parent+1]-maps.adjIndex[parent]);++adjCount)  {
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
						for (int fn=0;fn<face[f].nodeCount;++fn) {
							set<int> tempSet;
							nodeGhostSet.insert(pair<int,set<int> >(face[f].nodes[fn],tempSet) );
							for (int gn=0;gn<cellNodeCount;++gn) {
								if (raw.cellConnectivity[raw.cellConnIndex[gg]+gn]==faceNode(f,fn).globalId) {
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
							face[f].bc=GHOST_FACE;
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
		for (int c=0;c<cellCount;++c) {
			int n;
			for (int cn=0;cn<cell[c].nodeCount;++cn) {
				n=cell[c].nodes[cn];
				for (int ng=0;ng<node[n].ghosts.size();++ng) {
					g=node[n].ghosts[ng];
					flag=false;
					for (int cg=0;cg<cell[c].ghosts.size();++cg) {
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
//	cout << "[I Rank=" << Rank << "] Number of Inter-Partition Ghost Cells= " << ghostCount << endl;

	globalNumFaceNodes=0;
	globalFaceCount=0;
	bool include;
	for (int f=0;f<faceCount;++f) {
		include=true;
		if (face[f].bc==GHOST_FACE) {
			int g=-1*face[f].neighbor-1;
			if (ghost[g].partition<Rank) include=false;
		}
		if (include) {
			globalNumFaceNodes+=face[f].nodeCount;
			globalFaceCount++;
		}
	}

	int temp=0;
	MPI_Allreduce (&globalNumFaceNodes,&temp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	globalNumFaceNodes=temp;
	temp=0;
	MPI_Allreduce (&globalFaceCount,&temp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	globalFaceCount=temp;

	return 0;
	
} // end int Grid::create_ghosts

int Grid::get_volume_output_ids() {
	
	vector<int> output_node_counts;
	output_node_counts.resize(np);
	int count=0;
	for (int n=0;n<nodeCount;++n) {
		node[n].output_id=0;
		for (int ng=0;ng<node[n].ghosts.size();++ng) {
			int g=node[n].ghosts[ng];
			if (ghost[g].partition<Rank) {
				node[n].output_id=-1;
				break;
			}
		}
		if (node[n].output_id==0) {
			node[n].output_id=count;
			count++;
		}
	}
	
	output_node_counts[Rank]=count;
	MPI_Allgather(MPI_IN_PLACE,0,MPI_INT, 
				  &output_node_counts[0],1,MPI_INT, 
				  MPI_COMM_WORLD);
	
	// A sanity check here:
	int sum=0;
	for (int p=0;p<np;++p) sum+=output_node_counts[p];
	if (sum!=globalNodeCount) {
		cerr << "[E] Output node counts sum doesn't match global node count" << endl;
		exit(1);
	}
	
	node_output_offset=0;
	for (int p=0;p<Rank;++p) node_output_offset+=output_node_counts[p];
	for (int n=0;n<nodeCount;++n) if (node[n].output_id>=0) node[n].output_id+=node_output_offset;
	
	for (int p=1;p<np;++p) {
		int size;
		for (int pr=0;pr<p;++pr) {
			if (Rank==p) {
				vector<int> exchange_nodes;
				for (int n=0;n<nodeCount;++n) if (node[n].output_id==-1) exchange_nodes.push_back(node[n].globalId);
				// Send the size of the list
				size=exchange_nodes.size();
				MPI_Send(&size,1,MPI_INT,pr,p,MPI_COMM_WORLD);
				// Send the list
				MPI_Send(&exchange_nodes[0],size,MPI_INT,pr,p,MPI_COMM_WORLD);
				MPI_Recv(&exchange_nodes[0],size,MPI_INT,pr,p,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				// Fill in the node output id's
				int count=0;
				for (int n=0;n<nodeCount;++n) {
					if (node[n].output_id==-1) {
						node[n].output_id=exchange_nodes[count];
						count++;
					}
				}
				exchange_nodes.clear();
			}
			if (Rank==pr) {
				// Receive the list size
				MPI_Recv(&size,1,MPI_INT,p,p,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				vector<int> exchange_nodes;
				exchange_nodes.resize(size);
				// Receive the list
				MPI_Recv(&exchange_nodes[0],size,MPI_INT,p,p,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				// Fill in the output_id's of exchange nodes
				for (int i=0;i<size;++i) {
					if (maps.nodeGlobal2Local.find(exchange_nodes[i])!=maps.nodeGlobal2Local.end()) {
						exchange_nodes[i]=node[maps.nodeGlobal2Local[exchange_nodes[i]]].output_id;	
					} else { exchange_nodes[i]=-1; }
				}
				// Send back the list
				MPI_Send(&exchange_nodes[0],size,MPI_INT,p,p,MPI_COMM_WORLD);
				exchange_nodes.clear();
			}
		}
	}
	
	// Populate face and node lists for each boundary condition region
	bcCount=raw.bocoNodes.size();
	boundaryFaces.resize(bcCount);
	boundaryNodes.resize(bcCount);
	vector<set<int> > bcnodeset;
	set<int>::iterator it;
	bcnodeset.resize(bcCount);
	for (int f=0;f<faceCount;++f) {
		if (face[f].bc>=0) {
			boundaryFaces[face[f].bc].push_back(f);
			for (int fn=0;fn<face[f].nodeCount;++fn) {
				bcnodeset[face[f].bc].insert(face[f].nodes[fn]);
			}
		}
	}
	
	for (int b=0;b<bcCount;++b) {
		int counter=0;
		for (it=bcnodeset[b].begin();it!=bcnodeset[b].end();++it) {
			boundaryNodes[b].push_back(*it);
			counter++;
		}
	}
	bcnodeset.clear();

	return 0;
}

int Grid::get_bc_output_ids() {
	
	// Collect all the nodes lying on any bc region into one set
	set<int> bc_nodes;
	set<int>::iterator sit;
	int n;
	
	for (int f=0;f<faceCount;++f) {
		if (face[f].bc>=0) {
			for (int fn=0;fn<face[f].nodeCount;++fn) {
				bc_nodes.insert(face[f].nodes[fn]);
			}
		}
	}
	
	// Set all node output id's to -2 by default
	for (n=0;n<nodeCount;++n) node[n].bc_output_id=-2;
	
	vector<int> output_node_counts;
	output_node_counts.resize(np);
	int count=0;
	for (sit=bc_nodes.begin();sit!=bc_nodes.end();sit++) {
		n=*sit;
		node[n].bc_output_id=0;
		for (int ng=0;ng<node[n].ghosts.size();++ng) {
			int g=node[n].ghosts[ng];
			if (ghost[g].partition<Rank) {
				node[n].bc_output_id=-1;
				break;
			}
		}
		if (node[n].bc_output_id==0) {
			node[n].bc_output_id=count;
			count++;
		}
	}
	
	output_node_counts[Rank]=count;
	MPI_Allgather(MPI_IN_PLACE,0,MPI_INT, 
				  &output_node_counts[0],1,MPI_INT, 
				  MPI_COMM_WORLD);
	
	node_bc_output_offset=0;
	for (int p=0;p<Rank;++p) node_bc_output_offset+=output_node_counts[p];

	global_bc_nodeCount=0;
	for (int p=0;p<np;++p) global_bc_nodeCount+=output_node_counts[p];

	for (sit=bc_nodes.begin();sit!=bc_nodes.end();sit++) if (node[*sit].bc_output_id>=0) node[*sit].bc_output_id+=node_bc_output_offset;
	
	for (int p=1;p<np;++p) {
		int size;
		for (int pr=0;pr<p;++pr) {
			if (Rank==p) {
				vector<int> exchange_nodes;
				for (sit=bc_nodes.begin();sit!=bc_nodes.end();sit++) if (node[*sit].bc_output_id==-1) exchange_nodes.push_back(node[*sit].globalId);
				// Send the size of the list
				size=exchange_nodes.size();
				MPI_Send(&size,1,MPI_INT,pr,p,MPI_COMM_WORLD);
				// Send the list
				MPI_Send(&exchange_nodes[0],size,MPI_INT,pr,p,MPI_COMM_WORLD);
				MPI_Recv(&exchange_nodes[0],size,MPI_INT,pr,p,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				// Fill in the node output id's
				int count=0;
				for (sit=bc_nodes.begin();sit!=bc_nodes.end();sit++) {
					if (node[*sit].bc_output_id==-1) {
						node[*sit].bc_output_id=exchange_nodes[count];
						count++;
					}
				}
				exchange_nodes.clear();
			}
			if (Rank==pr) {
				// Receive the list size
				MPI_Recv(&size,1,MPI_INT,p,p,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				vector<int> exchange_nodes;
				exchange_nodes.resize(size);
				// Receive the list
				MPI_Recv(&exchange_nodes[0],size,MPI_INT,p,p,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				// Fill in the bc_output_id's of exchange nodes
				for (int i=0;i<size;++i) {
					if (maps.nodeGlobal2Local.find(exchange_nodes[i])!=maps.nodeGlobal2Local.end()) {
						exchange_nodes[i]=node[maps.nodeGlobal2Local[exchange_nodes[i]]].bc_output_id;	
					} else { exchange_nodes[i]=-1; }
				}
				// Send back the list
				MPI_Send(&exchange_nodes[0],size,MPI_INT,p,p,MPI_COMM_WORLD);
				exchange_nodes.clear();
			}
		}
	}
	
	return 0;
}
	
bool Cell::HaveNodes(int &nodelistsize, int nodelist []) {	
	bool match;
	for (int i=0;i<nodelistsize;++i) {
		match=false;
		for (int j=0;j<nodeCount;++j) if (nodelist[i]==nodes[j]) {match=true; break;}
		if (!match) return false;
	}
	return true;
}
