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

#ifndef VARIABLE
#define VARIABLE

#include "grid.h"

extern vector<Grid> grid;

template <class TYPE> 
class Variable {
public:
	int gid; // grid ID
	vector<bool> fixedonBC; // whether this variable is fixed on a certain BC of this grid
	vector<vector<TYPE> > bcValue; // Stores the bc data if specified on a certain bc
	vector<TYPE> cellData, faceData, nodeData, ghostData;
	bool cellStore, faceStore, nodeStore, ghostStore;
	TYPE temp;
	map<int,double>::iterator it;
	set<int>::iterator sit;
	// Function pointers
	// Store addresses of functions to be used when data is requested
	// Can be simple fetch from array if the variable is stored
	// Or can be an on-demand evaluation function
	// That Variable:: is important here
	TYPE & (Variable::*get_cell)(int); // get_cell is a pointer to a function that returns address of TYPE
	TYPE & (Variable::*get_ghost)(int);
	TYPE & (Variable::*get_face)(int);
	TYPE & (Variable::*get_node)(int);
	
	// Functions
	Variable (void) { 
		/* Set default choices */ 
		cellStore=true; 
		ghostStore=true; 
		faceStore=false; 
		nodeStore=false; 
		return;
	}
	void allocate (int g);
	TYPE &cell (int c); // This will invoke either fetch or the calculate function depending on what get_cell above points to
	TYPE &cell_fetch (int c);
	TYPE &cell_calculate (int c);
	
	TYPE &ghost (int g); // This will invoke either fetch or the calculate function depending on what get_ghost above points to
	TYPE &ghost_fetch (int g);
	TYPE &ghost_calculate (int g);
	
	TYPE &face (int f); // This will invoke either fetch or the calculate function depending on what get_face above points to
	TYPE &face_fetch (int f);
	TYPE &face_calculate (int f);
	
	TYPE &node (int n); // This will invoke either fetch or the calculate function depending on what get_node above points to
	TYPE &node_fetch (int n);
	TYPE &node_calculate (int n);
	
	TYPE &bc (int b,int f=-1);
	TYPE cell2node (int c, int n);
	
	vector<TYPE> cell_gradient (int c);
	void mpi_update(void);
	void dump_cell_data(string fileName);
	void read_cell_data(string fileName,vector<vector<int> > &partitionMap);
	
};

// There can't be a implementation .cc file for a template class. So, all the implementation is here in the header
template <class TYPE> 
void Variable<TYPE>::allocate (int g) {
	gid=g;
	if (cellStore) {
		cellData.resize(grid[gid].cellCount);
		get_cell=&Variable::cell_fetch;
	}
	if (ghostStore) {
		ghostData.resize(grid[gid].ghostCount);
		get_ghost=&Variable::ghost_fetch;
	} 
	if (faceStore) {
		faceData.resize(grid[gid].faceCount);		
		get_face=&Variable::face_fetch;
	} else {
		get_face=&Variable::face_calculate;
	}
	if (nodeStore) {
		nodeData.resize(grid[gid].nodeCount);
		get_node=&Variable::node_fetch;
	} else {
		get_node=&Variable::node_calculate;
	}
	
	fixedonBC.resize(grid[gid].bcCount);
	bcValue.resize(grid[gid].bcCount);
	for (int i=0;i<fixedonBC.size();++i) fixedonBC[i]=false;
	
	return;
}

template <class TYPE>
TYPE &Variable<TYPE>::cell (int c) { 
	// Call whatever get_cell is pointing to
	return (this->*get_cell)(c);
}

template <class TYPE>
TYPE &Variable<TYPE>::cell_fetch (int c) { 
	// Just return the data in cellData
	return cellData[c];
}

template <class TYPE>
TYPE &Variable<TYPE>::ghost (int g) { 
	// Call whatever get_ghost is pointing to
	return (this->*get_ghost)(g);
}

template <class TYPE>
TYPE &Variable<TYPE>::ghost_fetch (int g) { 
	// Just return the data in ghostData
	return ghostData[g];
}

template <class TYPE>
TYPE &Variable<TYPE>::face (int f) { 
	// Call whatever get_face is pointing to
	return (this->*get_face)(f);
}

template <class TYPE>
TYPE &Variable<TYPE>::face_fetch (int f) { 
	// Just return the data in faceData
	return faceData[f];
}

template <class TYPE>
TYPE &Variable<TYPE>::face_calculate (int f) { 
	if ( (grid[gid].face[f].bc>=0 && fixedonBC[grid[gid].face[f].bc]) || !cellStore) {		
		return bc(grid[gid].face[f].bc,f);
	}
	// Run the face averaging map from the grid class
	std::map<int,double>::iterator it;
	temp=0.;
	it=grid[gid].face[f].average.begin();
	for ( it=grid[gid].face[f].average.begin() ; it != grid[gid].face[f].average.end(); it++ ) {
		// Note: cell((*it).first) could be used but the get_cell method should be slightly faster
		if ((*it).first>=0) { // if contribution is coming from a real cell
			temp+=(*it).second*(this->*get_cell)((*it).first);
		} else { // if contribution is coming from a ghost cell
			temp+=(*it).second*(this->*get_ghost)(-1*((*it).first+1));
		}
	}
	return temp;
}

template <class TYPE>
TYPE &Variable<TYPE>::node (int n) { 
	// Call whatever get_node is pointing to
	return (this->*get_node)(n);
}

template <class TYPE>
TYPE &Variable<TYPE>::node_fetch (int n) { 
	// Just return the data in nodeData
	return nodeData[n];
}

template <class TYPE>
TYPE &Variable<TYPE>::node_calculate (int n) { 

	temp=0.;
	int count=0;
	// Loop node neighbor faces
	for (vector<int>::iterator itr=grid[gid].node[n].faces.begin(); itr!=grid[gid].node[n].faces.end(); itr++) {
		int fbc=grid[gid].face[*itr].bc;
		if (fbc>=0) {
			if (bcValue[fbc].size()==1) { // Uniform bc
				temp+=bcValue[fbc][0];
				count++;
			} else if (bcValue[fbc].size()>1) { // Distributed bc
				// Use cell gradient to get to the node
				temp+=cell2node(grid[gid].face[*itr].parent,n);
				count++;
			}
		}
	}
	if (count>0) {
		temp/=double(count);
		return temp;
	}
	// TODO: Wait until a general interpolation routine is available
	// Then check neighbor face bc's and interpolate to this node
	//if (!grid.node[n].bcs.empty()) {
	//	for (sit=grid.node[n].bcs.begin();sit!=raw.grid.node[n].bcs.end();sit++) {
	//	}
	//} else {
		// Run the node averaging map from the grid class
		temp=0.;
		it=grid[gid].node[n].average.begin();
		for ( it=grid[gid].node[n].average.begin() ; it != grid[gid].node[n].average.end(); it++ ) {
			// Note: cell((*it).first) could be used but the get_cell method should be slightly faster
			if ((*it).first>=0) { // if contribution is coming from a real cell
				temp+=(*it).second*(this->*get_cell)((*it).first);
			} else { // if contribution is coming from a ghost cell
				temp+=(*it).second*(this->*get_ghost)(-1*((*it).first+1));
			}
		}
		return temp;
	//}
}

template <class TYPE>
vector<TYPE> Variable<TYPE>::cell_gradient (int c) {
	
	map<int,Vec3D>::iterator it;
	map<int,double>::iterator fit;
	int f;
	Vec3D areaVec;
	vector<TYPE> grad (3,0.);
	// TODO: Calculate gradient using face average map
	//       Incorporate BC values (maybe optional)
	//       Get rid of the gradMap in the grid class
	
	/*
	// Add internal and interpartition face contributions
	for (it=grid[gid].cell[c].gradMap.begin();it!=grid[gid].cell[c].gradMap.end(); it++ ) {
		if ((*it).first>=0) { // if contribution is coming from a real cell
			for (int i=0;i<3;++i) grad[i]+=(*it).second[i]*(this->*get_cell)((*it).first);
		} else { // if contribution is coming from a ghost cell
			for (int i=0;i<3;++i) grad[i]+=(*it).second[i]*(this->*get_ghost)(-1*((*it).first+1));
		}
	} // end gradMap loop
	*/
	
	// The grad map loop above doesn't count the boundary faces
	// Add boundary face contributions
	for (int cf=0;cf<grid[gid].cell[c].faceCount;++cf) {
		f=grid[gid].cell[c].faces[cf];		
		//if (grid[gid].face[f].bc>=0) { // if a boundary face
			areaVec=grid[gid].face[f].normal*grid[gid].face[f].area/grid[gid].cell[c].volume;
			if (grid[gid].face[f].parent!=c) areaVec*=-1.;
			for (int i=0;i<3;++i) grad[i]+=(this->*get_face)(f)*areaVec[i];				
		//} // end if a boundary face
	} // end cell face loop
	
	return grad;
}

template <class TYPE>
TYPE &Variable<TYPE>::bc (int b,int f) {
	// TODO: For now
	if (bcValue[b].size()==1) {
		return bcValue[b][0];
	} else if (f>=0) {
		return bcValue[b][grid[gid].maps.face2bc[f]]; 
	}
}

template <class TYPE>
void Variable<TYPE>::dump_cell_data(string fileName) {

	ofstream file;
	// Write variables
	int size=sizeof(TYPE);
	for (int p=0;p<grid[gid].np;++p) {
		if (grid[gid].Rank==p) {
			if (grid[gid].Rank==0) file.open(fileName.c_str(),ios::binary); 
			else file.open(fileName.c_str(), ios::app | ios::binary);
			
			for (int c=0;c<grid[gid].cellCount;++c) file.write((char*) &cell(c),size);
			file.close();
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	return;
}

template <class TYPE>
void Variable<TYPE>::read_cell_data(string fileName,vector<vector<int> > &partitionMap) {
	
	string data="";
	TYPE dummy;
	ifstream file;
	file.open(fileName.c_str(),ios::in | ios::binary);

	int id;
	int size=sizeof(TYPE);
	for (int p=0;p<partitionMap.size();++p) {
		for (int c=0;c<partitionMap[p].size();++c) {
			// Get local cell id
			id=-1;
			if (grid[gid].maps.cellGlobal2Local.find(partitionMap[p][c])!=grid[gid].maps.cellGlobal2Local.end()) {
				id=grid[gid].maps.cellGlobal2Local[partitionMap[p][c]];
			}
			// If id is negative, that means the cell currently lies on another partition
			if (id>=0) { file.read((char*) &cell(id),size); } else { file.read((char*) &dummy,size); }
		}
	}
	file.close();

	return;
}

#endif













