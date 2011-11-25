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
#include "kdtree.h"

void strip_white_spaces(string &data);
void strip_quotes(string &data);
bool is_first_numeric(ifstream &file);
extern void StringExplode(string str, string separator, vector<string>* results);
extern GridRawData raw;
/*
 * Reads a tecplot grid file in FEPolhedron format
 *
 */

int Grid::readTEC() {

	raw.type=FACE;

	int TotalNumFaceNodes,NumConnectedBoundaryFaces,TotalNumBoundaryConnections;
	string line;

	// Open the grid file
	ifstream file;
	file.open(fileName.c_str(),ios::in | ios::binary);	

	while (true)  {
		// Break the loop if first non-whitespace character of the line is a number or "-"
		if (is_first_numeric(file)) break;
		
		// Grab a line from the file
		if (!getline(file,line)) {
			cerr << "[E] Error reading grid file: " << fileName << endl;
			exit(1);
		}
		// Strip white spaces in the line
		strip_white_spaces(line);
		// Split the line into expressions separated with commas
		vector<string> expression;
		StringExplode(line,",",&expression);
		// Each expressions is in a=b or a="b" format. Now get a and b for each
		vector<string> temp;
		for (int i=0;i<expression.size();++i) {
			StringExplode(expression[i],"=",&temp);
			if (temp.size()==2) { // If a=b is  found
				if (temp[0]=="Nodes") {
					globalNodeCount=atoi(temp[1].c_str());
					if (Rank==0) cout << "[I] Global node count = " << globalNodeCount << endl;
				} else if (temp[0]=="Elements") {
					globalCellCount=atoi(temp[1].c_str());
					if (Rank==0) cout << "[I] Global cell count = " << globalCellCount << endl;
				} else if (temp[0]=="Faces") {
					globalFaceCount=atoi(temp[1].c_str());
					if (Rank==0) cout << "[I] Global face count = " << globalFaceCount << endl;
				}  else if (temp[0]=="ZONETYPE") {
					if (temp[1]!="FEPolyhedron") {
						cerr << "[E] ZONETYPE=" << temp[1] << " is not supported" << endl;
						exit(1);
					}
				}  else if (temp[0]=="DATAPACKING") {
					if (temp[1]!="BLOCK") {
						cerr << "[E] DATAPACKING=" << temp[1] << " is not supported" << endl;
						exit(1);
					}
				} else if (temp[0]=="TotalNumFaceNodes") {
					TotalNumFaceNodes=atoi(temp[1].c_str());
					if (Rank==0) cout << "[I] TotalNumFaceNodes = " << TotalNumFaceNodes << endl;
				} else if (temp[0]=="NumConnectedBoundaryFaces") {
					NumConnectedBoundaryFaces=atoi(temp[1].c_str());
					if (NumConnectedBoundaryFaces!=0) {
						cerr << "[E] NumConnectedBoundaryFaces>0 is not supported" << endl;
						exit(1);
					}
				} else if (temp[0]=="TotalNumBoundaryConnections") {
					TotalNumBoundaryConnections=atoi(temp[1].c_str());
					if (TotalNumBoundaryConnections!=0) {
						cerr << "[E] TotalNumBoundaryConnections>0 is not supported" << endl;
						exit(1);
					}
				}
			}
			temp.clear();
		}
		expression.clear();
	}
	// Now the tecplot file header is read
	// Read the node coordinates
	if (Rank==0) cout << "[I] Reading node coordinates" << endl;
	raw.node.resize(globalNodeCount);
	for (int i=0;i<3;++i) for (int n=0;n<globalNodeCount;++n) file >> raw.node[n][i];

	// Skip the comment line
	if (!is_first_numeric(file)) getline(file,line);
	// Read number of nodes for each face
	if (Rank==0) cout << "[I] Reading node counts for each face" << endl;
	raw.faceNodeCount.resize(globalFaceCount);
	raw.faceConnIndex.resize(globalFaceCount);
	int faceNodeListSize=0;
	for (int f=0;f<globalFaceCount;++f) {
		file >> raw.faceNodeCount[f];
		raw.faceConnIndex[f]=faceNodeListSize;
		faceNodeListSize+=raw.faceNodeCount[f];
	}


	// Skip the comment line
	if (!is_first_numeric(file)) getline(file,line);
	// Read the node list for each face 
	if (Rank==0) cout << "[I] Reading node lists for each face" << endl;
	raw.faceConnectivity.resize(faceNodeListSize);
	for (int i=0;i<faceNodeListSize;++i) file >> raw.faceConnectivity[i];
	// Reduce indexing to start from 0
	for (int i=0;i<faceNodeListSize;++i) raw.faceConnectivity[i]--;

	// Read the left and right neighbors for each face
	// Skip the comment line
	if (!is_first_numeric(file)) getline(file,line);
	if (Rank==0) cout << "[I] Reading the left neighbors for each face" << endl;
	raw.left.resize(globalFaceCount);
	for (int f=0;f<globalFaceCount;++f) file >> raw.left[f]; 
	if (!is_first_numeric(file)) getline(file,line);
	if (Rank==0) cout << "[I] Reading the right neighbors for each face" << endl;
	raw.right.resize(globalFaceCount);
	for (int f=0;f<globalFaceCount;++f) file >> raw.right[f]; 
	// Reduce indexing to start from 0
	for (int f=0;f<globalFaceCount;++f) {raw.right[f]--; raw.left[f]--;}

	// Create a kdtree and insert all the nodes that lie in the boundaries
	kdtree *kd;
	kdres *kdres; // results
	int *data;
	kd = kd_create(3);	

	if (Rank==0) cout << "[I] Inserting boundary nodes to kdtree" << endl;	
	for (int f=0;f<globalFaceCount;++f) {
		if (raw.right[f]<0 || raw.left[f]<0) { // Boundary face
			// Insert all nodes of the face to kdtree
			int n;
			for (int fn=0;fn<raw.faceNodeCount[f];++fn) {
				n=raw.faceConnectivity[raw.faceConnIndex[f]+fn];
				data=&raw.faceConnectivity[raw.faceConnIndex[f]+fn];
				kd_insert3(kd,raw.node[n][0],raw.node[n][1],raw.node[n][2],data);	
			}
		}
	}

	// Read Each Boundary Zone
	int bc_count=0;
	bool line_grabbed=false;
	while (true) {
		int bcNodeCount=0;
		// Read header
		while (true)  {
			if (is_first_numeric(file)) break;
			// Grab a line from the file
			if (!line_grabbed) {
				if (!getline(file,line)) {
					cerr << "[E] Error reading grid file: " << fileName << endl;
					exit(1);
				}
			}
			line_grabbed=false;
			// Strip white spaces in the line
			strip_white_spaces(line);
			// Split the line into expressions separated with commas
			vector<string> expression;
			StringExplode(line,",",&expression);
			// Each expressions is in a=b or a="b" format. Now get a and b for each
			vector<string> temp;
			for (int i=0;i<expression.size();++i) {
				StringExplode(expression[i],"=",&temp);
				if (temp.size()==2) { // If a=b is  found
					if (temp[0]=="ZONET") {
						strip_quotes(temp[1]);
						raw.bocoNameMap.insert(pair<string,int>(temp[1],bc_count));
						bc_count++;
						if (Rank==0) cout << "[I] Found BC_" << bc_count << " :" << temp[1] << endl;
					} else if (temp[0]=="Nodes") {
						bcNodeCount=atoi(temp[1].c_str());
						if (Rank==0) cout << "[I] BC_" << bc_count << " node count = " << bcNodeCount << endl;
					}
				}
				temp.clear();
			}
			expression.clear();
		}
	
		// Now the BC zone header is read
		// Read the BC node coordinates
		if (Rank==0) cout << "[I] Reading BC_" << bc_count << " node coordinates" << endl;
		vector<Vec3D> bcnodes;
		bcnodes.resize(bcNodeCount);
		for (int i=0;i<3;++i) for (int n=0;n<bcNodeCount;++n) file >> bcnodes[n][i];
		
		// Search kdtree to find global node number
		if (Rank==0) cout << "[I] Searching kdtree to find BC_" << bc_count << " global node numbers" << endl;
		raw.bocoNodes.resize(bc_count);
		for (int n=0;n<bcNodeCount;++n) {
			kdres = kd_nearest3(kd,bcnodes[n][0],bcnodes[n][1],bcnodes[n][2]);
			data=(int *)kd_res_item_data(kdres);
			raw.bocoNodes[bc_count-1].insert(*data);	
		}
		// Empty bcnodes array
		bcnodes.clear();
		// Seek until another ZONE statement or EOF is reached
		while (true) {
			getline(file,line);
			strip_white_spaces(line);
			line_grabbed=true;
			if (line.substr(0,4)=="ZONE") break;
			if (file.eof()) break;
		}
		
		if (file.eof()) break;
	}

	// Destroy the kdtree
	kd_res_free(kdres);
	kd_free(kd);
	// Close the grid file
	file.close();

	globalNumFaceNodes=raw.faceConnectivity.size();

	// Now fill in cell connectivity (needed for partitioning)
	// Declare a temporary set of nodes for each cell
	vector<set<int> > cellnodes;
	cellnodes.resize(globalCellCount);
	for (int f=0;f<globalFaceCount;++f) {
		int ibegin=raw.faceConnIndex[f];
		int iend;
		if (f<globalFaceCount) iend=raw.faceConnIndex[f+1]-1;
		else iend=raw.faceConnectivity.size();
		for (int i=ibegin;i<=iend;++i) {
			if (raw.left[f]>=0)  cellnodes[raw.left[f]].insert(raw.faceConnectivity[i]);
			if (raw.right[f]>=0) cellnodes[raw.right[f]].insert(raw.faceConnectivity[i]);
		}
	}
	
	// Count the entries needed for cell connectivities
	int cell_conn_count=0;
	for (int c=0;c<globalCellCount;++c) cell_conn_count+=cellnodes[c].size();

	raw.cellConnectivity.resize(cell_conn_count);
	raw.cellConnIndex.resize(globalCellCount);
	// Fill in the cell connectivity array
	set<int>::iterator sit;
	int count=0;
	for (int c=0;c<globalCellCount;++c) {
		raw.cellConnIndex[c]=count;
		for (sit=cellnodes[c].begin();sit!=cellnodes[c].end();sit++) {
			raw.cellConnectivity[count]=*sit;
			count++;
		}
	}

	return 1;
}

void strip_white_spaces(string &data) {
        string whitespaces(" \t\f\v\n\r");
        size_t found=0;
        while (found!=string::npos) {
                found=data.find_first_of(whitespaces);
                if (found!=string::npos) data.erase(found,1);
        }
	return;
}

void strip_quotes(string &data) {
        string quotes("\"\'");
        size_t found=0;
        while (found!=string::npos) {
                found=data.find_first_of(quotes);
                if (found!=string::npos) data.erase(found,1);
        }
	return;
}

bool is_first_numeric(ifstream &file) {

	char c='\0';
	int pos=0;
	pos=file.tellg();
	c=file.peek();
	while(isspace(c)) {
		file.seekg(pos++);
		c=file.peek();
	}
	if (c=='-' || isdigit(c)) {
		return true;
	}

	return false;
}

