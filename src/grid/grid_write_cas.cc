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
#include "bc.h" 
#include "inputs.h"

extern vector<vector<BCregion> > bc;
extern InputFile input;

void Grid::write_cas(void) {
	// Open file
	ofstream file;
	file.open("grid.cas",ios::out);
	file << "(1 \"FreeCFD\")" << endl;
	file << "(2 3)" << endl;
	// Write number of nodes
	file << "(10 (0 1 " << hex << globalNodeCount << " 1 3))" << endl;
	// Write nodes
	file << "(10 (1 1 " << hex << globalNodeCount << " 1 3)(" << endl;
	for (int n=0;n<nodeCount;++n) {
		file << node[n][0] << "\t" << node[n][1] << "\t" << node[n][2] << endl;
	}
	file << "))" << endl;
	// Write number of cells
	file << "(12 (0 1 " << hex << globalCellCount << " 0))" << endl;
	// Write cells
	file << "(12 (1 1 " << hex << globalCellCount << " 1 0)(" << endl;
	for (int c=0;c<cellCount;++c) {
		if (cell[c].nodeCount==4) { // TETRA
			file << hex << 2  << endl;
		} else if (cell[c].nodeCount==5) { // PYRAMID
			file << hex << 5 << endl;
		} else if (cell[c].nodeCount==6) { // WEDGE
			file << hex << 6 << endl;
		} else if (cell[c].nodeCount==8) { // HEXA
			file << hex << 4 << endl;
		} else {
			cerr << "Unsupported cell type" << endl;
			exit(1);
		}
	}
	file << "))" << endl;
	// Write number of faces
	file << "(13 (0 1 " << hex << face.size() << " 0))" << endl;
	// Write faces
	// Loop faces to get interior faces
	vector<int> indices;
	for (int f=0;f<faceCount;++f) {
		if (face[f].bc==INTERNAL_FACE) {
			indices.push_back(f);
		}
	}
	int zone=0;
	if (indices.size()!=0) {
		zone++;
		file << "(13 (" << zone << " 1 " << hex << indices.size() << " " << hex << 2 << " 0)(" << endl;
		for (int i=0;i<indices.size();++i) {
				file << hex << face[indices[i]].nodeCount << "\t";
			for (int fn=face[indices[i]].nodeCount-1;fn>=0;--fn) {
				file << hex << face[indices[i]].nodes[fn]+1 << "\t";
			}
			// Parent and neighbor cells
			if (face[indices[i]].parent<0 || face[indices[i]].neighbor<0) {
				cerr << "ERRORRR" << endl;
				exit(1);
			}
			file << hex << face[indices[i]].parent+1 << "\t" << hex << face[indices[i]].neighbor+1 << endl;
		}
		file << "))" << endl;
	}
	
	// Loop faces to get wall faces
	indices.clear();
	for (int f=0;f<faceCount;++f) {
		if (bc[gid][face[f].bc].type==WALL) {
			indices.push_back(f);
		}
	}
	if (indices.size()!=0) {
		zone++;
		file << "(13 (" << zone << " 1 " << hex << indices.size() << " " << hex << 3 << " 0)(" << endl;
		for (int i=0;i<indices.size();++i) {
			file << hex << face[indices[i]].nodeCount << "\t";
			for (int fn=face[indices[i]].nodeCount-1;fn>=0;--fn) {
				file << hex << face[indices[i]].nodes[fn]+1 << "\t";
			}
			// Parent and neighbor cells
			file << hex << face[indices[i]].parent+1 << "\t" << 0 << endl;
		}
		file << "))" << endl;
	}
	
	// Loop faces to get symmetry faces
	indices.clear();
	for (int f=0;f<faceCount;++f) {
		if (bc[gid][face[f].bc].type==SYMMETRY) {
			indices.push_back(f);
		}
	}
	if (indices.size()!=0) {
		zone++;
		file << "(13 (" << zone << " 1 " << hex << indices.size() << " " << hex << 7 << " 0)(" << endl;
		for (int i=0;i<indices.size();++i) {
			file << hex << face[indices[i]].nodeCount << "\t";
			for (int fn=face[indices[i]].nodeCount-1;fn>=0;--fn) {
				file << hex << face[indices[i]].nodes[fn]+1 << "\t";
			}
			// Parent and neighbor cells
			file << hex << face[indices[i]].parent+1 << "\t" << 0 << endl;
		}
		file << "))" << endl;
	}
	
	// Loop faces to get velocity inlet faces
	indices.clear();
	for (int f=0;f<faceCount;++f) {
		if (bc[gid][face[f].bc].type==INLET && input.section("grid",gid).subsection("BC",face[f].bc).get_Vec3D("V").is_found) {
			indices.push_back(f);
		}
	}
	if (indices.size()!=0) {
		zone++;
		file << "(13 (" << zone << " 1 " << hex << indices.size() << " " << hex << 10 << " 0)(" << endl;
		for (int i=0;i<indices.size();++i) {
			file << hex << face[indices[i]].nodeCount << "\t";
			for (int fn=face[indices[i]].nodeCount-1;fn>=0;--fn) {
				file << hex << face[indices[i]].nodes[fn]+1 << "\t";
			}
			// Parent and neighbor cells
			file << hex << face[indices[i]].parent+1 << "\t" << 0 << endl;
		}
		file << "))" << endl;
	}
	
	// Loop faces to get mass flow inlet faces
	indices.clear();
	for (int f=0;f<faceCount;++f) {
		if (bc[gid][face[f].bc].type==INLET && !(input.section("grid",gid).subsection("BC",face[f].bc).get_Vec3D("V").is_found) ) {
			indices.push_back(f);
		}
	}
	if (indices.size()!=0) {
		zone++;
		file << "(13 (" << zone << " 1 " << hex << indices.size() << " " << hex << 20 << " 0)(" << endl;
		for (int i=0;i<indices.size();++i) {
			file << hex << face[indices[i]].nodeCount << "\t";
			for (int fn=face[indices[i]].nodeCount-1;fn>=0;--fn) {
				file << hex << face[indices[i]].nodes[fn]+1 << "\t";
			}
			// Parent and neighbor cells
			file << hex << face[indices[i]].parent+1 << "\t" << 0 << endl;
		}
		file << "))" << endl;
	}
	
	 // Loop faces to get outlet faces
	 indices.clear();
	 for (int f=0;f<faceCount;++f) {
	 if (bc[gid][face[f].bc].type==OUTLET) {
	 indices.push_back(f);
	 }
	 }
	 if (indices.size()!=0) {
	 zone++;
	 file << "(13 (" << zone << " 1 " << hex << indices.size() << " " << hex << 36 << " 0)(" << endl;
	 for (int i=0;i<indices.size();++i) {
	 file << hex << face[indices[i]].nodeCount << "\t";
	 for (int fn=face[indices[i]].nodeCount-1;fn>=0;--fn) {
	 file << hex << face[indices[i]].nodes[fn]+1 << "\t";
	 }
	 // Parent and neighbor cells
	 file << hex << face[indices[i]].parent+1 << "\t" << 0 << endl;
	 }
	 file << "))" << endl;
	 }
	 
	file.close();
	return;
}
