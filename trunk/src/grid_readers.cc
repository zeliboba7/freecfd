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
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <iomanip>
#include<cmath>
#include <algorithm>
using namespace std;
#include <cgnslib.h>

#include "grid.h"
string int2str(int number) ;
extern Grid grid;
extern int np, Rank;

GridRawData raw;
double block_stitch_tolerance=1.e-7;

int Grid::readCGNS() {
// parallel:OK

	int fileIndex,baseIndex,zoneIndex,sectionIndex,nBases,nZones,nSections,nBocos;
	char zoneName[20],sectionName[20]; //baseName[20]
	//  int nBases,cellDim,physDim;
	int size[3];
	globalNodeCount=0;
	globalCellCount=0;
	globalFaceCount=0;

	// Open the grid file for reading
	cg_open(fileName.c_str(),MODE_READ,&fileIndex);

	// Read number of bases
	cg_nbases(fileIndex,&nBases);

	// Number of bases is typically 1
	if (Rank==0) cout << "[I] Number of Bases= " << nBases << endl;

	//for (int baseIndex=1;baseIndex<=nBases;++baseIndex) {
	// TODO For now assuming there is only one base as I don't know in which cases there would be more
	baseIndex=1;
	// Read number of zones (number of blocks in the grid)
	cg_nzones(fileIndex,baseIndex,&nZones);
	int zoneNodeCount[nZones],zoneCellCount[nZones];
	if (Rank==0) cout << "[I] Number of Zones= " << nZones << endl;
	
	std::vector<double> coordX[nZones],coordY[nZones],coordZ[nZones];
	// zoneCoordMap returns the global is of a given node in a given zone
	std::vector<int> zoneCoordMap[nZones];

	int bocoCount=0;
	map<string,int>::iterator mit;
	for (int zoneIndex=1;zoneIndex<=nZones;++zoneIndex) { // For each zone
		// Read the zone
		cg_zone_read(fileIndex,baseIndex,zoneIndex,zoneName,size);
		// These are the number of cells and nodes in that zone
		zoneNodeCount[zoneIndex-1]=size[0];
		zoneCellCount[zoneIndex-1]=size[1];
		// Read number of sections
		cg_nsections(fileIndex,baseIndex,zoneIndex,&nSections);
		cg_nbocos(fileIndex,baseIndex,zoneIndex,&nBocos);
		if (Rank==0) cout << "[I] In Zone " << zoneName << endl;
		if (Rank==0) cout << "[I] ...Number of Nodes= " << size[0] << endl;
		if (Rank==0) cout << "[I] ...Number of Cells= " << size[1] << endl;
		if (Rank==0) cout << "[I] ...Number of Sections= " << nSections << endl;
		if (Rank==0) cout << "[I] ...Number of Boundary Conditions= " << nBocos << endl;

		// Read the node coordinates
		int nodeStart[3],nodeEnd[3];
		nodeStart[0]=nodeStart[1]=nodeStart[2]=1;
		nodeEnd[0]=nodeEnd[1]=nodeEnd[2]=size[0];

		coordX[zoneIndex-1].resize(size[0]);
		coordY[zoneIndex-1].resize(size[0]);
		coordZ[zoneIndex-1].resize(size[0]);
		cg_coord_read(fileIndex,baseIndex,zoneIndex,"CoordinateX",RealDouble,nodeStart,nodeEnd,&coordX[zoneIndex-1][0]);
		cg_coord_read(fileIndex,baseIndex,zoneIndex,"CoordinateY",RealDouble,nodeStart,nodeEnd,&coordY[zoneIndex-1][0]);
		cg_coord_read(fileIndex,baseIndex,zoneIndex,"CoordinateZ",RealDouble,nodeStart,nodeEnd,&coordZ[zoneIndex-1][0]);
		for (int i=0;i<size[0];++i) {
			// This map takes in zone index and node number within that zone
			// Returns the global index of that node
			// It is initialized to -1 for now
			zoneCoordMap[zoneIndex-1].push_back(-1);
		}

		// In case there are multiple connected zones, collapse the repeated nodes and fix the node numbering
		if (zoneIndex==1) { // If the first zone
			for (int c=0;c<coordX[0].size();++c) {
				// Global node count is incremented as new nodes are found.
				// When in the first zone, every node is new.
				zoneCoordMap[0][c]=globalNodeCount;
				globalNodeCount++;
			}
		} else {
			// Scan the coordinates of all the other zones before this one for duplicates	
			for (int c=0;c<coordX[zoneIndex-1].size();++c) {
				bool foundFlag=false;
				for (int z=0;z<zoneIndex-1;++z) {
					for (int c2=0;c2<coordX[z].size();++c2) {
						if (fabs(coordX[zoneIndex-1][c]-coordX[z][c2])<block_stitch_tolerance && fabs(coordY[zoneIndex-1][c]-coordY[z][c2])<block_stitch_tolerance && fabs(coordZ[zoneIndex-1][c]-coordZ[z][c2])<block_stitch_tolerance) {
							zoneCoordMap[zoneIndex-1][c]=zoneCoordMap[z][c2];
							foundFlag=true;
							break;
						}
					}
					if (foundFlag) break;
				}
				if (!foundFlag) {
					zoneCoordMap[zoneIndex-1][c]=globalNodeCount;
					globalNodeCount++;
				}
			}
		}

		// Read element ranges for all the boundary condition regions within the current zone
		int bc_range[nBocos][2];
		for (int bocoIndex=1;bocoIndex<=nBocos;++bocoIndex) {
			int dummy;
			cg_boco_read(fileIndex,baseIndex,zoneIndex,bocoIndex,bc_range[bocoIndex-1],&dummy);
		} // for boco

		// Loop sections within the zone
		// These include connectivities of cells and bonudary faces
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
				case PYRA_5:
					elemNodeCount=5; break;
				case PENTA_6:
					elemNodeCount=6; break;
				case HEXA_8:
					elemNodeCount=8; break;
			} //switch
			int elemNodes[elemEnd-elemStart+1][elemNodeCount];

			// Read element node connectivities
			cg_elements_read(fileIndex,baseIndex,zoneIndex,sectionIndex,*elemNodes,0);
			
			// Only pick the volume elements
			if (elemType==TETRA_4 || elemType==PYRA_5 || elemType==PENTA_6 || elemType==HEXA_8 ) {
				if (Rank==0) cout << "[I]    ...Found Volume Section " << sectionName << endl;
				// elements array serves as a start index for connectivity list elemConnectivity
				for (int elem=0;elem<=(elemEnd-elemStart);++elem) {
					raw.cellConnIndex.push_back(raw.cellConnectivity.size());
					for (int n=0;n<elemNodeCount;++n) {
						raw.cellConnectivity.push_back(zoneCoordMap[zoneIndex-1][elemNodes[elem][n]-1]);
					}
				}
				
				globalCellCount+=(elemEnd-elemStart+1);
			} else { // If not a volume element
				// Check if a boundary condition section
				
				bool bcFlag=false;
				string bocoName(sectionName);
				for (int nbc=0;nbc<nBocos;++nbc) { // for each bc in the current zone
					// Check if bc element range matches the section element range
					if (elemStart==bc_range[nbc][0] && elemEnd==bc_range[nbc][1]) {
						bcFlag=true;
						// Store the name and indices of boundary conditions in a map
						// Loop the current content of the map
						int bcIndex=-1;
						for (mit=raw.bocoNameMap.begin();mit!=raw.bocoNameMap.end();mit++) {
							// This will find the max index number inserted so far
							bcIndex=max(bcIndex,(*mit).second);
						}
						// Increment for new bc
						bcIndex++;
						// Insert the current bc section
						// (note: if the bc name is already in the map, this will not change anything, which is good)
						raw.bocoNameMap.insert(pair<string,int>(bocoName,bcIndex));
						break;
					}
				}
				if (bcFlag) {
					string bocoName(sectionName);
					if (Rank==0) cout << "[I]    ...Found BC Section (BC_" << raw.bocoNameMap[bocoName]+1 << ") :" << sectionName << endl;
					int bcIndex=raw.bocoNameMap[bocoName];
					// If the current bc region was not found before on a different zone
					if (bcIndex>=raw.bocoConnIndex.size()) {
						vector<int> bocoConnIndex;
						vector<int> bocoConnectivity;
						for (int elem=0;elem<=(elemEnd-elemStart);++elem) {
							bocoConnIndex.push_back(bocoConnectivity.size());
							for (int n=0;n<elemNodeCount;++n) bocoConnectivity.push_back(zoneCoordMap[zoneIndex-1][elemNodes[elem][n]-1]);
						}
						raw.bocoConnIndex.push_back(bocoConnIndex);
						raw.bocoConnectivity.push_back(bocoConnectivity);
					} else { // If the current bc region was found before on a different zone
						for (int elem=0;elem<=(elemEnd-elemStart);++elem) {
							raw.bocoConnIndex[bcIndex].push_back(raw.bocoConnectivity[bcIndex].size());
							for (int n=0;n<elemNodeCount;++n) {
								raw.bocoConnectivity[bcIndex].push_back(zoneCoordMap[zoneIndex-1][elemNodes[elem][n]-1]);
							}
						}
					}
				} // end if bc
			}// if
		} // for section
	} // for zone
//} // for base

	if (Rank==0) cout << "[I] Total Node Count= " << globalNodeCount << endl;
	// Merge coordinates of the zones
	raw.x.reserve(globalNodeCount);
	raw.y.reserve(globalNodeCount);
	raw.z.reserve(globalNodeCount);
	// for zone 0
	for (int n=0;n<coordX[0].size();++n) {
		raw.x.push_back(coordX[0][n]);
		raw.y.push_back(coordY[0][n]);
		raw.z.push_back(coordZ[0][n]);
	}
	for (int zone=1;zone<nZones;++zone) {
		for (int n=0;n<coordX[zone].size();++n) {
			bool unique=true;
			for (int i=0;i<zoneCoordMap[zone-1].size();++i) {
				if (zoneCoordMap[zone][n]<=zoneCoordMap[zone-1][i]) {
					unique=false;
					break;
				}
			}
			if (unique) {
				raw.x.push_back(coordX[zone][n]);
				raw.y.push_back(coordY[zone][n]);
				raw.z.push_back(coordZ[zone][n]);
			}
		}
	}
	if (raw.x.size()!=globalNodeCount) {
		cerr << "[E] raw.x.size() (=" << raw.x.size() << ") is different from globalNodeCount (=" << globalNodeCount << ")" << endl;
		exit(1);
	}
	if (Rank==0) cout << "[I] Total Cell Count= " << globalCellCount << endl;

} // end Grid::ReadCGNS

int Grid::readTEC() {

	globalNodeCount=0;
	globalCellCount=0;
	globalFaceCount=0;

	ifstream file;
	string data,chunk;
	
	file.open(fileName.c_str(),ios::in);

	string::size_type fileLoc=0;

	vector<string> varNames;
	// Get variable names
	getline(file,data,'\n');
	fileLoc+=data.size();
	int num_vars = count(data.begin(), data.end(), '\"');
	for (int i=0;i<(num_vars/2);++i) {
		data=data.substr(data.find("\"")+1,data.size());
		chunk=data.substr(0,data.find("\""));
		if (i>2) varNames.push_back(chunk); // don't include x y z
		data=data.substr(chunk.size()+1,data.size());
	}
	
	// Get number of nodes
	file.seekg(fileLoc);
	getline(file,data,'n');
	file.seekg(fileLoc+data.size()+4);
	file >> globalNodeCount;

	// Get number of cells	
	file.seekg(fileLoc);
	getline(file,data,'e');
	file.seekg(fileLoc+data.size()+4);
	file >> globalCellCount;
	
	cout << "[I] ...Number of nodes = " << globalNodeCount << endl;
	cout << "[I] ...Number of cells = " << globalCellCount << endl;
	

	getline(file,data,'\n');

	raw.x.reserve(globalNodeCount);
	raw.y.reserve(globalNodeCount);
	raw.z.reserve(globalNodeCount);
	
	
	double Ddummy;
	// Read node coordinates
	for (unsigned int n=0;n<globalNodeCount;++n) { file >> Ddummy; raw.x.push_back(Ddummy); }
	for (unsigned int n=0;n<globalNodeCount;++n) { file >> Ddummy; raw.y.push_back(Ddummy); }
	for (unsigned int n=0;n<globalNodeCount;++n) { file >> Ddummy; raw.z.push_back(Ddummy); }

	cout << "[I] ...Read node coordinates" << endl;		
	
	// Read variables
	for (unsigned int var=0;var<varNames.size();++var) {
		for (unsigned int n=0;n<globalNodeCount;++n) {
			file >> Ddummy;
			//node[n].vars.push_back(Ddummy);
		}	
	}
	cout << "[I] ...Read variables" << endl;
	
	// Fill in cell connectivity information
	for (unsigned int c=0;c<globalCellCount;++c) {
		int conn,connPrev;
		raw.cellConnIndex.push_back(raw.cellConnectivity.size());
		for (int i=0;i<8;++i) {
			connPrev=conn;
			file >> conn; conn-=1;
			if (i==0 || connPrev!=conn) {
				raw.cellConnectivity.push_back(conn);	
			}
			
		}
	}

	// Fill in face connectivity information
	bool bcFound;
	int bcCount,elemCount;
	bcCount=0;
	while (bcFound==false) {
		bcFound=false;
		if (file.ignore(int(1e20),'B').eof()) break;
		data=file.peek();
		if (data=="C") {
			bcFound==true;
			bcCount++;
			raw.bocoNameMap.insert(pair<string,int>("BC_"+int2str(bcCount),bcCount));
			file.ignore(100,'='); file.ignore(100,'=');
			file >> elemCount;
			getline(file,data,'\n');
			cout << "[I] ......BC_" << bcCount << " has " << elemCount << " faces" << endl;
			vector<int> bocoConnIndex;
			vector<int> bocoConnectivity;
			for (int elem=0;elem<elemCount;++elem) {
				int conn,connPrev;
				bocoConnIndex.push_back(bocoConnectivity.size());
				for (int i=0;i<4;++i) {
					connPrev=conn;
					file >> conn; conn-=1;
					if (i==0 || connPrev!=conn) {
						bocoConnectivity.push_back(conn);	
					}
				}		
			}
			raw.bocoConnIndex.push_back(bocoConnIndex);
			raw.bocoConnectivity.push_back(bocoConnectivity);
		}
	}
	
	cout << raw.cellConnIndex.size() << endl;
	for (int c=0;c<raw.cellConnIndex.size();++c) {
	  for (int i=0;i<8;++i) {
	    cout << raw.cellConnectivity[raw.cellConnIndex[c]+i] << "\t" << endl;
	  }
	  cout << endl;
	}

		
} // end Grid::ReadTEC
