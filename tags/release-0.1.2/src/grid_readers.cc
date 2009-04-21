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

	int fileIndex,baseIndex,nBases,nZones,nSections,nBocos;
	char zoneName[20],sectionName[20]; //baseName[20]
	//  int nBases,cellDim,physDim;
	int size[3];
	
	int POINT_LIST=1;
	int ELEMENT_LIST=2;
	
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
	// zoneCoordMap returns the global id of a given node in a given zone
	std::vector<int> zoneCoordMap[nZones];

	map<string,int>::iterator mit;
	
	for (int zoneIndex=1;zoneIndex<=nZones;++zoneIndex) { // For each zone
		vector<int> bc_method;
		vector<set<int> > bc_element_list;
		
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

		// Boundary condition regions may be given as element or point ranges, or element or point lists
		// For each zone and each boundary condition region, store a point list, convert if another method is used
		for (int bocoIndex=1;bocoIndex<=nBocos;++bocoIndex) {
			int dummy;
			// Find out the bc name and specification method 
			char bocoName[20];
			BCType_t bocotype;
			PointSetType_t ptset_type;
			int npnts;
			DataType_t NormalDataType;
 			cg_boco_info(fileIndex,baseIndex,zoneIndex,bocoIndex,bocoName,
 				      &bocotype,&ptset_type,&npnts,&dummy,&dummy,&NormalDataType,&dummy);

			// Check if this bc name matches to those found before
			bool new_bc=true;
			int bcIndex=-1;
			string bcName(bocoName);
			for (mit=raw.bocoNameMap.begin();mit!=raw.bocoNameMap.end();mit++) {
				if( (*mit).first == bcName) {
					new_bc=false;
					bcIndex=(*mit).second;
					break;
				}
			}
			
			// If no match, create new bc
    		        if (new_bc) {
				bcIndex=raw.bocoNameMap.size();
				raw.bocoNameMap.insert(pair<string,int>(bcName,bcIndex));
				raw.bocoNodes.resize(bcIndex+1);
			}
			
			if (Rank==0) cout << "[I] ...Found boundary condition BC_" << bcIndex+1 << " : " << bcName << endl;
			
			if (bc_method.size()<bcIndex+1) {
				bc_method.resize(bcIndex+1);
				bc_element_list.resize(bcIndex+1);
			}
			
			vector<int> list; list.resize(npnts);
			cg_boco_read(fileIndex,baseIndex,zoneIndex,bocoIndex,&list[0],&dummy);
			
			// Check the bc specification method
			if (ptset_type==PointList) {
				bc_method[bcIndex]=POINT_LIST;
				for (int i=0;i<list.size();++i) raw.bocoNodes[bcIndex].insert(zoneCoordMap[zoneIndex-1][list[i]-1]);
			} else if (ptset_type==ElementList) {
				bc_method[bcIndex]=ELEMENT_LIST;
				for (int i=0;i<list.size();++i) bc_element_list[bcIndex].insert(list[i]);
			} else if (ptset_type==PointRange) {
				bc_method[bcIndex]=POINT_LIST;
				for (int i=list[0];i<=list[1];++i) raw.bocoNodes[bcIndex].insert(zoneCoordMap[zoneIndex-1][i-1]);
			} else if (ptset_type==ElementRange) {
				// Convert element range to element list
				bc_method[bcIndex]=ELEMENT_LIST;
				for (int i=list[0];i<=list[1];++i) bc_element_list[bcIndex].insert(i);
			} else {
				if (Rank==0) cerr << "[E] Boundary condition specification is not recognized" << endl;
				exit(1);
			}
			
		} // for boco
		nBocos=raw.bocoNameMap.size();
		
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

			if (elemType==MIXED) {
				for (int elem=elemStart;elem<=elemEnd;++elem) {
					cg_ElementPartialSize(fileIndex,baseIndex,zoneIndex,sectionIndex,elem,elem,&elemNodeCount);
					int elemNodes[1][elemNodeCount];
					cg_elements_partial_read(fileIndex,baseIndex,zoneIndex,sectionIndex,elem,elem,*elemNodes,0);
					raw.cellConnIndex.push_back(raw.cellConnectivity.size());
					for (int n=1;n<elemNodeCount;++n) { // First entry is the cell type
						raw.cellConnectivity.push_back(zoneCoordMap[zoneIndex-1][elemNodes[0][n]-1]);
					}
				}
				globalCellCount+=(elemEnd-elemStart+1);
				if (Rank==0) cout << "[E]    ...Found Mixed Section " << sectionName << endl;
			} else {
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
					// Scan all the boundary condition regions
					for (int nbc=0;nbc<nBocos;++nbc) {
						if (bc_method[nbc]==ELEMENT_LIST) {
							for (int elem=0;elem<=(elemEnd-elemStart);++elem) {
								if (bc_element_list[nbc].find(elemStart+elem)!=bc_element_list[nbc].end()) {
									for (int n=0;n<elemNodeCount;++n) {
										raw.bocoNodes[nbc].insert(zoneCoordMap[zoneIndex-1][elemNodes[elem][n]-1]);
									}
								}
							}
						}
					}
				}// if
			} // if
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

	return 0;
	
} // end Grid::ReadCGNS

int Grid::readTEC() {

// 	globalNodeCount=0;
// 	globalCellCount=0;
// 	globalFaceCount=0;
// 
// 	ifstream file;
// 	string data,chunk;
// 	
// 	file.open(fileName.c_str(),ios::in);
// 
// 	string::size_type fileLoc=0;
// 
// 	vector<string> varNames;
// 	// Get variable names
// 	getline(file,data,'\n');
// 	fileLoc+=data.size();
// 	int num_vars = count(data.begin(), data.end(), '\"');
// 	for (int i=0;i<(num_vars/2);++i) {
// 		data=data.substr(data.find("\"")+1,data.size());
// 		chunk=data.substr(0,data.find("\""));
// 		if (i>2) varNames.push_back(chunk); // don't include x y z
// 		data=data.substr(chunk.size()+1,data.size());
// 	}
// 	
// 	// Get number of nodes
// 	file.seekg(fileLoc);
// 	getline(file,data,'n');
// 	file.seekg(fileLoc+data.size()+4);
// 	file >> globalNodeCount;
// 
// 	// Get number of cells	
// 	file.seekg(fileLoc);
// 	getline(file,data,'e');
// 	file.seekg(fileLoc+data.size()+4);
// 	file >> globalCellCount;
// 	
// 	cout << "[I] ...Number of nodes = " << globalNodeCount << endl;
// 	cout << "[I] ...Number of cells = " << globalCellCount << endl;
// 	
// 
// 	getline(file,data,'\n');
// 
// 	raw.x.reserve(globalNodeCount);
// 	raw.y.reserve(globalNodeCount);
// 	raw.z.reserve(globalNodeCount);
// 	
// 	
// 	double Ddummy;
// 	// Read node coordinates
// 	for (int n=0;n<globalNodeCount;++n) { file >> Ddummy; raw.x.push_back(Ddummy); }
// 	for (int n=0;n<globalNodeCount;++n) { file >> Ddummy; raw.y.push_back(Ddummy); }
// 	for (int n=0;n<globalNodeCount;++n) { file >> Ddummy; raw.z.push_back(Ddummy); }
// 
// 	cout << "[I] ...Read node coordinates" << endl;		
// 	
// 	// Read variables
// 	for (int var=0;var<varNames.size();++var) {
// 		for (int n=0;n<globalNodeCount;++n) {
// 			file >> Ddummy;
// 			//node[n].vars.push_back(Ddummy);
// 		}	
// 	}
// 	cout << "[I] ...Read variables" << endl;
// 	
// 	// Fill in cell connectivity information
// 	for (int c=0;c<globalCellCount;++c) {
// 		int conn,connPrev;
// 		raw.cellConnIndex.push_back(raw.cellConnectivity.size());
// 		for (int i=0;i<8;++i) {
// 			connPrev=conn;
// 			file >> conn; conn-=1;
// 			if (i==0 || connPrev!=conn) {
// 				raw.cellConnectivity.push_back(conn);	
// 			}
// 			
// 		}
// 	}
// 
// 	// Fill in face connectivity information
// 	bool bcFound;
// 	int bcCount,elemCount;
// 	bcCount=0;
// 	while (bcFound==false) {
// 		bcFound=false;
// 		if (file.ignore(int(1e20),'B').eof()) break;
// 		data=file.peek();
// 		if (data=="C") {
// 			bcFound==true;
// 			bcCount++;
// 			raw.bocoNameMap.insert(pair<string,int>("BC_"+int2str(bcCount),bcCount));
// 			file.ignore(100,'='); file.ignore(100,'=');
// 			file >> elemCount;
// 			getline(file,data,'\n');
// 			cout << "[I] ......BC_" << bcCount << " has " << elemCount << " faces" << endl;
// 			vector<int> bocoConnIndex;
// 			vector<int> bocoConnectivity;
// 			for (int elem=0;elem<elemCount;++elem) {
// 				int conn,connPrev;
// 				bocoConnIndex.push_back(bocoConnectivity.size());
// 				for (int i=0;i<4;++i) {
// 					connPrev=conn;
// 					file >> conn; conn-=1;
// 					if (i==0 || connPrev!=conn) {
// 						bocoConnectivity.push_back(conn);	
// 					}
// 				}		
// 			}
// 			raw.bocoConnIndex.push_back(bocoConnIndex);
// 			raw.bocoConnectivity.push_back(bocoConnectivity);
// 		}
// 	}
// 	
// 	cout << raw.cellConnIndex.size() << endl;
// 	for (int c=0;c<raw.cellConnIndex.size();++c) {
// 	  for (int i=0;i<8;++i) {
// 	    cout << raw.cellConnectivity[raw.cellConnIndex[c]+i] << "\t" << endl;
// 	  }
// 	  cout << endl;
// 	}

	return 0;
	
} // end Grid::ReadTEC
