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

GridRawData raw;
double block_stitch_tolerance=1.e-8;

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

	vector<double> min_x,max_x,min_y,max_y,min_z,max_z;
	min_x.resize(nZones);
	min_y.resize(nZones);
	min_z.resize(nZones);
	max_x.resize(nZones);
	max_y.resize(nZones);
	max_z.resize(nZones);
	
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
		
		min_x[zoneIndex-1]=coordX[zoneIndex-1][0];
		min_y[zoneIndex-1]=coordY[zoneIndex-1][0];
		min_z[zoneIndex-1]=coordZ[zoneIndex-1][0];
		max_x[zoneIndex-1]=coordX[zoneIndex-1][0];
		max_y[zoneIndex-1]=coordY[zoneIndex-1][0];
		max_z[zoneIndex-1]=coordZ[zoneIndex-1][0];

		for (int c=1;c<coordX[zoneIndex-1].size();++c) {
			min_x[zoneIndex-1]=min(coordX[zoneIndex-1][c],min_x[zoneIndex-1]);
			min_y[zoneIndex-1]=min(coordY[zoneIndex-1][c],min_y[zoneIndex-1]);
			min_z[zoneIndex-1]=min(coordZ[zoneIndex-1][c],min_z[zoneIndex-1]);
			max_x[zoneIndex-1]=max(coordX[zoneIndex-1][c],max_x[zoneIndex-1]);
			max_y[zoneIndex-1]=max(coordY[zoneIndex-1][c],max_y[zoneIndex-1]);
			max_z[zoneIndex-1]=max(coordZ[zoneIndex-1][c],max_z[zoneIndex-1]);
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

				for (int z=0;z<zoneIndex-1;++z) { // Loop other zones
					if (coordX[zoneIndex-1][c]<(min_x[z]-block_stitch_tolerance)) break;
					if (coordY[zoneIndex-1][c]<(min_y[z]-block_stitch_tolerance)) break;
					if (coordZ[zoneIndex-1][c]<(min_z[z]-block_stitch_tolerance)) break;
					
					if (coordX[zoneIndex-1][c]>(max_x[z]+block_stitch_tolerance)) break;
					if (coordY[zoneIndex-1][c]>(max_y[z]+block_stitch_tolerance)) break;
					if (coordZ[zoneIndex-1][c]>(max_z[z]+block_stitch_tolerance)) break;

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
				int connDataSize;
				cg_ElementDataSize(fileIndex,baseIndex,zoneIndex,sectionIndex,&connDataSize);
				vector<int> elemNodes;
				elemNodes.resize(connDataSize);
				cg_elements_read(fileIndex,baseIndex,zoneIndex,sectionIndex,&elemNodes[0],0);
				int connIndex=0;
				for (int elem=elemStart;elem<=elemEnd;++elem) {
					cg_npe(ElementType_t (elemNodes[connIndex]),&elemNodeCount);
					raw.cellConnIndex.push_back(raw.cellConnectivity.size());
					for (int n=1;n<=elemNodeCount;++n) { // First entry is the cell type
 						raw.cellConnectivity.push_back(zoneCoordMap[zoneIndex-1][elemNodes[connIndex+n]-1]);
 					}
					connIndex+=(elemNodeCount+1);
				}
				elemNodes.clear();
 				globalCellCount+=(elemEnd-elemStart+1);
				if (Rank==0) cout << "[I]    ...Found Mixed Section " << sectionName << endl;
			} else {
				int connDataSize;
				cg_ElementDataSize(fileIndex,baseIndex,zoneIndex,sectionIndex,&connDataSize);
				vector<int> elemNodes;
				elemNodes.resize(connDataSize);
				cg_elements_read(fileIndex,baseIndex,zoneIndex,sectionIndex,&elemNodes[0],0);
				int connIndex=0;
				if (elemType==TETRA_4 || elemType==PYRA_5 || elemType==PENTA_6 || elemType==HEXA_8 ) {
					if (Rank==0) cout << "[I]    ...Found Volume Section " << sectionName << endl;
					// elements array serves as a start index for connectivity list elemConnectivity
					for (int elem=0;elem<=(elemEnd-elemStart);++elem) {
						raw.cellConnIndex.push_back(raw.cellConnectivity.size());
						for (int n=0;n<elemNodeCount;++n) {
							raw.cellConnectivity.push_back(zoneCoordMap[zoneIndex-1][elemNodes[connIndex+n]-1]);
						}
						connIndex+=elemNodeCount;
					}
					globalCellCount+=(elemEnd-elemStart+1);
				} else { // If not a volume element
					// Scan all the boundary condition regions
					for (int nbc=0;nbc<nBocos;++nbc) {
						if (bc_method[nbc]==ELEMENT_LIST) {
							for (int elem=0;elem<=(elemEnd-elemStart);++elem) {
								if (bc_element_list[nbc].find(elemStart+elem)!=bc_element_list[nbc].end()) {
									for (int n=0;n<elemNodeCount;++n) {
										raw.bocoNodes[nbc].insert(zoneCoordMap[zoneIndex-1][elemNodes[connIndex+n]-1]);
									}
									connIndex+=elemNodeCount;
								}
							}
						}
					}
				}// if
				elemNodes.clear();
			} // if
		} // for section
	} // for zone
	nBocos=raw.bocoNameMap.size();
//} // for base

	if (Rank==0) cout << "[I] Total Node Count= " << globalNodeCount << endl;
	// Merge coordinates of the zones

	// for zone 0
	raw.node.resize(globalNodeCount);
	Vec3D temp;
	for (int z=0;z<nZones;++z) {
		for (int n=0;n<coordX[z].size();++n) {
			if (zoneCoordMap[z][n]>=0) {
				temp[0]=coordX[z][n];
				temp[1]=coordY[z][n];
				temp[2]=coordZ[z][n];
				raw.node[zoneCoordMap[z][n]]=temp;
			}
		}
		
	}
	
	/*
	raw.node.reserve(globalNodeCount);
	for (int n=0;n<coordX[0].size();++n) {
		temp[0]=coordX[0][n];
		temp[1]=coordY[0][n];
		temp[2]=coordZ[0][n];
		raw.node.push_back(temp);
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
				temp[0]=coordX[zone][n];
				temp[1]=coordY[zone][n];
				temp[2]=coordZ[zone][n];
				raw.node.push_back(temp);
			}
		}
	}
	
	if (raw.node.size()!=globalNodeCount) {
		cerr << "[E] raw.node.size() (=" << raw.node.size() << ") is different from globalNodeCount (=" << globalNodeCount << ")" << endl;
		exit(1);
	}
	 
	 */
	 
	if (Rank==0) cout << "[I] Total Cell Count= " << globalCellCount << endl;

	return 0;
	
} // end Grid::ReadCGNS
