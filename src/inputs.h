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
#ifndef INPUTS_H
#define INPUTS_H

#include "commons.h"
#include <map>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cctype>
using namespace std;
#include "vec3d.h"
extern int Rank;
extern ofstream infoFile;

template <typename T> class entry {
public:
	T value;
	T default_value;
	bool is_found;
	bool is_required;
	operator T() { return value; }
	bool operator== (const T &right) { return (value==right); }
	bool operator!= (const T &right) { return (value!=right); }
};


class InputBaseContainer {
public:	
	string rawData;
	string name,parentName;
	bool is_required, is_found;
	map<string,entry<int> > ints;
	map<string,entry<double> > doubles;
	map<string,entry<string> > strings;
	map<string,entry<Vec3D> > Vec3Ds;
	
	void register_int(string varName,bool is_required=true,int default_value=0) { 
		entry<int>  temp; temp.is_required=is_required; temp.default_value=default_value; ints.insert(pair<string,entry<int> > (varName,temp));
	}
	void register_double(string varName,bool is_required=true,double default_value=0) { 
		entry<double> temp; temp.is_required=is_required; temp.default_value=default_value; doubles.insert(pair<string,entry<double> > (varName,temp));
	}
	void register_string(string varName,bool is_required=true,string default_value="") { 
		entry<string> temp; temp.is_required=is_required; temp.default_value=default_value; strings.insert(pair<string,entry<string> > (varName,temp));
	}
	void register_Vec3D(string varName,bool is_required=true,Vec3D default_value=0) { 
		entry<Vec3D> temp;  temp.is_required=is_required; 
		default_value=default_value[0];
		temp.default_value=default_value; 
		Vec3Ds.insert(pair<string,entry<Vec3D> > (varName,temp));
	}
	
	void readEntries(void);
	void entry_not_found(string varName) {
		string path="";
		if (name!="") path=name+" -> ";
		if (parentName!="") path=parentName+" -> ";
		if (Rank==0) cerr << "[E] Required input entry " << path  <<  varName << " could not be found!!" << endl;
		exit(1);
	}
	
	entry<int> get_int(string varName) {return ints[varName];}
	entry<double> get_double(string varName) {return doubles[varName];}
	entry<string> get_string(string varName) {return strings[varName];}
	entry<Vec3D> get_Vec3D(string varName) {return Vec3Ds[varName];}
};


class Subsection: public InputBaseContainer {
	public:
		bool numbered;
		int count;
};

class Section: public InputBaseContainer {
public:
	map<string,Subsection> subsections;
	map<string,vector<Subsection> > numberedSubsections;
	void registerSubsection(string subsectionName, bool is_numbered=false, bool is_required=true) {
		if (is_numbered) {
			Subsection temp; temp.numbered=true; temp.count=0; temp.is_required=is_required;
			temp.name=subsectionName; temp.parentName=name;
			vector<Subsection> temp2;
			temp2.push_back(temp);
			numberedSubsections.insert(pair<string,vector<Subsection> > (subsectionName,temp2));
		} else {
			Subsection temp; temp.numbered=false; temp.is_required=is_required; 
			temp.name=subsectionName; temp.parentName=name;
			subsections.insert(pair<string,Subsection> (subsectionName,temp));
		}
	}

	Subsection &subsection(string subsectionName, int number=-1) {
		if (number>=0) { return numberedSubsections[subsectionName][number]; } 
		else { return subsections[subsectionName]; }
	}
};


class InputFile: public InputBaseContainer {
public:
	string fileName;
	map<string,Section> sections;
	
	InputFile(void);
	void setFile(string);
	void read(void);
	void readSection(string sectionName);
	void readSubsection(Subsection &sub);
	
	void registerSection(string sectionName, bool is_required=true) {
		Section temp; temp.is_required=is_required;  temp.name=sectionName;
		sections.insert(pair<string,Section> (sectionName,temp));
	}
	Section &section(string sectionName) {return sections[sectionName];}
	
	void stripComments(string &data);
	void strip_white_spaces(string &data);
	
};


/*
	
	
void read_doubles(map<string,double> &doubles, string sectionData);
void read_ints(map<string,int> &ints, string sectionData) ;
void read_strings(map<string,string> &strings, string sectionData) ;
void read_Vec3Ds(map<string,Vec3D> &Vec3Ds, string sectionData) ; */
bool extract_in_between(string &data, string begin, string end, string &result,bool check_char_before=false, string acceptList="");
int number_of_occurances(string haystack, string needle);
		
#endif
