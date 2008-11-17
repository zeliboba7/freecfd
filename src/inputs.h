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

#include <string>
#include <map>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
using namespace std;
#include "vec3d.h"

class Subsection {
public:
	map<string,double> doubles;
	map<string,int> ints;
	map<string,string> strings;
	map<string,Vec3D> Vec3Ds;
	void register_double(string varName);
	void register_int(string varName);
	void register_string(string varName);
	void register_Vec3D(string varName);
};

class numberedSubsection {
public:
	int count;
	vector<map<string,double> >doubles;
	vector<map<string,int> >ints;
	vector<map<string,string> > strings;
	vector<map<string,Vec3D> >Vec3Ds;
	void register_double(string varName);
	void register_int(string varName);
	void register_string(string varName);
	void register_Vec3D(string varName);
};

class Section {
public:
	map<string,double> doubles;
	map<string,int> ints;
	map<string,string> strings;
	map<string,Vec3D> Vec3Ds;
	map<string,Subsection> subsections;
	map<string,numberedSubsection> numberedSubsections;
	void register_double(string varName);
	void register_int(string varName);
	void register_string(string varName);
	void register_Vec3D(string varName);
	void register_subsection(string subsectionName);
	void register_numberedSubsection(string subsectionName);
};

class InputFile {
public:
	string fileName;
	map<string,Section> section;
	InputFile(void);
	void setFile(string);
	void read(void);
	int search(fstream&,string);
	void register_section(string sectionName);
	void read_section(string sectionName);
};

void read_doubles(map<string,double> &doubles, string sectionData);
void read_ints(map<string,int> &ints, string sectionData) ;
void read_strings(map<string,string> &strings, string sectionData) ;
void read_Vec3Ds(map<string,Vec3D> &Vec3Ds, string sectionData) ;


#endif
