#ifndef INPUTS_H
#define INPUTS_H

#include <string>
#include <map>
#include <vector>
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
	InputFile(string);
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
