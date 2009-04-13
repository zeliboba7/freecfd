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
#include "commons.h"
#include <iostream>
#include <fstream>

using namespace std;
#include "vec3d.h"
#include "inputs.h"

	 
extern string int2str(int number) ;

InputFile::InputFile(void) {
	;
}

void InputFile::setFile(string fName) {

	// Check if the input file exists
	fileName=fName;
	fstream file;
	file.open(fileName.c_str());
	//if (file.is_open() && Rank==0) { 
	if (file.is_open()) { 
		if (Rank==0) cout << "[I] Found input file " << fileName << endl; 
	} else {
		if (Rank==0) cerr << "[E] Input file " << fileName << " could not be found!!" << endl;
		exit(1);
	}
	// Read the whole input file to rawData
	string line;
	rawData="";
	while(getline(file, line)) rawData += line + "\n";
	file.close();
	
	// Strip all the inline or block comments (C++ style) from the rawData
	stripComments(rawData);
	// Strip all the white spaces from the rawData
	strip_white_spaces(rawData);

}

void InputFile::read(void) {
	// Read all the sections
	// Define iterator for section names
	map<string,Section>::iterator sectionIter;
	// For each section
	for (sectionIter=sections.begin();sectionIter!=sections.end();sectionIter++) {
		string sectionName=sectionIter->first;
		readSection(sectionName);
	}
	
	// Read variables
	readEntries();
	cout << rawData << endl;
	
}

void InputFile::readSection(string sectionName) {
	// Get section data
	section(sectionName).is_found=extract_in_between(rawData,sectionName+"{","}",section(sectionName).rawData,true,"};");
	// Exit if the section is required but couldn't be found
	if (section(sectionName).is_required && !section(sectionName).is_found) {
		if (Rank==0) cerr << "[E] Required input section " << sectionName << " could not be found!!" << endl;
		exit(1);
	} 
	
	// Read subsections
	map<string,Subsection>::iterator subsectionIter;
	for (subsectionIter=section(sectionName).subsections.begin();subsectionIter!=section(sectionName).subsections.end();subsectionIter++) {
		readSubsection(subsectionIter->second);
	}
	
	// Read numbered subsections
	map<string,vector<Subsection> >::iterator nsubsectionIter;
	for (nsubsectionIter=section(sectionName).numberedSubsections.begin();nsubsectionIter!=section(sectionName).numberedSubsections.end();nsubsectionIter++) {
		// Find how many subsections there are
		int count=number_of_occurances(section(sectionName).rawData,nsubsectionIter->first+"_");
		// First occurance was already initialized when numbered subsection was registered 
		// Copy duplicate the first entry count-1 times 
		nsubsectionIter->second[0].count=count;
		for (int i=1;i<count;++i) nsubsectionIter->second.push_back(nsubsectionIter->second[0]);
		// Fill in each one of them
		for (int i=0;i<count;++i) {
			// Append number to name
			nsubsectionIter->second[i].name+="_"+int2str(i+1);
			//string subsectionName=nsubsectionName+"_"+int2str(i+1);
			readSubsection(nsubsectionIter->second[i]);
		}
		
	}
	
	// Read variables
	section(sectionName).readEntries();
	
	if (Rank==0) cout << "[I] Read input section: " << sectionName << endl;
	
}

void InputFile::readSubsection(Subsection &sub) {
	// Get subsection data
	sub.is_found=extract_in_between(section(sub.parentName).rawData,sub.name+"(",");",sub.rawData,true,"};");
	// Exit if the subsection is required but couldn't be found
	if (sub.is_required && !sub.is_found) {
		if (Rank==0) cerr << "[E] Required input subsection " << sub.parentName << " -> " << sub.name << " could not be found!!" << endl;
		exit(1);
	} 
	sub.readEntries();
	return;
}

void InputBaseContainer::readEntries(void) {
	// Loop ints
	map<string,entry<int> >::iterator intIter;
	for (intIter=ints.begin();intIter!=ints.end();intIter++) {
		string varName=intIter->first;
		string varValue;
 		// Find the variable
		intIter->second.is_found=extract_in_between(rawData,varName+"=",";",varValue,true,"({;");
		if (intIter->second.is_found) {
			intIter->second.value=atoi(varValue.c_str());
		} else { // if not found
			if (intIter->second.is_required) {
				entry_not_found(varName);
			} else {
				intIter->second.value=intIter->second.default_value;	
			}
		}	
	} // end loop ints
	// Loop doubles
 	map<string,entry<double> >::iterator doubleIter;
 	for (doubleIter=doubles.begin();doubleIter!=doubles.end();doubleIter++) {
 		string varName=doubleIter->first;
		string varValue;
 		// Find the variable
		doubleIter->second.is_found=extract_in_between(rawData,varName+"=",";",varValue,true,"({;");
		if (doubleIter->second.is_found) {
			doubleIter->second.value=strtod(varValue.c_str(),NULL);
		} else { // if not found
			if (doubleIter->second.is_required) {
				entry_not_found(varName);
			} else {
				doubleIter->second.value=doubleIter->second.default_value;	
			}
		}
		//cout << varName << "\t" << doubleIter->second.value << endl; // DEBUG
		
	} // end loop doubles
	// Loop strings
	map<string,entry<string> >::iterator stringIter;
	for (stringIter=strings.begin();stringIter!=strings.end();stringIter++) {
		string varName=stringIter->first;
		string varValue;
 		// Find the variable
		stringIter->second.is_found=extract_in_between(rawData,varName+"=",";",varValue,true,"({;");
		if (stringIter->second.is_found) {
			stringIter->second.value=varValue;
		} else { // if not found
			if (stringIter->second.is_required) {
				entry_not_found(varName);
			} else {
				stringIter->second.value=stringIter->second.default_value;	
			}
		}
	} // end loop strings
	// Loop Vec3Ds
	map<string,entry<Vec3D> >::iterator Vec3DIter;
	for (Vec3DIter=Vec3Ds.begin();Vec3DIter!=Vec3Ds.end();Vec3DIter++) {
		string varName=Vec3DIter->first;
		string varValue;
 		// Find the variable
		Vec3DIter->second.is_found=extract_in_between(rawData,varName+"=",";",varValue,true,"({;");
		if (Vec3DIter->second.is_found) {
			// The data is read as string. Now convert it to double components of a vector
			char *pEnd;
			Vec3DIter->second.value.comp[0]=strtod(varValue.c_str()+1,&pEnd);
			Vec3DIter->second.value.comp[1]=strtod(pEnd+1,&pEnd);
			Vec3DIter->second.value.comp[2]=strtod(pEnd+1,NULL);
		} else { // if not found
			if (Vec3DIter->second.is_required) {
				entry_not_found(varName);
			} else {
				Vec3DIter->second.value=Vec3DIter->second.default_value;	
			}
		}
	} // end loop Vec3Ds
	
} // end InputBaseContainer::readEntries

void InputFile::stripComments(string &data) {
	string dummy;
	while (extract_in_between(data,"/*","*/",dummy));
	while (extract_in_between(data,"//","\n",dummy));
}

void InputFile::strip_white_spaces(string &data) {
	string whitespaces(" \t\f\v\n\r");
	size_t found;
	while (found!=string::npos) {
		found=data.find_first_of(whitespaces);
		if (found!=string::npos) data.erase(found,1);
	}
}
	
bool extract_in_between(string &data, string begin, string end, string &result,bool check_char_before, string acceptList) {
	size_t begin_pos, end_pos;
	begin_pos=0; end_pos=0;
	string pre;
	bool found=false;
	while (!found) {
		// Find the first occurance position of the beginning sequence
		begin_pos=data.find(begin,end_pos);
		if (begin_pos==string::npos) return false;
		// From where the first search left off, find the first occurance position of the ending sequence
		end_pos=data.find(end,begin_pos+begin.length());
		if (end_pos==string::npos) return false;
		// Check the character just before the beginning delimiter (if asked for)
		if (check_char_before) {
			if (begin_pos==0) {
				found=true;
			} else {
				pre=data.substr(begin_pos-1,1);
				if (pre.find_first_of(acceptList)!=string::npos) found=true;
				
			}
		} else {
			found=true;
		}
		if (found) {
			// Extract the what's in between
			result=data.substr(begin_pos+begin.length(),end_pos-begin_pos-begin.length());
			// Remove that chunk from the originial data 
			data.replace(begin_pos,end_pos+end.length()-begin_pos,"");
		}
	}
	return true;
}

int number_of_occurances(string haystack, string needle) {
	int count=-1;
	size_t found,start;
	while (found!=string::npos) {
		count++;
		start=(count==0)? 0 : found+needle.length();
		found=haystack.find(needle,start);
	}
	return count;
}