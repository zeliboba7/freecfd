#include <iostream>
#include <fstream>
#include <string>

using namespace std;
#include "vec3d.h"
#include "inputs.h"

InputFile::InputFile (string fName) {
	// Check if the input file exists
	fileName=fName;
	fstream file;
	file.open(fileName.c_str());
	if(file.is_open()) cout << "* Found input file " << fileName << endl;
	file.close();
}

void InputFile::read_section (string sectionName) {
	
	// Open the input file
	fstream file;
	file.open(fileName.c_str());
	
	string::size_type loc;
	string sectionData="";
	string subsectionData;
	string subsectionName;
	string::size_type endLoc;
	string nameWithIndex;
	string whitespaces (" \t\f\v\n\r");
	
	// Find the section
	if(search(file,sectionName)) {
		cout << "* \t Found section " << sectionName << endl;
	}
	// Find what is in between { and } after the section name appears in file
	getline(file,sectionData,'{'); getline(file,sectionData,'}');
	
	// Define iterator for double variables inside the section
	map<string,double>::iterator doubleIter;
	
	read_doubles(section[sectionName].doubles,sectionData);
	read_ints(section[sectionName].ints,sectionData);
	read_strings(section[sectionName].strings,sectionData);
	read_Vec3Ds(section[sectionName].Vec3Ds,sectionData);	
	
	// For each subsection
	
	map<string,Subsection>::iterator subsectionIter;
        for (subsectionIter=section[sectionName].subsections.begin();subsectionIter!=section[sectionName].subsections.end();subsectionIter++) {
            string subsectionName=subsectionIter->first;
	    // Find the section
	   if(sectionData.find("(",0)) {
		cout << "* \t\t Found subsection " << subsectionName << endl;
	    }
	    // Find what is in between ( and ) after the subsection name appears
            loc=sectionData.find(subsectionName,0);
            loc=sectionData.find("(",loc);
            endLoc=sectionData.find(")",loc);
            subsectionData=sectionData.substr(loc+1,endLoc-loc-1);
	    read_doubles(section[sectionName].subsections[subsectionName].doubles,subsectionData);
	    read_ints(section[sectionName].subsections[subsectionName].ints,subsectionData);
	    read_strings(section[sectionName].subsections[subsectionName].strings,subsectionData);
	    read_Vec3Ds(section[sectionName].subsections[subsectionName].Vec3Ds,subsectionData);
        }
	
	// For each numbered subsection
	
	map<string,numberedSubsection>::iterator numberedSubsectionIter;
        for (numberedSubsectionIter=section[sectionName].numberedSubsections.begin();numberedSubsectionIter!=section[sectionName].numberedSubsections.end();numberedSubsectionIter++) {
		
		subsectionName=numberedSubsectionIter->first;
		
		// Find how many of them are there
		for (int i=1;i<=20;++i) {
			char dummy[4];
			sprintf(dummy, "%4d", i);
			nameWithIndex=subsectionName+"_"+dummy;
			size_t found;
			int flag=0;
			while (flag==0) {
				found=nameWithIndex.find_first_of(whitespaces);
				if (found!=string::npos) {				
					nameWithIndex.erase(found,1);
				}else {
					flag=1;
				}
				
			}
			
			// Find the subsection
			if(sectionData.find(nameWithIndex,0)!=string::npos) {
				cout << "* \t\t Found subsection " << nameWithIndex << endl;
				if (i>1) {
					section[sectionName].numberedSubsections[subsectionName].doubles.push_back(section[sectionName].numberedSubsections[subsectionName].doubles[i-2]);
					section[sectionName].numberedSubsections[subsectionName].ints.push_back(section[sectionName].numberedSubsections[subsectionName].ints[i-2]);
					section[sectionName].numberedSubsections[subsectionName].strings.push_back(section[sectionName].numberedSubsections[subsectionName].strings[i-2]);
					section[sectionName].numberedSubsections[subsectionName].Vec3Ds.push_back(section[sectionName].numberedSubsections[subsectionName].Vec3Ds[i-2]);	
				}
				// Find what is in between ( and ) after the subsection name appears
				loc=sectionData.find(nameWithIndex,0);
				loc=sectionData.find("(",loc);
				endLoc=sectionData.find(")",loc);				
				subsectionData=sectionData.substr(loc+1,endLoc-loc-1);	
				read_doubles(section[sectionName].numberedSubsections[subsectionName].doubles[i-1],subsectionData);	
				read_ints(section[sectionName].numberedSubsections[subsectionName].ints[i-1],subsectionData);
				read_strings(section[sectionName].numberedSubsections[subsectionName].strings[i-1],subsectionData);
				read_Vec3Ds(section[sectionName].numberedSubsections[subsectionName].Vec3Ds[i-1],subsectionData);					
				section[sectionName].numberedSubsections[subsectionName].count=i;

			} else {
				break;
			}
			
		}		
		
        }

	file.close();
	cout << "* \t Finished reading section "<< sectionName << "\n" << endl;	

}


void InputFile::read (void) {
	// Open the input file
	fstream file;
	file.open(fileName.c_str());
    
	// Define iterator for section names
	map<string,Section>::iterator sectionIter;
	// For each section
	for (sectionIter=section.begin();sectionIter!=section.end();sectionIter++) {
		string sectionName=sectionIter->first;
		read_section(sectionName);
    }
    file.close();
    cout << "* Finished reading input file\n" << endl;
}

void InputFile::register_section(string sectionName) {
    Section temp; 
    section.insert(pair<string,Section>(sectionName,temp));
}

int InputFile::search (fstream& file, string searchString) {
    file.seekg(0);
    string line,dummy;
    string::size_type fileLoc=0;
    string::size_type lineLoc=0;
    while (getline(file,line,'{')) {        
        lineLoc=line.find(searchString);
        if (lineLoc!=line.npos) {
            file.seekg(fileLoc+lineLoc);
        return 1;
        }
        fileLoc+=line.length();
    }
    return 0;
}

void Section::register_double(string varName) {
    double temp; 
    doubles.insert(pair<string,double>(varName,temp));
}

void Section::register_int(string varName) {
    int temp; 
    ints.insert(pair<string,int>(varName,temp));
}

void Section::register_string(string varName) {
    string temp; 
    strings.insert(pair<string,string>(varName,temp));
}

void Section::register_Vec3D(string varName) {
    Vec3D temp; 
    Vec3Ds.insert(pair<string,Vec3D>(varName,temp));
}

void Section::register_numberedSubsection(string subsectionName) {
    numberedSubsection temp; 
    numberedSubsections.insert(pair<string,numberedSubsection>(subsectionName,temp));
    map<string,double> temp2;
    numberedSubsections[subsectionName].doubles.push_back(temp2);
    map<string,int> temp3;
    numberedSubsections[subsectionName].ints.push_back(temp3);
    map<string,string> temp4;
    numberedSubsections[subsectionName].strings.push_back(temp4);
    map<string,Vec3D> temp5;
    numberedSubsections[subsectionName].Vec3Ds.push_back(temp5);	
}

void Section::register_subsection(string subsectionName) {
    Subsection temp; 
    subsections.insert(pair<string,Subsection>(subsectionName,temp));
}

void Subsection::register_double(string varName) {
    double temp; 
    doubles.insert(pair<string,double>(varName,temp));
}

void Subsection::register_int(string varName) {
    int temp; 
    ints.insert(pair<string,int>(varName,temp));
}

void Subsection::register_string(string varName) {
    string temp; 
    strings.insert(pair<string,string>(varName,temp));
}

void Subsection::register_Vec3D(string varName) {
    Vec3D temp; 
    Vec3Ds.insert(pair<string,Vec3D>(varName,temp));
}

void numberedSubsection::register_double(string varName) {
    double temp;
    doubles[0].insert(pair<string,double>(varName,temp));
}

void numberedSubsection::register_int(string varName) {
    int temp;
    ints[0].insert(pair<string,int>(varName,temp));
}

void numberedSubsection::register_string(string varName) {
    string temp;
    strings[0].insert(pair<string,string>(varName,temp));
}

void numberedSubsection::register_Vec3D(string varName) {
    Vec3D temp;
    Vec3Ds[0].insert(pair<string,Vec3D>(varName,temp));
}

void read_doubles (map<string,double> &doubles, string sectionData) {
	// For each double variable
	map<string,double>::iterator doubleIter;
	for (doubleIter=doubles.begin();doubleIter!=doubles.end();doubleIter++) {
		string varName=doubleIter->first;
		// Find the variable
		string::size_type loc=0;
		int flag=0;
		while (flag==0) {
			loc=sectionData.find(varName,loc);
			if (loc==string::npos) break;
			string after=sectionData.substr(loc+varName.length(),1);
			string before;
			if (loc!=0) {
				before=sectionData.substr(loc-1,1);
			}else {
				before=" ";
			}
			if (after == " " | after == "=") {
				if (before==" " | before==";") {
					flag=1;
				}
			}
			loc=loc+varName.length();
		}
		// Find what is in between = and ; after the variable name appears
		loc=sectionData.find("=",loc);
		string::size_type endLoc=sectionData.find(";",loc);
		string temp=sectionData.substr(loc+1,endLoc-loc-1);
		char *pEnd;
		// The data is read as string. Now convert it to double
		doubleIter->second=strtod(temp.c_str(),&pEnd);
		cout << "* \t\t " << doubleIter->first << "=" << doubleIter->second << endl;
	}
}

void read_ints (map<string,int> &ints, string sectionData) {
	// For each integer variable
	map<string,int>::iterator intIter;
        for (intIter=ints.begin();intIter!=ints.end();intIter++) {
		string varName=intIter->first;
		// Find the variable
		string::size_type loc=0;
		int flag=0;
		while (flag==0) {
			loc=sectionData.find(varName,loc);
			if (loc==string::npos) break;
			string after=sectionData.substr(loc+varName.length(),1);
			string before;
			if (loc!=0) {
				before=sectionData.substr(loc-1,1);
			}else {
				before=" ";
			}
			if (after == " " | after == "=") {
				if (before==" " | before==";") {
					flag=1;
				}
			}
			loc=loc+varName.length();
		}
	    // Find what is in between = and ; after the variable name appears		
            loc=sectionData.find("=",loc);
            string::size_type endLoc=sectionData.find(";",loc);
            string temp=sectionData.substr(loc+1,endLoc-loc-1);
            char *pEnd;
	    // The data is read as string. Now convert it to integer
	    intIter->second=strtol(temp.c_str(),&pEnd,10);
            cout << "* \t\t " << intIter->first << "=" << intIter->second << endl;
        }
}

void read_strings (map<string,string> &strings, string sectionData) {
	// For each string variable
        map<string,string>::iterator stringIter;
        for (stringIter=strings.begin();stringIter!=strings.end();stringIter++) {
            string varName=stringIter->first;
		// Find the variable
		string::size_type loc=0;
		int flag=0;
		while (flag==0) {
			loc=sectionData.find(varName,loc);
			if (loc==string::npos) break;
			string after=sectionData.substr(loc+varName.length(),1);
			string before;
			if (loc!=0) {
				before=sectionData.substr(loc-1,1);
			}else {
				before=" ";
			}
			if (after == " " | after == "=") {
				if (before==" " | before==";") {
					flag=1;
				}
			}
			loc=loc+varName.length();
		}
	    // Find what is in between = and ; after the variable name appears		
            loc=sectionData.find("=",loc);
            string::size_type endLoc=sectionData.find(";",loc);
            string temp=sectionData.substr(loc+1,endLoc-loc-1);
            stringIter->second=temp;
            cout << "* \t\t " << stringIter->first << "=" << stringIter->second << endl;
        }
}

void read_Vec3Ds (map<string,Vec3D> &Vec3Ds, string sectionData) {
	// For each Vec3D variable
        map<string,Vec3D>::iterator vecIter;
        for (vecIter=Vec3Ds.begin();vecIter!=Vec3Ds.end();vecIter++) {
            string varName=vecIter->first;
		// Find the variable
		string::size_type loc=0;
		int flag=0;
		while (flag==0) {
			loc=sectionData.find(varName,loc);
			if (loc==string::npos) break;
			string after=sectionData.substr(loc+varName.length(),1);
			string before;
			if (loc!=0) {
				before=sectionData.substr(loc-1,1);
			}else {
				before=" ";
			}
			if (after == " " | after == "=") {
				if (before==" " | before==";") {
					flag=1;
				}
			}
			loc=loc+varName.length();
		}		
	    // Find what is in between [ and ] after the variable name appears
            loc=sectionData.find("[",loc);
            string::size_type endLoc=sectionData.find("]",loc);
            string temp=sectionData.substr(loc+1,endLoc-loc-1);
            char *pEnd;
	    // The data is read as string. Now convert it to double components of a vector
            vecIter->second.comp[0]=strtod(temp.c_str(),&pEnd);
            vecIter->second.comp[1]=strtod(pEnd+1,&pEnd);
            vecIter->second.comp[2]=strtod(pEnd+1,NULL);
            cout << "* \t\t " << vecIter->first << "="  << vecIter->second << endl;            
        }	
}
