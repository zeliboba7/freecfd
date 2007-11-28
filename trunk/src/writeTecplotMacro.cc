#include <iostream>
#include <fstream>
#include <string>
using namespace std;

extern double dt;
extern int timeStepMax, outFreq;

void writeTecplotMacro(int restart, int timeStepMax, int outFreq) {
ofstream file;
int timeStep;
char dummy[12];
string fileName;

file.open("tec.mcr", ios::out);

file << "#!MC 1100" << endl;
file << "$!VarSet |MFBD| = './'" << endl;
file << "$!READDATASET  '" ;

for (timeStep=0;timeStep<timeStepMax+restart;timeStep++) {

	if ((timeStep+1)%outFreq==0) {
		// Print timeStep integer to character
		sprintf(dummy,"%12d",timeStep+1);
		// Convert character to string and erase leading whitespaces
		string fileName=dummy;
		fileName.erase(0,fileName.rfind(" ",fileName.length())+1);
	    fileName="out"+fileName+".dat";
		file << "\"|MFBD|/" << fileName << "\" "  ;
	}
}

file << "'" << endl;
file << "READDATAOPTION = NEW" << endl;
file << "RESETSTYLE = YES" << endl;
file << "INCLUDETEXT = NO" << endl;
file << "INCLUDEGEOM = NO" << endl;
file << "INCLUDECUSTOMLABELS = NO" << endl;
file << "VARLOADMODE = BYNAME" << endl;
file << "ASSIGNSTRANDIDS = YES" << endl;
file << "INITIALPLOTTYPE = CARTESIAN3D" << endl;
file << "VARNAMELIST = '\"x\" \"y\" \"z\" \"rho\" \"u\" \"v\" \"w\" \"p\" '" << endl;
file << " $!FIELDLAYERS SHOWMESH = NO" << endl;
file << " $!GLOBALCONTOUR 1  VAR = 4" << endl;
file << " $!CONTOURLEVELS RESETTONICE" << endl;
file << "CONTOURGROUP = 1" << endl;
file << "APPROXNUMVALUES = 15" << endl;
file << "$!FIELDLAYERS SHOWCONTOUR = YES" << endl;

/*
int count=1;
for (timeStep=0;timeStep<timeStepMax;timeStep++) {
	if ((timeStep+1)%outFreq==0) {
file << "$!ACTIVEFIELDMAPS = [" << count << "]" << endl;
file << "$!EXTRACTFROMPOLYLINE" << endl;
file << "EXTRACTTHROUGHVOLUME = NO" << endl;
file << "EXTRACTLINEPOINTSONLY = NO" << endl;
file << "INCLUDEDISTANCEVAR = NO" << endl;
file << "NUMPTS = 100" << endl;
file << "EXTRACTTOFILE = NO" << endl;
file << "RAWDATA" << endl;
file << "2" << endl;
file << "0.306697096959 0.637188995808 0" << endl;
file << "-0.382811704963 0.27820927986 0" << endl;
		++count;
	}
}

int total=count-1;
file << "$!PLOTTYPE = XYLINE" << endl;
file << "$!DELETELINEMAPS  [1-6]" << endl;
count=1;
for (timeStep=0;timeStep<timeStepMax;timeStep++) {
	if ((timeStep+1)%outFreq==0) {
file << "$!CREATELINEMAP" << endl;
file << "$!LINEMAP [" << count << "]  NAME = 'Map " << count << "'" << endl;
file << "$!LINEMAP [" << count << "]  ASSIGN{YAXISVAR = 4}" << endl;
file << "$!LINEMAP [" << count << "]  ASSIGN{ZONE = "<< count+total << "}" << endl;
file << "$!ACTIVELINEMAPS += [" << count << "]" << endl;
		++count;
  }
}
*/
file.close();
return;
}
