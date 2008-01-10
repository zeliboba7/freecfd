#include <vector>
#include <cmath>
#include <iostream>
using namespace std;

#include "sparse.h"

VecSparse::VecSparse(void) {
	;
}

int VecSparse::position(int index) {
	for (int position=0;position<indices.size();++position) {
		if (index==indices[position]) return position;
	}
	return -1;
}

void VecSparse::put(int index, double a) {
	int pos=position(index);
	if (pos<0) {
		data.push_back(a);
		indices.push_back(index);
	} else {
		data[pos]+=a;
	}
}

double VecSparse::get(int index) {
	int pos=position(index);
	if (pos<0) {
		return 0.;
	} else {
		return data[pos];
	}
}

void VecSparse::flush(void) {
	data.clear();
	indices.clear();
}

VecSparse &VecSparse::operator*= (const double &right) {
	for (unsigned int i=0;i<data.size();++i) {
		data[i]*=right;
	}
	return *this;
}

VecSparse VecSparse::operator*(const double &right) {
	VecSparse temp;
	for (unsigned int i=0;i<data.size();++i) {
		temp.put(indices[i],data[i]*right);
	}
	return temp;
}

VecSparse &VecSparse::operator/= (const double &right) {
	for (unsigned int i=0;i<data.size();++i) {
		data[i]/=right;
	}
	return *this;
}

VecSparse VecSparse::operator/ (const double &right) {
	VecSparse temp;
	for (unsigned int i=0;i<data.size();++i) {
		temp.put(indices[i],data[i]/right);
	}
	return temp;
}

VecSparse &VecSparse::operator+= (const double &right) {
	for (unsigned int i=0;i<data.size();++i) {
		data[i]+=right;
	}
	return *this;
}

VecSparse VecSparse::operator+ (const double &right) {
	VecSparse temp;
	for (unsigned int i=0;i<data.size();++i) {
		temp.put(indices[i],data[i]+right);
	}
	return temp;
}

VecSparse &VecSparse::operator-= (const double &right) {
	for (unsigned int i=0;i<data.size();++i) {
		data[i]*=right;
	}
	return *this;
}

VecSparse VecSparse::operator- (const double &right) {
	VecSparse temp;
	for (unsigned int i=0;i<data.size();++i) {
		temp.put(indices[i],data[i]-right);
	}
	return temp;
}

VecSparse operator*(const double &left, const VecSparse &right) {
	VecSparse temp;
	for (unsigned int i=0;i<right.data.size();++i) {
		temp.put(right.indices[i],left*right.data[i]);
	}
	return temp;
}

VecSparse operator/ (const double &left, const VecSparse &right) {
	VecSparse temp;
	for (unsigned int i=0;i<right.data.size();++i) {
		temp.put(right.indices[i],left/right.data[i]);
	}
	return temp;
}

VecSparse operator+ (const double &left, const VecSparse &right) {
	VecSparse temp;
	for (unsigned int i=0;i<right.data.size();++i) {
		temp.put(right.indices[i],left+right.data[i]);
	}
	return temp;
}

VecSparse operator- (const double &left, const VecSparse &right) {
	VecSparse temp;
	for (unsigned int i=0;i<right.data.size();++i) {
		temp.put(right.indices[i],left-right.data[i]);
	}
	return temp;
}

VecSparse &VecSparse::operator+= (const VecSparse &right) {
	for (unsigned int i=0;i<right.data.size();++i) {
		put(right.indices[i],right.data[i]);
	}
	return *this;
}

VecSparse &VecSparse::operator-= (const VecSparse &right) {
	for (unsigned int i=0;i<right.data.size();++i) {
		put(right.indices[i],-1.*right.data[i]);
	}
	return *this;
}




