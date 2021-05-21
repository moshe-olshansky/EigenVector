#include <cstring>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <cmath>
#include <unistd.h>
using namespace std;

int getHiCInfo(string fname, vector<std::string> &CH,vector<int> &CHL, string &genome) {
	int version;
//	string norm("NONE");
	string unit("BP");
	ifstream fin;

	fin.open(fname, fstream::in);
	if (!fin.is_open()) return -1;

	string str;
	getline(fin, str, '\0' );
	fin.read((char*)&version, sizeof(int));
	if (version < 6) {
		cerr << "Version " << version << " no longer supported" << endl;
		 exit(1);
	}
	long master;
	fin.read((char*)&master, sizeof(long));
//	string genome;
	getline(fin, genome, '\0' );
	int nattributes;
	fin.read((char*)&nattributes, sizeof(int));

// reading and ignoring attribute-value dictionary
	for (int i=0; i<nattributes; i++) {
		string key, value;
		getline(fin, key, '\0');
		getline(fin, value, '\0');
	}
	int nChrs;
	fin.read((char*)&nChrs, sizeof(int));
	vector<std::string> chroms;
	vector<int> chrLen;
// chromosome map for finding matrix
	for (int i=0; i<nChrs; i++) {
		string name;
		int length;
		getline(fin, name, '\0');
		fin.read((char*)&length, sizeof(int));
		if (name != "ALL" && name != "All") {
				chroms.insert(chroms.end(),name);
				chrLen.insert(chrLen.end(),length);
			}
    	}
	fin.close();

	CH = chroms;
	CHL = chrLen;
	return nChrs;
}
