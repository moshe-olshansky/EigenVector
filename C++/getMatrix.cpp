#include <cstring>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <cmath>
#include <unistd.h>
using namespace std;

int getMatrix(string fname, int binsize, string norm, int **i, int **j, float **x, long *m, vector<std::string> &CH,vector<int> &CHL) {
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
	string genome;
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
	int offset[chroms.size()];
	int current_offset = 0;
	for (int i=0; i<chroms.size(); i++) {
		offset[i] = current_offset;
		current_offset += (int) ceil(chrLen.at(i)/((double) binsize));
	}
	fin.close();
	long nonZer = 0;
	for (int i=0; i<chroms.size(); i++) {
		for (int j=i; j<chroms.size();j++) {
			int num = getSize(norm, fname, chroms.at(i), chroms.at(j), unit, binsize);
			nonZer += num;
		}
	}
	int *ii = (int *) malloc(nonZer*sizeof(int));
	int *jj = (int *) malloc(nonZer*sizeof(int));
	float *xx = (float *) malloc(nonZer*sizeof(float));
	long pos = 0;
	vector<contactRecord> records;
	for (int i=0; i<chroms.size(); i++) {
		for (int j=i; j<chroms.size();j++) {
			records = straw(norm, fname, chroms.at(i), chroms.at(j), unit, binsize);
			long length=records.size();
			for (long k=0; k<length; k++) {
				ii[pos] = offset[i] + records[k].binX/binsize; 
				jj[pos] = offset[j] + records[k].binY/binsize; 
				xx[pos] = (float) records[k].counts;
				if (!isnan(xx[pos])) pos++;
			}
		}
	}
	fin.close();
	records.clear();
	records.shrink_to_fit();

	*i = ii;
	*j = jj;
	*x = xx;
	*m = pos;
	CH = chroms;
	CHL = chrLen;
	return current_offset;
}

int getMatrix(string fname, int binsize, string norm, int **i, int **j, double **x, long *m, vector<std::string> &CH,vector<int> &CHL) {
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
	string genome;
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
	int offset[chroms.size()];
	int current_offset = 0;
	for (int i=0; i<chroms.size(); i++) {
		offset[i] = current_offset;
		current_offset += (int) ceil(chrLen.at(i)/((double) binsize));
	}
	fin.close();
	long nonZer = 0;
	for (int i=0; i<chroms.size(); i++) {
		for (int j=i; j<chroms.size();j++) {
			int num = getSize(norm, fname, chroms.at(i), chroms.at(j), unit, binsize);
			nonZer += num;
		}
	}
	int *ii = (int *) malloc(nonZer*sizeof(int));
	int *jj = (int *) malloc(nonZer*sizeof(int));
	double *xx = (double *) malloc(nonZer*sizeof(double));
	long pos = 0;
	vector<contactRecord> records;
	for (int i=0; i<chroms.size(); i++) {
		for (int j=i; j<chroms.size();j++) {
			records = straw(norm, fname, chroms.at(i), chroms.at(j), unit, binsize);
			long length=records.size();
			for (long k=0; k<length; k++) {
				ii[pos] = offset[i] + records[k].binX/binsize; 
				jj[pos] = offset[j] + records[k].binY/binsize; 
				xx[pos] = (double) records[k].counts;
				if (!isnan(xx[pos])) pos++;
			}
		}
	}
	fin.close();
	records.clear();
	records.shrink_to_fit();

	*i = ii;
	*j = jj;
	*x = xx;
	*m = pos;
	CH = chroms;
	CHL = chrLen;
	return current_offset;
}

int getMatrix(string fname, int binsize, string norm, int **i, int **j, double **x, long *m) {
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
	string genome;
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
	int offset[chroms.size()];
	int current_offset = 0;
	for (int i=0; i<chroms.size(); i++) {
		offset[i] = current_offset;
		current_offset += (int) ceil(chrLen.at(i)/((double) binsize));
	}
	fin.close();
	long nonZer = 0;
	for (int i=0; i<chroms.size(); i++) {
		for (int j=i; j<chroms.size();j++) {
			int num = getSize(norm, fname, chroms.at(i), chroms.at(j), unit, binsize);
			nonZer += num;
		}
	}
	int *ii = (int *) malloc(nonZer*sizeof(int));
	int *jj = (int *) malloc(nonZer*sizeof(int));
	double *xx = (double *) malloc(nonZer*sizeof(double));
	long pos = 0;
	vector<contactRecord> records;
	for (int i=0; i<chroms.size(); i++) {
		for (int j=i; j<chroms.size();j++) {
			records = straw(norm, fname, chroms.at(i), chroms.at(j), unit, binsize);
			long length=records.size();
			for (long k=0; k<length; k++) {
				ii[pos] = offset[i] + records[k].binX/binsize; 
				jj[pos] = offset[j] + records[k].binY/binsize; 
				xx[pos] = (double) records[k].counts;
				if (!isnan(xx[pos])) pos++;
			}
		}
	}
	fin.close();
	records.clear();
	records.shrink_to_fit();

	*i = ii;
	*j = jj;
	*x = xx;
	*m = pos;
	return current_offset;
}

int getMatrix(string fname, int binsize, string norm, int **i, int **j, float **x, long *m, string chr) {
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
	string genome;
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
	int len;
	bool found = false;
	for (int i=0; i<nChrs; i++) {
		string name;
		getline(fin, name, '\0');
		fin.read((char*)&len, sizeof(int));
		if (name == chr) {
			found = true;
			break;
		}
    	}
	fin.close();
	if (!found) return(0);
	int nonZer = getSize(norm, fname, chr, chr, unit, binsize);
	int *ii = (int *) malloc(nonZer*sizeof(int));
	int *jj = (int *) malloc(nonZer*sizeof(int));
	float *xx = (float *) malloc(nonZer*sizeof(float));
	vector<contactRecord> records;
	long pos=0;
	records = straw(norm, fname, chr, chr, unit, binsize);
	long length=records.size();
	for (long k=0; k<length; k++) {
		ii[pos] = records[k].binX/binsize; 
		jj[pos] = records[k].binY/binsize; 
		xx[pos] = (float) records[k].counts;
		if (!isnan(xx[pos])) pos++;
	}
	fin.close();
	records.clear();
	records.shrink_to_fit();

	*i = ii;
	*j = jj;
	*x = xx;
	*m = pos;
	return (int) ceil(len/((double) binsize));
}
