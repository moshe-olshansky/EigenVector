#include <cstring>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <cmath>
#include <unistd.h>
using namespace std;

map<string, chromosome> readHeader(istream &fin, long &masterIndexPosition, string &genomeID, int &numChromosomes,
                                   int &version, long &nviPosition, long &nviLength);

unsigned int getMatrix(string fname, int binsize, string norm, string ob, bool interOnly, unsigned int **i0, unsigned int **j0, float **x0, long *m, vector<std::string> &CH,vector<int> &CHL) {
	string unit("BP");
	ifstream fin;

    long master = 0L;
    map<string, chromosome> chromosomeMap;
    string genomeID;
    int version = 0;
    long nviPosition = 0;
    long nviLength = 0;
    long totalFileSize;

	int nChrs;
	vector<std::string> chroms;
	vector<int> chrLen;
	fin.open(fname, fstream::in);
    chromosomeMap = readHeader(fin, master, genomeID, nChrs, version, nviPosition, nviLength);
    fin.close();
    for (map<string,chromosome>::iterator itr =  chromosomeMap.begin(); itr !=  chromosomeMap.end(); ++itr) {
         string chrom = itr->second.name;
         int length = itr->second.length;
	 if (chrom != "ALL" && chrom != "All") {
		chroms.insert(chroms.end(),chrom);
		chrLen.insert(chrLen.end(),length);
	}
    }

	unsigned int offset[chroms.size()];
	unsigned int current_offset = 0;
	for (int i=0; i<chroms.size(); i++) {
		offset[i] = current_offset;
		current_offset += (int) ceil(chrLen.at(i)/((double) binsize));
	}
	long nonZer = getNumRecordsForFile(fname,binsize,interOnly);
	unsigned int *ii = (unsigned int *) malloc(nonZer*sizeof(int));
	unsigned int *jj = (unsigned int *) malloc(nonZer*sizeof(int));
	float *xx = (float *) malloc(nonZer*sizeof(double));
	long pos = 0;
	int swap;
	int b1,b2;
	vector<contactRecord> records;
	int off = 0;
	if (interOnly) off = 1;
	for (int i=0; i<chroms.size(); i++) {
		for (int j=i+off; j<chroms.size();j++) {
			records = straw(ob,norm, fname, chroms.at(i), chroms.at(j), unit, binsize);
			long length=records.size();
			b1 = 0;
			b2 = 0;
			swap = 0;
			for (long k=0; k<length; k++) {
				if (records[k].binX > b1) b1 = records[k].binX;
				if (records[k].binY > b2) b2 = records[k].binY;
			}
			double junk = ((double) (chrLen.at(i)-chrLen.at(j)))*(b1-b2);
			if (junk < 0) swap = 1;
			for (long k=0; k<length; k++) {
				if (swap == 0) {
					ii[pos] = offset[i] + records[k].binX/binsize; 
					jj[pos] = offset[j] + records[k].binY/binsize; 
				}
				else {
					ii[pos] = offset[i] + records[k].binY/binsize; 
					jj[pos] = offset[j] + records[k].binX/binsize; 
				}
				xx[pos] = (float) records[k].counts;
				if (isnan(xx[pos])) xx[pos] = 0;
				pos++;
			}
			records.clear();
			records.shrink_to_fit();
		}
	}
	fin.close();

	*i0 = ii;
	*j0 = jj;
	*x0 = xx;
	*m = pos;
	CH = chroms;
	CHL = chrLen;
	return current_offset;
}
