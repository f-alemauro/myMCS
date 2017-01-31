/**
 * wrapper for the MCS class
 */

#include <cstdlib>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include "../include/MCS.h"
#include "../include/config.h"

using namespace std;
using namespace FMCS;

extern "C" {

/**
 * Write results to files
 */
void writeToFile(string data, string fileName)
{
	ofstream myFile;
	myFile.open(fileName.c_str());
	myFile << data;
	myFile.close();
	cout<< "File "<<fileName.c_str()<<"written!"<<endl;
}

/**
 * runMcs: creates the two molecules and launches the MCS computation.
 */
void runMcs(const char* structureStringOne, const char* structureStringTwo,
		int *atomMismatchLowerBound, int *atomMismatchUpperBound,
		int *bondMismatchLowerBound, int *bondMismatchUpperBound,
		int *matchTypeInt, int* runningModeInt,
		int *timeout,
		const char** resultIdxOne, const char** resultIdxTwo,
		const char** sdfOneSize, const char** sdfTwoSize, const char** mcsSize) {

	int substructureNumLimit = 1;
	int userDefinedLowerBound = 0;
	MCS::MatchType matchType;
	MCS::RunningMode runningMode;
	MCSCompound compoundOne, compoundTwo;
	int cmpOneSize, cmpTwoSize, mSize;
	string sdfOut;

	if (structureStringOne == NULL) {
		cout << "input structure one cannot be NULL...\n";
		return;
	}
	if (structureStringTwo == NULL) {
		cout << "input structure two cannot be NULL...\n";
		return;
	}

	//selecting the type of match; use RING_SENSETIVE for
	//ring closure
	switch(*matchTypeInt) {
	case 0: matchType = MCS::DEFAULT; break;
	case 1: matchType = MCS::AROMATICITY_SENSETIVE; break;
	case 2: matchType = MCS::RING_SENSETIVE; break;
	default:
		;
	}

	//selecting the type of running, fast or detail; use DETAIL for
	//ring closure
	switch(*runningModeInt) {
	case 0: runningMode = MCS::FAST; break;
	case 1: runningMode = MCS::DETAIL; break;
	default:
		;
	}

#ifdef HAVE_LIBOPENBABEL
	compoundOne.read(string(*structureStringOne), MCSCompound::SDF);
	compoundTwo.read(string(*structureStringTwo), MCSCompound::SDF);
#else
	compoundOne.read(string(structureStringOne));
	compoundTwo.read(string(structureStringTwo));
#endif

	//initialize and setup the MCS
	MCS mcs(compoundOne, compoundTwo,
			userDefinedLowerBound, substructureNumLimit,
			*atomMismatchLowerBound, *atomMismatchUpperBound,
			*bondMismatchLowerBound, *bondMismatchUpperBound,
			matchType, runningMode, *timeout);

	//retrieve the size of the two compounds
	cmpOneSize = mcs.getCompoundOne().size();
	cmpTwoSize = mcs.getCompoundTwo().size();

	//start the computation of the MCS
	mcs.calculate();

	mSize = mcs.size();


	//if (runningMode == MCS::DETAIL) {
	list<vector<size_t> > index1 = mcs.getFirstOriginalIndice();
	list<vector<size_t> > index2 = mcs.getSecondOriginalIndice();

	sdfOut = compoundOne.createDissimilarSDFs(index1.front());
	writeToFile(sdfOut, "out1.sdf");

	sdfOut = compoundTwo.createDissimilarSDFs(index2.front());
	writeToFile(sdfOut, "out2.sdf");

	sdfOut = compoundOne.createMCSSDFs(index1.front());
	writeToFile(sdfOut, "outMCS1.sdf");

	sdfOut = compoundTwo.createMCSSDFs(index2.front());
	writeToFile(sdfOut, "outMCS2.sdf");



	//	stringstream indexOneStringStream, indexTwoStringStream;
	//	for (list<vector<size_t> >::const_iterator i = index1.begin(); i != index1.end(); ++i) {
	//		for (vector<size_t>::const_iterator j = i->begin(); j != i->end(); ++j) {
	//			indexOneStringStream << *j << " ";
	//		}
	//		indexOneStringStream << "\n";
	//	}
	//	for (list<vector<size_t> >::const_iterator i = index2.begin(); i != index2.end(); ++i) {
	//		for (vector<size_t>::const_iterator j = i->begin(); j != i->end(); ++j) {
	//			indexTwoStringStream << *j << " ";
	//		}
	//		indexTwoStringStream << "\n";
	//	}
	//	static string indexOneString, indexTwoString;
	//	indexOneString = indexOneStringStream.str();
	//	indexTwoString = indexTwoStringStream.str();
	//
	//	cout << "indexOneString: " << endl;
	//	cout << indexOneString <<endl ;
	//
	//	cout << "indexTwoString: " << endl;
	//	cout << indexTwoString <<endl ;
	//
	//	*resultIdxOne = indexOneString.c_str();
	//	*resultIdxTwo = indexTwoString.c_str();
	//	}

	//	stringstream sizeStringStream;
	//
	//	sizeStringStream << cmpOneSize;
	//	static string cmpOneSizeString;
	//	cmpOneSizeString= sizeStringStream.str();
	//
	//	sizeStringStream.str("");
	//	sizeStringStream << cmpTwoSize;
	//	static string cmpTwoSizeString;
	//	cmpTwoSizeString = sizeStringStream.str();
	//
	//	sizeStringStream.str("");
	//	sizeStringStream << mSize;
	//	static string mSizeString;
	//	mSizeString = sizeStringStream.str();
	//
	//	*sdfOneSize = cmpOneSizeString.c_str();
	//	*sdfTwoSize = cmpTwoSizeString.c_str();
	//	*mcsSize = mSizeString.c_str();

}

/*
 * Entry point for the application
 */
int main(int argc, char *argv[]){

	vector<string> sdfSet;
	ifstream sdfFileStream;
	stringstream sdfBuffer;
	size_t last = 0;
	size_t next = 0;

	//check if the number of input argument is correct, else return.
	if(argc!=14){
		cout<<"Missing parameters; usage is ./mcswrap sdfFileName atomMismatchLowerBound atomMismatchUpperBound "
				"bondMismatchLowerBound bondMismatchUpperBound matchTypeInt runningModeInt timeout resultIdxOne "
				"resultIdxTwo sdfOneSize sdfTwoSize mcsSize"<<endl;
		return -1;
	}
	//std::string s = std::to_string(42);
	//cout<<s<<endl;
	//reading input parameters
	const char* fileName = (const char*)argv[1];
	int atomMismatchLowerBound = atoi(argv[2]);
	int atomMismatchUpperBound = atoi(argv[3]);
	int bondMismatchLowerBound = atoi(argv[4]);
	int bondMismatchUpperBound= atoi(argv[5]);
	int matchTypeInt= atoi(argv[6]);
	int runningModeInt= atoi(argv[7]);
	int timeout= atoi(argv[8]);
	const char** resultIdxOne = (const char**)argv[9];
	const char** resultIdxTwo = (const char**)argv[10];
	const char** sdfOneSize = (const char**)argv[11];
	const char** sdfTwoSize = (const char**)argv[12];
	const char** mcsSize = (const char**)argv[13];

	//Open the .sdf file stream
	sdfFileStream.open(fileName);
	sdfBuffer << sdfFileStream.rdbuf();
	string contents(sdfBuffer.str());

	//read each molecule in the .sdf file and fill the sdfSet variable
	//every element in sdfSet is a molecule
	while ((next = contents.find("$$$$", last)) != string::npos)
	{
		sdfSet.push_back(contents.substr(last, next-last-2));//subtract 2 for removing empty row at
		//the end of molecule
		last = next+6;//add 6 to delete "$$$$" and an empty row at the begin of the molecule

	}
	sdfFileStream.close();
	cout<<"Read "<< sdfSet.size()<<" molecules."<<endl;

	runMcs(sdfSet[0].c_str(), sdfSet[1].c_str(), &atomMismatchLowerBound, &atomMismatchUpperBound,
			&bondMismatchLowerBound, &bondMismatchUpperBound, &matchTypeInt,
			&runningModeInt, &timeout, resultIdxOne,
			resultIdxTwo, sdfOneSize, sdfTwoSize, mcsSize);
	return 0;
}
}
