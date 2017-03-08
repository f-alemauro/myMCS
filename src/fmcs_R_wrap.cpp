#include "../include/config.h"
#include "../include/MCS.h"
#include <cstdlib>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>

using namespace std;
using namespace FMCS;

extern "C" {
//	std::vector<string> sdfSet;
//	int userDefinedLowerBound = 0, substructureNumLimit = 1;
//	int a = 0;
//	int b = 0;
//	int c = 0;
//	int d = 0;
//	int e = 0;
//	int f = 0;
//	
//	int *atomMismatchLowerBound = &a;
//	int *atomMismatchUpperBound = &b;
//	int *bondMismatchLowerBound = &d;
//	int *bondMismatchUpperBound = &e;
//	int *timeout = &f;
//	MCS::MatchType matchType = MCS::RING_SENSETIVE;
//	MCS::RunningMode runningMode = MCS::DETAIL;
//
//	__declspec(dllexport) void test(){
//		cout << "test done!" << endl;
//	}
//	__declspec(dllexport) void setupMCS(const char* fileName, const char** stringOne, const char** stringTwo){
//		ifstream myReadFile;
//		myReadFile.open(fileName);
//		std::stringstream buffer;
//		buffer << myReadFile.rdbuf();
//		std::string contents(buffer.str());
//		size_t last = 0;
//		size_t next = 0;
//		while ((next = contents.find("$$$$", last)) != string::npos) {
//			sdfSet.push_back(contents.substr(last, next-last+5));
//			last = next +6;
//		}
//		myReadFile.close();
//		cout<<"Read "<< sdfSet.size()<<" molecules."<<endl;
//	
//		*stringOne =sdfSet[0].c_str();
//		*stringTwo =sdfSet[1].c_str();
//	}
//	//
//	__declspec(dllexport) void computeMCS(const char* structureStringOne, const char* structureStringTwo,
//			char** sdfMCS1, char** sdfMCS2, char** sdfOne, char** sdfTwo){
//	
//		if (structureStringOne == NULL) {
//			cout << "input structure one cannot be NULL...\n";
//			return;
//		}
//		if (structureStringTwo == NULL) {
//			cout << "input structure two cannot be NULL...\n";
//			return;
//		}
//		MCSCompound compoundOne, compoundTwo;
//	
//		compoundOne.read(string(structureStringOne));
//		compoundTwo.read(string(structureStringTwo));
//		MCS mcs(compoundOne, compoundTwo, userDefinedLowerBound, substructureNumLimit, *atomMismatchLowerBound, *atomMismatchUpperBound,
//				*bondMismatchLowerBound, *bondMismatchUpperBound,
//				matchType, runningMode, *timeout);
//	
//		mcs.calculate();
//	if (runningMode == MCS::DETAIL) {
//		list<vector<size_t> > index1 = mcs.getFirstOriginalIndice();
//		list<vector<size_t> > index2 = mcs.getSecondOriginalIndice();
//
//		string s = compoundOne.createDissimilarSDFs(index1.front());
//		*sdfOne = (char*)malloc(sizeof(char)*s.length()+1);
//		strcpy(*sdfOne, s.c_str());
//
//		string s1 = compoundTwo.createDissimilarSDFs(index2.front());
//		*sdfTwo = (char*)malloc(sizeof(char) * s1.length()+1);
//		strcpy(*sdfTwo, s1.c_str());
//
//		string s2 = compoundOne.createMCSSDFs(index1.front());
//		*sdfMCS1 = (char*)malloc(sizeof(char) * s2.length()+1);
//		strcpy(*sdfMCS1, s2.c_str());
//
//		string s3 = compoundTwo.createMCSSDFs(index2.front());
//		*sdfMCS2 = (char*)malloc(sizeof(char) * s3.length()+1);
//		strcpy(*sdfMCS2, s3.c_str());
//		}
//}

//void fmcs_R_wrap(const char** structureStringOne, const char** structureStringTwo,
//		int *atomMismatchLowerBound, int *atomMismatchUpperBound,
//		int *bondMismatchLowerBound, int *bondMismatchUpperBound,
//		int *matchTypeInt, int* runningModeInt,
//		int *timeout,
//		const char** resultIdxOne, const char** resultIdxTwo,
//		const char** sdfOneSize, const char** sdfTwoSize, const char** mcsSize) {
//	if (*structureStringOne == NULL) {
//		cout << "input structure one cannot be NULL...\n";
//		return;
//	}
//	if (*structureStringTwo == NULL) {
//		cout << "input structure two cannot be NULL...\n";
//		return;
//	}
//
//	int substructureNumLimit = 1;
//	int userDefinedLowerBound = 0;
//	MCS::MatchType matchType;
//	switch(*matchTypeInt) {
//	case 0: matchType = MCS::DEFAULT; break;
//	case 1: matchType = MCS::AROMATICITY_SENSETIVE; break;
//	case 2: matchType = MCS::RING_SENSETIVE; break;
//	default:
//		;
//	}
//	MCS::RunningMode runningMode;
//	switch(*runningModeInt) {
//	case 0: runningMode = MCS::FAST; break;
//	case 1: runningMode = MCS::DETAIL; break;
//	default:
//		;
//	}
//	MCSCompound compoundOne, compoundTwo;
//
//#ifdef HAVE_LIBOPENBABEL
//	compoundOne.read(string(*structureStringOne), MCSCompound::SDF);
//	compoundTwo.read(string(*structureStringTwo), MCSCompound::SDF);
//#else
//
//	compoundOne.read(string(*structureStringOne));
//
//	compoundTwo.read(string(*structureStringTwo));
//
//#endif
//	MCS mcs(compoundOne, compoundTwo,
//			userDefinedLowerBound, substructureNumLimit,
//			*atomMismatchLowerBound, *atomMismatchUpperBound,
//			*bondMismatchLowerBound, *bondMismatchUpperBound,
//			matchType, runningMode, *timeout);
//
//	mcs.calculate();
//
//	static int cmpOneSize, cmpTwoSize, mSize;
//	cmpOneSize = mcs.getCompoundOne().size();
//	cmpTwoSize = mcs.getCompoundTwo().size();
//	mSize = mcs.size();
//
//	if (runningMode == MCS::DETAIL) {
//
//		list<vector<size_t> > index1 = mcs.getFirstOriginalIndice();
//		list<vector<size_t> > index2 = mcs.getSecondOriginalIndice();
//
//		cout << mcs.getFirstSdfResultStringList().size() << " solution(s) found..." << endl;
//		stringstream indexOneStringStream, indexTwoStringStream;
//		for (list<vector<size_t> >::const_iterator i = index1.begin(); i != index1.end(); ++i) {
//			for (vector<size_t>::const_iterator j = i->begin(); j != i->end(); ++j) {
//				indexOneStringStream << *j << " ";
//			}
//			indexOneStringStream << "\n";
//		}
//
//		for (list<vector<size_t> >::const_iterator i = index2.begin(); i != index2.end(); ++i) {
//			for (vector<size_t>::const_iterator j = i->begin(); j != i->end(); ++j) {
//				indexTwoStringStream << *j << " ";
//			}
//			indexTwoStringStream << "\n";
//		}
//
//		static string indexOneString, indexTwoString;
//		indexOneString = indexOneStringStream.str();
//		indexTwoString = indexTwoStringStream.str();
//		cout << "indexOneString: " << endl;
//		cout << indexOneString <<endl ;
//
//		cout << "indexTwoString: " << endl;
//		cout << indexTwoString <<endl ;
//		*resultIdxOne = indexOneString.c_str();
//		*resultIdxTwo = indexTwoString.c_str();
//	}
//
//	stringstream sizeStringStream;
//	sizeStringStream << cmpOneSize;
//	static string cmpOneSizeString;
//	cmpOneSizeString= sizeStringStream.str();
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
//
//}
//
void fmcs_R_wrap_mod(const char* structureStringOne, const char* structureStringTwo,
	int *atomMismatchLowerBound, int *atomMismatchUpperBound,
	int *bondMismatchLowerBound, int *bondMismatchUpperBound,
	int *matchTypeInt, int* runningModeInt,
	int *timeout) {

	if (structureStringOne == NULL) {
		cout << "input structure one cannot be NULL...\n";
		return;
	}
	if (structureStringTwo == NULL) {
		cout << "input structure two cannot be NULL...\n";
		return;
	}

	int substructureNumLimit = 1;
	int userDefinedLowerBound = 0;
	MCS::MatchType matchType;
	switch (*matchTypeInt) {
	case 0: matchType = MCS::DEFAULT; break;
	case 1: matchType = MCS::AROMATICITY_SENSETIVE; break;
	case 2: matchType = MCS::RING_SENSETIVE; break;
	default:
		;
	}


	MCS::RunningMode runningMode;
	switch (*runningModeInt) {
	case 0: runningMode = MCS::FAST; break;
	case 1: runningMode = MCS::DETAIL; break;
	default:
		;
	}

	MCSCompound compoundOne, compoundTwo;


#ifdef HAVE_LIBOPENBABEL
	compoundOne.read(string(*structureStringOne), MCSCompound::SDF);
	compoundTwo.read(string(*structureStringTwo), MCSCompound::SDF);
#else
	compoundOne.read(string(structureStringOne));

	compoundTwo.read(string(structureStringTwo));

#endif
	MCS mcs(compoundOne, compoundTwo,
		userDefinedLowerBound, substructureNumLimit,
		*atomMismatchLowerBound, *atomMismatchUpperBound,
		*bondMismatchLowerBound, *bondMismatchUpperBound,
		matchType, runningMode, *timeout);


	mcs.calculate();
	static int cmpOneSize, cmpTwoSize, mSize;
	cmpOneSize = mcs.getCompoundOne().size();
	cmpTwoSize = mcs.getCompoundTwo().size();
	mSize = mcs.size();
	if (runningMode == MCS::DETAIL) {
		list<vector<size_t> > index1 = mcs.getFirstOriginalIndice();
		list<vector<size_t> > index2 = mcs.getSecondOriginalIndice();

		string sdfOut = compoundOne.createDissimilarSDFs(index1.front());
		ofstream myfile;
		myfile.open("out1.sdf");
		myfile << sdfOut;
		myfile.close();

		sdfOut = compoundTwo.createDissimilarSDFs(index2.front());
		myfile.open("out2.sdf");
		myfile << sdfOut;
		myfile.close();

		sdfOut = compoundOne.createMCSSDFs(index1.front());
		myfile.open("outMCS1.sdf");
		myfile << sdfOut;
		myfile.close();
		sdfOut = compoundTwo.createMCSSDFs(index2.front());
		myfile.open("outMCS2.sdf");
		myfile << sdfOut;
		myfile.close();

		stringstream indexOneStringStream, indexTwoStringStream;
		for (list<vector<size_t> >::const_iterator i = index1.begin(); i != index1.end(); ++i) {
			for (vector<size_t>::const_iterator j = i->begin(); j != i->end(); ++j) {
				indexOneStringStream << *j << " ";
			}
			indexOneStringStream << "\n";
		}
		for (list<vector<size_t> >::const_iterator i = index2.begin(); i != index2.end(); ++i) {
			for (vector<size_t>::const_iterator j = i->begin(); j != i->end(); ++j) {
				indexTwoStringStream << *j << " ";
			}
			indexTwoStringStream << "\n";
		}
		static string indexOneString, indexTwoString;
		indexOneString = indexOneStringStream.str();
		indexTwoString = indexTwoStringStream.str();

		cout << "indexOneString: " << endl;
		cout << indexOneString << endl;

		cout << "indexTwoString: " << endl;
		cout << indexTwoString << endl;
	}

	stringstream sizeStringStream;

	sizeStringStream << cmpOneSize;
	static string cmpOneSizeString;
	cmpOneSizeString = sizeStringStream.str();

	sizeStringStream.str("");
	sizeStringStream << cmpTwoSize;
	static string cmpTwoSizeString;
	cmpTwoSizeString = sizeStringStream.str();

	sizeStringStream.str("");
	sizeStringStream << mSize;
	static string mSizeString;
	mSizeString = sizeStringStream.str();

}

int main(int argc, char *argv[]){
	std::vector<string> sdfSet;
	string temp, actualSDF;
	int i = 0;

	if (argc != 9){
		cout << "Missing parameters; usage is ./mcswrap sdfFileName atomMismatchLowerBound atomMismatchUpperBound "
			"bondMismatchLowerBound bondMismatchUpperBound matchTypeInt runningModeInt timeout" << endl;
		return -1;
	}

	const char* fileName = (const char*)argv[1];
	int atomMismatchLowerBound = atoi(argv[2]);
	int atomMismatchUpperBound = atoi(argv[3]);
	int bondMismatchLowerBound = atoi(argv[4]);
	int bondMismatchUpperBound = atoi(argv[5]);
	int matchTypeInt = atoi(argv[6]);
	int runningModeInt = atoi(argv[7]);
	int timeout = atoi(argv[8]);


	ifstream myReadFile;
	myReadFile.open(fileName);

	std::stringstream buffer;
	buffer << myReadFile.rdbuf();
	std::string contents(buffer.str());

	size_t last = 0;
	size_t next = 0;
	while ((next = contents.find("$$$$", last)) != string::npos) {
		sdfSet.push_back(contents.substr(last, next - last + 5));
		last = next + 5;
	}
	myReadFile.close();
	fmcs_R_wrap_mod(sdfSet[0].c_str(), sdfSet[1].c_str(), &atomMismatchLowerBound, &atomMismatchUpperBound,
		&bondMismatchLowerBound, &bondMismatchUpperBound, &matchTypeInt,
		&runningModeInt, &timeout);
	return 0;
}
}//   extern "C"
