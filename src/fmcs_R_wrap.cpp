#include "../config.h"

#include <cstdlib>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>

#include "MCS.h"

using namespace std;
using namespace FMCS;

extern "C" {

void fmcs_R_wrap(const char** structureStringOne, const char** structureStringTwo,
		int *atomMismatchLowerBound, int *atomMismatchUpperBound,
		int *bondMismatchLowerBound, int *bondMismatchUpperBound,
		int *matchTypeInt, int* runningModeInt,
		int *timeout,
		const char** resultIdxOne, const char** resultIdxTwo,
		const char** sdfOneSize, const char** sdfTwoSize, const char** mcsSize) {
	if (*structureStringOne == NULL) {
		cout << "input structure one cannot be NULL...\n";
		return;
	}
	if (*structureStringTwo == NULL) {
		cout << "input structure two cannot be NULL...\n";
		return;
	}

	int substructureNumLimit = 1;
	int userDefinedLowerBound = 0;
	MCS::MatchType matchType;
	switch(*matchTypeInt) {
	case 0: matchType = MCS::DEFAULT; break;
	case 1: matchType = MCS::AROMATICITY_SENSETIVE; break;
	case 2: matchType = MCS::RING_SENSETIVE; break;
	default:
		;
	}
	MCS::RunningMode runningMode;
	switch(*runningModeInt) {
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

	compoundOne.read(string(*structureStringOne));

	compoundTwo.read(string(*structureStringTwo));

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

		cout << mcs.getFirstSdfResultStringList().size() << " solution(s) found..." << endl;
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
		cout << indexOneString <<endl ;

		cout << "indexTwoString: " << endl;
		cout << indexTwoString <<endl ;
		*resultIdxOne = indexOneString.c_str();
		*resultIdxTwo = indexTwoString.c_str();
	}

	stringstream sizeStringStream;
	sizeStringStream << cmpOneSize;
	static string cmpOneSizeString;
	cmpOneSizeString= sizeStringStream.str();
	sizeStringStream.str("");
	sizeStringStream << cmpTwoSize;
	static string cmpTwoSizeString;
	cmpTwoSizeString = sizeStringStream.str();

	sizeStringStream.str("");
	sizeStringStream << mSize;
	static string mSizeString;
	mSizeString = sizeStringStream.str();

	*sdfOneSize = cmpOneSizeString.c_str();
	*sdfTwoSize = cmpTwoSizeString.c_str();
	*mcsSize = mSizeString.c_str();

}

void fmcs_R_wrap_mod(const char* structureStringOne, const char* structureStringTwo,
		int *atomMismatchLowerBound, int *atomMismatchUpperBound,
		int *bondMismatchLowerBound, int *bondMismatchUpperBound,
		int *matchTypeInt, int* runningModeInt,
		int *timeout,
		const char** resultIdxOne, const char** resultIdxTwo,
		const char** sdfOneSize, const char** sdfTwoSize, const char** mcsSize) {


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
	switch(*matchTypeInt) {
	case 0: matchType = MCS::DEFAULT; break;
	case 1: matchType = MCS::AROMATICITY_SENSETIVE; break;
	case 2: matchType = MCS::RING_SENSETIVE; break;
	default:
		;
	}


	MCS::RunningMode runningMode;
	switch(*runningModeInt) {
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
		myfile.open ("out1.sdf");
		myfile << sdfOut;
		myfile.close();

		sdfOut = compoundTwo.createDissimilarSDFs(index2.front());
		myfile.open ("out2.sdf");
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
		cout << indexOneString <<endl ;

		cout << "indexTwoString: " << endl;
		cout << indexTwoString <<endl ;

		*resultIdxOne = indexOneString.c_str();

		*resultIdxTwo = indexTwoString.c_str();
	}

	stringstream sizeStringStream;

	sizeStringStream << cmpOneSize;
	static string cmpOneSizeString;
	cmpOneSizeString= sizeStringStream.str();

	sizeStringStream.str("");
	sizeStringStream << cmpTwoSize;
	static string cmpTwoSizeString;
	cmpTwoSizeString = sizeStringStream.str();

	sizeStringStream.str("");
	sizeStringStream << mSize;
	static string mSizeString;
	mSizeString = sizeStringStream.str();

	*sdfOneSize = cmpOneSizeString.c_str();
	*sdfTwoSize = cmpTwoSizeString.c_str();
	*mcsSize = mSizeString.c_str();

}

int main(int argc, char *argv[]){
	std::vector<string> sdfSet;
	string temp,actualSDF;
	int i = 0;

	if(argc!=14){
		cout<<"Missing parameters; usage is ./mcswrap sdfFileName atomMismatchLowerBound atomMismatchUpperBound "
				"bondMismatchLowerBound bondMismatchUpperBound matchTypeInt runningModeInt timeout resultIdxOne "
				"resultIdxTwo sdfOneSize sdfTwoSize mcsSize"<<endl;
		return -1;
	}

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

	ifstream myReadFile;
	myReadFile.open(fileName);


	std::stringstream buffer;
	buffer << myReadFile.rdbuf();
	std::string contents(buffer.str());

	size_t last = 0;
	size_t next = 0;
	while ((next = contents.find("$$$$", last)) != string::npos) {
		sdfSet.push_back(contents.substr(last, next-last+5));
		last = next +6;
	}

	myReadFile.close();
	cout<<"Read "<< sdfSet.size()<<" molecules."<<endl;

	fmcs_R_wrap_mod(sdfSet[0].c_str(), sdfSet[1].c_str(), &atomMismatchLowerBound,&atomMismatchUpperBound,
			&bondMismatchLowerBound,&bondMismatchUpperBound,&matchTypeInt,
			&runningModeInt,&timeout,resultIdxOne,
			resultIdxTwo,sdfOneSize, sdfTwoSize, mcsSize);

	stringstream strValue;
	strValue << *mcsSize;

	unsigned int mcssize;
	strValue >> mcssize;


	//for(i=0; i< mcssize ;i++)
	//{
	//    cout<<"resultIdxTwo: "<< *resultIdxTwo [i] <<" "<<endl;
	//
	//}

	return 0;
}
}  // extern "C"
