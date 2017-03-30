/**
    @author Alessio M. Franchi, Azadi Golbamaki
    @version 0.8, 2017/03/01
*/

#include "../include/log.h"
#include "../include/MCS.h"
#include <cstdlib>
#include <fstream>

using namespace std;
using namespace FMCS;

extern "C" {
	list<string> originalMcs1, originalMcs2, closedMcs1, closedMcs2, originaDelta1, originalDelta2, closedDelta1, closedDelta2;
	/**
	 * splitSDF (function exported from the library): given an SDF file name, it extracts the first two molecule's SDF.
	 * It returns void.
	 */
#ifdef _WIN32
	__declspec(dllexport) void splitSDF(const char* fileName, char** stringOne, char** stringTwo){
#elif __linux__
	void splitSDF(const char* fileName, const char** stringOne, const char** stringTwo){
#endif
		std::vector<string> sdfSet;
		ifstream myReadFile;
		myReadFile.open(fileName); //open the SDF file
		std::stringstream buffer;
		buffer << myReadFile.rdbuf();
		std::string contents(buffer.str());
		size_t last = 0, next = 0;
		//this loop search for the "$$$$" simbols to split each molecule in the SDF file
		while ((next = contents.find("$$$$", last)) != string::npos) { //repeat until end of file
			sdfSet.push_back(contents.substr(last, next - last + 5)); //take the substring from (last) position (next-last+5) long. +5 is for the final "$$$$"
			last = next + 5; //update the (last) position to skip the "%%%%"
		}
		myReadFile.close(); //close the SDF file
		LOG(logINFO) << "Extracted " << sdfSet.size() << " molecules from the input SDF file. Considering the first two.";
		*stringOne = (char*)malloc(sizeof(char) * sdfSet[0].length() + 1);
		strcpy(*stringOne, sdfSet[0].c_str());
		*stringTwo = (char*)malloc(sizeof(char) * sdfSet[1].length() + 1);
		strcpy(*stringTwo, sdfSet[1].c_str());

	}

	/**
	* getClosedDelta (function exported from the library): get the delta SDFs of the closed MCS
	* Return void.
	*/
#ifdef _WIN32
	__declspec(dllexport) void getClosedDelta(char** sdfMCS1, int* numOfMCS1, char** sdfMCS2, int* numOfMCS2){
#elif __linux__
	void getClosedDelta(char** sdfMCS1, int* numOfMCS1, char** sdfMCS2, int* numOfMCS2){
#endif

		*numOfMCS1 = closedDelta1.size();

		int lenght = 0;
		for (list<string>::iterator it = closedDelta1.begin(); it != closedDelta1.end(); it++){
			lenght = lenght + (*it).length();
		}
		lenght = lenght + *numOfMCS1;

		*sdfMCS1 = (char*)malloc(sizeof(char) * lenght);
		memset(*sdfMCS1, 0, sizeof(char) * lenght);
		int init = 0;
		for (list<string>::iterator it = closedDelta1.begin(); it != closedDelta1.end(); it++){
			strcpy(*sdfMCS1 + init, it->c_str());
			init = it->size() + 1;
		}

		*numOfMCS2 = closedDelta2.size();
		lenght = 0;
		for (list<string>::iterator it = closedDelta2.begin(); it != closedDelta2.end(); it++){
			lenght = lenght + it->size();
		}
		lenght = lenght + *numOfMCS2;

		*sdfMCS2 = (char*)malloc(sizeof(char) * lenght);
		memset(*sdfMCS2, 0, sizeof(char) * lenght);
		init = 0;
		for (list<string>::iterator it = closedDelta2.begin(); it != closedDelta2.end(); it++){
			strcpy(*sdfMCS2 + init, it->c_str());
			init = it->size() + 1;
		}
	}


	/**
	* getOriginalDelta (function exported from the library): get the delta SDFs of the original MCS
	* Return void.
	*/
#ifdef _WIN32
	__declspec(dllexport) void getOriginalDelta(char** sdfMCS1, int* numOfMCS1, char** sdfMCS2, int* numOfMCS2){
#elif __linux__
	void getOriginalDelta(char** sdfMCS1, int* numOfMCS1, char** sdfMCS2, int* numOfMCS2){
#endif

		*numOfMCS1 = originaDelta1.size();

		int lenght = 0;
		for (list<string>::iterator it = originaDelta1.begin(); it != originaDelta1.end(); it++){
			lenght = lenght + (*it).length();
		}
		lenght = lenght + *numOfMCS1;

		*sdfMCS1 = (char*)malloc(sizeof(char) * lenght);
		memset(*sdfMCS1, 0, sizeof(char) * lenght);
		int init = 0;
		for (list<string>::iterator it = originaDelta1.begin(); it != originaDelta1.end(); it++){
			strcpy(*sdfMCS1 + init, it->c_str());
			init = it->size() + 1;
		}

		*numOfMCS2 = originalDelta2.size();
		lenght = 0;
		for (list<string>::iterator it = originalDelta2.begin(); it != originalDelta2.end(); it++){
			lenght = lenght + it->size();
		}
		lenght = lenght + *numOfMCS2;

		*sdfMCS2 = (char*)malloc(sizeof(char) * lenght);
		memset(*sdfMCS2, 0, sizeof(char) * lenght);
		init = 0;
		for (list<string>::iterator it = originalDelta2.begin(); it != originalDelta2.end(); it++){
			strcpy(*sdfMCS2 + init, it->c_str());
			init = it->size() + 1;
		}
	}

	/**
	* getClosedMCS (function exported from the library): get the two SDFs of the closed MCS
	* Return void.
	*/
#ifdef _WIN32
	__declspec(dllexport) void getClosedMCS(char** sdfMCS1, int* numOfMCS1, char** sdfMCS2, int* numOfMCS2){
#elif __linux__
	void getClosedMCS(char** sdfMCS1, int* numOfMCS1, char** sdfMCS2, int* numOfMCS2){
#endif

		*numOfMCS1 = closedMcs1.size();
			
		int lenght = 0;
		for (list<string>::iterator it = closedMcs1.begin(); it != closedMcs1.end(); it++){
			lenght = lenght + (*it).length();
		}
		lenght = lenght + *numOfMCS1;

		*sdfMCS1 = (char*)malloc(sizeof(char) * lenght);
		memset(*sdfMCS1, 0, sizeof(char) * lenght);
		int init = 0;
		for (list<string>::iterator it = closedMcs1.begin(); it != closedMcs1.end(); it++){
			strcpy(*sdfMCS1 + init, it->c_str());
			init = it->size() + 1;
		}

		*numOfMCS2 = closedMcs2.size();
		lenght = 0;
		for (list<string>::iterator it = closedMcs2.begin(); it != closedMcs2.end(); it++){
			lenght = lenght + it->size();
		}
		lenght = lenght + *numOfMCS2;

		*sdfMCS2 = (char*)malloc(sizeof(char) * lenght);
		memset(*sdfMCS2, 0, sizeof(char) * lenght);
		init = 0;
		for (list<string>::iterator it = closedMcs2.begin(); it != closedMcs2.end(); it++){
			strcpy(*sdfMCS2 + init, it->c_str());
			init = it->size() + 1;
		}
	}

	/**
	* getOriginalMCS (function exported from the library): get the two SDFs of the original MCS
	* Return void.
	*/
#ifdef _WIN32
	__declspec(dllexport) void getOriginalMCS(char** sdfMCS1, int* numOfMCS1, char** sdfMCS2, int* numOfMCS2){
#elif __linux__
	void getOriginalMCS(char** sdfMCS1, int* numOfMCS1, char** sdfMCS2, int* numOfMCS2){
#endif
		
		*numOfMCS1 = originalMcs1.size();
		int lenght = 0;
		for (list<string>::iterator it = originalMcs1.begin(); it != originalMcs1.end(); it++){
			lenght = lenght + (*it).length();
		}
		lenght = lenght + *numOfMCS1;
		
		*sdfMCS1 = (char*)malloc(sizeof(char) * lenght);
		memset(*sdfMCS1, 0, sizeof(char) * lenght);
		int init = 0;
		for (list<string>::iterator it = originalMcs1.begin(); it != originalMcs1.end(); it++){
			strcpy(*sdfMCS1 + init, it->c_str());
			init = it->size()+1;
		}

		*numOfMCS2 = originalMcs2.size();
		lenght = 0;
		for (list<string>::iterator it = originalMcs2.begin(); it != originalMcs2.end(); it++){
			lenght = lenght + it->size();
		}
		lenght = lenght + *numOfMCS2;
		*sdfMCS2 = (char*)malloc(sizeof(char) * lenght);
		memset(*sdfMCS2, 0, sizeof(char) * lenght);
		init = 0;
		for (list<string>::iterator it = originalMcs2.begin(); it != originalMcs2.end(); it++){
			strcpy(*sdfMCS2 + init, it->c_str());
			init = it->size()+1;
		}
	}



	/**
	 * computeMCS (function exported from the library): given the two SDF to be compared, it computes the MCS, and the "extra" part of the two molecules.
	 * Return void.
	 */
#ifdef _WIN32
	__declspec(dllexport) void computeMCS(const char* structureStringOne, const char* structureStringTwo){
#elif __linux__
	void computeMCS(const char* structureStringOne, const char* structureStringTwo){
#endif

		MCS::MatchType matchType = MCS::RING_SENSETIVE;
		MCS::RunningMode runningMode = MCS::DETAIL;
		

		if (structureStringOne == NULL) {
			LOG(logERROR) << "Input structure one cannot be empty!";
			return;
		}

		if (structureStringTwo == NULL) {
			LOG(logERROR) << "Input structure tow cannot be empty!";
			return;
		}

		MCSCompound compoundOne, compoundTwo;
		//parsing the two SDFs to fill the compound structure
		compoundOne.read(string(structureStringOne));
		compoundTwo.read(string(structureStringTwo));

		int a = 0;
		int *atomMismatchUpperBound = &a;
		int e = 0;
		int *atomMismatchLowerBound = &a;
		int b = 0;
		int *bondMismatchLowerBound = &b;
		int c = 0;
		int *bondMismatchUpperBound = &c;
		int d = 0;
		int *timeout = &d;

		MCS mcs(compoundOne, compoundTwo, 0, 1, *atomMismatchLowerBound, *atomMismatchUpperBound, *bondMismatchLowerBound, *bondMismatchUpperBound, matchType, runningMode, *timeout);
		mcs.calculate();

		static int cmpOneSize, cmpTwoSize, mSize;
		cmpOneSize = mcs.getCompoundOne().size();
		cmpTwoSize = mcs.getCompoundTwo().size();
		mSize = mcs.size();
		if (runningMode == MCS::DETAIL) {

			cout << "Saving output..." << endl;
			string sdfOut;
			//ofstream myfile;

			list<vector<size_t> > originalIndex1 = mcs.getFirstOriginalIndice();
			for (int i = 0; i < originalIndex1.size(); i++){
				sdfOut = compoundOne.createDissimilarSDFs(originalIndex1.front());
				originaDelta1.push_back(sdfOut);
				//nameDiss = to_string(i) + "_originalDeltaFirst.sdf";
				//myfile.open(nameDiss);
				//myfile << sdfOut;
				//myfile.close();

				sdfOut = compoundOne.createMCSSDFs(originalIndex1.front());
				originalMcs1.push_back(sdfOut);
				
				//nameDiss = to_string(i) + "_originalMCSFirst.sdf";
				//myfile.open(nameDiss);
				//myfile << sdfOut;
				//myfile.close();
				originalIndex1.pop_front();
			}

			list<vector<size_t> > originalIndex2 = mcs.getSecondOriginalIndice();
			for (int i = 0; i < originalIndex2.size(); i++){
				sdfOut = compoundTwo.createDissimilarSDFs(originalIndex2.front());
				originalDelta2.push_back(sdfOut);
				//nameDiss = to_string(i) + "_originalDeltaSecond.sdf";
				//myfile.open(nameDiss);
				//myfile << sdfOut;
				//myfile.close();

				sdfOut = compoundTwo.createMCSSDFs(originalIndex2.front());
				originalMcs2.push_back(sdfOut);
				//nameDiss = to_string(i) + "_originalMCSSecond.sdf";
				//myfile.open(nameDiss);
				//myfile << sdfOut;
				//myfile.close();

				originalIndex2.pop_front();
			}

			list<vector<size_t> > closedIndex1 = mcs.getFirstClosedIndice();
			for (int i = 0; i < closedIndex1.size(); i++){
				sdfOut = compoundOne.createDissimilarSDFs(closedIndex1.front());
				closedDelta1.push_back(sdfOut);
				//nameDiss = to_string(i) + "_closedDeltaFirst.sdf";
				//myfile.open(nameDiss);
				//myfile << sdfOut;
				//myfile.close();

				sdfOut = compoundOne.createMCSSDFs(closedIndex1.front());
				closedMcs1.push_back(sdfOut);
				//nameDiss = to_string(i) + "_closedMCSFirst.sdf";
				//myfile.open(nameDiss);
				//myfile << sdfOut;
				//myfile.close();

				closedIndex1.pop_front();
			}

			list<vector<size_t> > closedIndex2 = mcs.getSecondClosedIndice();
			for (int i = 0; i < closedIndex2.size(); i++){
				sdfOut = compoundTwo.createDissimilarSDFs(closedIndex2.front());
				closedDelta2.push_back(sdfOut);
				//nameDiss = to_string(i) + "_closedDeltaSecond.sdf";
				//myfile.open(nameDiss);
				//myfile << sdfOut;
				//myfile.close();

				sdfOut = compoundTwo.createMCSSDFs(closedIndex2.front());
				closedMcs2.push_back(sdfOut);
				//nameDiss = to_string(i) + "_closedMCSSecond.sdf";
				//myfile.open(nameDiss);
				//myfile << sdfOut;
				//myfile.close();

				closedIndex2.pop_front();
			}

			/*
			//Retrieve all the MCS (they may be one or more)
			list<vector<size_t> > index1 = mcs.getFirstOriginalIndice();
			list<vector<size_t> > index2 = mcs.getSecondOriginalIndice();


			/*string s = compoundOne.createDissimilarSDFs(index1.front());
			*sdfOne = (char*)malloc(sizeof(char)*s.length()+1);
			strcpy(*sdfOne, s.c_str());

			string s1 = compoundTwo.createDissimilarSDFs(index2.front());
			*sdfTwo = (char*)malloc(sizeof(char) * s1.length()+1);
			strcpy(*sdfTwo, s1.c_str());

			string s2 = compoundOne.createMCSSDFs(index1.front());
			*sdfMCS1 = (char*)malloc(sizeof(char) * s2.length()+1);
			strcpy(*sdfMCS1, s2.c_str());

			string s3 = compoundTwo.createMCSSDFs(index2.front());
			*sdfMCS2 = (char*)malloc(sizeof(char) * s3.length()+1);
			strcpy(*sdfMCS2, s3.c_str());*/
		}
	}


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
		compoundOne.read(string(structureStringOne));
		compoundTwo.read(string(structureStringTwo));
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
			ofstream myfile;
			string nameDiss, nameMCS, sdfOut;
			list<vector<size_t> > originalIndex1 = mcs.getFirstOriginalIndice();
			list<vector<size_t> > originalIndex2 = mcs.getSecondOriginalIndice();
			cout << "Writing to file..." << endl;
			for (int i = 0; i < originalIndex1.size(); i++){
				sdfOut = compoundOne.createDissimilarSDFs(originalIndex1.front());
				nameDiss = to_string(i) + "_originalDeltaFirst.sdf";
				myfile.open(nameDiss);
				myfile << sdfOut;
				myfile.close();
				cout << "Dissimilar SDFs for molecule 1 written" << endl;

				sdfOut = compoundOne.createMCSSDFs(originalIndex1.front());
				nameMCS = to_string(i) + "_originalMCSFirst.sdf";
				myfile.open(nameMCS);
				myfile << sdfOut;
				myfile.close();
				originalIndex1.pop_front();
				cout << "Original MCS SDFs for molecule 1 written" << endl;
			}


			for (int i = 0; i < originalIndex2.size(); i++){
				sdfOut = compoundTwo.createDissimilarSDFs(originalIndex2.front());
				nameDiss = to_string(i) + "_originalDeltaSecond.sdf";
				myfile.open(nameDiss);
				myfile << sdfOut;
				myfile.close();
				cout << "Dissimilar SDFs for molecule 2 written" << endl;
				sdfOut = compoundTwo.createMCSSDFs(originalIndex2.front());
				nameMCS = to_string(i) + "_originalMCSSecond.sdf";
				myfile.open(nameMCS);
				myfile << sdfOut;
				myfile.close();
				originalIndex2.pop_front();
				cout << "Original MCS SDFs for molecule 2 written" << endl;
			}

			list<vector<size_t> > closedIndex1 = mcs.getFirstClosedIndice();
			for (int i = 0; i < closedIndex1.size(); i++){
				sdfOut = compoundOne.createDissimilarSDFs(closedIndex1.front());
				nameDiss = to_string(i) + "_deltaFirst.sdf";
				myfile.open(nameDiss);
				myfile << sdfOut;
				myfile.close();
				cout << "Closed dissimilar SDFs for molecule 1 written" << endl;
				sdfOut = compoundOne.createMCSSDFs(closedIndex1.front());
				nameMCS = to_string(i) + "_closedMCSFirst.sdf";
				myfile.open(nameMCS);
				myfile << sdfOut;
				myfile.close();
				closedIndex1.pop_front();
				cout << "Closed MCS SDFs for molecule 1 written" << endl;
			}

			list<vector<size_t> > closedIndex2 = mcs.getSecondClosedIndice();
			for (int i = 0; i < closedIndex2.size(); i++){
				sdfOut = compoundTwo.createDissimilarSDFs(closedIndex2.front());
				nameDiss = to_string(i) + "_deltaSecond.sdf";
				myfile.open(nameDiss);
				myfile << sdfOut;
				myfile.close();
				cout << "Closed dissimilar SDFs for molecule 2 written" << endl;
				sdfOut = compoundTwo.createMCSSDFs(closedIndex2.front());
				nameMCS = to_string(i) + "_closedMCSSecond.sdf";
				myfile.open(nameMCS);
				myfile << sdfOut;
				myfile.close();
				closedIndex2.pop_front();
				cout << "Closed MCS SDFs for molecule 2 written" << endl;
			}


			/*stringstream indexOneStringStream, indexTwoStringStream;
			for (list<vector<size_t> >::const_iterator i = index1.begin(); i != index1.end(); ++i) {
			for (vector<size_t>::const_iterator j = i->begin(); j != i->end(); ++j)
			indexOneStringStream << *j << " ";
			indexOneStringStream << "\n";
			}
			for (list<vector<size_t> >::const_iterator i = index2.begin(); i != index2.end(); ++i) {
			for (vector<size_t>::const_iterator j = i->begin(); j != i->end(); ++j)
			indexTwoStringStream << *j << " ";
			indexTwoStringStream << "\n";
			}
			static string indexOneString, indexTwoString;
			indexOneString = indexOneStringStream.str();
			indexTwoString = indexTwoStringStream.str();
			cout << "indexOneString: " << endl;
			cout << indexOneString << endl;
			cout << "indexTwoString: " << endl;
			cout << indexTwoString << endl;*/
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
		Log::ReportingLevel() = logDEBUG; //set the log level to DEBUG
		if (argc != 9){
			LOG(logERROR) << "Missing parameters; usage is ./runMCS sdfFileName atomMismatchLowerBound atomMismatchUpperBound "
				"bondMismatchLowerBound bondMismatchUpperBound matchTypeInt runningModeInt timeout";
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

		std::vector<string> sdfSet;
		ifstream myReadFile;
		myReadFile.open(fileName); //open the SDF file
		std::stringstream buffer;
		buffer << myReadFile.rdbuf();
		std::string contents(buffer.str());
		size_t last = 0, next = 0;
		//this loop search for the "$$$$" simbols to split each molecule in the SDF file
		while ((next = contents.find("$$$$", last)) != string::npos) { //repeat until end of file
			sdfSet.push_back(contents.substr(last, next - last + 5)); //take the substring from (last) position (next-last+5) long. +5 is for the final "$$$$"
			last = next + 5; //update the (last) position to skip the "%%%%"
		}
		myReadFile.close(); //close the SDF file
		LOG(logINFO) << "Extracted " << sdfSet.size() << " molecules from the input SDF file. Considering the first two.";
		fmcs_R_wrap_mod(sdfSet[0].c_str(), sdfSet[1].c_str(), &atomMismatchLowerBound, &atomMismatchUpperBound,
			&bondMismatchLowerBound, &bondMismatchUpperBound, &matchTypeInt,
			&runningModeInt, &timeout);
		return 0;
	}
}