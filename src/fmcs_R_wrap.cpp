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
	list<string> originalMcs1, originalMcs2, closedMcs1, closedMcs2, originaDelta1, originalDelta2, closedDelta1, closedDelta2, prunedMcs1, prunedMcs2, prunedDelta1,prunedDelta2;
	
	//Parameter variables
	int atomMismatchUpperBound = 0;
	int atomMismatchLowerBound = 0;
	int bondMismatchLowerBound = 0;
	int bondMismatchUpperBound = 0;
	MCS::MatchType matchType = MCS::RING_SENSETIVE;
	int timeout = 0;
	int substructureNumLimit = 1;
	int userDefinedLowerBound = 0;
	MCS::RunningMode runningMode = MCS::DETAIL;
	int loggingLevel = 2; //set the log level to Warning
	
	/**
	* setLogginLevel (function exported from the library): set the level for logging
	* Return void
	*/
#ifdef _WIN32
	__declspec(dllexport) void setLogginLevel(const int* level){
#elif __linux__
	void setLogginLevel(const int* level){
#endif
		loggingLevel = *level;
	}

	/**
	* setAtomMismatchBounds (function exported from the library): set the lower and upper bounds for the atom mismatch
	* Return void
	*/
#ifdef _WIN32
	__declspec(dllexport) void setAtomMismatchBounds(const int* lower, const int* upper){
#elif __linux__
	void setAtomMismatchBounds(const int* lower, const int* upper){
#endif
		atomMismatchUpperBound = *upper;
		atomMismatchLowerBound = *lower;
	}

	/**
	* setBondMismatchBounds (function exported from the library): set the lower and upper bounds for the bond mismatch
	* Return void
	*/
#ifdef _WIN32
	__declspec(dllexport) void setBondMismatchBounds(const int* lower, const int* upper){
#elif __linux__
	void setBondMismatchBounds(const int* lower, const int* upper){
#endif
		bondMismatchUpperBound = *upper;
		bondMismatchLowerBound = *lower;
	}

	/**
	* setBondMismatchBounds (function exported from the library): set the lower and upper bounds for the bond mismatch
	* Return void
	*/
#ifdef _WIN32
	__declspec(dllexport) void setMatchAndRunningMode(const int* mType, const int* rType){
#elif __linux__
	void setMatchType(const int* mType, const int* rType){
#endif

		
		switch (*mType) {
		case 0: matchType = MCS::DEFAULT; break;
		case 1: matchType = MCS::AROMATICITY_SENSETIVE; break;
		case 2: matchType = MCS::RING_SENSETIVE; break;
		default:
			;
		}

		switch (*rType) {
		case 0: runningMode = MCS::FAST; break;
		case 1: runningMode = MCS::DETAIL; break;
		default:
			;
		}
	}

	/**
	 * splitSDF (function exported from the library): given an SDF file name, it extracts the first two molecule's SDF.
	 * It returns void.
	 */
#ifdef _WIN32
	__declspec(dllexport) void splitSDF(const char* fileName, char** stringOne, char** stringTwo){
#elif __linux__
	void splitSDF(const char* fileName, char** stringOne, char** stringTwo){
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
	* getPrunedDelta (function exported from the library): get the delta SDFs of the pruned MCS
	* Return void.
	*/
#ifdef _WIN32
	__declspec(dllexport) void getPrunedDelta(char** sdfMCS1, int* numOfMCS1, char** sdfMCS2, int* numOfMCS2){
#elif __linux__
	void getPrunedDelta(char** sdfMCS1, int* numOfMCS1, char** sdfMCS2, int* numOfMCS2){
#endif

		*numOfMCS1 = prunedDelta1.size();

		int lenght = 0;
		for (list<string>::iterator it = prunedDelta1.begin(); it != prunedDelta1.end(); it++){
			lenght = lenght + (*it).length();
		}
		lenght = lenght + *numOfMCS1;

		*sdfMCS1 = (char*)malloc(sizeof(char) * lenght);
		memset(*sdfMCS1, 0, sizeof(char) * lenght);
		int init = 0;
		for (list<string>::iterator it = prunedDelta1.begin(); it != prunedDelta1.end(); it++){
			strcpy(*sdfMCS1 + init, it->c_str());
			init = it->size() + 1;
		}

		*numOfMCS2 = prunedDelta2.size();
		lenght = 0;
		for (list<string>::iterator it = prunedDelta2.begin(); it != prunedDelta2.end(); it++){
			lenght = lenght + it->size();
		}
		lenght = lenght + *numOfMCS2;

		*sdfMCS2 = (char*)malloc(sizeof(char) * lenght);
		memset(*sdfMCS2, 0, sizeof(char) * lenght);
		init = 0;
		for (list<string>::iterator it = prunedDelta2.begin(); it != prunedDelta2.end(); it++){
			strcpy(*sdfMCS2 + init, it->c_str());
			init = it->size() + 1;
		}
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
	__declspec(dllexport) void getPrunedMCS(char** sdfMCS1, int* numOfMCS1, char** sdfMCS2, int* numOfMCS2){
#elif __linux__
	void getPrunedMCS(char** sdfMCS1, int* numOfMCS1, char** sdfMCS2, int* numOfMCS2){
#endif

		*numOfMCS1 = prunedMcs1.size();
			
		int lenght = 0;
		for (list<string>::iterator it = prunedMcs1.begin(); it != prunedMcs1.end(); it++){
			lenght = lenght + (*it).length();
		}
		lenght = lenght + *numOfMCS1;

		*sdfMCS1 = (char*)malloc(sizeof(char) * lenght);
		memset(*sdfMCS1, 0, sizeof(char) * lenght);
		int init = 0;
		for (list<string>::iterator it = prunedMcs1.begin(); it != prunedMcs1.end(); it++){
			strcpy(*sdfMCS1 + init, it->c_str());
			init = it->size() + 1;
		}

		*numOfMCS2 = prunedMcs2.size();
		lenght = 0;
		for (list<string>::iterator it = prunedMcs2.begin(); it != prunedMcs2.end(); it++){
			lenght = lenght + it->size();
		}
		lenght = lenght + *numOfMCS2;

		*sdfMCS2 = (char*)malloc(sizeof(char) * lenght);
		memset(*sdfMCS2, 0, sizeof(char) * lenght);
		init = 0;
		for (list<string>::iterator it = prunedMcs2.begin(); it != prunedMcs2.end(); it++){
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
		switch (loggingLevel) {
		case 0: Log::ReportingLevel() = logDEBUG; break;
		case 1: Log::ReportingLevel() = logINFO; break;
		case 2: Log::ReportingLevel() = logWARNING; break;
		case 3: Log::ReportingLevel() = logERROR; break;
		default:
			;
		}

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

		MCS mcs(compoundOne, compoundTwo, userDefinedLowerBound, substructureNumLimit, atomMismatchLowerBound, atomMismatchUpperBound, bondMismatchLowerBound, bondMismatchUpperBound, matchType, runningMode, timeout);
		mcs.calculate();

		//static int cmpOneSize, cmpTwoSize, mSize;
		//cmpOneSize = mcs.getCompoundOne().size();
		//cmpTwoSize = mcs.getCompoundTwo().size();
		//mSize = mcs.size();
		if (runningMode == MCS::DETAIL) {
                    	cout << "Saving output..." << endl;
			string sdfOut;
		        list<vector<size_t> > originalIndex1 = mcs.getFirstOriginalIndice();
		        originalMcs1.clear();
		        originaDelta1.clear();
			for(list<vector<size_t> >::iterator it = originalIndex1.begin();it!=originalIndex1.end();it++){
                                sdfOut = compoundOne.createDissimilarSDFs(*it);
				originaDelta1.push_back(sdfOut);
				sdfOut = compoundOne.createMCSSDFs(*it);
				originalMcs1.push_back(sdfOut);
			}

			list<vector<size_t> > originalIndex2 = mcs.getSecondOriginalIndice();
			originalMcs2.clear();
			originalDelta2.clear();
			for(list<vector<size_t> >::iterator it = originalIndex2.begin();it!=originalIndex2.end();it++){				
                                sdfOut = compoundTwo.createDissimilarSDFs(*it);
				originalDelta2.push_back(sdfOut);
				sdfOut = compoundTwo.createMCSSDFs(*it);
				originalMcs2.push_back(sdfOut);
			}

			list<vector<size_t> > closedIndex1 = mcs.getFirstClosedIndice();
			closedDelta1.clear();
			closedMcs1.clear();
			for(list<vector<size_t> >::iterator it = closedIndex1.begin();it!=closedIndex1.end();it++){
				sdfOut = compoundOne.createDissimilarSDFs(*it);
				closedDelta1.push_back(sdfOut);
				sdfOut = compoundOne.createMCSSDFs(*it);
				closedMcs1.push_back(sdfOut);
			}

			list<vector<size_t> > closedIndex2 = mcs.getSecondClosedIndice();
			closedDelta2.clear();
			closedMcs2.clear();
			for(list<vector<size_t> >::iterator it = closedIndex2.begin();it!=closedIndex2.end();it++){
				sdfOut = compoundTwo.createDissimilarSDFs(*it);
				closedDelta2.push_back(sdfOut);
				sdfOut = compoundTwo.createMCSSDFs(*it);
				closedMcs2.push_back(sdfOut);
				
			}
                        
                        list<vector<size_t> > prunedIndex1 = mcs.getFirstPrunedIndice();
                        prunedDelta1.clear();
                        prunedMcs1.clear();
			for(list<vector<size_t> >::iterator it = prunedIndex1.begin();it!=prunedIndex1.end();it++){
				sdfOut = compoundOne.createDissimilarSDFs(*it);
				prunedDelta1.push_back(sdfOut);
				sdfOut = compoundOne.createMCSSDFs(*it);
				prunedMcs1.push_back(sdfOut);
			}

			list<vector<size_t> > prunedIndex2 = mcs.getSecondPrunedIndice();
			prunedDelta2.clear();
			                        prunedMcs2.clear();
			for(list<vector<size_t> >::iterator it = closedIndex2.begin();it!=closedIndex2.end();it++){
				sdfOut = compoundTwo.createDissimilarSDFs(*it);
				prunedDelta2.push_back(sdfOut);
				sdfOut = compoundTwo.createMCSSDFs(*it);
				prunedMcs2.push_back(sdfOut);
			}
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

		
		MCS::MatchType matchType;
		switch (*matchTypeInt) {
		case 0: matchType = MCS::DEFAULT; break;
		case 1: matchType = MCS::AROMATICITY_SENSETIVE; break;
		case 2: matchType = MCS::RING_SENSETIVE; break;
		default:
			matchType = MCS::RING_SENSETIVE;
		}

		MCS::RunningMode runningMode;
		switch (*runningModeInt) {
		case 0: runningMode = MCS::FAST; break;
		case 1: runningMode = MCS::DETAIL; break;
		default:
			runningMode = MCS::DETAIL;
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
		//static int cmpOneSize, cmpTwoSize, mSize;
		//cmpOneSize = mcs.getCompoundOne().size();
		//cmpTwoSize = mcs.getCompoundTwo().size();
		//mSize = mcs.size();
		if (runningMode == MCS::DETAIL) {
			ofstream myfile;
			string nameDiss, nameMCS, sdfOut;
			cout << "Writing to file..." << endl;
			
                        int i = 0;
                        
			list<vector<size_t> > originalIndex1 = mcs.getFirstOriginalIndice();
                        for(list<vector<size_t> >::iterator it = originalIndex1.begin();it!=originalIndex1.end();it++){
				sdfOut = compoundOne.createDissimilarSDFs(*it);
				nameDiss = "first_mol_original_delta_" +to_string(i) + ".sdf";
				myfile.open(nameDiss, ios_base::out);
				myfile << sdfOut;
				myfile.close();
				cout << "Dissimilar SDFs for molecule 1 written" << endl;

				sdfOut = compoundOne.createMCSSDFs(*it);
				nameMCS = "first_mol_original_MCS_" +to_string(i) + ".sdf";
				myfile.open(nameMCS, ios_base::out);
				myfile << sdfOut;
				myfile.close();
				
                                i++;
				cout << "Original MCS SDFs for molecule 1 written" << endl;
			}
                        list<vector<size_t> > originalIndex2 = mcs.getSecondOriginalIndice();      
                        i = 0;  
                        for(list<vector<size_t> >::iterator it = originalIndex2.begin();it!=originalIndex2.end();it++){
				sdfOut = compoundTwo.createDissimilarSDFs(*it);
				nameDiss = "second_mol_original_delta_" +to_string(i) + ".sdf";
				myfile.open(nameDiss, ios_base::out);
				myfile << sdfOut;
				myfile.close();
				cout << "Dissimilar SDFs for molecule 2 written" << endl;
				sdfOut = compoundTwo.createMCSSDFs(*it);
				nameMCS = "second_mol_original_MCS_" +to_string(i) + ".sdf";
				myfile.open(nameMCS, ios_base::out);
				myfile << sdfOut;
				myfile.close();
				i++;
				cout << "Original MCS SDFs for molecule 2 written" << endl;
			}

			list<vector<size_t> > closedIndex1 = mcs.getFirstClosedIndice();
                        i = 0;  
                        for(list<vector<size_t> >::iterator it = closedIndex1.begin();it!=closedIndex1.end();it++){
			
				sdfOut = compoundOne.createDissimilarSDFs(*it);
                                nameDiss = "first_mol_closed_delta_" +to_string(i) + ".sdf";
				myfile.open(nameDiss, ios_base::out);
				myfile << sdfOut;
				myfile.close();
				cout << "Closed dissimilar SDFs for molecule 1 written" << endl;
				sdfOut = compoundOne.createMCSSDFs(*it);
                                nameMCS = "first_mol_closed_mcs_" +to_string(i) + ".sdf";
				myfile.open(nameMCS, ios_base::out);
				myfile << sdfOut;
				myfile.close();
                                i++;
				cout << "Closed MCS SDFs for molecule 1 written" << endl;
			}

			list<vector<size_t> > closedIndex2 = mcs.getSecondClosedIndice();
                        i = 0;  
                        for(list<vector<size_t> >::iterator it = closedIndex2.begin();it!=closedIndex2.end();it++){
			
				sdfOut = compoundTwo.createDissimilarSDFs(*it);
				nameDiss = "second_mol_closed_delta_" +to_string(i) + ".sdf";
				myfile.open(nameDiss, ios_base::out);
				myfile << sdfOut;
				myfile.close();
				cout << "Closed dissimilar SDFs for molecule 2 written" << endl;
				sdfOut = compoundTwo.createMCSSDFs(*it);
				nameMCS = "second_mol_closed_mcs_" +to_string(i) + ".sdf";
				myfile.open(nameMCS, ios_base::out);
				myfile << sdfOut;
				myfile.close();
				i++;
				cout << "Closed MCS SDFs for molecule 2 written" << endl;
			}
                        
                        list<vector<size_t> > prunedIndex1 = mcs.getFirstPrunedIndice();
                        i = 0;  
                        for(list<vector<size_t> >::iterator it = prunedIndex1.begin();it!=prunedIndex1.end();it++){
				sdfOut = compoundOne.createDissimilarSDFs(*it);
				nameDiss = "first_mol_pruned_delta_" +to_string(i) + ".sdf";
				myfile.open(nameDiss, ios_base::out);
				myfile << sdfOut;
				myfile.close();
				cout << "Pruned dissimilar SDFs for molecule 1 written" << endl;
				sdfOut = compoundOne.createMCSSDFs(*it);
				nameMCS = "first_mol_pruned_mcs_" +to_string(i) + ".sdf";
				myfile.open(nameMCS, ios_base::out);
				myfile << sdfOut;
				myfile.close();
				i++;
				cout << "Pruned MCS SDFs for molecule 1 written" << endl;
			}
                        
                        list<vector<size_t> > prunedIndex2 = mcs.getSecondPrunedIndice();
                         i = 0;  
                        for(list<vector<size_t> >::iterator it = prunedIndex2.begin();it!=prunedIndex2.end();it++){
				sdfOut = compoundTwo.createDissimilarSDFs(*it);
				nameDiss = "second_mol_pruned_delta_" +to_string(i) + ".sdf";
				myfile.open(nameDiss, ios_base::out);
				myfile << sdfOut;
				myfile.close();
				cout << "Pruned dissimilar SDFs for molecule 2 written" << endl;
				sdfOut = compoundTwo.createMCSSDFs(*it);
				nameMCS = "second_mol_pruned_mcs_" +to_string(i) + ".sdf";
				myfile.open(nameMCS, ios_base::out);
				myfile << sdfOut;
				myfile.close();
				i++;
				cout << "Pruned MCS SDFs for molecule 2 written" << endl;
			}
		}
        

/*		stringstream sizeStringStream;
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
		mSizeString = sizeStringStream.str();*/
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
