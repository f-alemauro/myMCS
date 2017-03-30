#include "../include/config.h"
#include "../include/MCS.h"
#include "../include/MCSMap.h"
#include "../include/util.h"
#include <ctime>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;
namespace FMCS {


/**
 * Opening and reading the rules file.
 */
void MCS::readRuleFiles(string fileName)
{
	ifstream ruleFile(fileName.c_str());
	string line;
	stringstream ss;
	cout << "Starting reading rules file!" << endl;
	while (getline(ruleFile, line)) {
		ss << line;
		string atom1 = "", atom2 = "";
		ss >> atom1 >> atom2;
		//if both atoms are not null create a rule
		if (atom1 != "" && atom2 != "") {
			//get the corresponding indexes of the two atoms
			int atomType1 = MCSCompound::Atom::atomTypeMap[getUpper(atom1)];
			int atomType2 = MCSCompound::Atom::atomTypeMap[getUpper(atom2)];
			cout << "New rule: " << atom1 << " -> " << atom2 << endl;
			rules[atomType1][atomType2] = true;
		}
	}
	cout << "Rules file read!" << endl;
}
/**
 * Constructor of the MCS class
 */
MCS::MCS(const MCSCompound& compoundOne, const MCSCompound& compoundTwo,
		size_t userDefinedLowerBound, size_t substructureNumLimit,
		size_t atomMishmatchLower, size_t atomMismatchUpper,
		size_t bondMismatchLower, size_t bondMismatchUpper,
		MatchType matchType, RunningMode runningMode, int timeout)

: compoundOne(compoundOne.size() > compoundTwo.size() ? compoundTwo : compoundOne),
  compoundTwo(compoundOne.size() > compoundTwo.size() ? compoundOne : compoundTwo),
  userDefinedLowerBound(userDefinedLowerBound), substructureNumLimit(substructureNumLimit),
  atomMismatchLowerBound(atomMishmatchLower), atomMismatchUpperBound(atomMismatchUpper),
  bondMismatchLowerBound(bondMismatchLower), bondMismatchUpperBound(bondMismatchUpper),
  matchType(matchType), runningMode(runningMode),
  atomMismatchCurr(0), bondMismatchCurr(0), currSubstructureNum(0),
   bestSize(0), identicalGraph(false),
  haveBeenSwapped(compoundOne.size() > compoundTwo.size() ? true : false) {

	readRuleFiles("rules");
}

/**
 * The core of the MCS computation
 */
void MCS::calculate() {
	vector<size_t> tmpIdVector;
	const MCSCompound::Atom* atomsOne, * atomsTwo;
	size_t atomCountOne, atomCountTwo;
	int resultCount;
	clearResult();

		if (compoundOne.getSdfString() == compoundTwo.getSdfString()) {

			identicalGraph = true;
			if (runningMode == DETAIL) {


				//get the atoms of compound one and their number
				atomsOne = compoundOne.getAtoms();
				atomCountOne = compoundOne.getAtomCount();
				//for each atom, get its corresponding original index ID
				for (int i = 0; i < atomCountOne; ++i)
					tmpIdVector.push_back(atomsOne[i].originalId);
				originalMCS1.push_back(tmpIdVector);

				tmpIdVector.clear();

				//get the atoms of compound one and their number
				atomsTwo = compoundTwo.getAtoms();
				atomCountTwo = compoundTwo.getAtomCount();
				//for each atom, get its corresponding original index ID
				for (int i = 0; i < atomCountTwo; ++i)
					tmpIdVector.push_back(atomsTwo[i].originalId);
				originalMCS2.push_back(tmpIdVector);
			}
	//if the two compound are different calculate the MCS
		}
		else{
		
			//generating and compiling the matrix of ring similarity
			int nRingOne = compoundOne.ringAtomsMap.size();
			int nRingTwo = compoundTwo.ringAtomsMap.size();
			ringsCompareMap = (int **)malloc(nRingOne*sizeof(int**));
			if (ringsCompareMap == NULL)
				return;
			for (int i = 0; i < nRingOne; i++){
				ringsCompareMap[i] = (int *)malloc(nRingTwo*sizeof(int*));
				if (ringsCompareMap[i] == NULL)
					return;
			}


			for (int r = 0; r < nRingOne; r++)
				for (int c = 0; c < nRingTwo; c++)
					if (compoundOne.ringSmartMap.find(r)->second == compoundTwo.ringSmartMap.find(c)->second)
						ringsCompareMap[r][c] = 1;
					else 
						ringsCompareMap[r][c] = 0;
		
			massimo();
		}
	if (runningMode == DETAIL) {
		resultCount = 0;
    	//for each MCS found (there may be more than one MCS, all with the same number of atoms), close the rings
		for (std::list<MCSMap>::const_iterator iMap = bestList.begin(); iMap != bestList.end(); ++iMap) {
			++resultCount;
			//stringstream resultCountStringStream;
			//resultCountStringStream << resultCount;
			//string resultCountString = resultCountStringStream.str();
			//stringstream resultStringStreamOne(
			//compoundOne.subgraph(iMap->getKeyList(), size(), string("fmcs_") + resultCountString));
			//stringstream resultStringStreamTwo(
			//compoundTwo.subgraph(iMap->getValueList(), size(), string("fmcs_") + resultCountString));
			
			const size_t* idArrayOnePtr = iMap->getKeyList();
			const size_t* idArrayTwoPtr = iMap->getValueList();
			int length = size();
			vector<size_t> idxOne, idxTwo;
			for (int i = 0; i < length; ++i) {
				idxOne.push_back(compoundOne.atoms[idArrayOnePtr[i]].originalId);
				idxTwo.push_back(compoundTwo.atoms[idArrayTwoPtr[i]].originalId);
			}
			originalMCS1.push_back(idxOne);
			originalMCS2.push_back(idxTwo);

			if (matchType == RING_SENSETIVE) {
				//serching for idxOne
				for(map<size_t, vector<size_t> >::const_iterator mappa = compoundOne.ringAtomsMap.begin(); mappa!= compoundOne.ringAtomsMap.end(); mappa++){
					if (mappa->second.size()<7){
						int counter = 0;
						for(vector<size_t>::const_iterator atoms= mappa->second.begin();atoms!=mappa->second.end()&& counter<3 ;atoms++)
							if (std::find(idxOne.begin(), idxOne.end(), compoundOne.atoms[*atoms].originalId) != idxOne.end())
								counter ++;
							if (counter>2)
								for(vector<size_t>::const_iterator atoms= mappa->second.begin();atoms!=mappa->second.end();atoms++)
									if (std::find(idxOne.begin(), idxOne.end(), compoundOne.atoms[*atoms].originalId) == idxOne.end())
										idxOne.push_back(compoundOne.atoms[*atoms].originalId);
					}
				}
				//serching for idxTwo
				for(map<size_t, vector<size_t> >::const_iterator mappa = compoundTwo.ringAtomsMap.begin(); mappa!= compoundTwo.ringAtomsMap.end(); mappa++){
					if (mappa->second.size()<7){
						int counter = 0;
						for(vector<size_t>::const_iterator atoms= mappa->second.begin();atoms!=mappa->second.end() && counter<3;atoms++)
							if (std::find(idxTwo.begin(), idxTwo.end(), compoundTwo.atoms[*atoms].originalId) != idxTwo.end())
								counter ++;
							if (counter>2)
								for(vector<size_t>::const_iterator atoms= mappa->second.begin();atoms!=mappa->second.end();atoms++)
									if (std::find(idxTwo.begin(), idxTwo.end(), compoundTwo.atoms[*atoms].originalId) == idxTwo.end())
										idxTwo.push_back(compoundTwo.atoms[*atoms].originalId);
					}
				}
			}            
			closedMCS1.push_back(idxOne);                      
			closedMCS2.push_back(idxTwo);
		}
                
	}
}

	void MCS::massimo() {
		MCSList<size_t> atomListOne = compoundOne.getAtomList();
		MCSList<size_t> atomListTwo = compoundTwo.getAtomList();
		grow(atomListOne, atomListTwo);
	}

	bool MCS::compatible(size_t atomOne, size_t atomTwo,
			size_t& bondMisCount, bool& introducedNewComponent) const {

		MCSList<size_t> targetNeighborMapping;

		const MCSList<size_t>& atomOneNeighborList = compoundOne[atomOne];
		const size_t* atomOneNeighborsPtr = atomOneNeighborList.get();
		size_t atomOneNeighborSize = atomOneNeighborList.size();

		for (size_t i = 0; i < atomOneNeighborSize; ++i) {
			if (currentMapping.containsKey(atomOneNeighborsPtr[i])) {
				targetNeighborMapping.push_back(atomOneNeighborsPtr[i]);
			}
		}

		MCSList<size_t> currNeighborMapping;
		
		const MCSList<size_t>& atomTwoNeighborList = compoundTwo[atomTwo];
		const size_t* atomTwoNeighborsPtr = atomTwoNeighborList.get();
		size_t atomTwoNeighborSize = atomTwoNeighborList.size();
		for (size_t i = 0; i < atomTwoNeighborSize; ++i) {
			size_t k = currentMapping.getKey(atomTwoNeighborsPtr[i]);
			if (k != MCSMap::npos) {
				currNeighborMapping.push_back(k);
			}
		}
		if (!targetNeighborMapping.equals(currNeighborMapping)) {
			return false;
		}
		if (targetNeighborMapping.size() == 0) {
			introducedNewComponent = true;
		}
		size_t targetNeighborMappingSize = targetNeighborMapping.size();
		const size_t* targetNeighborMappingPtr = targetNeighborMapping.get();

		if (matchType == DEFAULT) {
			for (size_t i = 0; i < targetNeighborMappingSize; ++i) {
				size_t n = currentMapping.getValue(targetNeighborMappingPtr[i]);
				const MCSCompound::Bond bondOne = compoundOne.atoms[atomOne].getBond(targetNeighborMappingPtr[i]);
				const MCSCompound::Bond bondTwo = compoundTwo.atoms[atomTwo].getBond(n);
				if (bondOne.bondType != bondTwo.bondType) {
					++bondMisCount;
				}
			}
		} else if (matchType == AROMATICITY_SENSETIVE) {
			for (size_t i = 0; i < targetNeighborMappingSize; ++i) {
				size_t n = currentMapping.getValue(targetNeighborMappingPtr[i]);
				const MCSCompound::Bond bondOne = compoundOne.getAtom(atomOne).getBond(targetNeighborMappingPtr[i]);
				const MCSCompound::Bond bondTwo = compoundTwo.getAtom(atomTwo).getBond(n);
				if (bondOne.isAromatic != bondTwo.isAromatic) {
					++bondMisCount;
				} else if (!bondOne.isAromatic && bondOne.bondType != bondTwo.bondType) {
					++bondMisCount;
				}
			}
		} else { //ring sensetive
			int flagMis = 0;
			for (size_t i = 0; i < targetNeighborMappingSize; ++i) {
				size_t n = currentMapping.getValue(targetNeighborMappingPtr[i]);
				const MCSCompound::Bond bondOne = compoundOne.atoms[atomOne].getBond(targetNeighborMappingPtr[i]);
				const MCSCompound::Bond bondTwo = compoundTwo.atoms[atomTwo].getBond(n);
				
				if (bondOne.isInARing != bondTwo.isInARing)
					++bondMisCount;

				else if (bondOne.isInARing){ //here the two bonds are both in a ring
                                        
					if (bondOne.isAromatic != bondTwo.isAromatic){ //if the aromaticity is different -> miscount++
						//cout << "Aromaticity error!" << endl;						
						//getchar();
						++bondMisCount;
					}
                                        
                                        
					
                                        else { // if both are in an aromatic or both in an aliphatic ring                                            
                                            
                                            if (bondOne.bondType != bondTwo.bondType) // this check is only for aliphatics                                                     
                                                ++bondMisCount;
                                            
                                            else {    //else bondtype check
                                                //if the two bonds have the same aromaticity, then check the rings...
                                            //check if the two rings have the same atoms.
                                                
                                                size_t oneAtomA = bondOne.firstAtom;
						size_t oneAtomB = bondOne.secondAtom;
						size_t secondAtomA = bondTwo.firstAtom;
						size_t secondAtomB = bondTwo.secondAtom;
						vector<size_t> ringIntersectionOne, ringIntersectionSecond;
                                                
	
						//extract the id of the ring the atoms are in
						
						vector<size_t> ringsIDoneA = compoundOne.atoms[oneAtomA].ringId;
                                                vector<size_t> ringsIDoneB = compoundOne.atoms[oneAtomB].ringId;
						vector<size_t> ringsIDsecondA = compoundTwo.atoms[secondAtomA].ringId;
						vector<size_t> ringsIDsecondB = compoundTwo.atoms[secondAtomB].ringId;
					
						//computing the intersection between the two lists of ring IDs (we extract the two lists for the two bonds of rings that contain both the atoms of a bond)
						sort(ringsIDoneA.begin(), ringsIDoneA.end());
						sort(ringsIDoneB.begin(), ringsIDoneB.end());
						sort(ringsIDsecondA.begin(), ringsIDsecondA.end());
						sort(ringsIDsecondB.begin(), ringsIDsecondB.end());

						set_intersection(ringsIDoneA.begin(), ringsIDoneA.end(), ringsIDoneB.begin(), ringsIDoneB.end(), back_inserter(ringIntersectionOne));
						set_intersection(ringsIDsecondA.begin(), ringsIDsecondA.end(), ringsIDsecondB.begin(), ringsIDsecondB.end(), back_inserter(ringIntersectionSecond));
						
						
						int i, j, k, z;
						vector<size_t> listOfAtomsInRingOne, listOfAtomsInRingSecond;
	
						for (i = 0; i < ringIntersectionOne.size(); i++) {
							size_t ID = ringIntersectionOne[i];
							for (j = 0; j < ringIntersectionSecond.size(); j++) {
								size_t IDSecond = ringIntersectionSecond[j];
								if (ringsCompareMap[ID][IDSecond] == 1){
									flagMis = 1;
								}
							}
						}
						if (flagMis == 0){
							++bondMisCount;
						}
                                                
                                            } //else bondtype check    
                                            
					} //END if "both rings are aromatic"
				} //END if "both are in a ring!"
				else if (bondOne.bondType != bondTwo.bondType) {
					++bondMisCount;
				}
			}
		}

		return true;
	}

	size_t MCS::top(MCSList<size_t>& atomList) {

		size_t bestCandidateAtom = atomList.front();
		size_t candidateAtom = static_cast<size_t> (-1);
		size_t i, bestIdx = 0, candidateIdx;
		size_t atomListSize = atomList.size();
		const size_t* atomPtr = atomList.get();
		for (i = 0; i < atomListSize; ++i) {
			if (compoundOne[atomPtr[i]].size() > compoundOne[bestCandidateAtom].size()) {
				bestCandidateAtom = atomPtr[i];
				bestIdx = i;
			}
			const MCSList<size_t>& neighborAtomList = compoundOne[atomPtr[i]];
			size_t neighborAtomListSize = neighborAtomList.size();
			const size_t* neighborAtomPtr = neighborAtomList.get();
			for (size_t j = 0; j < neighborAtomListSize; ++j) {
				if (currentMapping.containsKey(neighborAtomPtr[j])) {
					if (candidateAtom == static_cast<size_t> (-1) || compoundOne[atomPtr[i]].size() > compoundOne[candidateAtom].size()) {
						candidateAtom = atomPtr[i];
						candidateIdx = i;
						break;
					}
				}
			}
		}

		if (candidateAtom == static_cast<size_t> (-1)) {
			atomList.eraseIdx(bestIdx);
			return bestCandidateAtom;
		}

		if (candidateAtom != static_cast<size_t> (-1)) {
			atomList.eraseIdx(candidateIdx);
		}
		return candidateAtom;
	}

	void MCS::boundary() {
		if (runningMode == FAST) {
			if (currentMapping.size() > bestSize) {
				if (atomMismatchCurr < atomMismatchLowerBound || bondMismatchCurr < bondMismatchLowerBound) {
					return;
				}
				bestSize = currentMapping.size();
			}
		} else {
			if (currentMapping.size() == size()) {
				if (atomMismatchCurr < atomMismatchLowerBound || bondMismatchCurr < bondMismatchLowerBound) {
					return;
				}
				bestList.push_back(currentMapping);
			} else if (currentMapping.size() > size()) {

				if (atomMismatchCurr < atomMismatchLowerBound || bondMismatchCurr < bondMismatchLowerBound) {
					return;
				}
				bestList.clear();
				bestList.push_back(currentMapping);
			}
		}
	}

	void MCS::grow(MCSList<size_t>& atomListOne, MCSList<size_t>& atomListTwo) {
		MCSList<size_t> atomListOneCopy = atomListOne;
		MCSList<size_t> atomListTwoCopy = atomListTwo;
		MCSList<size_t> atomListOneDegrees;
		MCSList<size_t> atomListTwoDegrees;

		size_t atomListOneSize = atomListOne.size();
		const size_t* atomListOnePtr = atomListOne.get();
		for (size_t i = 0; i < atomListOneSize; ++i) {
			if (!currentMapping.containsKey(atomListOnePtr[i])) {
				int degree = 0;
				const MCSList<size_t>& neighborAtomList = compoundOne.atoms[atomListOnePtr[i]].neighborAtoms;
				size_t neighborAtomListSize = neighborAtomList.size();
				const size_t* neighborAtomPtr = neighborAtomList.get();
				for (size_t j = 0; j < neighborAtomListSize; ++j) {
					if (currentMapping.containsKey(neighborAtomPtr[j])) {
						++degree;
					}
				}
				atomListOneDegrees.push_back(degree);
			}
		}

		size_t atomListTwoSize = atomListTwo.size();
		const size_t* atomListTwoPtr = atomListTwo.get();
		for (size_t i = 0; i < atomListTwoSize; ++i) {
			if (!currentMapping.containsValue(atomListTwoPtr[i])) {
				int degree = 0;
				const MCSList<size_t>& neighborAtomList = compoundTwo.atoms[atomListTwoPtr[i]].neighborAtoms;
				size_t neighborAtomListSize = neighborAtomList.size();
				const size_t* neighborAtomPtr = neighborAtomList.get();
				for (size_t j = 0; j < neighborAtomListSize; ++j) {
					if (currentMapping.containsValue(neighborAtomPtr[j])) {
						++degree;
					}
				}
				atomListTwoDegrees.push_back(degree);
			}
		}

		size_t currentBound = currentMapping.size();
		size_t atomListOneDegreesSize = atomListOneDegrees.size();
		const size_t* atomListOneDegreesPtr = atomListOneDegrees.get();
		for (size_t i = 0; i < atomListOneDegreesSize; ++i) {
			if (atomListTwoDegrees.contains(atomListOneDegreesPtr[i])) {
				++currentBound;
				atomListTwoDegrees.erase(atomListOneDegreesPtr[i]);
			}
		}

		if (runningMode == FAST) {
			if (currentBound < userDefinedLowerBound || currentBound <= bestSize) {
				return;
			}
		} else {

			if (currentBound < userDefinedLowerBound || currentBound < size()) {
				return;
			}
		}


		while (true) {

			if (atomListOneCopy.empty() || atomListTwoCopy.empty()) { // atomListTwoCopy
				//cout<<"Empty atomListCopy"<<endl;
				boundary();
				return;
			}

			size_t topCandidateAtom = top(atomListOneCopy);
			size_t atomListTwoSize = atomListTwoCopy.size(); // atomListTwoCopy
			const size_t* atomListTwoPtr = atomListTwoCopy.get(); // atomListTwoCopy
			for (size_t i = 0; i < atomListTwoSize; ++i) {

				bool atomMismatched = false;
				bool atomMismatchAllowed = true;
				int atom1 = compoundOne.getAtom(topCandidateAtom).atomType;
				int atom2 = compoundTwo.getAtom(atomListTwoPtr[i]).atomType;
				if (atom1 != atom2) {

					if (rules.count(atom1) > 0) {
						if (!rules[atom1][atom2]) {
							atomMismatchAllowed = false;
						}
					} else if (rules.count(atom2) > 0) {
						if (!rules[atom2][atom1]) {
							atomMismatchAllowed = false;
						}
					}
					++atomMismatchCurr;
					atomMismatched = true;
				}

				if (!(atomMismatchCurr > atomMismatchUpperBound) && atomMismatchAllowed) {

					size_t bondMisCount = 0;
					bool introducedNewComponent = false;
					if (compatible(topCandidateAtom, atomListTwoPtr[i], bondMisCount, introducedNewComponent)) {

						if (!(bondMismatchCurr + bondMisCount > bondMismatchUpperBound)) {

							bondMismatchCurr = bondMismatchCurr + bondMisCount;

							if (introducedNewComponent) {
								++currSubstructureNum;
							}

							if (!(currSubstructureNum > substructureNumLimit)) {

								currentMapping.push_back(topCandidateAtom, atomListTwoPtr[i]);

								atomListTwo.erase(atomListTwoPtr[i]);

								grow(atomListOneCopy, atomListTwo);


								atomListTwo.push_back(atomListTwoPtr[i]);
								currentMapping.pop_back();

							}
							if (introducedNewComponent) {
								--currSubstructureNum;
							}

							bondMismatchCurr = bondMismatchCurr - bondMisCount;
						}
					}

				}
				if (atomMismatched) {
					--atomMismatchCurr;
				}
			}
		}
	}
	/**
	 * Reset all the result from the previous computation.
	 */
	void MCS::clearResult() {
		bestSize = 0;
		bestList.clear();
		identicalGraph = false;
		currentMapping.clear();
		sdfSet1.clear();
		sdfSet2.clear();
//		timeoutStop = false;
		//_isTimeout = false;
	}
}
