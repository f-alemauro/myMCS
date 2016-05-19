//#include "../config.h"

#include "MCSMap.h"
#include "MCSCompound.h"
#include "MCSRingDetector.h"
#include "util.h"

#ifdef HAVE_LIBOPENBABEL

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>

using namespace OpenBabel;

#endif

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <map>
#include <string>
#include <stdio.h>
#include <list>
#include <algorithm>
#include <iterator>



using namespace std;

namespace FMCS {

const char MCSCompound::Atom::elements[111][3] = {"Ru", "Re", "Rf", "Rg", "Ra", "Rb", "Rn", "Rh", "Be", "Ba", "Bh", "Bi", "Bk", "Br", "H", "P", "Os", "Ge", "Gd", "Ga", "Pr", "Pt", "Pu", "C", "Pb", "Pa", "Pd", "Cd", "Po", "Pm", "Hs", "Ho", "Hf", "Hg", "He", "Md", "Mg", "K", "Mn", "O", "Mt", "S", "W", "Zn", "Eu", "Zr", "Er", "Ni", "No", "Na", "Nb", "Nd", "Ne", "Np", "Fr", "Fe", "Fm", "B", "F", "Sr", "N", "Kr", "Si", "Sn", "Sm", "V", "Sc", "Sb", "Sg", "Se", "Co", "Cm", "Cl", "Ca", "Cf", "Ce", "Xe", "Tm", "Cs", "Cr", "Cu", "La", "Li", "Tl", "Lu", "Lr", "Th", "Ti", "Te", "Tb", "Tc", "Ta", "Yb", "Db", "Dy", "Ds", "At", "I", "U", "Y", "Ac", "Ag", "Ir", "Am", "Al", "As", "Ar", "Au", "Es", "In", "Mo"};

map<string, int> MCSCompound::Atom::atomTypeMap;

bool MCSCompound::Atom::atomTypeMapInitialized = atomTypeMapInit();

bool MCSCompound::Atom::atomTypeMapInit() {
	for (int i = 0; i < 111; ++i) {
		stringstream symbolStringStream;
		//symbolStringStream.width(3);
		//symbolStringStream.setf(ios::left, ios::adjustfield);
		symbolStringStream << elements[i];
		string atomSymbol = getUpper(symbolStringStream.str());
		atomTypeMap[atomSymbol] = i + 1;
	}
	return true;
}

MCSCompound::MCSCompound()
: bondCount(0), atomCount(0), atoms(NULL), bonds(NULL) {
}

MCSCompound::~MCSCompound() {
	if (atoms != NULL) {
		delete[] atoms;
		atoms = NULL;
	}
	if (bonds != NULL) {
		delete[] bonds;
		atoms = NULL;
	}
}

#ifdef HAVE_LIBOPENBABEL

MCSCompound::MCSCompound(const MCSCompound& other)
: bondCount(0), atomCount(0), atoms(NULL), bonds(NULL), SdfContentString(other.SdfContentString), SmiContentString(other.SmiContentString) {

	if (other.atoms != NULL) {
		atoms = new Atom[other.atomCount];
		memcpy(atoms, other.atoms, sizeof (Atom) * other.atomCount);
		atomCount = other.atomCount;
	}
	if (other.bonds != NULL) {
		bonds = new Bond[other.bondCount];
		memcpy(bonds, other.bonds, sizeof (Bond) * other.bondCount);
		bondCount = other.bondCount;
	}
}

#else

MCSCompound::MCSCompound(const MCSCompound& other)
: bondCount(0), atomCount(0), atoms(NULL), bonds(NULL), SdfContentString(other.SdfContentString) {

	if (other.atoms != NULL) {
		atoms = new Atom[other.atomCount];
		memcpy(atoms, other.atoms, sizeof (Atom) * other.atomCount);
		atomCount = other.atomCount;
	}
	if (other.bonds != NULL) {
		bonds = new Bond[other.bondCount];
		memcpy(bonds, other.bonds, sizeof (Bond) * other.bondCount);
		bondCount = other.bondCount;
	}
}

#endif

const MCSCompound& MCSCompound::operator=(const MCSCompound& that) {
	if (this == &that) {
		return *this;
	}

	if (atoms != NULL) {
		delete[] atoms;
		atoms = NULL;
	}
	if (bonds != NULL) {
		delete[] bonds;
		bonds = NULL;
	}

	bondCount = 0;
	atomCount = 0;
	SdfContentString = that.SdfContentString;

#ifdef HAVE_LIBOPENBABEL
	SmiContentString = that.SmiContentString;
#endif

	if (that.atoms != NULL) {
		atoms = new Atom[that.atomCount];
		memcpy(atoms, that.atoms, sizeof (Atom) * that.atomCount);
		atomCount = that.atomCount;
	}
	if (that.bonds != NULL) {
		bonds = new Bond[that.bondCount];
		memcpy(bonds, that.bonds, sizeof (Bond) * that.bondCount);
		bondCount = that.bondCount;
	}

	return *this;
}
#ifdef HAVE_LIBOPENBABEL

void MCSCompound::read(const string& sdfString, ReadType type) {

	switch (type) {
	case SMI:
		parseSMI(sdfString);
		break;
	case SDF:
		parseSDF(sdfString);
		break;
	}

	for (int i = 0; i < bondCount; ++i) {
		atoms[bonds[i].firstAtom].neighborAtoms.push_back(bonds[i].secondAtom);
		atoms[bonds[i].firstAtom].neighborBonds.push_back(&bonds[i]);
		atoms[bonds[i].secondAtom].neighborAtoms.push_back(bonds[i].firstAtom);
		atoms[bonds[i].secondAtom].neighborBonds.push_back(&bonds[i]);
	}

}
#else

void MCSCompound::read(const std::string& sdfString) {
	molB = parseSDF(sdfString);
	size_t newAtomCount = 0;
	size_t newBondCount = 0;

	for (int i = 0; i < bondCount; ++i) {
		atoms[bonds[i].firstAtom].neighborAtoms.push_back(bonds[i].secondAtom);
		atoms[bonds[i].firstAtom].neighborBonds.push_back(&bonds[i]);
		atoms[bonds[i].secondAtom].neighborAtoms.push_back(bonds[i].firstAtom);
		atoms[bonds[i].secondAtom].neighborBonds.push_back(&bonds[i]);
	}

	MCSRingDetector ringDector(*this);
	ringDector.detect();
}

#endif

string MCSCompound::subgraph(const size_t* index, size_t indexLength, const string& newCompoundName) const {

	stringstream content(this->SdfContentString);
	string res;
	string line;
	vector<string> lines, subLines;

	while (getline(content, line)) {
		lines.push_back(line);
	}

	subLines.push_back(lines[0] + "_" + newCompoundName);
	subLines.push_back("FMCS substructure");

	subLines.push_back("Auto Generated from FMCS");
	subLines.push_back("counts line");

	stringstream atomNumss(lines[3].substr(0, 3));
	stringstream bondNumss(lines[3].substr(3, 3));
	int numAtoms, numbonds;
	atomNumss >> numAtoms;
	bondNumss >> numbonds;

	MCSMap atomMap;
	for (int i = 0; i < indexLength; ++i) {
		subLines.push_back(lines[4 + index[i]]);
		atomMap.push_back(index[i], i);
	}
	int numResBonds = 0;
	for (int i = 0; i < numbonds; ++i) {
		string bondline = lines[4 + numAtoms + i];
		stringstream beginAtomss(bondline.substr(0, 3)), endAtomss(bondline.substr(3, 3));
		int beginAtom, endAtom;
		beginAtomss >> beginAtom;
		endAtomss >> endAtom;
		if (atomMap.containsKey(beginAtom - 1) && atomMap.containsKey(endAtom - 1)) {
			size_t newBeginAtom = atomMap.getValue(beginAtom - 1) + 1;
			size_t newEndAtom = atomMap.getValue(endAtom - 1) + 1;
			stringstream resultStringStream;
			resultStringStream.width(3);
			resultStringStream << newBeginAtom;
			resultStringStream.width(3);
			resultStringStream << newEndAtom << bondline.substr(6);
			subLines.push_back(resultStringStream.str());
			++numResBonds;
		}
	}

	stringstream resLine3ss;
	resLine3ss.width(3);
	resLine3ss << indexLength;
	resLine3ss.width(3);
	resLine3ss << numResBonds << lines[3].substr(6);
	subLines[3] = resLine3ss.str();

	for (vector<string>::iterator i = subLines.begin(); i != subLines.end(); ++i) {
		res += *i;
		res += "\n";
	}
	res += "M END\n";
	res += "$$$$";
	return res;
}

const MCSCompound::Bond& MCSCompound::Atom::getBond(int anotherAtom) const {
	return *neighborBonds[neighborAtoms.where(anotherAtom)];
}

size_t MCSCompound::getNeighborID(size_t e, size_t me) const {
	size_t other;
	if (bonds[e].firstAtom == me) {
		other = bonds[e].secondAtom;
	} else if (bonds[e].secondAtom == me) {
		other = bonds[e].firstAtom;
	} else {
		other = static_cast<size_t> (-1);
	}
	return other;
}

size_t MCSCompound::addNewRingAtom(std::string name) {
	newRingAtoms.push_back(Atom(atomCount + newRingAtoms.size() + 1, 0, newRingAtoms.size(), name, true, std::vector<size_t>()));
}

MCSList<size_t> MCSCompound::getAtomList() const {
	MCSList<size_t> l;
	for (size_t i = 0; i < atomCount; ++i) {
		l.push_back(i);
	}
	return l;
}

#ifdef HAVE_LIBOPENBABEL

void MCSCompound::parseSMI(const string& smiString) {

	stringstream smiss, outss;
	smiss << smiString;
	OBConversion conv(&smiss, &outss);
	if (conv.SetInAndOutFormats("SMI", "SDF")) {
		OBMol mol;
		conv.Read(&mol);
		conv.Write(&mol);
	}
	parseSDF(outss.str().c_str());
}
#endif

#ifdef HAVE_LIBOPENBABEL

void MCSCompound::parseSDF(const string& sdfString) {

	stringstream ss;
	stringstream ssSMI;
	stringstream ssSDF;

	ss << sdfString;
	ss >> compoundName;

	OBConversion conv(&ss, &ssSDF);

	if (conv.SetInAndOutFormats("SDF", "SDF")) {
		OBMol mol;
		if (conv.Read(&mol)) {
			mol.DeleteHydrogens();
			atoms = new Atom[mol.NumAtoms()];
			bonds = new Bond[mol.NumBonds()];
			int i = 0;

			FOR_ATOMS_OF_MOL(atom, mol) {
				atoms[i] = Atom(atom->GetIdx() - 1, atom->GetAtomicNum());
				++i;
			}
			i = 0;

			FOR_BONDS_OF_MOL(bond, mol) {
				int bondType = 0;
				if (bond->IsKSingle()) {
					bondType = 1;
				} else if (bond->IsKDouble()) {
					bondType = 2;
				} else if (bond->IsKTriple()) {
					bondType = 3;
				} else {
					delete[] atoms;
					atoms = NULL;
					delete[] bonds;
					atoms = NULL;
					throw InvalidBondTypeException();
				}
				bonds[i] = Bond(bond->GetIdx(), bond->GetBeginAtomIdx() - 1, bond->GetEndAtomIdx() - 1, bondType, bond->IsAromatic(), bond->IsInRing());
				++i;
			}

			this->bondCount = mol.NumBonds();
			this->atomCount = mol.NumAtoms();

			conv.Write(&mol);
			SdfContentString = ssSDF.str();
		}
	}

	OBConversion conv2(&ssSDF, &ssSMI);
	if (conv2.SetInAndOutFormats("SDF", "SMI")) {
		OBMol mol;
		conv2.Read(&mol);
		conv2.Write(&mol);
		SmiContentString = ssSMI.str();
	}
}
#else

string MCSCompound::deleteHydrogens(const string& sdf, vector<size_t>& originalIds, molBlocks& molB) {
	stringstream originalStringStream;
	originalStringStream << sdf;
	string compoundNameLine;
	string informationLine;
	string commentLine;
	getline(originalStringStream, compoundNameLine);
	getline(originalStringStream, informationLine);
	getline(originalStringStream, commentLine);

	string oldCountsLine;
	getline(originalStringStream, oldCountsLine);

	string atomCountString = oldCountsLine.substr(0, 3);
	string bondCountString = oldCountsLine.substr(3, 3);
	int oldAtomCount = atoi(atomCountString.c_str());
	int oldBondCount = atoi(bondCountString.c_str());

	string newAtomBlock;
	int *newAtomIndex = new int[oldAtomCount];
	int newAtomCount = 0;
	for (int i = 0; i < oldAtomCount; ++i) {
		string atomBlockLine;
		getline(originalStringStream, atomBlockLine);
		string atomSymbolRawString = atomBlockLine.substr(31, 3);
		stringstream rawStringStream(atomSymbolRawString);
		string atomSymbolString;
		rawStringStream >> atomSymbolString;
		if (atomSymbolString != "H") {
			newAtomBlock += atomBlockLine;
			newAtomBlock += "\n";
			newAtomIndex[i] = newAtomCount + 1;
			++newAtomCount;
			originalIds.push_back(i + 1);
		} else {
			newAtomIndex[i] = -1;
		}
	}
	string newBondBlock;
	int newBondCount = 0;
	for (int i = 0; i < oldBondCount; ++i) {
		string oldBondBlockLine;
		getline(originalStringStream, oldBondBlockLine);
		int oldFirstAtomIndex = atoi(oldBondBlockLine.substr(0, 3).c_str());
		int oldSecondAtomIndex = atoi(oldBondBlockLine.substr(3, 3).c_str());
		int newFirstAtomIndex = newAtomIndex[oldFirstAtomIndex - 1];
		int newSecondAtomIndex = newAtomIndex[oldSecondAtomIndex - 1];
		if (newFirstAtomIndex != -1 && newSecondAtomIndex != -1) {
			stringstream newBondBlockLineStringStream;
			newBondBlockLineStringStream.width(3);
			newBondBlockLineStringStream << newFirstAtomIndex;
			newBondBlockLineStringStream.width(3);
			newBondBlockLineStringStream << newSecondAtomIndex;
			newBondBlockLineStringStream << oldBondBlockLine.substr(6);
			newBondBlock += newBondBlockLineStringStream.str();
			newBondBlock += "\n";
			++newBondCount;
		}
	}

	string chgISOlinesTMP;
	string chgISOlines;
	size_t posCHG;
	size_t posISO;

	while (originalStringStream.good()) {
		getline(originalStringStream, chgISOlinesTMP);
		posCHG = chgISOlinesTMP.find("M  CHG");
		if (posCHG != string::npos) {
			chgISOlines += chgISOlinesTMP + "\n";
		}
		posISO = chgISOlinesTMP.find("M  ISO");
		if (posISO != string::npos) {
			chgISOlines += chgISOlinesTMP + "\n";
		}

	}

	stringstream newCountLineStringStream;
	newCountLineStringStream.width(3);
	newCountLineStringStream << newAtomCount;
	newCountLineStringStream.width(3);
	newCountLineStringStream << newBondCount;
	newCountLineStringStream << oldCountsLine.substr(6);
	string newCountLine = newCountLineStringStream.str();
	delete newAtomIndex;
	newAtomIndex = NULL;

	string newSDFString = compoundNameLine + "\n"
			+ informationLine + "\n"
			+ commentLine + "\n"
			+ newCountLine + "\n"
			+ newAtomBlock
			+ newBondBlock
			+ "M END\n"
			+ "$$$$";

	string infoString = compoundNameLine + "\n"
			+ informationLine + "\n"
			+ commentLine;
	molB.infoBlock = infoString;
	molB.atomBlock = newAtomBlock;
	molB.bondBlock = newBondBlock;
	molB.chgISO = chgISOlines;

	return newSDFString;
}

MCSCompound::molBlocks MCSCompound::parseSDF(const string& sdf) {
	stringstream sdfStringStream;
	vector<size_t> originalIds;
	molBlocks molB;
	SdfContentString = deleteHydrogens(sdf, originalIds, molB);

	sdfStringStream << SdfContentString;
	string compoundNameLine;
	string informationLine;
	string commentLine;
	getline(sdfStringStream, compoundNameLine);
	getline(sdfStringStream, informationLine);
	getline(sdfStringStream, commentLine);
	string countsLine;
	getline(sdfStringStream, countsLine);
	string atomCountString = countsLine.substr(0, 3);
	string bondCountString = countsLine.substr(3, 3);

	atomCount = atoi(atomCountString.c_str());
	bondCount = atoi(bondCountString.c_str());

	atoms = new Atom[atomCount];
	bonds = new Bond[bondCount];
	cout << "Reading atoms..." << endl;

	for (size_t i = 0; i < atomCount; ++i) {
		string atomBlockLine;
		getline(sdfStringStream, atomBlockLine);

		string atomSymbolRawString = atomBlockLine.substr(31, 3);
		stringstream rawStringStream(atomSymbolRawString);
		string atomSymbol;
		rawStringStream >> atomSymbol;
		std::vector<size_t> id;
		atoms[i] = Atom(i, originalIds[i], MCSCompound::Atom::atomTypeMap[getUpper(atomSymbol)], atomSymbol, false, id);
	}

	for (size_t i = 0; i < bondCount; ++i) {
		string bondBlockLine;
		getline(sdfStringStream, bondBlockLine);
		int firstAtom = -1, secondAtom = -1, bondType = -1;
		firstAtom = atoi(bondBlockLine.substr(0, 3).c_str()) - 1;
		secondAtom = atoi(bondBlockLine.substr(3, 3).c_str()) - 1;
		bondType = atoi(bondBlockLine.substr(6, 3).c_str());
		bonds[i] = Bond(i, firstAtom, secondAtom, bondType, false, false);
	}
	return molB;
}

#endif

const MCSCompound::Bond* MCSCompound::getBond(size_t firstAtom, size_t secondAtom) const {

	for (int i = 0; i < bondCount; ++i) {
		if ((bonds[i].firstAtom == firstAtom && bonds[i].secondAtom == secondAtom)
				|| (bonds[i].firstAtom == secondAtom && bonds[i].secondAtom == firstAtom)) {
			return bonds + i;
		}
	}

	return NULL;
}

string MCSCompound::MCS2SDF(vector<size_t> mcs, bool isMCS){
	list<std::vector<size_t> > bondsList = strings2int(molB.bondBlock, true);
		list<std::vector<size_t> > subgraph;
		std::vector<size_t> listOfAtoms;
		string finalSDF;
		for (list<vector<size_t> >::iterator i = bondsList.begin(); i != bondsList.end(); ++i) {
			if (*(i->begin()) == 0) {
				subgraph.clear();
				listOfAtoms.clear();
				size_t atom1, atom2;
				atom1 = *(i->begin() + 1);
				atom2 = *(i->begin() + 2);
				if (((std::find(mcs.begin(), mcs.end(), atom1) != mcs.end()) && (std::find(mcs.begin(), mcs.end(), atom2) != mcs.end()))==isMCS) {
					ricerca(atom1, subgraph, listOfAtoms, bondsList, mcs, isMCS);

					subgraph = updateBondList(listOfAtoms,subgraph);

					string bondString = generateBondString(subgraph);
					string atomString = generateAtomString(listOfAtoms);
					stringstream infoLiness;
                                        string CHGstring = evaluateCHGs(listOfAtoms);
                                        
                                        cout<<CHGstring<<endl;
					infoLiness  << " "<<listOfAtoms.size()<<" "<<subgraph.size()<<"  0  0  0  0            999 V2000"<<endl;
                                        
					finalSDF += molB.infoBlock+"\n";
					finalSDF += infoLiness.str();;
					finalSDF += atomString;
					finalSDF += bondString;
					finalSDF += "M  END\n$$$$\n";

				} else
					*(i->begin()) = 1;
			}
		}
		return finalSDF;

}

string MCSCompound::createMCSSDFs(vector<size_t> mcs) {
	return MCS2SDF(mcs,true);

}

string MCSCompound::createDissimilarSDFs(vector<size_t> mcs) {
	return MCS2SDF(mcs,false);
}


list<vector<size_t> > MCSCompound::updateBondList(std::vector<size_t> listOfAtoms,list<std::vector<size_t> > listOfSubgraph){
	std::sort(listOfAtoms.begin(), listOfAtoms.end());
	std::string line;
	int atomLinesCounter = 1;
	for (vector<size_t>::iterator k = listOfAtoms.begin(); k != listOfAtoms.end(); ++k) {
		for (list<vector<size_t> >::iterator resultI = listOfSubgraph.begin(); resultI != listOfSubgraph.end(); ++resultI) {
			if ((*resultI)[1] == *k){
				cout<<"Change: "<<*k<<"-->"<<atomLinesCounter<<endl;
				(*resultI)[1] = atomLinesCounter;
			}
			if ((*resultI)[2] == *k){
				cout<<"Change: "<<*k<<"-->"<<atomLinesCounter<<endl;
				(*resultI)[2] = atomLinesCounter;
			}
		}
		atomLinesCounter++;

	}
	return listOfSubgraph;
}

string MCSCompound::evaluateCHGs(std::vector<size_t> mcs){
    std::stringstream chgLine_tmp;
    if (!molB.chgISO.empty()){
        cout<<"chgISO block is: "<<molB.chgISO.c_str()<<endl;
        std::istringstream charges (molB.chgISO);
        std::string line;    
        
        while (std::getline(charges, line)) {
            std::size_t atomNumber = atoi(line.substr(8,1).c_str());
            cout<<"the atomNumber in the chg or iso line is: "<< atomNumber<<endl; 
		
    std::istringstream iss(line.substr(12,line.length()));
    std::string token;
    int atom_found = 0;
    while (std::getline(iss, token, ' '))
    {
        if (token == ""){
            if (atom_found == 1){
                chgLine_tmp << " ";
            }
        } else if (token != ""){
            if (atom_found == 1){
                chgLine_tmp << token;
                atom_found = 0;
            }
            else if ((std::find(mcs.begin(), mcs.end(), atoi(token.c_str())) != mcs.end())){
                atom_found = 1;
                chgLine_tmp << token << " ";
            }
        }    
    }
           
    }
    }
    return chgLine_tmp.str();
}
string MCSCompound::generateBondString(list<std::vector<size_t> > listOfSubgraph){
	std::ostringstream bondStringSS;
	for (list<vector<size_t> >::iterator j = listOfSubgraph.begin(); j != listOfSubgraph.end(); ++j){
		if ((*((*j).begin()+1) < 10) && (*((*j).begin()+2) < 10)) {
			bondStringSS<<"  ";
			std::copy((*j).begin()+1, (*j).end()-1, std::ostream_iterator<int>(bondStringSS, "  "));
			bondStringSS <<(*j).back()<<"\n";
		}else if ((*((*j).begin()+1) < 10) && (*((*j).begin()+2) >= 10)){
			bondStringSS<<"  "<<*((*j).begin()+1)<<" "<<*((*j).begin()+2)<<"  ";
			std::copy((*j).begin()+3, (*j).end()-1, std::ostream_iterator<int>(bondStringSS, "  "));
			bondStringSS <<(*j).back()<<"\n";
		}else if ((*((*j).begin()+1) >= 10) && (*((*j).begin()+2) < 10)){
			bondStringSS<<" "<<*((*j).begin()+1)<<"  "<<*((*j).begin()+2)<<"  ";
			std::copy((*j).begin()+3, (*j).end()-1, std::ostream_iterator<int>(bondStringSS, "  "));
			bondStringSS <<(*j).back()<<"\n";
		}else{
			bondStringSS<<" "<<*((*j).begin()+1)<<" "<<*((*j).begin()+2)<<"  ";
			std::copy((*j).begin()+3, (*j).end()-1, std::ostream_iterator<int>(bondStringSS, "  "));
			bondStringSS <<(*j).back()<<"\n";
		}
	}
	return bondStringSS.str();
}

string MCSCompound::generateAtomString(std::vector<size_t> listOfAtoms){
	string atomString;
	std::sort(listOfAtoms.begin(), listOfAtoms.end());
	std::istringstream fB(molB.atomBlock);
	std::string line;
	vector<size_t>::iterator k = listOfAtoms.begin();
	int lineCounter = 1;
	while(std::getline(fB, line)){
		if ((*k)==lineCounter++){
			atomString+=line+"\n";
			k++;
		}
	}
	return atomString;
}

void MCSCompound::ricerca(size_t atomTarget, std::list<std::vector<size_t> >& result, std::vector<size_t>& resultAtomList, std::list<std::vector<size_t> >& bondsList, vector<size_t> mcs, bool isMCS) {
	if (std::find(resultAtomList.begin(), resultAtomList.end(), atomTarget) == resultAtomList.end()){
		resultAtomList.push_back(atomTarget);
	}
	list<vector<size_t> >::iterator i = bondsList.begin();
	while (i != bondsList.end()) {
		if (*(i->begin()) == 0) {
			size_t atom1, atom2;
			atom1 = *(i->begin() + 1);
			atom2 = *(i->begin() + 2);
			size_t bondMate;
			if (((std::find(mcs.begin(), mcs.end(), atom1) != mcs.end()) && (std::find(mcs.begin(), mcs.end(), atom2) != mcs.end()))==isMCS){
				if (atom1 == atomTarget) {
					bondMate = atom2;
					*(i->begin()) = 1;
					result.push_back(*i);
					ricerca(bondMate, result, resultAtomList,bondsList, mcs, isMCS);
					i++;
				} else if (atom2 == atomTarget) {
					bondMate = atom1;
					*(i->begin()) = 1;
					result.push_back(*i);
					ricerca(bondMate, result,resultAtomList, bondsList, mcs, isMCS);
					i++;
				} else
					i++;
			} else {//la riga non va considerata
				*(i->begin()) = 1;
				i++;
			}
		} else
			i++;
	}
}
}
