
#include "util.h"

#include <cctype>
#include <string>

using namespace std;

string getUpper(const string& s) {
	string upper(s);
	for (int i = 0; i < upper.length(); ++i) {
		upper[i] = toupper(upper[i]);
	}
	return upper;
}

list<std::vector<int> > strings2int (string block, bool addBit){
	list<std::vector<int> > outList;
	std::vector<int> tmpVector;
	std::istringstream fB(block);
	std::string tmpLine;
	int number;
	while (std::getline(fB, tmpLine)) {
		tmpVector.clear();
		if (addBit)
			tmpVector.push_back(0);
		std::stringstream stream(tmpLine);
		while (stream >> number) {
			tmpVector.push_back(number);
		}
		outList.push_back(tmpVector);
	}
	return outList;
}
list<std::vector<size_t> > strings2size_t (string block, bool addBit){
	list<std::vector<size_t> > outList;
	std::vector<size_t> tmpVector;
	std::istringstream fB(block);
	std::string tmpLine;
	int number;
	while (std::getline(fB, tmpLine)) {
		tmpVector.clear();
		if (addBit)
			tmpVector.push_back(0);
		std::stringstream stream(tmpLine);
		while (stream >> number) {
			tmpVector.push_back(number);
		}
		outList.push_back(tmpVector);
	}
	return outList;
}

