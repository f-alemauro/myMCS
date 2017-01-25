

#ifndef _UTIL_H
#define _UTIL_H

#include <string>
#include <list>
#include <vector>
#include <sstream>
std::string getUpper(const std::string& s);
std::list<std::vector<int> > strings2int(std::string block, bool addBit);
std::list<std::vector<size_t> > strings2size_t(std::string block, bool addBit);
int factorial(int n);
#endif // _UTIL_H
