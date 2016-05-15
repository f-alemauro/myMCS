

#ifndef _UTIL_H
#define _UTIL_H

#include <string>
#include <list>
#include <vector>
#include <sstream>
std::string getUpper(const std::string& s);


std::list<std::vector<size_t> > strings2int(std::string block, bool addBit);
#endif // _UTIL_H
