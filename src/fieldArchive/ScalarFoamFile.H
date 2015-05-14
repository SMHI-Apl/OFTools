#ifndef SCALARFOAMFILE_H_
#define SCALARFOAMFILE_H_

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>

class ScalarFoamFile
{
public:
	std::ifstream ifStr;
	std::string filePath;
	std::string fileName;
	//streampos filePos;
	ScalarFoamFile(std::string InNameField);
	ScalarFoamFile(const ScalarFoamFile &sff);
	int getNextValue(double &outVal);
	int init();
	~ScalarFoamFile();
	ScalarFoamFile &operator=(const ScalarFoamFile &sff);
        std::string getDirName();
};
std::vector<std::string> split(std::string str,std::string sep);


#endif /*SCALARFOAMFILE_H_*/
