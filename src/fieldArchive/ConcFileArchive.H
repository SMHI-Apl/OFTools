#ifndef CONCFILEARCHIVE_H_
#define CONCFILEARCHIVE_H_

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <iostream>
#include <iomanip>
#include "ScalarFoamFile.H"

class ConcFileArchive
{
public:
	int nWdirs,nWspeeds,nConcNames;
	std::vector<ScalarFoamFile> fileVec;
	std::vector<double> wdirs,wspeeds,concs;
	std::vector<std::string> concNames;
	std::vector<std::string> fileNames;
	
	ConcFileArchive(std::vector<std::string> inFileNames);
	~ConcFileArchive();
		
	int nextCellValue();
	//double getConc(std::string concName, double wspeed, double wdir,double sigma);

	std::vector<double> getConcMetWeights(double wspeed,double wdir,double sigma);
	std::vector<double> getConcNameWeights(std::string concName);
	std::vector<double> getWspeedWeights(double wspeed);
	std::vector<double> getWdirWeights(double wdir,double sigma);
	
	int getConcNameInd(std::string concName);
	int getWspeedInd(double wspeed);
	int getWdirInd(double wdir);
	int getConcInd(std::string concName, double wspeed, double wdir);
};


int getLowInd(std::vector<double> vec,double val);
int getHighInd(std::vector<double> vec,double val);
double sum(std::vector<double> v1);
std::vector<double> multi(std::vector<double> v1, std::vector<double> v2);
std::vector<double> multi(std::vector<double> v1, double val);
std::vector<double> add(std::vector<double> v1, std::vector<double> v2);
std::vector<double> sort(std::vector<double> v);
std::vector<double> unique(std::vector<double> v);
std::vector<std::string> unique(std::vector<std::string> v);
double P(double z);


#endif /*CONCFILEARCHIVE_H_*/
