#ifndef TIMESERIES_H_
#define TIMESERIES_H_

#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include "tsRow.H"

class timeSeries
{
	friend std::ostream &operator<<(std::ostream &, timeSeries &);
public:
	int nrows;
	int ncols;
	std::list<tsRow> rows;
	std::vector<std::string> header;
			
	timeSeries();
	timeSeries(std::string fileName);
	timeSeries(const timeSeries &inTs);
	timeSeries & operator =(const timeSeries &rhs);
	
	int appendRow(tsRow inrow);
        void deleteRow(std::list<tsRow>::iterator row);
        int sort();
	int setCurrentCol(int ind);
	int getColInd(std::string colName);
        int matchDates(timeSeries &);
        void removeNodata(double nodata=-9999.0,int colInd=-1);
	double percentile(double perc);
        double condPerc(double perc, timeSeries ts, int nClasses, int dirCol, int wspeedCol, double minDir,double maxDir, double excludeVal=-9999);
        int maxIndex();
        int getYear(int index);
        int getMonth(int index);
        int getDay(int index);
        int getHour(int index);
        double getValue(int index);
	double average();
	timeSeries dailyAverage(int nValid);
        timeSeries dailyMax(int nValid);


	~timeSeries();
};

#endif /*TIMESERIES_H_*/
