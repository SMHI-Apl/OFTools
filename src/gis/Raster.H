#ifndef RASTER_H_
#define RASTER_H_

#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <iostream>
#include <stdlib.h>



class Raster {
public:
	int nrows;
	int ncols;
	double xll;
	double yll;
	double xur;
	double yur;
	double cellsize;
	double nodata;
	std::vector< std::vector<double> > data;
			
	Raster();
	Raster(int NROWS, int NCOLS,
	       double XLL, double YLL,
	       double CELLSIZE ,int NODATA=-9999,
	       double VALUE=0.0);

        Raster(const Raster &rast,double value);
        Raster(const Raster &rast);
	Raster & operator =(const Raster &rast);
  // const Raster operator*(const double s);
 
        bool conformal(const Raster &rast); 
        bool match(const Raster &rast);
 
        double interpValue(double x, double y);
        double getValue(double x, double y);	
	int getIndex(double x, double y, int& row, int& col);
	int write(const char fileName[]);
	int read(const char fileName[]);
        
        double getX(int col);
        double getY(int row);
        
        double max();
        double min();
        double sum();
        
        Raster where(const Raster &rast,const Raster &rast1,const Raster &rast2);
        Raster where(const Raster &rast,const double &val1,const double &val2);
       
        void scale(double value);
        void whereAdd(const Raster &rast,double value);
	
	~Raster();


};

Raster operator /(const Raster &rast1, const Raster &rast2);
Raster operator *(const Raster &rast1, const Raster &rast2);
Raster operator +(const Raster &rast1, const Raster &rast2);
Raster operator -(const Raster &rast1, const Raster &rast2);
Raster operator *(const double value, const Raster &rast);
Raster operator *(const Raster &rast, const double value);
Raster operator +(const Raster &rast, const double value);
Raster operator -(const Raster &rast, const double value);
Raster operator /(const Raster &rast, const double value);
Raster operator ==(const Raster &rast,const double value);
Raster operator <(const Raster &rast,const double value);
Raster operator >(const Raster &rast,const double value);
Raster operator <=(const Raster &rast,const double value);
Raster operator >=(const Raster &rast,const double value);
Raster intersect(const Raster &rast1,const Raster &rast2);

#endif /*RASTER_H_*/
