#ifndef TSROW_H_
#define TSROW_H_

#include <iostream>
#include <vector>

class tsRow
{
	friend std::ostream &operator<<(std::ostream &, tsRow &);
	friend int operator>>(std::istream &, tsRow &);
public:
	int year;
    int month;
    int day;
    int hour;
    int ncols;
    int currentCol;
    std::vector<double> data;
    
    tsRow();
    tsRow(int ncolsIn);
    tsRow(const tsRow &);
    tsRow(int inyear, int inmonth, int inday, int inhour, std::vector<double> indata);
    ~tsRow();
    tsRow &operator=(const tsRow &rhs);
    int operator==(const tsRow &rhs) const;
    int operator<(const tsRow &rhs) const;
    int setCurrentCol(int ind);
};

#endif /*TSROW_H_*/
