#include "tsRow.H"

tsRow::~tsRow()
{
}

tsRow::tsRow()   // Constructor
{
   year = 0;
   month = 0;
   day=0;
   hour=0;
   ncols=0;
   currentCol=0;
}

tsRow::tsRow(int ncolsIn)   // Constructor
{
   year = 0;
   month = 0;
   day=0;
   hour=0;
   ncols=ncolsIn;
   data.resize(ncols);
   currentCol=0;
}


tsRow::tsRow(const tsRow &copyin)   // Copmonth constructor to handle pass bmonth value.
{                             
   year = copyin.year;
   month = copyin.month;
   day = copyin.day;
   hour = copyin.hour;
   ncols = copyin.ncols;
   data = copyin.data;
   currentCol=copyin.currentCol;
}

tsRow::tsRow(int inyear,int inmonth, int inday, int inhour, std::vector<double> indata)
{
	year = inyear;
   	month = inmonth;
   	day = inday;
   	hour = inhour;
   	data = indata;
   	ncols=indata.size();
   	currentCol=0;
}  	


std::ostream &operator<<(std::ostream &output, tsRow &inRow)
{
   output << inRow.year << '\t' << inRow.month << '\t' << inRow.day<<'\t'<<inRow.hour << '\t';
   std::vector<double>::iterator i;
   i=inRow.data.begin();
   
   for(i=inRow.data.begin();i!= inRow.data.end()-1; i++)
   { 
   	output << *i << "\t";
   }
   output << *(inRow.data.end()-1)<<"\n";
   
   return output;
}

int operator>>(std::istream &inStr, tsRow &outRow)
{
	int checkCount=0;
	if(inStr >> outRow.year)
		checkCount++;
	if(inStr >> outRow.month)
		checkCount++;
	if(inStr >> outRow.day)
		checkCount++;
	if(inStr >> outRow.hour)
		checkCount++;
	std::vector<double>::iterator i;
	double val;
	for(i=outRow.data.begin();i!= outRow.data.end(); i++)
	{	
		if(inStr>>val)
			checkCount++;
	 	*i = val;
   	}
   	if(checkCount==(4+outRow.ncols))
   		return 1;
   	else
   		return 0;
}

tsRow& tsRow::operator=(const tsRow &rhs)
{
   this->year = rhs.year;
   this->month = rhs.month;
   this->day = rhs.day;
   this->hour = rhs.hour;
   this->ncols=rhs.ncols;
   this->data=rhs.data;
   this->currentCol=rhs.currentCol;
   return *this;
}

int tsRow::operator==(const tsRow &rhs) const
{
   if( this->data[this->currentCol] != rhs.data[this->currentCol]) return 0;
   return 1;
}

// This function is required for built-in STL list functions like sort
int tsRow::operator<(const tsRow &rhs) const
{
   if(this->data[this->currentCol] < rhs.data[this->currentCol]) return 1;
   return 0;
}

int tsRow::setCurrentCol(int ind)
{
	if(ind<0 || ind>=this->ncols)
		return 0;
	else
	{
		currentCol=ind;
		return 1;
	}
}
