#include "timeSeries.H"


timeSeries::timeSeries()
{
 nrows=0;
 ncols=0;
}
timeSeries::timeSeries(const timeSeries &inTs)
{
	rows=inTs.rows;
	header=inTs.header;
	ncols=inTs.ncols;
	nrows=inTs.nrows;
}

timeSeries::~timeSeries()
{
}


std::ostream &operator<<(std::ostream &output, timeSeries &inTs)
{
    if(int(inTs.header.size()) < inTs.ncols) {
            std::cout<<"Error: invalid header in timeseries"<<"\n";
            exit(1);
        }
    
	output<<"year\tmonth\tday\thour\t";        
	std::vector<std::string>::iterator headerIter;
	for(headerIter=inTs.header.begin();
		headerIter!= inTs.header.end()-1;
		headerIter++)
	{ 
   		output << *headerIter << "\t";
   	}
   	output << *(inTs.header.end()-1)<< "\n";
	std::list<tsRow>::iterator row;
	for(row=inTs.rows.begin();
		row!= inTs.rows.end();
		row++)
   	{ 
   		output << *row;
   	}
   	return output;
}
timeSeries::timeSeries(std::string fileName)
{
	nrows=0;
	ncols=0;
	std::ifstream istr(fileName.c_str());
	if(not istr)
	{
		std::cout<<"Error: timeseries file: "<<fileName<<" could not be opened"<<"\n";
		exit(1);
	}
		
	char headerLine[200];
	istr.getline(headerLine,200);
	std::string headerStr(headerLine);
	std::string buf;
	std::stringstream ss(headerStr);
	
	while(ss>>buf)
		if(buf!="year" && buf!="month" && buf!="day" && buf!="hour") 
			header.push_back(buf);
	ncols=header.size();
	
	char dataRow[200];
        char empty_row[] = "\n";
	istr.getline(dataRow,200);
	while(not istr.eof() && dataRow !=empty_row)
	{				
			std::string rowStr(dataRow);
			std::stringstream rowStream(rowStr);
			tsRow inRow(ncols);
			inRow.ncols=ncols;
			if(not(rowStream>>inRow))
			{
			  std::cout<<"Error in input file: "<<fileName<<"\n";
				break;
			}
			nrows++;
			rows.push_back(inRow);
			istr.getline(dataRow,200);
	}
}
timeSeries& timeSeries::operator=(const timeSeries &rhs)
{
	this->header=rhs.header;
	this->rows=rhs.rows;
	this->ncols=rhs.ncols;
	this->nrows=rhs.nrows;
	return *this;
}
	
int timeSeries::appendRow(tsRow inrow)
{
	if(inrow.ncols!=this->ncols && this->ncols!=0)
		return 0;
	this->nrows++;
	this->rows.push_back(inrow);
	this->ncols=inrow.ncols;
	return 1;
}

void timeSeries::deleteRow(std::list<tsRow>::iterator row)
{
  this->rows.erase(row);
  this->nrows-=1;
}

int timeSeries::getColInd(std::string colName)
{
	std::vector<std::string>::iterator i;
	int ind=0;
	for(i=this->header.begin();i!=this->header.end();i++)
	{
		if(*i==colName)
			return ind;
		else
			ind++;
	}
	return -1;
}

int timeSeries::setCurrentCol(int ind)
{
	std::list<tsRow>::iterator row;
	for(row=this->rows.begin();row!=this->rows.end();row++)
	{
		if(not ((*row).setCurrentCol(ind)))
			return 0;
	}
	return 1;
}

int timeSeries::sort()
{
	this->rows.sort();
	return 1;
}

double timeSeries::percentile(double perc)
{
        perc=perc/100.0;
	double nValsSmaller, rest,lowVal,highVal; 
	int nWholeValsSmaller;
	std::list<tsRow>::iterator row;
	
	nValsSmaller=perc * double(this->nrows);
	nWholeValsSmaller = int(nValsSmaller);
	rest=nValsSmaller-double(nWholeValsSmaller);

	row=this->rows.begin();
	for(int i=1;i<nWholeValsSmaller;i++)
		row++;
	lowVal=(*row).data[(*row).currentCol];
	if(rest>0)
		row++;
	highVal=(*row).data[(*row).currentCol];
	if(highVal!=lowVal)
		return (lowVal+rest*(highVal-lowVal));
	else
		return lowVal;
}





//Calculate percentile for all hours for which a condition on a secon ts is fulfilled
//The condition(min and max  values) is evaluated for the current column of the given ts
//example: Percentile of wind speed with condition on wind direction
double timeSeries::condPerc(double perc, timeSeries ts, int nClasses, int dirCol, int wspeedCol, double minDir,double maxDir, double excludeVal)
{
  timeSeries tmpTs(*this);        
  double nValsSmaller, rest,lowVal,highVal; 
  int nWholeValsSmaller;
  std::list<tsRow>::iterator tmpRow,row;

  tmpTs.matchDates(ts);
  ts.matchDates(tmpTs);
  
  int nExcluded=0;
  row = ts.rows.begin();
  
  bool deleteRow = false;
  std::list<tsRow>::iterator toBeErased;
  toBeErased = tmpTs.rows.begin();
 
  for(tmpRow=tmpTs.rows.begin();tmpRow!=tmpTs.rows.end();tmpRow++)
    {
      //if outside interval, row from last iteration is deleted 
      if(deleteRow)
	tmpTs.deleteRow(toBeErased);
      
      deleteRow=false;
      
      //Checks if the ts value is equal to the excludeVal (e.g. nodata value)
      if (row->data[dirCol]==excludeVal)
	{
	  nExcluded+=1;
	  deleteRow=true;
	}
      else 
	{
	  //Checks if value is outside of interval
	  //Fixa dir >360 eller <0
	  bool overlapsZero=false;
	  
	  if(minDir>maxDir)
	      overlapsZero=true;
    
	  if(overlapsZero)
	    {
	      if( not( (row->data[dirCol]<maxDir) || (row->data[dirCol]>=minDir)) ) 
		deleteRow=true;
	    }	      
	  else
	    {
	      if((row->data[dirCol]>=maxDir) || (row->data[dirCol]<minDir) ) 
		deleteRow=true;
	    }
	}   
	  //Stores the current row to be deleted in the next iteration
	  toBeErased=tmpRow;
	  row++;
     } 

  perc=perc/100.0;	
 
  //Devides the excluded values over the classes
  double excludedValsPerClass=nExcluded/double(nClasses);

  //Adds the excludedValsPerClass to get the percentile correct
  nValsSmaller=perc * double(tmpTs.nrows+excludedValsPerClass);

  nWholeValsSmaller = int(nValsSmaller);
  rest=nValsSmaller-double(nWholeValsSmaller);

  tmpTs.setCurrentCol(wspeedCol);
  tmpTs.sort();

  tmpRow=tmpTs.rows.begin();
  for(int i=1;i<nWholeValsSmaller;i++)
    tmpRow++;
  lowVal=tmpRow->data[wspeedCol];

  if(rest>0)
    tmpRow++;
  highVal=tmpRow->data[wspeedCol];
  if(highVal!=lowVal)
    return (lowVal+rest*(highVal-lowVal));
  else
    return lowVal;
}

int timeSeries::maxIndex()
{
	std::list<tsRow>::iterator row;
	double maxVal=-9e99;
	int ind=0,maxInd=-9999;
	
	for(row=this->rows.begin();row!=this->rows.end();row++)
	{
	  if(row->data[row->currentCol]>maxVal)
	    maxInd=ind;
	  ind++;
	}
	  
	return maxInd;
}


double timeSeries::average()
{
	double sum=0.0,n=0.0; 
	std::list<tsRow>::iterator row;
	for(row=this->rows.begin();row!=this->rows.end();row++)
	{
		sum+=(*row).data[(*row).currentCol];
		n+=1.0;
	}
	return sum/n;
}

int timeSeries::getYear(int index)
{
  std::list<tsRow>::iterator row;
  int ind=0;  
  for(row=this->rows.begin();row!=this->rows.end();row++)
    {
      if(ind<index)
	ind++;
      else
	return row->year;
    }
  return -9999; 
}


int timeSeries::getMonth(int index)
{
  std::list<tsRow>::iterator row;
  int ind=0;  
  for(row=this->rows.begin();row!=this->rows.end();row++)
    {
      if(ind<index)
	ind++;
      else
	return row->month;
    }
  return -9999; 
}

int timeSeries::getDay(int index)
{
  std::list<tsRow>::iterator row;
  int ind=0;  
  for(row=this->rows.begin();row!=this->rows.end();row++)
    {
      if(ind<index)
	ind++;
      else
	return row->day;
    }
  return -9999; 
}

int timeSeries::getHour(int index)
{
  std::list<tsRow>::iterator row;
  int ind=0;  
  for(row=this->rows.begin();row!=this->rows.end();row++)
    {
      if(ind<index)
	ind++;
      else
	return row->hour;
    }
  return -9999; 
}



double timeSeries::getValue(int index)
{
  std::list<tsRow>::iterator row;
  int ind=0;  
  for(row=this->rows.begin();row!=this->rows.end();row++)
    {
      if(ind<index)
	ind++;
      else
	return row->data[row->currentCol];
    }
  return -9999; 
}

//Function presumes that time series values are give in time order
timeSeries timeSeries::dailyAverage(int nValid)
{
  int y1,m1,d1,y2,m2,d2;
  int nvals=0;
  std::list<tsRow>::iterator row;
  timeSeries dailyTs;
  
  row=this->rows.begin();
  std::vector<double> meanVec((*row).ncols,0.0);
  
  y1=(*row).year;
  m1=(*row).month;
  d1=(*row).day;
  
  for(row=this->rows.begin();row!=this->rows.end();row++)
    {
      if(row==this->rows.begin())
	{
	  for(int col=0;col<int(meanVec.size());col++)
	    meanVec[col]+=(*row).data[col];
	  nvals++;
	}
      else
	{
	  y2=(*row).year;
	  m2=(*row).month;
	  d2=(*row).day;
	  if(y2==y1 && m2==m1 && d2==d1)
	    {
	      nvals++;
	      for(int col=0;col<int(meanVec.size());col++)
		meanVec[col]+=(*row).data[col];
	    }
	  else
	    {
	      if(nvals>=nValid)
		{
		  for(int col=0;col<int(meanVec.size());col++)
		    meanVec[col]/=nvals;
		  tsRow dailyRow(y1,m1,d1,-1,meanVec);
		  dailyTs.appendRow(dailyRow);
		  
		  nvals=1;
		  meanVec=std::vector<double>(meanVec.size(),0.0);
		  y1=y2;
		  m1=m2;
		  d1=d2;
		}
	      else
		{
		  nvals=1;
		  meanVec=std::vector<double>(meanVec.size(),0.0);
		  y1=y2;
		  m1=m2;
		  d1=d2;
		}
	    }
	}
    }
  return dailyTs;
}

//Function presumes that time series values are give in time order
timeSeries timeSeries::dailyMax(int nValid)
{
  int y1,m1,d1,y2,m2,d2;
  int nvals=0;
  std::list<tsRow>::iterator row;
  timeSeries dailyTs;
  
  row=this->rows.begin();
  std::vector<double> maxVec((*row).ncols,-999999999.0);
  
  y1=(*row).year;
  m1=(*row).month;
  d1=(*row).day;
  
  for(row=this->rows.begin();row!=this->rows.end();row++)
    {
      if(row==this->rows.begin())
	{
	  for(int col=0;col<int(maxVec.size());col++)
	    {
	      if(maxVec[col]<(*row).data[col])
		maxVec[col]=(*row).data[col];
	    }
	  nvals++;
	}
      else
	{
	  y2=(*row).year;
	  m2=(*row).month;
	  d2=(*row).day;
	  if(y2==y1 && m2==m1 && d2==d1)
	    {
	      nvals++;
	      std::vector<double>::iterator col;
	      for(int col=0;col<int(maxVec.size());col++)
		{
		  if(maxVec[col]<(*row).data[col])
		    maxVec[col]=(*row).data[col];
		}
	    }
	  else
	    {
	      if(nvals>=nValid)
		{
		  tsRow dailyRow(y1,m1,d1,-1,maxVec);
		  dailyTs.appendRow(dailyRow);
		  
		  nvals=1;
		  maxVec=std::vector<double>(maxVec.size(),-999999999.0);
		  y1=y2;
		  m1=m2;
		  d1=d2;
		}
	      else
		{
		  nvals=1;
		  maxVec=std::vector<double>(maxVec.size(),0.0);
		  y1=y2;
		  m1=m2;
		  d1=d2;
		}
	    }
	}
}
return dailyTs;
}

int timeSeries::matchDates(timeSeries &ts)
{
  int y1,m1,d1,h1,y2,m2,d2,h2;
  std::list<tsRow>::iterator row1;
  std::list<tsRow>::iterator row2;
  std::list<tsRow>::iterator toBeErased;
  
  int removedDates=0;
  bool rowToRemove=0;
  bool dateFound=0;
  for(row1=this->rows.begin();row1!=this->rows.end();row1++)
    {
      int removedDates=0;
  
      if(rowToRemove)
	{
	  this->rows.erase(toBeErased);
	  this->nrows-=1;
	  rowToRemove=0;
	}

      dateFound=0;


      y1=(*row1).year;
      m1=(*row1).month;
      d1=(*row1).day;
      h1=(*row1).hour;      
      
      row2= ts.rows.begin();
      while(!dateFound && row2!=ts.rows.end())
	{
	  y2=(*row2).year;
	  m2=(*row2).month;
	  d2=(*row2).day;
	  h2=(*row2).hour;
	  if(y2==y1 && m2==m1 && d2==d1 && h2==h1)
	      dateFound=1;
	  row2++;
	}
      
      if(!dateFound)
	{
	
	  removedDates++;
	  rowToRemove=1;
	  toBeErased=row1;

	}
    }
  
  //If the last row is missing it is handled here
  if(rowToRemove)
    {
      this->rows.erase(toBeErased);
      this->nrows-=1;
    }
      
  return removedDates;
}


void timeSeries::removeNodata(double nodata, int colInd)
  {
    if (colInd==-1)
      {
	for(int col=0;col<this->ncols;col++)
	  {
	    this->setCurrentCol(col);
	    tsRow nodataRow(this->ncols);
	    nodataRow.setCurrentCol(col);
	    nodataRow.data[col]=nodata;
	    //remove all rows which are equal to nodataRow
	    //equality is tested only on currentCol!
	    this->rows.remove(nodataRow);
	  }
      }
    else
      {
	this->setCurrentCol(colInd);
	tsRow nodataRow(this->ncols);
	nodataRow.setCurrentCol(colInd);
	nodataRow.data[colInd]=nodata;
	//remove all rows which are equal to nodataRow
	//equality is tested only on currentCol!
	this->rows.remove(nodataRow);
      }

    this->nrows=int(this->rows.size());
    return;
  }
