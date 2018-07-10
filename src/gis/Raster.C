#include "Raster.H"
Raster::Raster()
{
 nrows=0;
 ncols=0;
 xll=0.0;
 yll=0.0;
 xur=0.0;
 yur=0.0;
 nodata=-9999;
 cellsize=1000;
}
Raster::Raster(int NROWS, int NCOLS, double XLL, double YLL,double CELLSIZE ,int NODATA,double VALUE)
{
 nrows=NROWS;
 ncols=NCOLS;
 xll=XLL;
 yll=YLL;
 xur=XLL+(NCOLS-1)*CELLSIZE;
 yur=YLL+(NROWS-1)*CELLSIZE;
 nodata=NODATA;
 cellsize=CELLSIZE;
 std::vector<double> rasterRow(NCOLS,VALUE);
 std::vector< std::vector<double> > newRaster(NROWS, rasterRow);
 data=newRaster;
}


Raster::Raster(const Raster &rast, double value)
{
	ncols=rast.ncols;
	nrows=rast.nrows;
	nodata=rast.nodata;
	cellsize=rast.cellsize;
	xll=rast.xll;
	yll=rast.yll;
	yur=rast.yur;
	xur=rast.xur;
	std::vector<double> rasterRow(ncols,value);
	std::vector< std::vector<double> > newRaster(nrows, rasterRow);
	data=newRaster;
}

Raster::Raster(const Raster &rast)
{
	ncols=rast.ncols;
	nrows=rast.nrows;
	nodata=rast.nodata;
	cellsize=rast.cellsize;
	xll=rast.xll;
	yll=rast.yll;
	yur=rast.yur;
	xur=rast.xur;
	std::vector<double> rasterRow(ncols,nodata);
	std::vector< std::vector<double> > newRaster(nrows, rasterRow);
	data=newRaster;
}

Raster::~Raster()
{
}

int Raster::write(const char fileName[])
{
  std::ofstream fout;
  fout.open(fileName);
  fout << "nrows         "<< nrows <<"\n";
  fout << "ncols         "<< ncols << "\n";
  fout << "xllcorner     "<< xll << "\n";
  fout << "yllcorner     "<< yll << "\n";
  fout << "cellsize      "<< cellsize << "\n";
  fout << "NODATA_value  -9999"<<"\n";
  
  std::vector< std::vector<double> >::iterator row;
  std::vector<double>::iterator col;
  
  for(row=data.begin();row!=data.end();row++)
  {
  	for(col=row->begin();col!=row->end();col++)
  	{
  		if(col!=row->begin())
  			fout<<" ";
  		fout<<*col;
  	}
  	fout<<"\n";
  }
  fout.close();
  return 1;
}

int Raster::read(const char fileName[])
 {
   try
     {
       std::string tmp;
       std::ifstream fin(fileName);

       fin>>tmp; fin >> ncols;
       fin>>tmp; fin >> nrows;
       fin>>tmp; fin >> xll;
       fin>>tmp; fin >> yll;
       fin>>tmp; fin >> cellsize;       
       fin>>tmp; fin >> nodata;
       
       if(ncols==0 or nrows==0)
	 return 0;

       xur=xll+(ncols-1)*cellsize;
       yur=yll+(nrows-1)*cellsize;     

       double value;
       for (int row=0;row<nrows;row++)
	 {
	   std::vector<double> rowVector;
	   for(int col=0;col<ncols;col++)
	     {
	       fin >> value;
	       rowVector.push_back(value); 
	     }
	   data.push_back(rowVector);
	 }
       return 1;
     }
   catch(...)
     {
       return 0;
     }
 }
 
Raster& Raster::operator=(const Raster &rhs)
{
  nrows=rhs.nrows;
  ncols=rhs.ncols;
  xll=rhs.xll;
  yll=rhs.yll;
  nodata=rhs.nodata;
  cellsize=rhs.cellsize;
  yur=rhs.yur;
  xur=rhs.xur;
  std::vector<double> rasterRow(ncols,nodata);
  std::vector< std::vector<double> > newRaster(nrows, rasterRow);
  data=newRaster;
  for(int row=0;row<this->nrows;row++)
    for(int col=0;col<this->ncols;col++)
      data[row][col]=rhs.data[row][col];
  return *this;
}

int Raster::getIndex(double x, double y,int& row, int& col)
{

  if((x<xll+cellsize*0.5 or y<yll+cellsize*0.5) or (x > xur+cellsize*0.5 or y > yur+cellsize*0.5))
    return 0;   
  else
    {

      row=nrows-int((y-(yll+0.5*cellsize))/cellsize)-1;
      col=int((x-(xll+0.5*cellsize))/cellsize);
   
        if(row == nrows)
            row--;
        if(col == ncols)
            col--;

        return 1;
    }
}

bool Raster::inside(double x, double y)
{
  if((x<xll or y<yll) or (x > xur or y > yur))
    return false;
  else
    return true;
}  

double Raster::getValue(double x, double y)
{
	int row,col;
	if(!getIndex(x,y,row,col))
	  {
              std::cout<<"Value in "<<x<<","<<y<<" is outside of raster boundaries";
	    exit(1);
	  }
	else
	  return data[row][col];
	
}

double Raster::getX(int col)
{
  return xll+col*cellsize+0.5*cellsize;
}

double Raster::getY(int row)
{
  return yll+(nrows-(row+1))*cellsize+0.5*cellsize;
} 


double Raster::max()
{
  double maxValue=-9e99;
  for (int row=0;row<nrows;row++)
    {
      for(int col=0;col<ncols;col++)
	{
	  if(data[row][col]>maxValue and data[row][col]!=nodata)
	    maxValue=data[row][col];
	}
    }
  if(maxValue==-9e99)
    return nodata;
  else
    return maxValue;
}

double Raster::min()
{
  double minValue=9e99;
  for (int row=0;row<nrows;row++)
    {
      for(int col=0;col<ncols;col++)
	{
	  if(data[row][col]<minValue and data[row][col]!=nodata)
	    minValue=data[row][col];
	}
    }
  if(minValue==9e99)
    return nodata;
  else
    return minValue;
}

double Raster::sum()
{
  double sumValue=0;
  for (int row=0;row<nrows;row++)
    {
      for(int col=0;col<ncols;col++)
	{
	  if(data[row][col]!=nodata)
	    sumValue+=data[row][col];
	}
    }
  return sumValue;
}

void Raster::scale(double value)
{
  for(int row=0;row<nrows;row++)
    for(int col=0;col<ncols;col++)
      {
	if(data[row][col]!=nodata)
	data[row][col]=value*data[row][col];
      } 
}


double Raster::interpValue(double x, double y)
{
  int row,col;
  std::vector<double> vals,dists;
  double d=0;
  
  if(!getIndex(x,y,row,col))
    {
      std::cout<<"interpValue: Value is outside of raster boundaries";
      std::cout<<"\n: x="<<x;
      std::cout<<"\n: y="<<y;
      exit(1);
    }

  
  int kvadrant=0;

  if(x==this->getX(col) and y==this->getY(row))
    kvadrant=0;
  else if(x>this->getX(col) and y>this->getY(row))
    kvadrant=1;
  else if(x<this->getX(col) and y>this->getY(row))
    kvadrant=2;
  else if(x<this->getX(col) and y<this->getY(row))
    kvadrant=3;
  else if(x>this->getX(col) and y<this->getY(row))
    kvadrant=4;
  else if(x==this->getX(col) and y>this->getY(row))
    kvadrant=5;
  else if(x==this->getX(col) and y<this->getY(row))
    kvadrant=6;
  else if(x>this->getX(col) and y==this->getY(row))
    kvadrant=7;
  else if(x<this->getX(col) and y==this->getY(row))
    kvadrant=8;

  
  int row1,row2,row3,row4,col1,col2,col3,col4;

  // def of kvadrant# changed to mathematical std def
  //  2|1
  //  -+-
  //  3|4

  switch (kvadrant)
    { 
    case 0:

      row1=row;
      col1=col;
      row2=row;
      col2=col;
      row3=row;
      col3=col;
      row4=row;
      col4=col;

      break;

    case 1:
      row1=row-1;
      col1=col;
      row2=row;
      col2=col;
      row3=row;
      col3=col+1;
      row4=row-1;
      col4=col+1;

      if (row1<0 or col3>=this->ncols)
	{
	  std::cout<<"1: Value is outside of raster boundaries";
	  std::cout<<": row1="<<row1;
	  std::cout<<": col3="<<col3;
	  std::cout<<": x="<<x;
	  std::cout<<": y="<<y;
	  std::cout<<"\n";
	  exit(1);
	}
      break;
      
    case 2:
      row1=row-1;
      col1=col-1;
      row2=row;
      col2=col-1;
      row3=row;
      col3=col;
      row4=row-1;
      col4=col;

      if (row1<0 or col1<0)
	{
	  std::cout<<"2: Value is outside of raster boundaries";
	  std::cout<<": row1="<<row1;
	  std::cout<<": col1="<<col1;
	  std::cout<<": x="<<x;
	  std::cout<<": y="<<y;
	  std::cout<<"\n";
	  exit(1);
	}

      break;

    case 3:
      row1=row;
      col1=col-1;
      row2=row+1;
      col2=col-1;
      row3=row+1;
      col3=col;
      row4=row;
      col4=col;
      if (row2>=this->nrows or col2<0)
	{
	  std::cout<<"3: Value is outside of raster boundaries";
	  std::cout<<": row2="<<row2;
	  std::cout<<": col2="<<col2;
	  std::cout<<": x="<<x;
	  std::cout<<": y="<<y;
	  std::cout<<"\n";
	  exit(1);
	}

      break;
      
    case 4:
      row1=row;
      col1=col;
      row2=row+1;
      col2=col;
      row3=row+1;
      col3=col+1;
      row4=row;
      col4=col+1;
      if (row2>=this->nrows or col3>=this->ncols)
	{
	  std::cout<<"4: Value is outside of raster boundaries";
	  std::cout<<": row2="<<row2;
	  std::cout<<": col3="<<col3;
	  std::cout<<": x="<<x;
	  std::cout<<": y="<<y;
	  std::cout<<"\n";
	  exit(1);
	}

      break;
      
    case 5:
      row1=row-1;
      col1=col;
      row2=row;
      col2=col;
      row3=row;
      col3=col;
      row4=row-1;
      col4=col;

      if (row1<0)
	{
	  std::cout<<"5: Value is outside of raster boundaries";
	  std::cout<<": row1="<<row1;
	  std::cout<<": x="<<x;
	  std::cout<<": y="<<y;
	  std::cout<<"\n";
	  exit(1);
	}
      break;
      
    case 6:
      row1=row;
      col1=col;
      row2=row+1;
      col2=col;
      row3=row+1;
      col3=col;
      row4=row;
      col4=col;

      if (row2>=this->nrows)
	{
	  std::cout<<"6: Value is outside of raster boundaries";
	  std::cout<<": row2="<<row2;
	  std::cout<<": x="<<x;
	  std::cout<<": y="<<y;
	  std::cout<<"\n";
	  exit(1);
	}

      break;

    case 7:
      row1=row;
      col1=col;
      row2=row;
      col2=col;
      row3=row;
      col3=col+1;
      row4=row;
      col4=col+1;
      if (col3>=this->ncols)
	{
	  std::cout<<"7: Value is outside of raster boundaries";
	  std::cout<<": col3="<<col3;
	  std::cout<<": x="<<x;
	  std::cout<<": y="<<y;
	  std::cout<<"\n";
	  exit(1);
	}

      break;
      
    case 8:
      row1=row;
      col1=col-1;
      row2=row;
      col2=col-1;
      row3=row;
      col3=col;
      row4=row;
      col4=col;
      if (col1<0)
	{
	  std::cout<<"8: Value is outside of raster boundaries";
	  std::cout<<": col1="<<col1;
	  std::cout<<": x="<<x;
	  std::cout<<": y="<<y;
	  std::cout<<"\n";
	  exit(1);
	}

      break;
      
    }
	  	  
  // DEBUG OUTPUT
  /*
  std::cout<<"\n ";
  std::cout<<"\n x="<<x;
  std::cout<<"\n y="<<y;
  std::cout<<"\n row1="<<row1;
  std::cout<<"\n row2="<<row2;
  std::cout<<"\n row3="<<row3;
  std::cout<<"\n row4="<<row4;
  std::cout<<"\n col1="<<col1;
  std::cout<<"\n col2="<<col2;
  std::cout<<"\n col3="<<col3;
  std::cout<<"\n col4="<<col4;
  std::cout<<"\n x1="<<this->getX(col1);
  std::cout<<"\n x2="<<this->getX(col4);
  std::cout<<"\n y1="<<this->getY(row1);
  std::cout<<"\n y2="<<this->getY(row2);
  */

  d=sqrt(pow(this->getX(col1)-x,2)+pow(this->getY(row1)-y,2));		
  dists.push_back(d);
  vals.push_back(data[row1][col1]);
  
  d=sqrt(pow(this->getX(col2)-x,2)+pow(this->getY(row2)-y,2));		
  dists.push_back(d);
  vals.push_back(data[row2][col2]);
  
  d=sqrt(pow(this->getX(col3)-x,2)+pow(this->getY(row3)-y,2));		
  dists.push_back(d);
  vals.push_back(data[row3][col3]);
	
  d=sqrt(pow(this->getX(col4)-x,2)+pow(this->getY(row4)-y,2));		
  dists.push_back(d);
  vals.push_back(data[row4][col4]);
  
  double a=0;
  double b=0;
  for(int j=0;j<4;j++)
    {
      if (dists[j]>0)
	{
	  a=a+vals[j]/pow(dists[j],2);
	  b=b+1./pow(dists[j],2);
	}
      else
	return vals[j];
    }
  return a/b;

}

bool Raster::conformal(const Raster &rast)
{       
  if (cellsize!=rast.cellsize)
    return false;
  double xsteps=(xll-rast.xll)/double(cellsize);
  double ysteps=(yll-rast.yll)/double(cellsize);
  double xrest=xsteps-floor(xsteps);
  double yrest=ysteps-floor(ysteps);
  if(xrest!=0 or yrest!=0)
    return false;
  else
    return true;
}

bool Raster::match(const Raster &rast)
{
  if(this->conformal(rast) and xll==rast.xll 
     and yll==rast.yll and ncols==rast.ncols 
     and nrows==rast.nrows)
    return true;
  else
    return false;
}
    
void Raster::whereAdd(const Raster &rast,double value)
{
  for(int row=0;row<nrows;row++)
    {    
      for(int col=0;col<ncols;col++)
	{
	  if (rast.data[row][col]>0 and rast.data[row][col]!=rast.nodata)
	    data[row][col]=data[row][col]+value;
	  if (rast.data[row][col]==rast.nodata)
	    data[row][col]=rast.nodata;
	}
    }
}

void Raster::translate(const double &x, const double &y) {
  xll += x;
  xur += x;
  yll += y;
  yur += y;
}

Raster where(const Raster &rast,const double &val1, const double &val2)
{
  Raster temp(rast,0);
  for(int row=0;row<temp.nrows;row++)
    {    
      for(int col=0;col<temp.ncols;col++)
	{
	  if (rast.data[row][col]>0 and rast.data[row][col]!=rast.nodata)
	    temp.data[row][col]=val1;
	  else
	    temp.data[row][col]=val2;
	}
    }
  return temp;
}  

Raster where(const Raster &rast,const Raster &rast1, const Raster &rast2)
{
  Raster temp(rast,0);
  for(int row=0;row<temp.nrows;row++)
    {    
      for(int col=0;col<temp.ncols;col++)
	{
	  if (rast.data[row][col]>0 and rast.data[row][col]!=rast.nodata)
	    temp.data[row][col]= rast1.data[row][col];
	  else
	    temp.data[row][col]= rast2.data[row][col];
	}
    }
  return temp;
}  


Raster operator /(const Raster &rast1, const Raster &rast2)
{
  Raster temp(rast1);

  for(int row=0;row<temp.nrows;row++)
    for(int col=0;col<temp.ncols;col++)
      {
	if(rast1.data[row][col]!=rast1.nodata and rast2.data[row][col]!=rast2.nodata)
	  temp.data[row][col]=rast1.data[row][col]/rast2.data[row][col];
	else
	  temp.data[row][col]=temp.nodata;
      }
	return temp;
}

Raster operator *(const Raster &rast1, const Raster &rast2)
{
  Raster temp(rast1);

  for(int row=0;row<temp.nrows;row++)
    for(int col=0;col<temp.ncols;col++)
      {
	if(rast1.data[row][col]!=rast1.nodata and rast2.data[row][col]!=rast2.nodata)
	  temp.data[row][col]=rast1.data[row][col]*rast2.data[row][col];
	else
	  temp.data[row][col]=temp.nodata;
      }
  return temp;
}
Raster operator +(const Raster &rast1, const Raster &rast2)
{
  Raster temp(rast1);

  for(int row=0;row<temp.nrows;row++)
    for(int col=0;col<temp.ncols;col++)
      {
	if(rast1.data[row][col]!=rast1.nodata and rast2.data[row][col]!=rast2.nodata)
	  temp.data[row][col]=rast1.data[row][col]+rast2.data[row][col];
	else
	  temp.data[row][col]=temp.nodata;
      }
return temp;
}
Raster operator -(const Raster &rast1, const Raster &rast2)
{
  Raster temp(rast1);

  for(int row=0;row<temp.nrows;row++)
    for(int col=0;col<temp.ncols;col++)
      {
	if(rast1.data[row][col]!=rast1.nodata and rast2.data[row][col]!=rast2.nodata)
	  temp.data[row][col]=rast1.data[row][col]-rast2.data[row][col];
	else
	  temp.data[row][col]=temp.nodata;
      }
  return temp;
}

Raster operator *(const double value, const Raster &rast)
{
  Raster temp(rast);
  for(int row=0;row<temp.nrows;row++)
    for(int col=0;col<temp.ncols;col++)
      {
	if(rast.data[row][col]!=rast.nodata)
	  temp.data[row][col]=value*rast.data[row][col];
	else
	  temp.data[row][col]=temp.nodata;
      }
	return temp;
}

Raster operator *(const Raster &rast, const double value)
{
  Raster temp(rast);
  for(int row=0;row<temp.nrows;row++)
    for(int col=0;col<temp.ncols;col++)
      {
	if(rast.data[row][col]!=rast.nodata)
	  temp.data[row][col]=value*rast.data[row][col];
	else
	  temp.data[row][col]=temp.nodata;
      }
return temp;
}

Raster operator +(const Raster &rast, const double value)
{
  Raster temp(rast);
  for(int row=0;row<temp.nrows;row++)
    for(int col=0;col<temp.ncols;col++)
      { 
	if(rast.data[row][col]!=rast.nodata)
	  temp.data[row][col]=value+rast.data[row][col];
	else
	  temp.data[row][col]=temp.nodata;
      }
  return temp;
}

Raster operator -(const Raster &rast, const double value)
{
  Raster temp(rast);
  for(int row=0;row<temp.nrows;row++)
    for(int col=0;col<temp.ncols;col++)
      {
	if(rast.data[row][col]!=rast.nodata)
	  temp.data[row][col]=rast.data[row][col]-value;
	else
	  temp.data[row][col]=temp.nodata;
      }
  return temp;
}

Raster operator /(const Raster &rast, const double value)
{
  Raster temp(rast);
  for(int row=0;row<temp.nrows;row++)
    for(int col=0;col<temp.ncols;col++)
      {
      if(rast.data[row][col]!=rast.nodata)
	temp.data[row][col]=rast.data[row][col]/value;
      else
	temp.data[row][col]=temp.nodata;
      }
  return temp;
}

Raster operator ==(const Raster &rast,const double value)
{
  Raster temp(rast,0);
  for(int row=0;row<temp.nrows;row++)
    {    
      for(int col=0;col<temp.ncols;col++)
	{
	  if (rast.data[row][col]==value  and rast.data[row][col]!=rast.nodata)
	    temp.data[row][col]=1;
	}
    }
  return temp;
}
  
Raster operator <(const Raster &rast,const double value)
{
  Raster temp(rast,0);
  for(int row=0;row<temp.nrows;row++)
    {    
      for(int col=0;col<temp.ncols;col++)
	{
	  if (rast.data[row][col]<value  and rast.data[row][col]!=rast.nodata)
	    temp.data[row][col]=1;
	}
    }
  return temp;
}

Raster operator >(const Raster &rast,const double value)
{
  Raster temp(rast,0);
  for(int row=0;row<temp.nrows;row++)
    {    
      for(int col=0;col<temp.ncols;col++)
	{
	  if (rast.data[row][col]>value  and rast.data[row][col]!=rast.nodata)
	    temp.data[row][col]=1;
	}
    }
  return temp;
}

Raster operator <=(const Raster &rast,const double value)
{
  Raster temp(rast,0);
  for(int row=0;row<temp.nrows;row++)
    {    
      for(int col=0;col<temp.ncols;col++)
	{
	  if (rast.data[row][col]<=value and rast.data[row][col]!=rast.nodata)
	    temp.data[row][col]=1;
	}
    }
  return temp;
}

Raster operator >=(const Raster &rast,const double value)
{
  Raster temp(rast,0);
  for(int row=0;row<temp.nrows;row++)
    {    
      for(int col=0;col<temp.ncols;col++)
	{
	  if (rast.data[row][col]>=value and rast.data[row][col]!=rast.nodata)
	    temp.data[row][col]=1;
	}
    }
  return temp;
}

Raster intersect(const Raster &rast1,const Raster &rast2)
{
  Raster temp(rast1,0);
  for(int row=0;row<temp.nrows;row++)
    {    
      for(int col=0;col<temp.ncols;col++)
	{
	  if (rast1.data[row][col]>0 and rast2.data[row][col]>0)
	    temp.data[row][col]=1;
	}
    }
  return temp;
}
