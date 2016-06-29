#include "ConcFileArchive.H"

ConcFileArchive::ConcFileArchive(std::vector<std::string> inFileNames)
{
	std::vector<std::string>::iterator fn,pp;
	std::vector<std::string> allConcNames;
	std::vector<double> allWdirs, allWspeeds;
	std::string pathRoot, dirNamePart, fileNamePart;
	for(fn=inFileNames.begin();fn!=inFileNames.end();fn++)
	{
	  pathRoot="";
	  std::vector<std::string> pathParts,fileNameParts;
	  std::string separator("/");
	  pathParts=split(*fn,separator);
	  dirNamePart=pathParts[pathParts.size()-2];
	  fileNamePart=pathParts[pathParts.size()-1];
	  fileNameParts=split(fileNamePart,".");
	 	  
	  for( pp=pathParts.begin();pp!=pathParts.end()-2;pp++)
	    pathRoot+="/"+*pp;    
  
	  for(int i=0;i<int(dirNamePart.size());i++)
	    {
	      std::string oldSeparator="_";
	      std::string newSeparator=" ";
	      if(dirNamePart.substr(i,1)==oldSeparator)
		dirNamePart.replace(i,1,newSeparator);
	    }
	  
	  std::stringstream nameStream(dirNamePart);
	  double wdir,wspeed;
	  std::string concName,dummy;

	  nameStream>>dummy;//gets rid of the wspeed prefix
	  nameStream>>wspeed;
	  nameStream>>dummy;//gets rid of the wdir prefix
	  nameStream>>wdir;
	  
	  allWdirs.push_back(wdir);
	  allWspeeds.push_back(wspeed);
	  allConcNames.push_back(fileNamePart);
	}
		
	wdirs=unique(allWdirs);
	wdirs=sort(wdirs);
	wspeeds=unique(allWspeeds);
	wspeeds=sort(wspeeds);
	concNames=unique(allConcNames);
	nWdirs=int(wdirs.size());
	nWspeeds=int(wspeeds.size());
	nConcNames=int(concNames.size());
	

	for(int concInd=0;concInd<nConcNames;concInd++)
	{
		for(int speedInd=0;speedInd<nWspeeds;speedInd++)
		{
			for(int dirInd=0;dirInd<nWdirs;dirInd++)
			{
				std::stringstream fnStream;
				fnStream<<pathRoot<<"/"<<"wspeed_";
				fnStream<<std::setiosflags(std::ios::fixed)<<std::setprecision(1)<<wspeeds[speedInd];
				fnStream<<"_wdir_";
				fnStream<<std::setiosflags(std::ios::fixed)<<std::setprecision(1)<<wdirs[dirInd];
				fnStream<<"/"<<concNames[concInd];
				std::string fileName;
				fnStream >> fileName;
				ScalarFoamFile foamFile(fileName);
				fileVec.push_back(foamFile);
				concs.push_back(-9999.0);
				fileNames.push_back(fileName);
			}
		}
	}
}

ConcFileArchive::~ConcFileArchive()
{
}


int ConcFileArchive::nextCellValue()
{
	double cellConc=-9998.0;
	std::vector<ScalarFoamFile>::iterator fileIter;
	int ind=0;
	for(fileIter=this->fileVec.begin();fileIter!=this->fileVec.end();fileIter++)
	{
		if((*fileIter).getNextValue(cellConc))
		{
			this->concs[ind]=cellConc;
			ind++;
		}
		else
			return 0;
	}		
	return 1;
}

// double ConcFileArchive::getConc(std::string concName, double wspeed, double wdir, double sigma)
// {
	// std::vector<double> metWeights=this->getConcWeights(concName,wspeed,wdir, sigma);
	// std::vector<double> concVec;
	// concVec = multi(metWeights,this->concs);
	// double sum=0.0;
	// std::vector<double>::iterator i;
	
	// for(i=concVec.begin();i!=concVec.end();i++)
		// sum+=*i;
	// return sum;
// }

std::vector<double> ConcFileArchive::getConcMetWeights(double wspeed,double wdir,double sigma)
{
	std::vector<double> weights(this->concs.size(),1.0);
	std::vector<double> concNameWeights,wspeedWeights,wdirWeights,concWeights;
	
	wspeedWeights=this->getWspeedWeights(wspeed);
	wdirWeights=this->getWdirWeights(wdir,sigma);
	
	concWeights=multi(wspeedWeights,wdirWeights);
	
	return concWeights;
}


	
std::vector<double> ConcFileArchive::getConcNameWeights(std::string concName)
{
	std::vector<double> weights(this->concs.size(),0.0);
	std::vector<double>::iterator wdirIter,wspeedIter;
	
	for(wspeedIter=this->wspeeds.begin();
    wspeedIter!=this->wspeeds.end();
    wspeedIter++)
    {
		for(wdirIter=this->wdirs.begin();
		wdirIter!=this->wdirs.end();
		wdirIter++)
		{
			int concInd = this->getConcInd(concName,  *wspeedIter, *wdirIter);
			weights[concInd]=1.0;			
		}
	}
	return weights;		
}
	
std::vector<double> ConcFileArchive::getWspeedWeights(double wspeed)
{
	std::vector<double> w(this->concs.size(),0.0);
	int ind1= getLowInd(this->wspeeds,wspeed);
	int ind2= getHighInd(this->wspeeds,wspeed);	
	std::vector<std::string>::iterator cIter;
	std::vector<double>::iterator wdIter;	
	
	for(cIter=this->concNames.begin();
		cIter!=this->concNames.end();
		cIter++)
		{
			for(wdIter=this->wdirs.begin();
				wdIter!=this->wdirs.end();
				wdIter++)
				{
					if(ind1==ind2)	
					{
						int concInd = this->getConcInd( *cIter,  this->wspeeds[ind1], *wdIter);
						w[concInd]=1.0;
					}
					else
					{
						int concInd1 = this->getConcInd( *cIter,  this->wspeeds[ind1], *wdIter);
						int concInd2 = this->getConcInd( *cIter,  this->wspeeds[ind2], *wdIter);
						w[concInd2]=(wspeed - this->wspeeds[ind1])/(this->wspeeds[ind2]-this->wspeeds[ind1]);
						w[concInd1]= 1.0-w[concInd2];
					}
				}
		}
	return w;
}

std::vector<double> ConcFileArchive::getWdirWeights(double wdir,double sigma)
{
	//First calculates normal distributed weights for all -90<=wdirs<=90
	//The weights are normalized to sum up to 1.0
	//and then sorted into a std::vector with the same length as concs
	std::vector<double> weights;
	std::vector<double> concWeights(this->concs.size(),0.0);
	double relw,relw1,relw2,mw1,mw2;

	for(int i=0;i<int(this->wdirs.size());i++)
	{
	    relw=this->wdirs[i]-wdir;
	    if(relw<-180)
	    	relw=relw+360;
	    if(relw>180)
	    	relw=relw-360; 
	    if(relw <=90 && relw >=-90)
	    {
	        if( i==0)
	            relw1=this->wdirs[this->wdirs.size()-1]-wdir;
	        else
	            relw1=this->wdirs[i-1]-wdir;
	            
            if(relw1 < -180.0)
                relw1=relw1+360.0;
            if(relw1 > 180.0)
                relw1=relw1-360.0;
                
	        if(i==int(this->wdirs.size())-1)
	            relw2=this->wdirs[0]-wdir;
	        else
	            relw2=this->wdirs[i+1]-wdir;
	            
			if(relw2<-180.0)
				relw2=relw2+360;
			if(relw2>180.0)
				relw2=relw2-360;        
	        
	        mw1=relw1+(relw-relw1)/2.0;
	        mw2=relw+(relw2-relw)/2.0;
	        
	        weights.push_back(P(mw2/sigma)-P(mw1/sigma));
	    }
	    else
	        weights.push_back(0.0);
	}

    double sum=0.0;
    std::vector<double>::iterator wi;
    for(wi=weights.begin();wi!=weights.end();wi++)
    	sum+=*wi;
    for(wi=weights.begin();wi!=weights.end();wi++)
    	(*wi)=(*wi)/sum;
    
    std::vector<double>::iterator wspeedIter,wdirIter;
    std::vector<std::string>::iterator concNameIter;
    
    for(concNameIter=this->concNames.begin();
    concNameIter!=this->concNames.end();
    concNameIter++)
    {
		for(wspeedIter=this->wspeeds.begin();
		wspeedIter!=this->wspeeds.end();
		wspeedIter++)
		{
			for(int i =0;i<int(this->wdirs.size());i++)
			{
				int concInd = this->getConcInd( *concNameIter,  *wspeedIter, this->wdirs[i]);
				concWeights[concInd]=weights[i];
			}
		}
	}
	return concWeights;		
}

int ConcFileArchive::getConcInd(std::string concName, double wspeed, double wdir)
{
	std::vector<double>::iterator wdirIter,wspeedIter;
	std::vector<std::string>::iterator concNameIter;
	
	int wdirInd=this->getWdirInd(wdir);
	int wspeedInd=this->getWspeedInd(wspeed);
	int concNameInd=this->getConcNameInd(concName);
		
	int vecInd;
	vecInd=concNameInd*(nWdirs*nWspeeds)+wspeedInd*nWdirs +wdirInd;
	return vecInd;
}

int ConcFileArchive::getWspeedInd(double wspeed)
{
	std::vector<double>::iterator i;
	int ind=0;
	for(i=this->wspeeds.begin();i!=this->wspeeds.end();i++)
	{
		if(*i==wspeed)
			return ind;
		else
			ind++;
	}
	return -1;
}

int ConcFileArchive::getWdirInd(double wdir)
{
	std::vector<double>::iterator i;
	int ind=0;
	for(i=this->wdirs.begin();i!=this->wdirs.end();i++)
	{
		if(*i==wdir)
			return ind;
		else
			ind++;
	}
	return -1;
}

int ConcFileArchive::getConcNameInd(std::string concName)
{
	std::vector<std::string>::iterator i;
	int ind=0;
	for(i=this->concNames.begin();i!=this->concNames.end();i++)
	{
		if(*i==concName)
			return ind;
		else
			ind++;
	}
	return -1;
}

int getHighInd(std::vector<double> vec, double val)
{
	std::vector<double>::iterator i;
	int ind=0;
	for(i=vec.begin();i!=vec.end();i++)
	{
		if(*i>=val)
			return ind;
		ind++;
	}
	return int(vec.size())-1;
}

int getLowInd(std::vector<double> vec, double val)
{
	std::vector<double>::iterator i;
	int ind=0;
	for(i=vec.begin();i!=vec.end();i++)
	{
		if(*i==val)
			return ind;
		if(*i>val && ind==0)
			return 0;
		if(*i>val)
			return ind-1;
		ind++;
	}
	return int(vec.size())-1;
}

double sum(std::vector<double> v1)
{
	std::vector<double>::iterator i;
	double summa=0;
	for(i=v1.begin();i!=v1.end();i++)
		summa+=(*i);
	return summa;
}

std::vector<double> multi(std::vector<double> v1, std::vector<double> v2)
{
	if(v1.size()!=v2.size())
		std::cout<<" Trying to multiply std::vectors of different sizes"<<"\n";
	
	int size=int(v1.size());
	
	std::vector<double> res(size,0);
	
	for(int i=0;i<size;i++)
		res[i]=v1[i]*v2[i];
	return res;
}
	
std::vector<double> add(std::vector<double> v1, std::vector<double> v2)
{
	if(v1.size()!=v2.size())
		std::cout<<" Trying to add vectors of different sizes"<<"\n";
	
	int size=int(v1.size());
	
	std::vector<double> res(size,0);
	
	for(int i=0;i<size;i++)
		res[i]=v1[i]+v2[i];
	
	return res;
}

std::vector<double> multi(std::vector<double> v1, double val)
{
	int size=int(v1.size());	
	std::vector<double> res(size,0);
	for(int i=0;i<size;i++)
	{
		res[i]=v1[i] * val;
	}
	return res;
}

std::vector<double> unique(std::vector<double> v)
{
	std::vector<double>::iterator i1, i2;
	std::vector<double> uniqueVec;
	int foundBefore;
	for(i1=v.begin();i1!=v.end();i1++)
	{
		foundBefore=0;
		for(i2=uniqueVec.begin();i2!=uniqueVec.end();i2++)
		{
			if(*i1==*i2)
			{
				foundBefore=1;
			}
		}
		if(not foundBefore)
		{
			uniqueVec.push_back(*i1);
		}		
	}
	return uniqueVec;
}

std::vector<std::string> unique(std::vector<std::string> v)
{
	std::vector<std::string>::iterator i1, i2;
	std::vector<std::string> uniqueVec;
	
	for(i1=v.begin();i1!=v.end();i1++)
	{
		int foundBefore=0;
		for(i2=uniqueVec.begin();i2!=uniqueVec.end();i2++)
		{
			if(*i1==*i2)
				foundBefore=1;
		}
		if(not foundBefore)
			uniqueVec.push_back(*i1);
				
	}
	return uniqueVec;
}

std::vector<double> sort(std::vector<double> v)
{
	std::vector<double>::iterator i1, i2;
	std::vector<double> sortedVec;
	int i=0;
	for(i1=v.begin();i1!=v.end();i1++)
	{
		if(i==0)
		{
			sortedVec.push_back(*i1);
			i++;
		}
		else
		{
			for(i2=sortedVec.begin();i2!=sortedVec.end();i2++)
			{				
				if((*i1)==(*i2) or (*i1)<(*i2))
				{
					sortedVec.insert(i2,*i1);
					break;
				}
			}
			if(i2==sortedVec.end())
				sortedVec.push_back(*i1);
		}
	}		
	return sortedVec;
}

double P(double z)
{
    if(z>6.0)
    {
        return 1.0;
    }
    if(z< -6.0)
    {
        return 0.0;
    }
    double b1=0.31938153;
    double b2=-0.356563782;
    double b3=1.781477937;
    double b4=-1.821255978;
    double b5=1.330274429;
    double p=0.2316419;
    double c2=0.3989423;
    
    double a=fabs(z);
    double t=1.0/(1.0+a*p);
    double b=c2*exp((-z)*(z/2.0));
    double n=((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
    n=1-b*n;
    if(z<0.0)
    {
        n=1.0-n;
    }
    return n;
}

