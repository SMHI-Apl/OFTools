#include "ScalarFoamFile.H"

ScalarFoamFile::ScalarFoamFile(std::string infileName)
{
	this->filePath=infileName;
	std::vector<std::string> pathParts= split(infileName,"/");
	this->fileName=pathParts[int(pathParts.size())-1];
	if(not this->init())
		exit(1);
	//filePos=ifStr.tellg();
	//ifStr.close();	
}

ScalarFoamFile::ScalarFoamFile(const ScalarFoamFile &sff)
{
	filePath=sff.filePath;
	fileName=sff.fileName;
	if(not this->init())
		exit(1);
	//filePos=sff.filePos;	
}

int ScalarFoamFile::init()
{
	bool firstCellFound=0;
	if(not this->ifStr.is_open())
		this->ifStr.open(this->filePath.c_str());
		
	if(not this->ifStr.is_open())
	{
		std::cout<<"Error in ScalarFoamFile.cpp: Could not open file: "<<this->filePath<<"\n";
		return 0;
	}
	while(not this->ifStr.eof()&& not firstCellFound)
	{
		std::string buf;
		this->ifStr>>buf;
		if(buf=="internalField")
		{
			while(not this->ifStr.eof() && not firstCellFound)
			{
				std::string anyStr;
				this->ifStr>>anyStr;
				std::string stopStr="(";					
				if(anyStr==stopStr)
					firstCellFound=1;	
			}
		}
	}
	return 1;	
}

int ScalarFoamFile::getNextValue(double &outVal)
{
	std::string word;
	//this->ifStr.open(this->filePath.c_str());
	if(not this->ifStr.is_open())
	{
		std::cout<<"Error in ScalarFoamFile.cpp: File not opened/initialized"<<"\n";	
		exit(1);
	}

	//this->ifStr.seekg(filePos);

	this->ifStr>>word;
	//this->filePos=ifStr.tellg();
	//this->ifStr.close();
	std::stringstream lineStream(word);
	std::string stopStr=")";
	if(word!=stopStr)
	{
		lineStream>>outVal;
		return 1;
	}
	else
		return 0;
}		

ScalarFoamFile::~ScalarFoamFile()
{
	if(this->ifStr.is_open())
		this->ifStr.close();
}

ScalarFoamFile& ScalarFoamFile::operator=(const ScalarFoamFile &sff)
{
	this->filePath=sff.filePath;
	this->fileName=sff.fileName;
	if(not this->init())
		exit(1);
	//this->filePos=sff.filePos;	
   	return *this;
}


std::string ScalarFoamFile::getDirName()
{
  std::vector<std::string> pathParts;
  pathParts = split(filePath,"/");
  return pathParts[pathParts.size()-2];
}

std::vector<std::string> split(std::string str,std::string sep)
{
	int startPos=0;
	std::vector<std::string> res;
	for(int i=0;i<int(str.length());i++)
	{
		if(str.substr(i,1)==sep)
		{
			int len=i-startPos;
			if(len>0)
				res.push_back(str.substr(startPos,len));
			else
				res.push_back("");
			startPos=i+1;
		}

		if(i == int(str.length())-1)
		{
			int len=i-startPos+1;
			if(len>0)
				res.push_back(str.substr(startPos,len));
			else
				res.push_back("");
		}

	}
	return res;
}
