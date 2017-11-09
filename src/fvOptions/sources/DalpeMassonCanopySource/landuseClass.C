#include "landuseClass.H"


Foam::landuseClass::~landuseClass()
{
}

Foam::landuseClass::landuseClass() {
  
  code_ = -1;
  name_ = word("default");
  // Cd_.dimensions() = dimensionSet(0,0,0,0,0,0,0); 
  Cd_=0;

  //LAI_.dimensions() = dimensionSet(0,0,0,0,0,0,0); 
  LAI_=0;

  //z0_.dimensions() = dimensionSet(0,1,0,0,0,0,0);
  z0_=0;

  //height_.dimensions() = dimensionSet(0,1,0,0,0,0,0);
  height_=0;

  //LADmax_.dimensions() = dimensionSet(0,-1,0,0,0,0,0);
  LADmax_=0;
}


Foam::landuseClass::landuseClass(const dictionary& dict, word name) {

  
  dictionary landuseClassDict(dict.subDict(name));

  landuseClassDict.lookup("code") >> code_;

  name_ = name;
  // Cd_.dimensions() = dimensionSet(0,0,0,0,0,0,0); 
  // LAI_.dimensions() = dimensionSet(0,0,0,0,0,0,0);
  // z0_.dimensions() = dimensionSet(0,1,0,0,0,0,0);
  // height_.dimensions() = dimensionSet(0,1,0,0,0,0,0);
  // LADmax_.dimensions() = dimensionSet(0,-1,0,0,0,0,0);

  Cd_ = landuseClassDict.lookupOrDefault("Cd", 0.2, false, false);
  height_ = landuseClassDict.lookupOrDefault("height", 0.0, false, false);
  z0_ = landuseClassDict.lookupOrDefault("z0", 0.001, false, false);
  LAI_ = landuseClassDict.lookupOrDefault("LAI", 0.0, false, false);
  LADmax_ = landuseClassDict.lookupOrDefault("LADmax", -1.0, false, false);

  landuseClassDict.lookup("LADProfile", true) >> LADProfile_;

  if (LADmax_ == -1) {
    LADmaxFromLAI();
  }
}

Foam::scalar Foam::landuseClass::integrateLAD()
{
  //uses simpsons rule to integrate between the points
  scalar LAItmp=0;
  for(label fi=1;fi<LADProfile_.size();fi++)
    {
      if(LADProfile_[fi]==0 or LADProfile_[fi-1]==0)
	LAItmp += 1/2.0*abs(LADProfile_[1]-LADProfile_[0])*LADmax_*0.1*height_;
      else
	LAItmp += min(LADProfile_[fi],LADProfile_[fi-1])*LADmax_*0.1*height_
	  +1/2.0*abs(LADProfile_[fi]-LADProfile_[fi-1])*0.1*height_;
    }
  return LAItmp;
}

void Foam::landuseClass:: LADmaxFromLAI()
{
  LADmax_ = 0;
  scalar LAItmp = 0;
  while(LAItmp < LAI_)
    {
      LAItmp = Foam::landuseClass::integrateLAD();
      LADmax_ += 0.01;
    }
}
