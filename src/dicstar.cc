#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>



#include "dicstar.h"


DicStar::DicStar()
  : BaseStar(0.,0.,0.) {
  Set_to_Zero();
}

DicStar::DicStar(double xx, double yy, double ff)
  : BaseStar(xx,yy,ff) {
  Set_to_Zero();
}

DicStar::DicStar(const std::vector<string>& firstKeys, const std::vector<string>& newkeys)
  : BaseStar() {
  Set_to_Zero();
  firstkeys = firstKeys;
  for(unsigned int i=0; i<newkeys.size() ;i++) {
    key.push_back(newkeys[i]);
    val.push_back(0);
  }
}



void DicStar::Set_to_Zero() {
  
  
}

int DicStar::setval(const string &thekey,double newval) {
  for(unsigned int i=0; i<key.size() ;i++) {
    if(key[i]==thekey) {
      val[i]=newval;
      return 1;
    }
  } 
  cerr << "DicStar::setval unknown key " << thekey << endl; 
  return 0;
}

double DicStar::getval(const string &thekey) {
  for(unsigned int i=0; i<key.size() ;i++) {
    if(key[i]==thekey) {
      return val[i];
    }
  }
  cerr << "DicStar::getval unknown key " << thekey << endl; 
  return 0;
}



void DicStar::Read(istream& r, const char *Format) {
  
  BaseStar::read_it(r, Format);
  //int format = DecodeFormat(Format, "DicStar");
  for(unsigned int i=0;i<val.size();i++)
    r >> val[i];

}


DicStar* DicStar::read(const std::vector<string> &firstKeys, const std::vector<string>& newkeys, istream& r, const char *Format) {
  
  DicStar *pstar = new DicStar(firstKeys,newkeys);
  pstar->Read(r,Format);
  return(pstar);
}

DicStar* DicStar::read( istream& r, const char *Format) {
  
  DicStar *pstar = new DicStar();
  pstar->Read(r,Format);
  return(pstar);
}
void DicStar::dumpn(ostream& s) const {
  s << " x : " << x;
  s << " y : " << y;
  s << " flux : " << flux;
  
  for(unsigned int i=0;i<val.size();i++) {
    s << " " << key[i] << " : " << val[i];
  }
}

void DicStar::dump(ostream& s) const {
  dumpn(s);
  s << endl;
}

void DicStar::writen(ostream& s) const {
  s << x << " " ;
  s << y << " " ;
  s << flux << " " ;
  for(unsigned int i=0;i<val.size();i++) {
    s << val[i]  << " " ;
  }
}

void DicStar::write(ostream& s) const {
  writen(s);
  s << endl; 
}

string DicStar::WriteHeader_(ostream & pr, const char *i) const {
  
  if (i==NULL) i="";
  //  string baseStarFormat =  BaseStar::WriteHeader_(pr, i);
  

  for(unsigned int j=0;j<firstkeys.size();j++) 
    pr << "# " << firstkeys[j] << i << " : " << endl;
  

  for(unsigned int j=0;j<key.size();j++) {
    pr << "# " << key[j] << i << " : " << endl;
  }
  
  //static char format[256];
  string format = " DicStar 1";
  //sprintf(format,"%s DicStar %d",baseStarFormat, 1);
  return format;
}


/************************** FINDEFINITION DicStar ************************/

#include "starlist.h"
#include "starlist.cc" /* since starlist is a template class */

template class StarList<DicStar>;

DicStarList::DicStarList(const string &FileName) {
  
  
  
  ifstream rd(FileName.c_str());
  if (!rd)
    {
      cout << "DicList cannot open :" << FileName << endl;
      return ;
    }
  char c ;
  char buff[400];
  char onekey[400];
  char column[256];
  char *format = 0;
  ClearList();
  int line = 0;
  while( rd >> c ) // to test eof
    {
      rd.unget() ;
      if ( (c == '#') ) // we jump over the line  (not always ...)
        {
	  
	  rd.getline(buff,400);
	  line++;
	  /* hack something reading " format <StarType> <integer>" to drive the decoding (in Star::read) */
	  char *p = strstr(buff,"format");
	  if (p) /* this test is enough because the format is the last line of the header ... */
	      format = p + strlen("format");
          else
	    {
	    //	    if(line>3)
	    if(sscanf(buff,"# %s %s",onekey,column)==2 
	       && column[0] == ':') 
	      {
		if (firstKey.size() <3) firstKey.push_back(onekey);
		else  key.push_back(onekey);
		continue;
	      }
	    if(sscanf(buff,"# %s",onekey) == 1 && strcmp(onekey,"end")==0) 
	      {
	      for(unsigned int i = 0 ; i< key.size() ; i++)
		cout << key[i] << " ";
		cout << endl;
	      }
	    }
	}
      else // no '#'
	{
	  DicStar* s = DicStar::read(firstKey,key,rd, format); 
	  if (!s) 
	    {
	      return ;
	    }
	  push_back(s);
	}
    }
}

#ifdef USE_ROOT
template class StarListWithRoot<DicStar>;
ClassImpT(StarListWithRoot,DicStar);

/* comments to drive the Makefile part that runs rootcint
RUN_ROOTCINT

LINKDEF_CONTENT : #pragma link C++ class DicStar;
LINKDEF_CONTENT : #pragma link C++ class list<DicStar*>;
LINKDEF_CONTENT : #pragma link C++ class StarList<DicStar>-;
LINKDEF_CONTENT : #pragma link C++ function operator << (ostream&, const StarList<DicStar>&);
LINKDEF_CONTENT : #pragma link C++ class StarListWithRoot<DicStar>-;
LINKDEF_CONTENT : #pragma link C++ class StarList<DicStar>::iterator;
LINKDEF_CONTENT : #pragma link C++ typedef DicStarIterator;
*/
#include "root_dict/jimstardict.cc"
#endif /* USE_ROOT */
