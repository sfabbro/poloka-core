#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>



#include "dicstar.h"


#include "fastifstream.h"
#include "starlistexception.h"
#include <polokaexception.h>

DicStar::DicStar()
  : BaseStar(0.,0.,0.) {
  Set_to_Zero();
  format = "BaseStar 2 ";
  firstkeys.push_back("x");
  firstkeys.push_back("y");
  firstkeys.push_back("flux");
}

DicStar::DicStar(double xx, double yy, double ff)
  : BaseStar(xx,yy,ff) {
  Set_to_Zero();
  format = "BaseStar 2 ";
  firstkeys.push_back("x");
  firstkeys.push_back("y");
  firstkeys.push_back("flux");
}

DicStar::DicStar(const std::vector<string>& firstKeys, const std::vector<string>& newkeys, const char *Format) 
  : BaseStar(), format(Format)
{
  Set_to_Zero();
  firstkeys = firstKeys;
  for(unsigned int i=0; i<newkeys.size() ;i++) {
    key.push_back(newkeys[i]);
    val.push_back(0);
  }
}

bool DicStar::HasKey(const string &Key) const
{
  for(unsigned int i=0; i<key.size() ;i++)
    if(key[i]==Key) return true;
  return false;
}

void DicStar::AddKey(const string &KeyName, const double &Val)
{
  key.push_back(KeyName);
  val.push_back(Val);
}


unsigned DicStar::NKeys() const
{
  return (key.size()+3);
}


void DicStar::Set_to_Zero() 
{
  x = y = flux = 0;
  rank = 0;
  for(unsigned int i=0;i<val.size();i++) val[i] = 0;
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

double DicStar::getval(const string &thekey) const {
  for(unsigned int i=0; i<key.size() ;i++) {
    if(key[i]==thekey) {
      return val[i];
    }
  }
  if (thekey == firstkeys[0]) return x;
  if (thekey == firstkeys[1]) return y;
  if (thekey == firstkeys[2]) return flux;

  
  cerr << "DicStar::getval unknown key " << thekey << endl; 
  throw(PolokaException(string("DicStar::getval unknown key ")+thekey));
  return 0;
}


void DicStar::Read(fastifstream& r, const char *Format) {
  
  BaseStar::read_it(r, Format);
  //int format = DecodeFormat(Format, "DicStar");
  for(unsigned int i=0;i<val.size();i++)
    r >> val[i];

}


DicStar* DicStar::read(const std::vector<string> &firstKeys, const std::vector<string>& newkeys, fastifstream& r, const char *Format) {
  
  if (firstKeys.size() <3)
    throw(StarListException(" need at least 3 keys to read a DicStar "));
  DicStar *pstar = new DicStar(firstKeys,newkeys,Format);
  pstar->Read(r,Format);
  return(pstar);
}

DicStar* DicStar::read( fastifstream& r, const char *Format) {
  DicStar *pstar = new DicStar();
  abort(); // this routine should not be used : it does not know anything about the keys.
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
  if (firstkeys.size()==3)
    s << x << " "  << y << " " << flux << " " ;
  else BaseStar::writen(s);
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
  
  return format; // stored when reading or constructing 
}


/************************** FINDEFINITION DicStar ************************/

#include "starlist.h"
#include "starlist.cc" /* since starlist is a template class */

template class StarList<DicStar>;

DicStarList::DicStarList(const string &FileName) {
  init(FileName);
}
void DicStarList::init(const string &FileName) {
  fastifstream rd(FileName.c_str());
  vector<string> allKeys;

  if (!rd)
    {
      throw(StarListException("DicList cannot open :"+FileName));
      return ;
    }
  char c ;
#define LINE_LENGTH 4096
  char buff[LINE_LENGTH];
  char onekey[400];
  char column[256];
  const char *format = NULL;
  ClearList();
  int star_count = 0;
  while( rd >> c ) // to test eof
    {
      rd.unget() ;
      // nrl: we do not want to ignore the '@'s anymore 
      if (c == '@') 
	{ 
	  rd.getline(buff,LINE_LENGTH);
	  if (star_count)
	    {
	      cout << " file : " << FileName 
		   << ": ignoring header line found after actual data" << endl
		   << buff << endl;
	      continue;
	    }
	  GlobVal().ProcessLine(buff);
	  continue;
	}	  
      else if (c == '#')
        {	  
	  rd.getline(buff,LINE_LENGTH);
	  if (star_count)
	    {
	      cout << " file : " << FileName 
		   << ": ignoring header line found after actual data" << endl
		   << buff << endl;
	      continue;
	    }
	  /* hack something reading " format <StarType> <integer>" to drive the decoding (in Star::read) */
	  char *p =buff+1; /* skip '#' */
	  p += strspn(p," \t"); /* skip leading spaces */
	  if (strstr(p,"format") == p)  /* this test is enough because the format is the last line of the header ... */
	    {
	      char *toto = p + strlen("format");
	      size_t l = strlen(toto)+2;
	      // this s a memory leak , but it only happens once per list.
	      toto = (char *) malloc(l);
	      strcpy(toto,p + strlen("format"));
	      format = toto;
	    }
          else // it should be a regular key, or "end"
	    {
	    //	    if(line>3)
	    if(sscanf(buff+1," %[^: ] %s",onekey,column)==2 
	       && column[0] == ':' ) 
	      {
		allKeys.push_back(onekey);
		continue;
	      }
	    if(sscanf(buff,"# %s",onekey) == 1 && strcmp(onekey,"end")==0) 
	      {
		// DEBUG
	      }
	    }
	}
      else // no '#', no '@'
	{
	  if (firstKey.size() == 0) // handle the split between BaseStar keys and others
	    {	      
	      unsigned baseStarLength;
	      if ( format == NULL || (baseStarLength = NValsBaseStar(format)) < 3)
		{
		  cout << "*******************************************************************" << endl;
		  cout << " the file " << FileName << " does not contain a BaseStar format tag" << endl;
		  cout << " it is very likely that the decoding of this file is going nuts. " << endl;
		  cout << " If your file begins with x,y, flux, add \'#format BaseStar 2\' in its header. " << endl;
		  cout << " If your file begins with x,y, sx,sy,rhoxy, flux, add \'#format BaseStar 3\' in its header. " << endl;		    
		  cout << "*******************************************************************" << endl;
		  baseStarLength = 3;
		  format = "BaseStar 2 ";
		}
	      for (unsigned k=0; k<allKeys.size(); ++k)
		if (k<baseStarLength) firstKey.push_back(allKeys[k]); else key.push_back(allKeys[k]);
	      }	      
	  DicStar* s = DicStar::read(firstKey,key,rd, format); 
	  if (rd.fail())
	    {
	      if (s) delete s;
	      throw(StarListException("bad extraction in DicStarList reader, file="+FileName));
	    }
	  if (!s) 
	    {
	      return ;
	    }
	  s->rank = ++star_count; // call the first one 1 !
	  push_back(s);
	}
    }
}

/* return an "empty" star which has the same fields as the 
   o;er ones of the list */
DicStar *DicStarList::EmptyStar() const
{
  if (empty())
    {
      cerr << " requested a DicStarList::EmptyStar() on an empty list"
	   << " which is meaningless " << endl;
      return NULL;
    }
  DicStar *empty = new DicStar(*front());
  empty->Set_to_Zero();
  return empty;
}
  
bool DicStarList::HasKey(const string &Key) const
{
  return(front()->HasKey(Key));
}
  

BaseStarList* Dic2Base(DicStarList* This)
{
  return (BaseStarList*)This;
}

const BaseStarList* Dic2Base(const DicStarList* This)
{
  return (const BaseStarList*) This;
}

BaseStarList& Dic2Base(DicStarList& This)
{
  return (BaseStarList&)This;
}

const BaseStarList& Dic2Base(const DicStarList& This)
{
  return (const BaseStarList&) This;
}



