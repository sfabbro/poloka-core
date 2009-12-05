// 
// \file globalval.cc
// 
// Last modified: $Date: 2009/12/05 08:24:12 $
// By:            $Author: astier $
// 
#include <iostream>
#include <sstream>
#include "globalval.h"
#include "fileutils.h" // for DecomposeString 
#include <stdlib.h> // for atof
#include <string.h>// for strlen


bool GlobalVal::HasKey(const string &Key) const
{
  return (find(Key) != end());
}


unsigned GlobalVal::NKey() const
{
  return this->size();
}


template<typename T>
bool GlobalVal::GenericAddKey(const string& Key, const T& Values) {
  if (HasKey(Key))
    {
      cerr << " cannot have twice the same key val in GlobalVal " << Key 
	   << endl;
      return false;
    }
    

  typename T::const_iterator I;
  for(I=Values.begin();I!=Values.end();I++) {
    stringstream sstrm;    
    sstrm << *I;
    (*this)[Key].push_back(sstrm.str());
  }
  return true;
}


bool GlobalVal::AddKey(const string& Key, const list<string>& Values)
{
  return GenericAddKey(Key, Values);
}


bool GlobalVal::AddKey(const string &Key, const vector<string> &Values)
{
  if (HasKey(Key))
    {
      cerr << " cannot have twice the same key val in GlobalVal " << Key 
	   << endl;
      return false;
    }
  
  (*this)[Key] = Values;
  return true;
}

bool GlobalVal::AddKey(const string &Key, const string &Value)
{
  vector<string> tmp;
  tmp.push_back(Value);
  return AddKey(Key,tmp);
}


bool GlobalVal::AddKey(const string &Key, const vector<double> &Values)
{
  if (HasKey(Key))
    {
      cerr << " cannot have twice the same key val in GlobalVal " << Key 
	   << endl;
      return false;
    }
  
  unsigned int i;
  for(i=0;i<Values.size();i++) {
    stringstream sstrm;
    sstrm << Values[i];
    (*this)[Key].push_back(sstrm.str());
  }
  //  (*this)[Key] = Values;
  return true;
}


bool GlobalVal::AddKey(const string& Key, const list<double>& Values)
{
  return GenericAddKey(Key, Values);
}


bool GlobalVal::AddKey(const string &Key, const double &Value)
{
  vector<double> tmp;
  tmp.push_back(Value);
  return AddKey(Key,tmp);
}


string GlobalVal::getStringValue(const string& Key) const
{
  const_iterator i = find(Key);
  if (i == end())
    {
      cerr << " could not read key " << Key << endl;
      return "";
    }
  else return i->second[0];
}


vector<string> GlobalVal::getStringValues(const string& Key) const
{
  const_iterator i = find(Key);
  if(i == end())
    {
      cerr << " could not read key " << Key << endl;
      vector<string> ret;
      return ret;
    }
  return i->second;
}


double GlobalVal::getDoubleValue(const string &Key) const
{
  const_iterator i = find(Key);
  if (i == end())
    {
      cerr << " could not read key " << Key << endl;
      return 99;
    }
  else return atof(i->second[0].c_str());
}

void GlobalVal::setDoubleValue(const string &Key, double val) 
{
  iterator i = find(Key);
  if (i == end())
    {
      cerr << " could not read key " << Key << endl;
      return ;
    }
  else
    {
      char c[100] ;
      sprintf(c,"%f",val);
      string sc = c ;
      i->second[0] = c ;
    }
  return ;
}




vector<double> GlobalVal::getDoubleValues(const string &Key) const
{
  vector<double> ret;
  const_iterator i = find(Key);
  if (i == end())
    {
      cerr << " could not read key " << Key << endl;
      return ret;
    }
  unsigned int k;
  for(k=0;k<i->second.size();k++)
    ret.push_back(atof(i->second[k].c_str()));
  return ret;
}


#include <sstream>

vector<string> GlobalVal::OutputLines() const
{
  vector<string> out;
  for (const_iterator i = begin(); i != end(); ++i)
    {
      //      const vector<double> &values = i->second;
      const vector<string> &values = i->second;
      if (values.size() == 0) continue;
      ostringstream s;
      s << i->first;
      for (unsigned k=0; k < values.size(); ++k) s << ' ' << values[k];
      out.push_back(s.str());
    }
  return out;
}


#include <fstream>
GlobalVal::GlobalVal(const std::string &FileName)
{
  std::ifstream r(FileName.c_str());
  char c ;
  char buff[4096];
  while( r >> c ) // to test eof
    {
      r.unget() ;
      if ( (c == '@') ) 
	{
	  r.getline(buff,4096); 
	  ProcessLine(buff);
	  continue;
	}
      if ( (c=='#') || isdigit(c)) break;
    }
  r.close();
}


#include <stdlib.h>

#ifdef STORAGE
static void explode(vector<string>& v, string const& str, char tok=',')
{
  string s, sc=str;
  string::size_type /* p0=0, */ p1;

  if(str == "*") {
    stringstream sstrm;
    for(int i=0;i<5;i++) {
      sstrm << i;
      v.push_back(sstrm.str());
      sstrm.str("");
    }
    return;
  }

  do {
    p1 = sc.find_first_of(tok);
    s = sc.substr(0, p1);
    v.push_back(s);
    if(p1 != string::npos)
      sc = sc.substr(p1+1, string::npos);
  } while(p1 != string::npos);

}
#endif


bool GlobalVal::ProcessLine(const string &Line)
{
  // use standard C stuff because it enables error checking
  const char *s = Line.c_str();
  if (*s == '@') s++;

  char buf[256];
  if (sscanf(s,"%s",buf)!= 1)
    {
      cerr << " could no read a key in " << Line << endl;
      return false;
    }
  s += strlen(buf);
  string key(buf);
  if (HasKey(key))
    {
      cerr <<" already have a key labelled " << key << endl;
      return false;
    }

  vector<string> &values = (*this)[key];
  string strtmp = s;
  DecomposeString(values, strtmp, " ");
  return true;
}

