#ifndef DICTFILE__H
#define DICTFILE__H


#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <list>

using namespace std;

class DictFile;

class Dictionnary : public map<string,int> 
{ // just add a const [] to map
  public :
    int Locate(const string &Key) const 
    {
      const_iterator i = find(Key);
      return (i == end()) ? -1 : i->second;
    }
};

class DictFile;

//! only its Value(string Key) routine is useful
class DictFileEntry
{
  private :
    vector<string> elements;
  const DictFile &file;

  public :
    DictFileEntry(const char *Line, const DictFile &F);


  struct Val
    {
      const string & val;
      Val( const string &Sval) : val(Sval) {};
      operator double() const { return atof(val.c_str());};
      operator string() const { return val;};
  };

  //! this routine enable double toto = line.Value("STUFF");
  Val Value(const string &Key, const bool DiesIfAbsent = true) const;

  bool HasKey(const string &Key) const;

};

  



/*! a class to access files like
# Flux : B Flux
# Fluxerr : error on B Flux
# FluxPsf : B FluxPsf
# FluxPsferr : error on B FluxPsf
# Day : Day
# Dayerr : error on Day
# AirMass : Airmass when relevant
# Absorption : Absorption when relevant
# Band : Band of obs
# Instrument : Instrument of abs
# format LightCurvePoint 1
# end
@BAND B
@INSTRUMENT toto
-13.5392 0.0509903 0 1e+30 -3.60747 0.01 0 1 B STANDARD
-13.5291 0.0412349 0 1e+30 -3.47747 0.01 0 1 B STANDARD
-13.4792 0.0412312 0 1e+30 -2.55747 0.01 0 1 B STANDARD
-13.4792 0.0509836 0 1e+30 -1.67747 0.01 0 1 B STANDARD
-13.4292 0.0510007 0 1e+30 0.442531 0.01 0 1 B STANDARD
*/
class DictFile : public list <DictFileEntry>
{
  Dictionnary dict;
  typedef   map<string, string> GKeyMap;
  GKeyMap globalKeys;
  string fileName;

  public :
    DictFile(const string &FileName);

  const Dictionnary &Dict() const { return dict;}
  const string &FileName() const {return fileName;}

  bool HasKey(const string &Key) const { return dict.Locate(Key) != -1;}
  bool HasGlobalKey( const string &Key) const 
    { return (globalKeys.find(Key) != globalKeys.end()); };

  string GlobalValue(const string &Key, const bool FatalIfAbsent = true) const;

};


typedef DictFile::const_iterator DictFileCIterator;


#ifdef EXAMPLE
example of usage of classe in this file :

DictList file(FileName);
double redshift = file.FlobalValue("REDSHIFT"); // dies if absent.
for (DictFileCIterator line = file.begin(); line != file.end(); 
  ++line)
{
  MyStruct a;
  if (line->HasKey("Stuff")) a->stuff = line->Value("Stuff");
}

#endif
    




#endif /* DICTFILE__H */
