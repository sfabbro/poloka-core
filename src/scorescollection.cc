#include "scorescollection.h"


#include <vector>
#include <map>
#include <iostream>

using namespace std;



class ScoresCollection : public map<string, vector<double> >
{

 public :
  string context;

    ScoresCollection(const string &C) : context(C) {}

    void write(ostream &S,const bool WithContext) const;

};

typedef ScoresCollection::iterator ScoresIterator;
typedef ScoresCollection::const_iterator ScoresCIterator;


#include <iomanip>

void ScoresCollection::write(ostream &S, const bool WithContext) const
{
  int oldprec = S.precision();
  S << setprecision(10);
  string prefix = (WithContext)? "@"+context+"_" : "@";
  for (ScoresCIterator i=begin(); i!=end(); ++i)
    {
      const vector<double> &vals = i->second;
      if (vals.size() == 1)
	S << prefix << i->first << ' ' << vals[0] << endl;
      else
	for (unsigned k=0; k < vals.size(); ++k)
	  S << prefix << i->first << k << ' ' << vals[k] << endl;
    }
  S << setprecision(oldprec);
}

#include <list>


typedef list<ScoresCollection> ColCol;
typedef ColCol::const_iterator ColCIterator;
typedef ColCol::iterator ColIterator;

ColCol allCols;



ScoresCollection* FindCollection(const string &C)
{
  for (ColIterator i=allCols.begin(); i!= allCols.end(); ++i)
    if (i->context == C) return &(*i);
  return NULL;
}


void StoreScore(const string &Context, const string &Name, 
		const vector<double> & Vals)
{
  ScoresCollection *c = FindCollection(Context);
  if (c== NULL)
    {
      allCols.push_back(ScoresCollection(Context));
      c = &allCols.back();
    }
  (*c)[Name] = Vals;
}

void StoreScore(const string &Context, const string &Name, const double Value)
{
  vector<double> tmp(1,Value);
  StoreScore(Context, Name, tmp);
}

void WriteScores(const string &Context, ostream &S, const bool WithContext)
{
  const ScoresCollection* c= FindCollection(Context);
  if (c) c->write(S, WithContext);
}

#include <fstream>

void WriteScores(const string &Context, const string &FileName, 
		 const bool WithContext)
{
  ofstream file(FileName.c_str());
  WriteScores(Context,file, WithContext);
}


void WriteAllScores(ostream &S, const bool WithContext)
{
  for (ColIterator i=allCols.begin(); i!= allCols.end(); ++i)
    i->write(S, WithContext);
}


void WriteAllScores(const string &FileName, const bool WithContext,
		    const bool Append)
{
  ofstream file(FileName.c_str(), (Append) ? ios::app|ios::out : ios::out);
  WriteAllScores(file, WithContext);
  file.close();
}
