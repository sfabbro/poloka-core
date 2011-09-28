// -*- C++ -*-
// 
// \file globalval.h
// 
// Last modified: $Date: 2010/09/09 18:10:42 $
// By:            $Author: seb $
// 

#ifndef GLOBALVAL__H
#define GLOBALVAL__H


#include <string>
#include <list>
#include <vector>
#include <map>
#include <iostream>

using namespace std;

//! to store in files things like "Key value(s)" things.
// what exactly does this class bring except complicating a map?
// answer : use a map if you prefer.
class GlobalVal : private  map<string, vector<string> > {
public :

  GlobalVal() {};

  bool AddKey(const string &Key, const vector<string> &Values);
  
  bool AddKey(const string &Key, const string &Value);
  
  bool AddKey(const string& Key, const list<string>& Values);
  
  bool AddKey(const string &Key, const vector<double> &Values);
  
  bool AddKey(const string &Key, const double &Value);

  bool AddKey(const string& Key, const list<double>& Values);
  
  unsigned NKey() const;

  bool HasKey(const string &Key) const;
  
  string         getStringValue(const string& Key) const;
  
  vector<string> getStringValues(const string& Key) const;
  
  double         getDoubleValue(const string& Key) const;
  void           setDoubleValue(const string &Key, double val) ;
  vector<double> getDoubleValues(const string& Key) const;
  void           setDoubleValues(const string &Key, const vector<double>& vals);
  vector<string> OutputLines() const;
  void Dump() const;
  void AppendTo(GlobalVal & glob) const ;
  void AppendTo(map<string, string>& globalKeys) const ;

  bool ProcessLine(const string &Line);

  // this constructor is to be used when you only need to read '@' lines in a file
  GlobalVal(const std::string &FileName);


private:
  template<typename T>
  bool GenericAddKey(const string& Key, const T& Values);
};

ostream& operator << (ostream &, const GlobalVal &G);



#endif /* GLOBALVAL__H */
