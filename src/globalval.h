// -*- C++ -*-
// 
// \file globalval.h
// 

#ifndef GLOBALVAL__H
#define GLOBALVAL__H


#include <string>
#include <list>
#include <vector>
#include <map>


using namespace std;

//! to store in files things like "Key value(s)" things.
class GlobalVal : private  map<string, vector<string> > {
public :
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
  
  vector<double> getDoubleValues(const string& Key) const;

  vector<string> OutputLines() const;

  bool ProcessLine(const string &Line);


private:
  template<typename T>
  bool GenericAddKey(const string& Key, const T& Values);
};



#endif /* GLOBALVAL__H */
