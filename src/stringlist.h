#ifndef STRINGLIST__H
#define STRINGLIST__H
#include <list>
#include <string>
#include <iostream>

using namespace std;


typedef list<string>::iterator StringIterator;
typedef list<string>::const_iterator StringCIterator;

class StringList : public list<string> {
  public:
  //!
  StringCIterator Locate( const string &Needle) const;
  //!
  StringIterator Locate( const string &Needle);
  //!
  bool Contains(const string &Needle) const;
  //!
  void Remove(const string &Banned);
  //!
  void Substitute(const string &Original, const string &Substitution);
};

#ifndef SWIG
ostream & operator <<(ostream &stream, const StringList &List);
#endif

#endif /* STRINGLIST__H */
