#include <string>
#include "stringlist.h"
#include "algorithm"

StringCIterator StringList::Locate( const string &Needle) const
{
  StringCIterator i = begin();
  for ( ; i!= end() ; ++i) if (*i == Needle) break;
  return i;
}

StringIterator StringList::Locate( const string &Needle)
{
  return find(begin() , end(), Needle);
}

bool StringList::Contains(const string &Needle) const
{
  return (Locate(Needle) != end());
}

void StringList::Remove(const string &Banned)
{
  StringIterator i = find(begin(), end(), Banned);
  if (i!=end()) erase(i);
}

void StringList::Substitute(const string &Original, const string &Substitution)
{
  StringIterator i = Locate(Original);
  if (i != end()) *i = Substitution;
}

string StringList::AllEntries() const
{
  string result;
  for (StringCIterator i = begin(); i != end(); ++i) result += *i+" ";
  return result;
}

ostream & operator <<(ostream &stream, const StringList &List)
{
  for (StringCIterator i = List.begin(); i != List.end(); ++i) 
      stream << *i << ' '; return stream; 
}

