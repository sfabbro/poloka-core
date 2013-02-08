#include <string>
#include <algorithm>
#include <fstream>

#include <poloka/stringlist.h>

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

void StringList::Append(const string & toapp) 
{
  
  for (StringIterator i = begin(); i != end(); ++i) 
    *i = *i + toapp ;
  return ;
}


// tested: seems OK
bool StringList::ReadFromFile(const string &FileName)
{
  ifstream l(FileName.c_str());
  if (!l) return false;
  string s;
  while (true)
    {
      l >> s;
      if (!l) return true;
      push_back(s);
    }
  l.close();
  return true;
}

void StringList::WriteinColumn(ostream &stream)
{
for (StringCIterator i = begin(); i != end(); ++i) 
      stream << *i << endl ; 
 return ;
}

void StringList::WriteinColumn(const string &FileName)
{
  ofstream pr(FileName.c_str());
  WriteinColumn(pr); pr.close();
  return ;
}

ostream & operator <<(ostream &stream, const StringList &List)
{
  for (StringCIterator i = List.begin(); i != List.end(); ++i) 
      stream << *i << ' '; return stream; 
}

