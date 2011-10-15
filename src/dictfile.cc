#include "dictfile.h"
#include "fileutils.h"
//#include "fatalerror.h"
#include <stdio.h>
#include <iostream>
#include <fstream>


void FatalError(const std::string &Message, const bool Abort=true)
{
  std::cerr << Message << std::endl << " we stop here ... " << std::endl;
  if (Abort) abort();
}

DictFileEntry::DictFileEntry(DictFile& F)
  : file(F)
{
  // maybe we should update the elements here.
}

DictFileEntry::DictFileEntry(const char *Line, DictFile &F)
  : file(F) 
{
  DecomposeString(elements, Line);
}

DictFileEntry::Val DictFileEntry::Value(const string &Key, 
				       const bool DiesIfAbsent) const
{
  int offset = file.Dict().Locate(Key);
  if (offset >= 0 && offset < int(elements.size())) 
    return Val(elements[offset]);
  if (DiesIfAbsent) FatalError("did not find Key "+Key+"in file "
			       +file.FileName());
  cout << "return 'empty' string" << endl;
  return Val("empty");
}


bool DictFileEntry::HasKey(const string &Key) const
{
  return file.HasKey(Key);
}


void DictFileEntry::writen(ofstream & pr) const
{
  for (int ii = 0 ; ii < file.Dict().size() ; ii++)
    pr << elements[ii] << " " ;

}

      
void DictFileEntry::AddKey(const string &Key, const string &Val)
{
  if (!file.HasKey(Key)) file.AddKey(Key);
  int offset = file.Dict().Locate(Key);
  if (offset == int(elements.size()))
    elements.push_back(Val);
  else
    cerr << " DicFileEntry::AddKey problem with Key " << Key 
	 << " Val " << Val << endl;
}

void DictFileEntry::AddKey(const string &Key, const double &Val)
{
  char s[80];
  sprintf(s,"%12.12E",Val);
  AddKey(Key,string(s));
}


void DictFileEntry::ModKey(const string &Key, const string &Val)
{
  if (!file.HasKey(Key)) {
    cerr << " DicFileEntry::ModKey no such key " << Key << " use AddKey" << endl;
    return;
  }
  int offset = file.Dict().Locate(Key);
  elements[offset]=Val;
}

void DictFileEntry::ModKey(const string &Key, const double &Val)
{
  char s[80];
  sprintf(s,"%f",Val);
  ModKey(Key,string(s));
}



#include <cstring>

DictFile::DictFile()
  : fileName("")
{
  // no entries
}

DictFile::DictFile(const string &FileName) : fileName(FileName)
{
  FILE *f = fopen(FileName.c_str(),"r");
  if (!f) FatalError(" cannot open "+FileName);  
  char line[4096];
  bool foundEnd = false;
  while (fgets(line,4096,f))
    {
      // skip leading spaces
      char *start = line;
      while (*start == ' ') start++;
      // cut at \n
      char *cr = rindex(start,'\n'); if (cr) *cr = '\0';

      if (*start == '\0') continue;
      if (*start == '#')
	{
	  start ++; // skip '#'
	  while (*start == ' ') start++;
	  if (strstr(start,"end") == start) {foundEnd = true; continue;}
	  char *column = strchr(start,':');
	  if (column)
	    {
	      if (foundEnd)
		{
		  std::cerr << " bizarre file structure : found tags after end"
			    << std::endl;
		}
	      column = start + strcspn(start," \t:");
	      *column = '\0';
	      int presentSize = dict.size();
	      string tag(start);
	      dict[tag] = presentSize;
	    }
	}
      else if (*start == '@')
	{
	  vector<string> words;
	  DecomposeString(words, start+1);
	  // IGNORE ... 
	  //	  if (words.size() != 2)
	  //	    {
	  //	      std::cerr << " DictFile : lines starting with '@' should have 2 words" 
	  //			<< std::endl << " erroneous line (in  " << fileName 
	  //			<< ") :" << std::endl
	  //			<< start << std::endl;
	  //	      FatalError(" Giving up ");
	  //	    }
	  if (HasGlobalKey(words[0]))
	    FatalError(" DictFile : Key " +words[0]+
		       " appears more than once in "+fileName);
	  globalKeys[words[0]] = words[1];
	}	  
      else // no '@, no '#'
	{
	  push_back(DictFileEntry(start,*this));
	}
    }
  fclose(f);
}

void DictFile::DumpKeys() const 
{
  for(Dictionnary::const_iterator it = dict.begin(); it!=dict.end(); ++it) {
    cout << it->second << " " << it->first << endl;
  }
}





bool DictFile::Write(const string &FileName) const
{
  unsigned presentSize = dict.size();
  FILE * f = fopen(FileName.c_str(), "w");
  if (!f)
    {
      cerr << " DistFile::Write() : could not open " << FileName << endl;
      return false;
    }

  // write global keys and values
  for(GKeyMap::const_iterator it = globalKeys.begin() ; it != globalKeys.end(); it++) {
    fprintf(f,"@%s %s\n",it->first.c_str(),it->second.c_str());
  }


  // invert the dictionnary:
  map<int,string> tags;
  for (Dictionnary::const_iterator it = dict.begin(); it != dict.end(); ++it)
    tags[it->second] = it->first;

  //write the header
  for (unsigned i = 0; i < presentSize; ++i)
    fprintf(f,"# %s :\n",tags[i].c_str());

  fprintf(f,"# end\n");

  // write data
  for (const_iterator it = begin(); it != end(); ++it)
    {
      const DictFileEntry &entry = *it;
      for (unsigned int i =0; i < presentSize; ++i)
	{
	  const string& val(entry.Value(tags[i]));
	  fprintf(f,"%s ", val.c_str());
	}
      fprintf(f,"\n");
    }
  fclose(f);
  return true;
}


bool DictFile::Write_Data(ofstream & pr) const
{unsigned presentSize = dict.size();
// invert the dictionnary:
  map<int,string> tags;
  for (Dictionnary::const_iterator it = dict.begin(); it != dict.end(); ++it)
    {
      tags[it->second] = it->first; 
      //cerr  << "tag : " << tags[it->second] << endl ;
    }
  // write data
  for (const_iterator it = begin(); it != end(); ++it)
    {
      const DictFileEntry &entry = *it;
      for (unsigned int i =0; i < presentSize; ++i)
	{
	  string val = entry.Value(tags[i]);
	  pr << val << " " ;
	}
      pr << endl ;
    }
  return true;
}


void DictFile::AddGlobalKeys(const GKeyMap & Keys )
{
  for(GKeyMap::const_iterator it = Keys.begin() ; it != Keys.end(); it++) {
    AddGlobalKey(it->first.c_str(),it->second.c_str() );
  }
}

void DictFile::AddKey(const string &Key)
{
  if (!HasKey(Key))
    {
      int presentSize = dict.size();
      dict[Key] = presentSize;
    }
}

void DictFile::RmKey(const string &Key)
{ 
  if (!HasKey(Key))
    {
      cerr << "No Key " << Key << " found" << endl ;
      return ;
    }
  int num = dict.Locate(Key);  
  // modifies data
  for (iterator it = begin(); it != end(); ++it)
    {
      DictFileEntry &entry = *it;
      if (!entry.EraseElement(num)) //pb ds le erase
	return;
    }

// et on l'enleve du dictionnaire
  dict.erase(Key);
}


string DictFile::GlobalValue(const string &Key, const bool ExitIfAbsent) const
{			  
  GKeyMap::const_iterator i = globalKeys.find(Key);
  if (i == globalKeys.end())
    {
      std::cerr << " DictFile : trying to read global key '" 
		<< Key << "'" << std::endl
		<< " which was not found in file " << fileName << std::endl;
      if (ExitIfAbsent) FatalError(" giving up");
      return "";
      }
  return i->second;
}


void DictFile::RemoveGlobalKey( const string &Key) {
  if( HasGlobalKey(Key) ) {
    globalKeys.erase(globalKeys.find(Key));
  }else{
    cerr << "warning no such key " << Key << " in global values of " << fileName << endl;
  }
}

void DictFile::AddGlobalKey( const string &Key, const string &Val) {
  if( HasGlobalKey(Key) ) {
    cerr << "warning overwriting " << Key << "=" << GlobalValue(Key) << " in " << fileName << endl;
  }
  globalKeys[Key]=Val;
}

