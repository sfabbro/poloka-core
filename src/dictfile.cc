#include "dictfile.h"
#include "fileutils.h"
//#include "fatalerror.h"
#include <stdio.h>

void FatalError(const std::string &Message, const bool Abort=true)
{
  std::cerr << Message << std::endl << " we stop here ... " << std::endl;
  if (Abort) abort();
}


DictFileEntry::DictFileEntry(const char *Line, const DictFile &F)
  : file(F) 
{
  DecomposeString(elements, Line);
}

DictFileEntry::Val DictFileEntry::Value(const string &Key, 
				       const bool DiesIfAbsent) const
{
  int offset = file.Dict().Locate(Key);
  if (offset >= 0 && offset < int(elements.size())) return Val(elements[offset]);
  if (DiesIfAbsent) FatalError("did not find Key "+Key+"in file "
			       +file.FileName());
  return Val(" ");
}


bool DictFileEntry::HasKey(const string &Key) const
{
  return file.HasKey(Key);
}



#include <cstring>

DictFile::DictFile(const string &FileName) : fileName(FileName)
{
  FILE *f = fopen(FileName.c_str(),"r");
  if (!f) FatalError(" cannot open "+FileName);  
  char line[1024];
  bool foundEnd = false;
  while (fgets(line,1024,f))
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
	  if (words.size() != 2)
	    {
	      std::cerr << " DictFile : lines starting with '@' should have 2 words" 
			<< std::endl << " erroneous line (in  " << fileName 
			<< ") :" << std::endl
			<< start << std::endl;
	      FatalError(" Giving up ");
	    }
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
