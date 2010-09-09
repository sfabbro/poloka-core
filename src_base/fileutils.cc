#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "fileutils.h"
#include "stringlist.h"

#include <string>
#include <cctype>
#include <iostream>
#include <algorithm>
#include <iterator>
#include "polokaexception.h"


void DirName(const char *FullFileName, char *DirName)  /* DirName is assumed to point to an area at least as long as FullFileName */
{
char *end;
strcpy(DirName,FullFileName);
end = DirName+strlen(DirName) - 1;
if (*end == '/') *end = '\0'; /* ignore a trailing / as the dirname command */
end = strrchr(DirName,'/');
if (end) *end = '\0';
else strcpy(DirName,".");
}

string DirName(const string &FullFileName)
{
char toto[256];
DirName(FullFileName.c_str(), toto);
return string(toto);
}

void BaseName(const char *FullFileName, char *Name)
{
char a_path[256];
const char *begin;
strcpy(a_path,FullFileName);
if (a_path[strlen(a_path)-1] == '/') a_path[strlen(a_path)-1] = '\0';
begin = strrchr(a_path,'/');
if (!begin) strcpy(Name, a_path);
else strcpy(Name,begin + 1);
}

string BaseName(const string &Name)
{
char toto[512];
BaseName(Name.c_str() ,toto);
return string(toto);
}


/* 
   return the dot position where extension starts : for toto.fits.gz,
   returns 4. in dbimage.cc, we assume that extension of toto.fits.gz
   is fits.gz rather than gz. So, do NOT change here!
*/
static unsigned dotposition(const string &FileName)
{
  const char *slash =  strrchr(FileName.c_str(),'/');
  if (!slash) slash = FileName.c_str();
  const char *dot = strchr(slash,'.');
  if (!dot) return FileName.size();
  return dot-FileName.c_str();
}



string CutExtension(const string &FileName)
{
  return FileName.substr(0,dotposition(FileName));
}



std::string FileExtension(const std::string &FileName)
{
  unsigned index = dotposition(FileName);
  if (index == FileName.size()) return " ";
  return FileName.substr(dotposition(FileName)+1);
}



int MKDir(const char *path, const bool Warn)
{
  if (FileExists(path)) 
    {
      if (IsDirectory(path)) return 1;
      else 
	{
	  cerr << " MKDir : regular file exists at this path : " << path << endl;
	  throw(PolokaException(" MKDir : regular file exists, file="+string(path)));
	}
    }
char command[1024];
sprintf(command,"mkdir %s",path);
int status  = system(command);
// cout <<" debug :status for " << command << ' ' << status << endl;
if (status)
  {
    if (Warn) cerr << " MKDir : could not create " << path << endl;
    throw(PolokaException("(IO) MKDir could not create directory "+string(path)));

  }
 return (!status);
}


string FullFileName(const string &FileName)
{
if (strlen(FileName.c_str()) && FileName[0] == '/') return FileName;
//char cwd[512];
//getcwd(cwd,512); /* very nice : may not be equal to $cwd from shell (if there are links involved) */
return  AddSlash(getenv("PWD")) + FileName;
}

static string RemoveDoubleSlash(const string &FileName)
{
  std::string result = FileName;
    // string::find does not tell clearly when the pattern is not found
    while (strstr(result.c_str(),"//"))
      {
	int pos = result.find("//");
	result.erase(pos,1);
      }
    return result;
}


string RelativeLink(const char *RealFile, const char *LinkLocation)
  //  caveats : does not check if there are symbolic links in the given pathes.. 
  //            misses physically identical pathes if they do not have the same name
  //  this kind of utility routine should be stolen in shell or unix source code.
{
string realFile = FullFileName(RealFile);
string link     = FullFileName(LinkLocation);
//
 realFile = RemoveDoubleSlash(realFile);
 link     = RemoveDoubleSlash(link);

 if (realFile == link) return ""; 


/* find the common part */
int bound = min(realFile.length(), link.length());
int i;
for (i = 0; i< bound; ++i) 
   if (link[i] != realFile[i]) break;
while (link[i] != '/') i--;
i++;
if (i == 1) /* they only share the leading '/' */ return realFile;
  /* assemble the result */  
string result;
for (const char* p = link.c_str() + i; (p=strchr(p,'/')) ; )
          { result += "../"; ++p; }
result = result + (realFile.c_str()+i);
return result;
}

int MakeRelativeLink(const std::string &RealFile, 
		     const std::string &LinkLocation)
{
  string link_value = RelativeLink(StandardPath(RealFile).c_str(), 
				   StandardPath(LinkLocation).c_str());

  if (link_value =="") return 1; // trying to link a file on itself

  // if the link already exists, symlink does not overwrite it so
  remove(LinkLocation.c_str());
  if (0==symlink(link_value.c_str(), LinkLocation.c_str()))
    return 1;
  std::cerr << " could not create a link " << std::endl
	    << LinkLocation << "->" << link_value << std::endl;
  return 0;
  //string command = " ln -fs " + link_value + " " + LinkLocation;
  //if (system(command.c_str()) == 0) return 1;
  //cerr << " could not issue " << command << endl;
  //return 0;
}

string AddSlash(const string &path)
{
if (path[path.length()-1] == '/') return path ; else return path+"/";
}


int FileExists(const char *FileName)
{
return (access(FileName,F_OK) == 0);
}

int FileExists(const string &FileName)
{
return (access(FileName.c_str(),F_OK) == 0);
}


bool IsOfRightSize(const string &FileName, int rightsize)
{
  struct stat a_stat;
  stat(FileName.c_str(),&a_stat);
  cerr << "Size of file " << FileName << ": " << a_stat.st_size << endl ;
  return(a_stat.st_size == rightsize);
}
typedef   std::vector<std::string> SVect;


bool RemoveFiles(const std::string &FileNames)
{
  if (FileNames == "") return true;
  SVect files;
  DecomposeString(files,FileNames);
  bool ok = true;
  for (unsigned i = 0; i< files.size(); ++i)
    {
      const std::string &file = files[i];

      cerr << "Removing " << file.c_str() << "...." << endl ;
      if ( !FileExists(file.c_str()))
	{
	  cerr << " ... file " <<  file.c_str() << " doesn't exist" << endl ;
	}
      else
	{
	  if (0 != remove(file.c_str()))
	    {
	      std::cerr << " something weird happened when removing " 
			<< file << std::endl;
	      ok = false;
	    }
	}
    }
  return ok;
}

int FileIsWritable(const char *FileName)
{
struct stat a_stat;
if (stat(FileName, &a_stat) == 0) /* success ! */
  {
    return (a_stat.st_mode & (S_IWUSR|S_IWGRP|S_IWOTH));
  }
return 0;
}

/* borrowed in ximtool/fitsio.c */
int  IsFits (const string &Name)
{
        FILE *fp;
	int value = 0;
        char keyw[12];

        if ((fp = fopen (Name.c_str(), "r"))) {
            fscanf (fp, "%[^= ]=", keyw);
            if( (strcmp ("SIMPLE", keyw) == 0 /* && val == 'T' */) ||
		(strcmp ("XTENSION", keyw) == 0) )
                value = 1;
            fclose (fp);
        }
        return value;
}



#include <sys/types.h>
#include <dirent.h>



/* open a directory and returns the file names it contains */
int DirectoryContents(const char *DirName, StringList &FileNames)
{
DIR *this_dir = opendir(DirName);
FileNames.clear();
if (!this_dir) return 0;
struct dirent* entry;
while ((entry=readdir(this_dir)))
  {
  if (!entry->d_name) continue;
  if (strcmp(entry->d_name,".") == 0 || strcmp(entry->d_name,"..") == 0) continue;
  FileNames.push_back(entry->d_name);
  }
FileNames.sort();
return 1;
}


#ifdef I_DO_NOT_KNOW_HOW_TO_DO
#include <sys/stat.h>
#include <unistd.h>

int IsLink(const string &FileName)
{
struct stat stat_struct;
  /* lstat return 0 on success */
return (lstat(FileName.c_str(), &stat_struct) ==0  && 
      S_ISLNK(stat_struct.st_mode));

}

#endif /* I_DO_NOT_KNOW_HOW_TO_DO */


int IsDirectory(const string &FileName)
{
struct stat stat_struct;
  /* stat return 0 on success */
if (stat(FileName.c_str(), &stat_struct) != 0) return 0;
return S_ISDIR(stat_struct.st_mode);
}

//! return true for filename ending by ".Z" or ".gz"
bool IsZipped(const string &FileName)
{
  const char*p = strrchr(FileName.c_str(),'.');
  return (strstr(p,".Z") == p || strstr(p,".gz") == p );
}


#include <fstream>

static int shell_command_output(const string &command, StringList &output)
{
  string total_command;
  char tmpname[128];
  
  //char tmpname[128] = "/tmp/fileXXXXXX";
  //mkstemp(tmpname);
  //char toto[20] = "/fileXXXXXX";
  //char tata[20] = P_tmpdir;
  //char *tmpname = strcat(tata,toto);

  // It's the good way to use mkstemp but the temp directory is hard coded.
  // If you absolutely want to avoid warnings you should use these lines
  // Don't commit this, POR FAVOR!
  // tmpnam is the only function to create temp files 
  // using the correct temporary directory and to be POSIX.
  tmpnam(tmpname);
  // READ THE COMMENT BEFORE CHANGING THE ABOVE LINE

  if (tmpname[0] == '\0') // sometimes, I cannot get a temporary file name 
    strcat(tmpname,"/tmp/toto");
  total_command = command +" > "+string(tmpname);
  output.clear();
  if (system(total_command.c_str()) != 0) 
    {
      cerr << " could not execute : " << total_command << endl;
      return 0;
    }
  ifstream a_file(tmpname);
  while (!a_file.eof())
    {
      string a_line;
      a_file >> a_line;
      if (a_line.length() == 0) continue;
      output.push_back(a_line);
  }
  remove(tmpname);
  return 1;
}


int ExpandPath(const char *a_path, StringList &Expansions)
{
if (!a_path) return 0;
shell_command_output(" echo " + string(a_path), Expansions);
return Expansions.size();
}


StringList SplitString(const string &String, const char Sep)
{
StringList Result;
const char* next;
const char* start = String.c_str();
if (*start == Sep) ++start;
while ((next = strchr(start,Sep)))
   {
     char cp[256];
     strncpy(cp, start , next-start);
     cp[next-start] = '\0';
     Result.push_back(string(cp));
     start = next+1;
   }
Result.push_back(string(start));
if (getenv("SPLIT_DEBUG"))
  {
    for (StringCIterator p = Result.begin(); p != Result.end(); ++p) cout << *p << ' '; cout << endl;
  }
return Result;
}


string StandardPath(const string &Path)
  /*  Checks that a path does not contain things improperly handled 
      by RelativeLink or correct what can be corrected. */
{
string standard_path = Path;
//remove leading ./
 if (standard_path.find("./") == 0) standard_path.erase(0,2);

if (standard_path == ".")
  standard_path = getenv("PWD");
else if (strstr(standard_path.c_str(),".."))
  {
  string complete_path = FullFileName(standard_path);
  StringList dirs = SplitString(complete_path,'/');
  for (StringIterator p = dirs.begin(); p != dirs.end() ; )
    {
    if (p->length() == 0) { p = dirs.erase(p); if (p!=dirs.begin()) --p; continue;} 
    StringIterator next = p;
    ++next;
    if (next == dirs.end()) break;
    if (*p != ".." && *next == "..") // found a/..
      {
	p = dirs.erase(p); // erase a
	p = dirs.erase(p);   // erase ..
        if (p != dirs.end()) --p; // in case you had b/a/../.., restart at b.
      }
    else ++p;
    }
  standard_path = "";
  for (StringCIterator p = dirs.begin(); p != dirs.end() ; ++p)
    {
    
    standard_path += "/"+(*p);
    }
  }
return standard_path;
}


/***************************************************************************/
char* StringDuplicate (
 const char* This
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  char* string;
  int   length;
/*.........................................................................*/
  if(This==NULL)   return NULL;
  length           = strlen(This);
  string           = (char*)malloc((length+1)*sizeof(char));
  if(string==NULL) return NULL;
  string[length]   = '\0';
  return           strcpy(string,This);
}


/***************************************************************************/
int StringMatchPattern (
 const char* This   
,const char* a_pattern 
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
   int count;
  int      wcount;
  int      lpattern,lstring;
  char*    string;
  char*    token;
  const char* token2;
  int      ok,match;
  int      wnumber;
#define MAX_WORDS 100
  char*    words[MAX_WORDS]; 
/*.........................................................................*/
  if ( (a_pattern==NULL) && (This==NULL) ) return 1;
  if ( (a_pattern==NULL) && (This!=NULL) ) return 0;
  if ( (a_pattern!=NULL) && (This==NULL) ) return 0;
  lpattern  = strlen(a_pattern);
  lstring   = strlen(This);
  if ((lpattern==0)&&(lstring==0)) return 1;
  if ((lpattern==0)&&(lstring!=0)) return 0;
  if ((lpattern!=0)&&(lstring==0)) return 0;
/* pattern is * */
  if(strcmp(a_pattern,"*")==0) return 1;
  wcount = 0;
  for(count=0;count<lpattern;count++) {if(a_pattern[count]=='*') wcount++;}
/* no wildcard */
  if(wcount==0)
    {
      return (strcmp(a_pattern,This)==0 ? 1 : 0 );
    }

/* complex pattern */
  token    = string = StringDuplicate(a_pattern);
  wnumber  = 0;
  while(1)
    { char* pos;
      pos   = strstr (token,"*");
      if(pos!=NULL)
        {
          *pos = '\0';
          if(*token!='\0') 
            {
              if(wnumber>=MAX_WORDS) {fprintf(stderr, "StringMatchPattern : overflow of MAX_WORDS\n");} 
              else                   {words[wnumber] = token;wnumber++;}
            }
        token = pos + 1;
        }
      else /*last word*/
        {
          if(*token!='\0') 
            {
              if(wnumber>=MAX_WORDS) {fprintf(stderr, "StringMatchPattern : overflow of MAX_WORDS\n");} 
              else                   {words[wnumber] = token;wnumber++;}
            }
          break;
        }
    }
/* check that at least one word is not empty */
  ok = 0;
  for(count=0;count<wnumber;count++)
    { 
      if(*(words[count])!='\0') {ok = 1;break;}
    }
  if(ok==0) {free(string);return 1;} /* only wildcards */

/* loop on words */
  match    = 1;
  token2    = This;
  for(count=0;count<wnumber;count++)
    { int   lword;
      lword = strlen(words[count]); 
      if(lword>0) 
        { 
	  char* pos;
	  if(count==0)
	    {
	      if(a_pattern[0]!='*') /*Begin of pattern (words[0]) and This must match.*/
		{
		  if(strncmp(token2,words[count],lword)!=0) 
		    {
		      match = 0; /*Different.*/
		      break;      
		    }
		  token2 = token2 + lword;
		  continue;
		}
	    }
	  pos                = const_cast<char *>(strstr (token2,words[count]));
	  if(pos==NULL)      {match=0;break;}
	  if((count==(wnumber-1)) && (a_pattern[lpattern-1]!='*') ) /*Last word.*/
	    {
	      if(strcmp(This+lstring-lword,words[count])!=0) match = 0; /*Compare last word and end of This.*/
	      break;
	    }
	  else
	    {
	      token2 = pos + lword;
	    }
        }
    }
  free(string);
/*printf ("debug:match:%d:%s|%s|\n",match,This,a_pattern);*/
  return        match;
}
/***************************************************************************/
string StringToUpper(const string Source)
  //Converts a string to uppercase
{
  int n = Source.length();
  char* dest = new char[n+1];
  for(int i=0; i <= n  ;i++) dest[i] = toupper(Source[i]); 
  return string(dest);
}

string StringToLower(const string Source)
//Converts a string to lowercase
{
  int n = Source.length();
  char* dest = new char[n+1];
  for(int i=0; i <= n  ;i++) dest[i] = tolower(Source[i]); 
  return string(dest);
}

string AddPrefix(const string &Prefix, const string &FileName)
{
  return DirName(FileName)+"/"+Prefix+BaseName(FileName);
}

int RemovePattern(string &Source, const string &aPattern)
{
  // cout << " Removing '" << aPattern << "' from '"<< Source << "'\n";
  int iter = 0;
  string::size_type pos = Source.find(aPattern);
  while (pos != Source.npos)
    {
      Source.erase(pos, aPattern.length()); 
      iter++;
      pos = Source.find(aPattern);
    }
  return iter;
}

static bool is_space(char c) { return isspace(c) ; }

string DeleteWhiteSpaces(const string &Source)
{
  string dest;
  remove_copy_if(Source.begin(), Source.end(), back_inserter(dest), is_space);
  return dest;
}

void DecomposeString(vector<string> &SubStrings, const string &Source, 
		     const char *tokens)
{
  string::size_type start = 0;
  string::size_type tokpos = Source.find_first_of(tokens,start);

  while (tokpos != Source.npos) 
    {
      string sub = Source.substr(start,tokpos-start);
      RemovePattern(sub," ");
      RemovePattern(sub,"\t");
      if (sub.length() > 0) SubStrings.push_back(sub);
      start = tokpos+1;
      tokpos = Source.find_first_of(tokens,start);
    }
  // get the last element
  string sublast = Source.substr(start, Source.length());
  RemovePattern(sublast," ");
  RemovePattern(sublast,"\t");
  if (sublast.length() > 0) SubStrings.push_back(sublast);

  // for (unsigned i=0; i<SubStrings.size(); ++i) 
  //   cout << i << "='"<< SubStrings[i] << "'"<<endl;
}

std::string SubstituteExtension(const std::string &Original, 
				const std::string &NewExtension)
{
  unsigned int dotpos = dotposition(Original);
  if (dotpos > Original.size()) return Original;
  return Original.substr(0,dotpos)+NewExtension;
}


std::string SubstitutePattern(const std::string &InputString, 
			      const std::string &Pattern,
			      const std::string &Substitution)
{
  const char *pinput = InputString.c_str();
  int lpat = Pattern.size();
  int lpatout = Substitution.size();
  const char *where;

  int sizein = InputString.size();
  int npatt = 0 ;
  const char *pinputn = InputString.c_str();
  const char *wheren;
  while ((wheren = strstr(pinputn,Pattern.c_str())))
    {
      pinputn = wheren+lpat;
      npatt++;      
    }
  //  cerr << "Input String Length : " << sizein << endl ;
  //  cerr << "Pattern occurrence : " << npatt << endl ;
  int sizeout = sizein + npatt * (lpatout - lpat + 10 );
  //  cerr << "Estimated Substituted String Length : " << sizeout << endl ;

  char *output = new char[sizeout];
  char *pout = output;
  while ((where = strstr(pinput,Pattern.c_str())))
    {
      strncpy(pout,pinput,where-pinput);
      pout+= where-pinput;
      strcpy(pout, Substitution.c_str());
      pout += Substitution.size();
      pinput = where+lpat;
    }
  strcpy(pout, pinput);// copy the remainder
  std::string toto(output);
  delete [] output;
  return toto;
}

