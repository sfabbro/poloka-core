#ifndef FILEUTILS__H
#define FILEUTILS__H

#include <string>
#include <vector>
#include <stringlist.h>



//! directory name
void DirName(const char *FullFileName, char *DirName);

//! directory name
string DirName(const string &FullFileName);

//! drops the directory
void BaseName(const char *FullFileName, char *Name);

//! drops the directory
string BaseName(const std::string &Name);

//! remove everything after the last period(.)
string CutExtension(const std::string &Name);

//!
int   FileExists(const char *FileName);
//!
int   FileExists(const std::string &FileName);
//!
int   FileIsWritable(const char *FileName);

//! does the file begin as a fits file?
int   IsFits (const std::string &fname);

//! returns filenames in this directory.
/*! Drops files beginning by . and .. */
int   DirectoryContents(const char *DirName, StringList &FileNames);

//! return true for filename ending by ".Z" or ".gz"
bool IsCompressed(const string &FileName);


//!
int   IsDirectory(const std::string &FileName);

#ifdef DOES_NOT_WPRK
//! Is the file a symbolic link? 
int   IsLink(const std::string &FileName);
#endif


//! mkdir shell command
int MKDir(const char *path, const bool Warn=true);

//! expands a file name (prepends cwd if missing)
std::string FullFileName(const string  &FileName);

//! used to setup 'relative' symbolic links
/*!  Relative symbolic links do not refer to the common part of the path.
     RelativeLink("/a/b/c/d", /a/b/c/e") returns ../d */
std::string RelativeLink(const char *RealFile, const char *LinkLocation);

//! same as above but issue the command to the system.
int MakeRelativeLink(const std::string &RealFile, const std::string &LinkLocation);

//! AddPrefix("save","dir/foo.c") returns "dir/savefoo.c"
std::string AddPrefix(const std::string &Prefix, const std::string &FileName);

//! adds a / at the end of path is none.
std::string AddSlash(const std::string &path);

//! 'standardize' a path (used by RelativeLink).
std::string StandardPath(const std::string &Path);

//! Expands a path by invoking the shell.
 int ExpandPath(const char *a_path, StringList &Expansions);

//! I am not sure it does more than strdup.
char* StringDuplicate ( const char* This);

//! tests if a string matches a pattern with * in it.
int StringMatchPattern ( const char* This   ,const char* a_pattern );

//! translates to upper case.
std::string StringToUpper(const std::string Source);

//! translate to lower case.
std::string StringToLower(const std::string Source);

//! removes a given pattern
int RemovePattern(std::string &Source, const std::string &aPattern);

//! Decompose string into substring separated by a char token
void DecomposeString(std::vector<std::string> &SubStrings, const string &Source, const char token=' ');

//! substitutes ".xxx" by "<newExtension>" at the end of "Original"
/*! to change toto.fits to toto.head call 
  SubstituteExtension("toto.fits", ".head");
*/
std::string SubstituteExtension(std::string Original, std::string NewExtension);


#endif /* FILEUTILS__H */
