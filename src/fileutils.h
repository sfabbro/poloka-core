#ifndef FILEUTILS__H
#define FILEUTILS__H

#include <vector>
#include <stringlist.h>

//! directory name
void DirName(const char *FullFileName, char *DirName);

//! directory name
string DirName(const string &FullFileName);

//! drops the directory
void BaseName(const char *FullFileName, char *Name);

//! drops the directory
string BaseName(const string &Name);

//! remove everything after the last period(.)
string CutExtension(const string &Name);

//!
int   FileExists(const char *FileName);
//!
int   FileExists(const string &FileName);
//!
int   FileIsWritable(const char *FileName);

//! does the file begin as a fits file?
int   IsFits (const string &fname);

//! returns filenames in this directory.
/*! Drops files beginning by . and .. */
int   DirectoryContents(const char *DirName, StringList &FileNames);

//!
int   IsDirectory(const string &FileName);

#ifdef DOES_NOT_WPRK
//! Is the file a symbolic link? 
int   IsLink(const string &FileName);
#endif

//! mkdir shell command
int MKDir(const char *path, const bool Warn=true);

//! expands a file name (prepends cwd if missing)
string FullFileName(const char *FileName);

//! used to setup 'relative' symbolic links
/*!  Relative symbolic links do not refer to the common part of the path.
     RelativeLink("/a/b/c/d", /a/b/c/e") returns ../d */
string RelativeLink(const char *RealFile, const char *LinkLocation);

//! same as above but issue the command to the system.
int MakeRelativeLink(const char*RealFile, const char *LinkLocation);

//! AddPrefix("save","dir/foo.c") returns "dir/savefoo.c"
string AddPrefix(const string &Prefix, const string &FileName);

//! adds a / at the end of path is none.
string AddSlash(const string &path);

//! 'standardize' a path (used by RelativeLink).
string StandardPath(const string &Path);

//! Expands a path by invoking the shell.
 int ExpandPath(const char *a_path, StringList &Expansions);

//! I am not sure it does more than strdup.
char* StringDuplicate ( const char* This);

//! tests if a string matches a pattern with * in it.
int StringMatchPattern ( const char* This   ,const char* a_pattern );

//! translates to upper case.
string StringToUpper(const string Source);

//! translate to lower case.
string StringToLower(const string Source);

//! removes a given pattern
int RemovePattern(string &Source, const string &aPattern);

//! Decompose string into substring separated by a char token
void DecomposeString(vector<string> &SubStrings, const string &Source, const char token=' ');

#endif /* FILEUTILS__H */
