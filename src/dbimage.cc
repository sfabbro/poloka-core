#include <vector>
#include <iostream>
#include <cstdio>
#include <cstdlib> /* for atoi */
#include <unistd.h> // for link
#include "fileutils.h"
#include "dbimage.h"


struct Path
{
string path;
string symbolicName;
int date;
  Path(const string &Path, const string &SPath, const int Date) : 
           path(Path), symbolicName(SPath), date(Date) {};
};


typedef list<Path>                 PathList;
typedef list<Path>::iterator       PathIterator;
typedef list<Path>::const_iterator PathCIterator;


class DbConfigFile {
public :
  DbConfigFile(const char* FileName);
  void AddImagePath(const char *a_path, const char *a_path_name);
  /* given an image name, returns its location 
     (as given in the config file ) */

  const Path* GetImagePath(const string &ImageName) const;
  friend ostream& operator << (ostream& stream, const DbConfigFile & config);
  void dump(ostream &stream = cout) const { stream << *this;}
  PathList GetImagePathes(const string &PathName) const;
  string FileName() const { return fileName;};

private :
  string fileName;
  vector<Path> image_pathes;
  void add_simple_image_path(const string &a_path, const string &a_path_name);
};


static DbConfigFile* DbConfig();
static int DumpLevel = 2;

int DbConfigSetDumpLevel(const int level)
{
int old_dump_level = DumpLevel;
DumpLevel = level;
return old_dump_level;
}


string DbConfigFileName() { return DbConfig()->FileName();}

extern "C" { /* callable from the parser (in C) */

void DbConfigAddImagePath(const char * a_path, const char *a_path_name)
{
DbConfig()->AddImagePath(a_path,a_path_name);
}


	   }  /* end of extern C */


static string without_blanks(const string &a_string)
{
string result = a_string;
int j = 0;
for (unsigned int i=0; i< a_string.length() ; ++i)
  {
  if (a_string[i] == ' ') continue;
  result[j++] = a_string[i];
  }
return result;
}

static int extract_date(const char* a_path)
{
char dirname[512];
char upper_dirname[512];
int date = 0;
strcpy(dirname, a_path);

 do {
    char basename[512];
    BaseName(dirname,basename);
    date = atoi(basename);
    if (date > 19950000) break;
    DirName(dirname, upper_dirname);
    if (strcmp(dirname,upper_dirname) == 0) return 0;
    strcpy(dirname, upper_dirname);
 } while (strlen(dirname) > 1);
return date;
}



#include <stdarg.h>

static void FatalError(const char * Message,...)
{
cerr << " ******** a fatal error occured in the image database handling : *****" << endl;
va_list args;
va_start (args,Message);
vfprintf (stderr,Message,args);
va_end(args);
cerr << " ******* the job stops here *********** " << endl;
exit(2);
}


static char* locate_config_file()
{char *fileName;
// first environment variable
if ((fileName = getenv("DBCONFIG")))  return fileName;

// second .dbconfig in PWD : built the complete filename for printout purposes
StringList expansion;
static char filename[256];

if (ExpandPath("$PWD/.dbconfig", expansion))
  {
  strcpy(filename,(*expansion.begin()).c_str());
  if (FileExists(filename)) return filename;
  }
// third .dbconfig in home directory : 
if (ExpandPath("$HOME/.dbconfig", expansion))
  {
  strcpy(filename,(*expansion.begin()).c_str());
  if (FileExists(filename)) return filename;
  }
FatalError(" DbInit :: No config file found \n");



return NULL;
}




void DbConfigDump(ostream &stream)
{
DbConfig()->dump(stream);
}


DbConfigFile::DbConfigFile(const char *FileName)
{
fileName = FileName;
}


static DbConfigFile* DbConfig()
{
static DbConfigFile *the_config = NULL;
if (!the_config) 
  {
  char *db_config_file = locate_config_file();

  the_config =  new DbConfigFile(locate_config_file());
  /* there is a trick here : DbConfigFileParse will call DbConfig()
     to get the pointer on the unique instance of DbConfig, and 
     the returned value will be correct since the_config is assigned */
  if (!DbConfigFileParse(db_config_file)) 
    FatalError(" cannot parse the config file %s\n", db_config_file);
  }
return the_config;
}

void DbInit()
{
DbConfig();
}


void DbConfigFile::add_simple_image_path(const string &a_path, const string &a_path_name)
{
string path = without_blanks(a_path);
int date = extract_date(path.c_str());
if (!FileExists(a_path.c_str()))
  {
    cerr << " the image path " << path << "  does not exist " << endl;
    return;
  }
// if (date == 0) FatalError(" cannot isolate a date (yyyymmdd) in the path %s\n", path.c_str());
path = AddSlash(path);
image_pathes.push_back(Path(path,a_path_name,date));
}

void  DbConfigFile::AddImagePath(const char *a_path, const char *a_path_name)
{
  /* here we should  translate a_path if it contains a $ or a * (using ls) */
if (strchr(a_path,'*') || strchr(a_path,'$') || a_path[0] == '~' )
  {
  StringList pathes;
  ExpandPath(a_path, pathes);
  for (StringIterator pi = pathes.begin(); pi != pathes.end(); ) 
        if (!FileExists((*pi).c_str())) pi = pathes.erase(pi); else ++pi;
  if (pathes.size() == 0) 
    {
      cerr << " cannot expand " << a_path << endl;
      return;
    }
  for (StringCIterator pi = pathes.begin(); pi != pathes.end(); ++pi ) add_simple_image_path(*pi, a_path_name); 
  }
else add_simple_image_path(a_path, a_path_name);
}


static int file_is_in_path(const string& ImageName, const string &Path)
{
string the_image_directory_name = Path+ImageName;
return (FileExists(the_image_directory_name.c_str()) && IsDirectory(the_image_directory_name.c_str()) );
}

PathList DbConfigFile::GetImagePathes(const string &PathName) const
{
PathList out;
for (unsigned int i=0; i<image_pathes.size(); ++i)
  {
  if (StringMatchPattern(image_pathes[i].symbolicName.c_str(),PathName.c_str()))
    {
    out.push_back(image_pathes[i]);
    }
  }
return out;
}


const Path* DbConfigFile::GetImagePath(const string &ImageName) const
{
for (unsigned int i = 0; i< image_pathes.size(); ++i)
  {
  if (file_is_in_path(ImageName, image_pathes[i].path))
    {
    return &(image_pathes[i]);
    }
  }
return 0;
}

ostream& operator << (ostream& stream, const DbConfigFile & config)
{
stream << "Db configuration read from :" << config.fileName << endl;
stream << "Image Pathes : " << endl;
for (unsigned int i = 0; i < config.image_pathes.size(); i++)
stream << " symbolic tag : " << config.image_pathes[i].symbolicName << " , actual path : " << config.image_pathes[i].path << endl;
return stream;
} 


/**************************** implementation of DbImage ********************************/

#define DIR_IS_IMAGE ".dbstuff"


static string db_image_tag(const string &a_directory)
{
return AddSlash(a_directory) + DIR_IS_IMAGE;
}

bool is_image(const string& a_directory)
{
return (FileExists(db_image_tag(a_directory)));
}


DbImage::DbImage(const string &ImageName)
{
  imageName = without_blanks(ImageName);
  init_from_name();
}


DbImage::DbImage(const char *ImageName)
{
  imageName = without_blanks(ImageName);
  init_from_name();
}


void DbImage::init_from_name()
{
  if (imageName[0] == '/') /* given as an absolute path : do not need the config. */
    {
      directory = AddSlash(imageName);
      imageName = BaseName(imageName); 
      if (is_image(directory))
	{
	  path = new Path(directory, imageName, 0); /* build a path, but it will not be used in fact */
	}
      else path = 0;
    }
  else /* use the db config to locate the image */
    {
      path = DbConfig()->GetImagePath(imageName);
      if (path) directory = AddSlash(AddSlash(path->path) + imageName);
      if (!is_image(directory)) path = 0; /* kills the image as a valid image */
    }
  saveEverythingElse = false; // this image already exists
  //  if (!FileExists(EverythingElseFileName())) saveEverythingElse = true;
}




bool DbImage::create(const string &ActualPath) // the actual creator
{
  if (!FileExists(ActualPath)) return false;
  string dir = AddSlash(ActualPath)+Name();
  directory = AddSlash(dir);
  // create the actual image directory
  if (!MKDir(dir.c_str())) return false;
  // create the file saying that this directory is a DbImage (.dbstuff up to now)
  string tag =  db_image_tag(dir);
  //cerr << tag << endl;
  fclose(fopen(tag.c_str(), "w"));
  // get the pointer to the path (which is used at least to tell if the DbImage is OK)
  path = DbConfig()->GetImagePath(Name());
  if (!path) path = new Path(ActualPath,"adhoc",0);
  
  // check that when retrieved, the image is the one just created
  DbImage db2(Name());
  db2.saveEverythingElse = false; // the one we are creating now should save its own data (which may not exists yet)
#if 0
  // very uninteresting
  if (! (*this == db2))
    {
      cerr << " when creating DbImage " << Name() 
	   << ", we see that the newborn image is not " << endl 
	   << " the one that would be retreived with the same name " << endl;
      cerr << " just created : directory = " << directory;
      cerr << " retrieved : directory = " << directory;
    }
#endif
  saveEverythingElse = true; // it is a new image, so we have to write stuff down.
  return true;
}


bool DbImage::Create(const string &Where)
{
  // if Where is a genuine directory :
if (FileExists(Where)) return create(Where);
//! if not, check that after shell interpretation, it is unambiguous:
 PathList path_list=DbConfig()->GetImagePathes(Where);
 int count = path_list.size();
 if (count == 0)
   {
     cerr << " DbImage::Create could not create image " << Name() << " in " << Where << endl;
     return false;
   }
 if (count > 1)
   {
     cerr << " DbImage::Create : ambiguous tag given : " << Where << " several pathes match using first one ! " << endl;
   }
 return create(path_list.begin()->path);
}


DbImage::DbImage(const string &ImageName, const Path* APath)
{
imageName = ImageName;
path = APath;
if (path) directory = AddSlash(AddSlash(path->path) + imageName);
if (!is_image(directory) ) path = 0; /* it kills it */
// cout << " debug directory : " << directory << endl;
}


bool DbImage::operator == (const DbImage &Right) const
{
  return ((imageName == Right.imageName) && 
	  (directory == Right.directory) &&
	  (path == Right.path));
}

static string GzTestedName(string Name)
{
  // By default, if nothing, returns gzname
  string gzname = Name+".gz"; 
  //  if (FileExists(gzname)) return gzname;
  if (FileExists(Name)) return Name;
  else return gzname;
}

string DbImage::FitsImageName(const DbImageKind Kind) const 
{
  // commented out for "virtual" image installation
  // if (!IsValid()) return string(""); 
if (Kind == Raw)        return directory + "raw.fits";
if (Kind == Elixir)     return ElixirName();
if (Kind == Calibrated) return directory + "calibrated.fits";
return "";
}

string DbImage::ElixirName() const
{
  return directory + "elixir.fits";
}

string DbImage::ImageCatalogName(const DbImageCatalogKind Kind) const 
{
if (!IsValid()) return string("");
if (Kind == SExtractor)         return directory + "se.list";
if (Kind == Fitted_for_seeing) return directory + "se.fit_seeing.list";
if (Kind == DaophotAls) return directory + "calibrated.als";
if (Kind == DaophotNst) return directory + "calibrated.nst";
if (Kind == DaophotPk) return directory + "calibrated.pk";
if (Kind == DaophotLst) return directory + "calibrated.lst";
if (Kind == DaophotNei) return directory + "calibrated.nei";
if (Kind == DaophotAp) return directory + "calibrated.ap";
return "";
}

string DbImage::ImagePsfName(const DbImagePsfKind Kind) const 
{
if (!IsValid()) return string("");
if (Kind == DaophotPsf) return directory + "calibrated.psf";
return "";
}

void DbImage::dump(ostream &stream)  const
{
if (DumpLevel >=2)
  stream << path->symbolicName << " : " << imageName << endl;
else if (DumpLevel == 1)
  stream << imageName << endl;
}

string DbImage::ImageBackName(const DbImageKind Kind) const
{
return FitsImageName(Kind);  
/* The code that writes image_backs adds some stuff 
to the file name to separate the average and rms images */
}

string DbImage::FitsFlatName() const
{
if (!IsValid()) return string("");
return directory + "flat.fits";
}

string DbImage::FitsBiasName() const
{
if (!IsValid()) return string("");
return directory + "bias.fits";
}

string DbImage::FitsDarkName() const
{
if (!IsValid()) return string("");
return directory + "dark.fits";
}

// the dead.fits files are not gzipped because we don't want to test the existence of the file in this routine as it is used by relative_ln to create the link to dead mask even if the file doesn't exit
string DbImage::FitsDeadName() const
{
if (!IsValid()) return string("");
return directory + "dead.fits";
}

// bad.fits is readen in input by sextractor which doesn't decompress
string DbImage::FitsBadName() const
{
if (!IsValid()) return string("");
return directory + "bad.fits";
}

string DbImage::FitsWeightName() const
{
if (!IsValid()) return string("");
return directory + "weight.fits";
}

string DbImage::FitsCosmicName() const
{
if (!IsValid()) return string("");
return GzTestedName(directory + "cosmic.fits");
}


string DbImage::FitsSatelliteName() const
{
if (!IsValid()) return string("");
return GzTestedName(directory + "satellite.fits");
}

string DbImage::FitsFringeName() const
{
if (!IsValid()) return string("");
return directory + "fringe.fits";
}

string DbImage::FitsBackName() const
{
if (!IsValid()) {return string("");}
return directory + "back.fits";
}

string DbImage::FitsMiniBackName() const
{
if (!IsValid()) return string("");
return directory + "miniback.fits";
}


string DbImage::FitsSaturName() const
{
if (!IsValid()) return string("");
return GzTestedName(directory + "satur.fits");
}

string DbImage::ImageMatchUsnoName() const
{
if (!IsValid()) return string("");
return directory + "match_usno.dat";
}

string DbImage::GetFileName(const char* WhichFile) const
{
  if (strcmp(WhichFile,"raw")==0)    return FitsImageName(Raw);
  if (strcmp(WhichFile,"cal")==0)    return FitsImageName(Calibrated);
  if (strcmp(WhichFile,"flat")==0)   return FitsFlatName();
  if (strcmp(WhichFile,"sat")==0)    return FitsSaturName();
  if (strcmp(WhichFile,"dead")==0)   return FitsDeadName();
  if (strcmp(WhichFile,"back")==0)   return FitsBackName();
  if (strcmp(WhichFile,"cos")==0) return FitsCosmicName();
  if (strcmp(WhichFile,"satel")==0) {  return FitsSatelliteName();}
  if (strcmp(WhichFile,"miniback")==0)   return FitsMiniBackName();
  if (strcmp(WhichFile,"bias")==0)   return FitsBiasName();
  if (strcmp(WhichFile,"cat")==0)    return ImageCatalogName(SExtractor);
  if (strcmp(WhichFile,"fitcat")==0) return ImageCatalogName(Fitted_for_seeing);
  if (strcmp(WhichFile,"fringe")==0) return FitsFringeName();
  if (strcmp(WhichFile,"psf")==0) return ImagePsfName(DaophotPsf);
  if (strcmp(WhichFile,"dark")==0)    return FitsDarkName();
  if (strcmp(WhichFile,"dir")==0) return directory;
  if (strcmp(WhichFile,"usno")==0) return ImageMatchUsnoName();
  if (strcmp(WhichFile,"weight")==0) return FitsWeightName();
  return "";
}

DbImage::~DbImage()
{
  if (saveEverythingElse) writeEverythingElse();
}

#ifdef USE_ROOT
#include "TFile.h"
#endif /* USE_ROOT */

bool DbImage::writeEverythingElse()
{
  if (!saveEverythingElse) return false;
  string fileName = EverythingElseFileName();
  //  cout << " writing \"EverythingElse\" file : " << fileName << endl;
#ifdef USE_ROOT
  TFile tFile(fileName.c_str(),"RECREATE");
  if (!tFile.IsOpen())
    {
      cerr << " could not open " << fileName << endl;
      return false;
    }
  Write();
  // tFile.Close() called by destructor
  saveEverythingElse = false;
  return true;
#else
  //  cerr << " no way to write " << fileName << " without defining USE_ROOT " << endl;
  return false;
#endif /* USE_ROOT */
}

#ifdef USE_ROOT
ClassImp(DbImage);
#endif /* USE_ROOT */

/***************************   DbImageList ****************************/

DbImageList::DbImageList(const list<string> names)
{
  for (list<string>::const_iterator it= names.begin(); it != names.end() ; ++it)
    {
      string name = *it;
      DbImage dbimage(name);
      if (dbimage.IsValid())
	{
	  this->push_back(dbimage);
	  continue;
	}
      else
	{
	  if (Collect(name.c_str()) == 0)
	    {
	      cerr << " could not interpret " << name << " neither as a DbImage nor as a symbolic path " << endl;
	    }
	}
    }
}

DbImageList::DbImageList(const string &PathName)
{
Collect(PathName.c_str());
}


DbImageList::DbImageList(const char * PathName)
{
Collect(PathName);
}

int DbImageList::Collect(const char * PathName)
{
PathList path_list = DbConfig()->GetImagePathes(PathName);
unsigned int old_size = this->size();
if (path_list.size() == 0)
  {  /* interpret it directly as a genuine file path, I am not sure it is 
         a very good idea */
  Path a_path(PathName,"toto",extract_date(PathName));
  if (a_path.date != 0)   path_list.push_back(a_path); 
  }
for (PathCIterator ipath = path_list.begin(); ipath != path_list.end(); ++ipath)
  {
  StringList file_list;
  DirectoryContents((*ipath).path.c_str(), file_list);
  for (StringCIterator ifile = file_list.begin(); ifile != file_list.end(); ++ifile)
    {
    DbImage an_image(*ifile, &(*ipath));
    if (an_image.IsValid())     push_back(an_image);
    }
  }
return int(this->size() - old_size);
}



void DbImageList::dump(ostream &stream) const
{
for (DbImageCIterator ii = begin(); ii != end(); ++ii)
  {
  stream << *ii;
  }
}

void DbImageList::FilterByDate(const int a_date)
{
for (DbImageIterator dbi = begin(); dbi != end(); )
  {
  if ((*dbi).path->date != a_date)
    {
    dbi = erase(dbi);
    }
  else ++dbi;
  }
}

#ifdef STORAGE
void DbImageList::SelectedDump(const int WhichInfos) const 
{
for (DbImageCIterator ii = begin(); ii != end(); ++ii)
  {
  (*ii).SelectedDump(WhichInfos);
  }
}
#endif

/************* image installation *******************/
/* assumes that a_path comes from StandardPath (to make sure that RelativeLink works properly) */

int InstallImage(const char *a_path, const char *a_file, DbImageKind kind)
{
    /* figure out the symbolic name */
char symb_name[256];
BaseName(a_file,symb_name); /* strip leading path */
char *point = strrchr(symb_name,'.');
if (point) *point = '\0'; /* remove everything after the last dot */

/* build the directory name */
string base_dir = AddSlash(a_path) + string(symb_name);

DbImage dbimage(symb_name);
bool ok_dir = dbimage.Create(a_path); // create the directory
string fits_name = dbimage.FitsImageName(kind);
// cout << " debug fits_name : " << fits_name << endl;

/* ready to setup the link */
string link_value = RelativeLink(StandardPath(a_file).c_str(), 
				 fits_name.c_str());
 if (!ok_dir || symlink(link_value.c_str(), fits_name.c_str()) == 0) 
   {
     cout << " Installing " << symb_name << endl;
     return 1;
   }

#if 0 
string command = " ln -fs " + link_value + " " + fits_name;
if (getenv("INSTALL_DEBUG")) cerr << " debug : " << command << endl;
// no point in passing the "ln" command if the directory could not be created
// however, output the name of the would be dbimage (on stdout). 
// many error messages go to stderr.
if (!ok_dir || system(command.c_str()) == 0) 
  {
    cout << symb_name << endl;
    return 1;
  }
#endif
return 0;
}


int AssignInfo(const DbImage &Image, const string &FlatFitsFileName, const char *WhichInfo)
{
string flatName = Image.GetFileName(WhichInfo);
if (flatName.length() == 0)
  {
    cerr << " AssignInfo received " << WhichInfo << " as an info name" << endl;
    return 0;
  }

 if (FileExists(flatName)) unlink(flatName.c_str()); 

 // string link_value = RelativeLink(StandardPath(FlatFitsFileName).c_str(), flatName.c_str());
 cout << " Assigning " << flatName << " : " 
      << symlink(FlatFitsFileName.c_str(), flatName.c_str()) << endl;

#if 0
string command = "ln -fs " + link_value + " " + flatName;
 cout << " debug " << command << endl;
if (system(command.c_str()) == 0) return 1;
#endif

return 0;
}



#ifdef USE_ROOT
/*
RUN_ROOTCINT
LINKDEF_CONTENT : #pragma link C++ enum _DbImageKind;
LINKDEF_CONTENT : #pragma link C++ enum DbImageCatalogKind;
LINKDEF_CONTENT : #pragma link C++ enum _DbImagePsfKind;
LINKDEF_CONTENT : #pragma link C++ class DbImage;



*/


#include "root_dict/dbimagedict.cc"
#endif /* USE_ROOT */
