#include <map>
#include <vector>
#include <iostream>
#include <cstdio>
#include <cstdlib> // for atoi
#include <cstring> // for strcmp
#include <cstdarg>
#include <unistd.h> // for link

#include <poloka/fileutils.h>
#include <poloka/stringlist.h>
#include <poloka/dbimage.h>
#include <poloka/dbconfigexception.h>

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
  void AddCatalogPath(const char * a_path);
  /* given an image name, returns its location 
     (as given in the config file ) */

  const Path* GetImagePath(const string &ImageName) const;
  //  friend ostream& operator << (ostream& stream, const DbConfigFile & config);
  void dump(ostream &stream = cout) const;
  PathList GetImagePathes(const string &PathName) const;
  string FileName() const { return fileName;};
  string FindCatalog(const string &FileName, const bool Warn) const;
  string NewImageName(const string &TypeName, const string &DirName, 
		      const string &BaseName) const;

  void AddNewImageName(const string &TypeName,
		       const string &BaseName,
		       const string &NewName);



private :
  string fileName;
  vector<Path> image_pathes;
  void add_simple_image_path(const string &a_path, const string &a_path_name);
  vector<string> catalog_pathes;

  typedef map<string,string> MapSS;
  typedef map<string,string>::const_iterator MapSSCIterator;
  void init_image_names();

  MapSS NewImageNames;
/* we would need a map with 2 indices : the TypeName ("ReducedImage", "ImageSum"), and the filename ("calibrated", "weight"). To do this in a simple way ,we just address a map with the sum of both strings
 */
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

void DbConfigAddCatalogPath(const char * a_path)
{
DbConfig()->AddCatalogPath(a_path);
}

void DbConfigAddNewImageNames(const char *TypeName,
			     const char *NewNames)
{
  vector<string> names;
  DecomposeString(names, NewNames);
  for (unsigned int i=0; i < names.size(); ++i)
    {
      string basename = CutExtension(names[i]);
      DbConfig()->AddNewImageName(TypeName, basename, names[i]);
    }
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



static void FatalError(const char * Message,...)
{
cerr << " ******** a fatal error occured in the image database handling : *****" << endl;
 char mess[8192];
va_list args;
va_start (args,Message);
 vsnprintf (mess,8192,Message,args);
va_end(args);
 cerr << mess << endl;
 throw (DbConfigException(mess));
}


static char* locate_config_file()
{char *fileName;
// first environment variable
if ((fileName = getenv("POLOKA_DB_CONFIG")))  return fileName;

// second .dbconfig in PWD : built the complete filename for printout purposes
StringList expansion;
static char filename[256];

if (ExpandPath("$PWD/.poloka/dbconfig", expansion))
  {
  strcpy(filename,(*expansion.begin()).c_str());
  if (FileExists(filename)) return filename;
  }
// third .dbconfig in home directory : 
if (ExpandPath("$HOME/.poloka/dbconfig", expansion))
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

/********************* DbConfigFile ***********************/

DbConfigFile::DbConfigFile(const char *FileName)
{
  fileName = FileName;
  init_image_names(); 
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


void  DbConfigFile::AddCatalogPath(const char *a_path)
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
  for (StringCIterator pi = pathes.begin(); pi != pathes.end(); ++pi ) 
    catalog_pathes.push_back(AddSlash(string(*pi)));
  }
else catalog_pathes.push_back(AddSlash(string(a_path)));
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

void DbConfigFile::dump(ostream& stream) const
{
  stream << "Db configuration read from :" << fileName << endl;
  // dbimages
  stream << "#### Image Pathes : " << endl;
  for (unsigned int i = 0; i < image_pathes.size(); i++)
    stream << " symbolic tag : " << image_pathes[i].symbolicName 
	   << " , actual path : " << image_pathes[i].path << endl;
  // astrometric catalogs
  if (catalog_pathes.size())
    {
      stream << "##### (Reference)Catalog Pathes : " << endl;
      for (unsigned int i=0; i<catalog_pathes.size(); ++i)
	stream << "    " << catalog_pathes[i] << endl;
    }
      // compression handling
  cout << "#### ImageNames used when creating new images (compression handling):" << endl;
  for (  MapSSCIterator i = NewImageNames.begin();
	 i!= NewImageNames.end(); ++i)
    stream << i->first << " : " << i->second << endl;
} 


void DbConfigExample()
{
  cout  
<< "# this is a comment" << endl
<< "ImagePath" << endl
<< "{" << endl
<< "# 'here' is where non existing images are created. You'd better define it" << endl
<< "here : ." << endl
<< "cfht99 : /snovad15/cfht99/1999* ,   /somewhere_else/cfht99/" << endl
<< "vlt99 : /snovad1/vlt99/1999*" << endl
<< "newstuff : /snovad8/wiyn99" << endl
<< "}" << endl
<< "" << endl
<< "# where to find astrometric catalogs (non-USNO catalogs usually)" << endl
<< "CatalogPath" << endl
<< "{" << endl
<< "  ." << endl
<< "  /data/my_catalogs" << endl
<< "  /data/catalogs/D*" << endl
<< "}" << endl
<< endl
<< "#what are the image names for the various DbImage derived classes" << endl
<< "# .fits : regular fits image" << endl
<< "# .fits.fz : rice compressed fits image" << endl
<< "# .fits.gz gzip compression" << endl
<< "# mind the spaces around '{' and '}'" << endl
<< "ImageNames" << endl
<< "{" << endl
<< "#default if not overwritten:" << endl
<< "  { satur.fits.gz }" << endl
<< "#for the ImageSum class" << endl
<< "  ImageSum { calibrated.fits }" << endl
<< "  TransformedImage { calibrated.fits satur.fits.gz }" << endl
 << "}" << endl
 ;
}



string DbConfigFile::FindCatalog(const string &FileName, const bool Throw) const
{
  if (FileName == "") return "";
  if (FileName[0] == '/')
    {
      if (FileExists(FileName)) return FileName;
    }
  else for (unsigned i = 0; i < catalog_pathes.size(); ++i)
    {
      std::string cand = catalog_pathes[i]+FileName;
      if (FileExists(cand)) return cand;
    }
  if (Throw) 
    FatalError("DbConfigFindCatalog :  cannot locate %s",FileName.c_str());
  return "";
}

// locator for catalogs
string DbConfigFindCatalog(const string &FileName, const bool Throw)
{
  return DbConfig()->FindCatalog(FileName, Throw);
}



/************** implementation of DbImage *************************/



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

#define DIR_IS_IMAGE ".dbstuff"


static string db_image_tag(const string &a_directory)  {
  return AddSlash(a_directory) + DIR_IS_IMAGE;
}

static bool is_image(const string& a_directory) {
  return (FileExists(db_image_tag(a_directory)));
}

bool DbImage::IsValid() const 
{ 
  return ((directory != "" ) && is_image(directory));
}


void DbImage::init_from_name()
{
  if (imageName[0] == '/') /* given as an absolute path : do not need the config. */
    {
      directory = AddSlash(imageName);
      imageName = BaseName(imageName); 
    }
  else /* use the db config to locate the image */
    {
       const Path* path = DbConfig()->GetImagePath(imageName);
      if (path) directory = AddSlash(AddSlash(path->path) + imageName);
    }
  imageName = BaseName(imageName);
  //if(!IsValid()) {
  //std::cerr << "ERROR in  DbImage::init_from_name cannot find " << imageName << endl;
  //}
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
  
  // check that when retrieved, the image is the one just created
  DbImage db2(Name());
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
  return true;
}


bool DbImage::Create(const string &Where)
{
  // if it is already there, don't change.
  if (IsValid()) return true;
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
 return (create(path_list.begin()->path) && StoreTypeName());
}


DbImage::DbImage(const string &ImageName, const Path* APath)
{
imageName = ImageName;
if (APath) directory = AddSlash(AddSlash(APath->path) + imageName);
}


bool DbImage::operator == (const DbImage &Right) const
{
  return ((imageName == Right.imageName) && 
	  (directory == Right.directory));
}





/* The two following routines have to do with reloading saved 
   DbImage's and inheriters. with objio, they appear to be
   useless. But the "type.name" files are used by snlsdb, so
   they have to be written anyway. 
*/

#define TYPE_FILE_NAME "type.name"


bool DbImage::StoreTypeName()
{
  string fileName = Dir()+TYPE_FILE_NAME;
  FILE *file = fopen(fileName.c_str(),"w");
  if (!file) 
    {
      cerr << " could not open " <<  fileName << endl;
      return false;
    }
  fprintf(file,"%s",TypeName().c_str());
  fclose(file);
  return true;
}
			       

string DbImage::StoredTypeName() const
{
  char name[256];
  string fileName = Dir()+TYPE_FILE_NAME;
  FILE *file = fopen(fileName.c_str(),"r");
  if (!file)
    {
      cerr << "cannot open in read mode " << fileName << endl;
      return "NoType";
    }
  fscanf(file,"%s",name);
  fclose(file);

  return string(name);
}



void DbConfigFile::AddNewImageName(const string &TypeName,
				   const string &BaseName,
				   const string &NewName)
{
  NewImageNames[TypeName+BaseName] = NewName;
}


//defaults, may be overwritten in dbconfig file:

void DbConfigFile::init_image_names()
{
  AddNewImageName("","raw","raw.fits.fz");
  AddNewImageName("","calibrated","calibrated.fits.fz");
  AddNewImageName("","weight","weight.fits.fz");
  AddNewImageName("","sub","sub.fits.fz");
  AddNewImageName("","subweight","sub.weight.fits.fz");
  AddNewImageName("","elixir","elixir.fits");
  AddNewImageName("","back","back.fits");
  AddNewImageName("","miniback","miniback.fits");

  AddNewImageName("","dead","dead.fits");
  AddNewImageName("","flat","flat.fits");
  AddNewImageName("","dark","dark.fits");
  AddNewImageName("","fringe","fringe.fits");
  AddNewImageName("","bias","bias.fits");

  AddNewImageName("","cosmic","cosmic.fits.gz");
  AddNewImageName("","segmentation","segmentation.fits.gz");
  AddNewImageName("","satellite","satellite.fits.gz");
  AddNewImageName("","satur","satur.fits.gz");

  // swarp writes in 32 bits, so no Rice compression
  AddNewImageName("SwarpStack","calibrated","calibrated.fits");  
  AddNewImageName("SwarpStack","weight","weight.fits");
  // no rice for satur, gz does it faster and just as good, without header pbs
}

/* returns the filename, from stored instructions, that come
   from  init_image_names and the parsing of dbconfig
*/
string DbConfigFile::NewImageName(const string &TypeName, 
				  const string &Dir, 
				  const string &BaseName) const
{
  //DEBUG
#ifdef DEBUG
  static bool called = false;
  if (!called)
    {
      called = true;
      cout << " image names contents " << endl;
      for (  MapSSCIterator i = NewImageNames.begin();
	     i!= NewImageNames.end(); ++i)
	cout << i->first << " : " << i->second << endl;
    }
#endif



  MapSSCIterator i = NewImageNames.find(TypeName+BaseName);
  if (i!= NewImageNames.end()) return Dir+i->second;
  // no, go to default
  i = NewImageNames.find(BaseName);
  if (i!= NewImageNames.end()) return Dir+i->second;
  return Dir+BaseName+".fits";
}

/* the routine that provides names for existing or new images.
   may become a DbImage member function 
*/
static string image_name(const string &TypeName, const string &Dir, 
			 const string &BaseName)
{
  string withoutExtension = Dir+BaseName;
  string fileName = withoutExtension+".fits";
  if (FileExists(fileName))  return fileName;
  fileName = withoutExtension+".fits.gz";
  if (FileExists(fileName))  return fileName;
  fileName = withoutExtension+".fits.fz";
  if (FileExists(fileName))  return fileName;
  fileName = withoutExtension+".head";
  if(FileExists(fileName)) return fileName;
  // the file does not exists yet. so generate its
  // name according to defaults
  return DbConfig()->NewImageName(TypeName,Dir,BaseName);


  // first, check if there is an entry that corresponds
  // to this ReducedImage type and this basename.
}


string DbImage::FitsImageName(const DbImageKind Kind) const 
{
  // commented out for "virtual" image installation
  // if (!IsValid()) return string(""); 
if (Kind == Raw)        return image_name(TypeName(), directory, "raw");
if (Kind == Elixir)     return image_name(TypeName(), directory, "elixir");
if (Kind == Calibrated) return image_name(TypeName(), directory, "calibrated");
if (Kind == Subtracted) return image_name(TypeName(), directory, "sub");
return "";
}


string DbImage::FitsWeightImageName(const DbImageKind Kind) const 
{
if (Kind == Calibrated) return image_name(TypeName(), directory, "weight");
if (Kind == Subtracted) return image_name(TypeName(), directory, "sub.weight");
return "";
}

string DbImage::ElixirName() const
{
  return image_name(TypeName(), directory, "elixir");
}

string DbImage::ImageCatalogName(const DbImageCatalogKind Kind) const 
{
if (!IsValid()) return string("");
if (Kind == SExtractor)         return directory + "se.list";
if (Kind == Subtraction)         return directory + "det.list";
return "";
}

string DbImage::AperCatalogName() const
{
  return directory+"aperse.list";
}

string DbImage::FixedAperCatalogName() const
{
  return directory+"fixed_aperse.list";
}


string DbImage::StarCatalogName() const
{
  return directory+"standalone_stars.list";
}

string DbImage::ImagePsfName(const DbImagePsfKind Kind) const 
{
if (!IsValid()) return string("");
if (Kind == DaophotPsf) return directory + "calibrated.psf";
return "";
}

void DbImage::dump(ostream &stream)  const
{
  if (DumpLevel >=2) {
    stream << directory << " : " << imageName;
    if(IsValid())
      stream << " is valid" << endl;
    else
      stream << " is not valid" << endl;
  } else if (DumpLevel == 1)
  stream << imageName << endl;
}

string DbImage::FitsSubName() const
{
  return image_name(TypeName(), directory, "sub");
}

string DbImage::FitsFlatName() const
{
  return image_name(TypeName(), directory, "flat");
}

string DbImage::FitsBiasName() const
{
  return image_name(TypeName(), directory, "bias");
}

string DbImage::FitsDarkName() const
{
  return image_name(TypeName(),directory,"dark");
}


string DbImage::FitsDeadName() const
{
  return image_name(TypeName(),directory,"dead");
}

// bad.fits is readen in input by sextractor which doesn't decompress
string DbImage::FitsBadName() const
{
  return image_name(TypeName(),directory,"bad");
}

string DbImage::FitsWeightName() const
{
  return image_name(TypeName(),directory,"weight");
}

string DbImage::FitsSubWeightName() const
{
  return image_name(TypeName(),directory,"sub.weight");
}

string DbImage::FitsCosmicName() const
{
  return image_name(TypeName(),directory,"cosmic");
}


string DbImage::FitsSegmentationName() const
{
  return image_name(TypeName(),directory,"segmentation");
}


string DbImage::FitsSatelliteName() const
{
  return image_name(TypeName(),directory,"satellite");
}

string DbImage::FitsFringeName() const
{
  return image_name(TypeName(),directory,"fringe");
}

string DbImage::FitsBackName() const
{
  return image_name(TypeName(),directory,"back");
}

string DbImage::FitsMiniBackName() const
{
  return image_name(TypeName(),directory,"miniback");
}


string DbImage::FitsSaturName() const
{
  return image_name(TypeName(),directory,"satur");
}

string DbImage::ImageMatchUsnoName() const
{
if (!IsValid()) return string("");
return directory + "match_usno.dat";
}

string DbImage::GetFileName(const char* WhichFile) const
{
  if (strcmp(WhichFile,"raw")==0)    return FitsImageName(Raw);
  if (strcmp(WhichFile,"elixir")==0)  return FitsImageName(Elixir);
  if (strcmp(WhichFile,"cal")==0)    return FitsImageName(Calibrated);
  if (strcmp(WhichFile,"flat")==0)   return FitsFlatName();
  if (strcmp(WhichFile,"sat")==0)    return FitsSaturName();
  if (strcmp(WhichFile,"dead")==0)   return FitsDeadName();
  if (strcmp(WhichFile,"back")==0)   return FitsBackName();
  if (strcmp(WhichFile,"cos")==0) return FitsCosmicName();
  if (strcmp(WhichFile,"seg")==0) return FitsSegmentationName();
  if (strcmp(WhichFile,"satel")==0) {  return FitsSatelliteName();}
  if (strcmp(WhichFile,"miniback")==0)   return FitsMiniBackName();
  if (strcmp(WhichFile,"bias")==0)   return FitsBiasName();
  if (strcmp(WhichFile,"cat")==0)    return ImageCatalogName(SExtractor);
  if (strcmp(WhichFile,"fringe")==0) return FitsFringeName();
  if (strcmp(WhichFile,"psf")==0) return ImagePsfName(PolokaPsf);
  if (strcmp(WhichFile,"dark")==0)    return FitsDarkName();
  // directory names in standard unix tools assume no trailing slashes
  if (strcmp(WhichFile,"dir")==0) return directory.substr(0, directory.size()-1);
  if (strcmp(WhichFile,"usno")==0) return ImageMatchUsnoName();
  if (strcmp(WhichFile,"weight")==0) return FitsWeightImageName(Calibrated);
  if (strcmp(WhichFile,"sub")==0)    return FitsImageName(Subtracted);
  if (strcmp(WhichFile,"subweight")==0) return FitsWeightImageName(Subtracted);
  return "";
}


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
  path_list.push_back(a_path); 
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

// void DbImageList::FilterByDate(const int a_date)
// {
// for (DbImageIterator dbi = begin(); dbi != end(); )
//   {
//     if ((*dbi).path->date != a_date)
//     {
//     dbi = erase(dbi);
//     }
//   else ++dbi;
//   }
// }

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
    return EXIT_FAILURE;
  }
 
 if (FileExists(flatName)) unlink(flatName.c_str()); 

 // string link_value = RelativeLink(StandardPath(FlatFitsFileName).c_str(), flatName.c_str());
 cout << " Assigning " << flatName << " : " 
      << symlink(FlatFitsFileName.c_str(), flatName.c_str()) << endl;

#if 0
string command = "ln -fs " + link_value + " " + flatName;
 cout << " debug " << command << endl;
if (system(command.c_str()) == 0) return EXIT_SUCCESS;
#endif

return EXIT_SUCCESS;
}


