// This may look like C code, but it is really -*- C++ -*-
#ifndef DBIMAGE__H
#define DBIMAGE__H


/*! \file dbimage.h
\brief documentation for the DbImage class and the \ref dbconfig file, and more
generally, the \ref database_page.
*/


/*! \page database_page DataBase

\subsection intro_database What does our so-called data base?
The purpose of having some database code is to decouple most of the processing
code from the actual way data is stored. The data can then be accessed through
routines which transform abstract names into genuine file names. Implementing this
way provides the tremendeous advantage that many actual database organizations
can be accomadated through the same software interface. Most (if not all) code
does not require any modification in case the actual data base used changes.
  At the moment, the proposed implementation only deals with images,
not subtractions. What the data base software handles is file names,
not contents. 

\subsection dbconfig The data base configuration file.

The actual location where data is to be searched for is
given through a configuration file. This configuration file is searched
 -# if the environment variable DBCONFIG is defined, as the file name it provides,
 -# as .dbconfig in the current directory
 -# as $HOME/.dbconfig

\subsection dbconfig_example Example of Db configuration file
A Db config file consists in mnemonic tags followed by actual file pathes.
Those file pathes can use * [] but not {}.
Here is an example of a Db config file:
\code

#this is a comment
ImagePath
{
here : .
cfht99 : /snovad15/cfht99/1999*
vlt99 : /snovad1/vlt99/1999*
newstuff : /snovad8/wiyn99
}

\endcode
*/


#ifdef TODO

- the handling of saveEverythingElse in the copies
- why no copy constructor?

#endif



#ifdef __cplusplus /* this file is also included by the dbfileParse.c 
which is C code generated from dbfile.y */

#include <string>
#include <iostream> /* for prototyping */

using namespace std;


typedef enum DbImageKind { Raw = 1, Calibrated, Elixir};
typedef enum DbImageCatalogKind { SExtractor = 1, Fitted_for_seeing , DaophotAls, DaophotAp, DaophotNei, DaophotLst, DaophotPk, DaophotNst};
typedef enum DbImagePsfKind { DaophotPsf = 1};

class Path;

/*! A DbImage refers to one image as the telescope provides it (more precisely, one CCD),
  together with associated data used for the reduction (flat and bias frames) or produced 
  during the reduction (lists of stars). */


#include "rootstuff.h"

class DbImage 
#ifdef USE_ROOT
: public TObject
#endif /* USE_ROOT */
{

public:



private :
  string imageName;
  string directory;//!
  const Path *path;//!
  bool create(const string &ActualPath);

protected :
  bool saveEverythingElse; //! // wether we have to write the EverythingElse file in destructor

public :

  
  //! a constructor: its argument is a unique image identifier (eg r124280). 
  explicit DbImage(const string &ImageName);
  
  //! a constructor: its argument is a unique image identifier (eg r124280). 
  explicit DbImage(const char *ImageName);
  
  DbImage():imageName(""),directory(""),path(NULL), saveEverythingElse(false){};
  //! for images obtained from the above constructor, checks that the image could be located.
  bool IsValid() const { return (path != NULL /*directory != "" */);};

  //! returns the image name.
  string Name() const {return imageName;}
    
  string Dir() const {return directory;}

  DbImage(const string &ImageName, const Path* APath); /* used when collecting DbImages that match a name (with wildcards) */

  //! to tell if we have twice the same image.
  bool operator == (const DbImage &Right) const;


  //! returns the FitsImage file name (a file name for the file system). 
  /*!   The Kind argument  can be Raw or FlatFielded. Nothing
    ensures that the file exists. One may use the FileExists routine 
    to check. */
  string FitsImageName(const DbImageKind Kind) const ;

  //! out of elixir name
  string ElixirName() const;

  //!  same for the ImageBack (background value and rms interpolated map) 
			    /*! Same remark as above concerning the existence. */
  string ImageBackName(const DbImageKind Kind) const ;

  //! the name of the fits file containig the flatfield used for flatfielding.
  string FitsFlatName() const;

  //! same for bias.
  string FitsBiasName() const;


  //! same for dark.
  string FitsDarkName() const;

  //! name of the weight imag
  string FitsWeightName() const;

  //! same for dead pixel map. 
  string FitsDeadName() const;
  
  //! same for bad pixel map (built from weight map). 
  string FitsBadName() const;
  
  //! same for cosmic pixel map. 
  string FitsCosmicName() const;

  
  //! same for satellite pixel map. 
  string FitsSatelliteName() const;
  
  //! same for fringe pattern map map. 
  string FitsFringeName() const;

  //! background image
  string FitsBackName() const;


  //! min background image
  string FitsMiniBackName() const;
  //! same for saturated stars pixels map. 
  string FitsSaturName() const;

  //! return the results of the usno match
  string ImageMatchUsnoName() const;

  //! returns the list of stars detected and measured on the Image (a file name for the file system). 
  /*!  The Kind argument  can be SExtractor or Fitted_for_seeing. 
    See FitsImageName for caution instructions*/
  string ImageCatalogName(const DbImageCatalogKind Kind=SExtractor) const ;

  //! returns the name where the psf parameters and look-up table for residuals is stored
  string ImagePsfName(const DbImagePsfKind Kind=DaophotPsf) const; 


  string EverythingElseFileName() const {return Dir()+"reddata.root";}

  //!To create the directories where the fits images, catalogues
  //! will be put: ex: ~/FakeDb/test: DbImage dbim("test"); dbim.Create("~/FakeDb/");
  bool Create(const string &Path);
    
  /* a somehow generic entry for all the above ones */
  string GetFileName(const char* WhichFile) const;

  friend ostream& operator <<(ostream &stream, const DbImage &I) { I.dump(stream); return stream;}

  void dump(ostream &stream = cout) const;


  ~DbImage();

  friend class DbImageList;

#ifdef USE_ROOT
  ClassDef(DbImage,1);
#endif /* USE_ROOT */

  // routines which have to do with IOs
  void init_from_name();
  bool writeEverythingElse();

};

#include <list>

//! DbImages can be globally located and stored into image lists. 

class DbImageList : public list<DbImage> {


  public :
    //! locates all DbImage's in a given symbolic path name (see db configuration file section).
    /*!   This routine handles wildcards : * will match any path name, 'new*' will match any path
       name beginning with 'new'. */
    DbImageList(const char * PathName);
    DbImageList(const string &PathName);
  //! tries to interpret ecah argument either as a DbImage name or as a symbolic path
    DbImageList(const list<string> names);

  DbImageList() {};

  //! discards from a DbImageList all images not matching a date 
    void FilterByDate(const int a_date);


    //! DOCF collects  and appends DbImages inside PathName. 
    int Collect(const char * PathName);

    void dump(ostream &stream = cout) const;
    //! example of usage : cout << DbImageList("*"); . It works as the first statement of your program !
    friend ostream & operator << (ostream& stream, const DbImageList& L) {L.dump(stream); return stream;};
    // SelectedDump(const int WhichInfos) const;
   
};


typedef list<DbImage>::iterator       DbImageIterator;
typedef list<DbImage>::const_iterator DbImageCIterator;

//! reads the config file. Called automatically if not user called.
/*! If you programm does a lot of things before actually calling any DbImage constructor, 
it may be a good idea to call it early in your main to check if the
configuration file is correct. */

void DbInit();

//! dumps to cout the (interpreted) contents of the configuration file. 
/*! It triggers DbInit if not already done. */
void DbConfigDump(ostream &stream = cout);

/* called by image_install.cc ;  this one knows the db structure */
int InstallImage(const char *a_path, const char *a_file, DbImageKind kind);

int AssignInfo  (const DbImage &Image, const string &FlatFitsFileName, const char * WhichInfo);

int DbConfigSetDumpLevel(const int level);

string DbConfigFileName();


#endif /* __cplusplus */



#ifndef __CINT__
#ifdef __cplusplus
extern "C" {
#endif

  /* these ones are here because the parser is in C (it is yacc/lex generated) and is
     hence separated from dbimage.cc. */

void DbConfigAddImagePath(const char * a_path, const char *a_path_name);
int  DbConfigFileParse(const char*ConfigFileName);

#ifdef __cplusplus
}
#endif
#endif /* __CINT__ */


#endif /* DBIMAGE__H */
