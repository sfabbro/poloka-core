#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <fitsio.h>
#include "fitsimage.h"
#include "fitstoad.h" 
#include "fileutils.h"
#include "frame.h"
#include "stringlist.h"


#include "fitsexception.h"

extern "C" {
#include <fitsio2.h> // for ffr4fi{1,2}
}


/******** FitsKey *********/

//#define DEBUG

/* there are 2 kinds of keys : genuine fits keys and "virtual" ones. 
The first kind ones have a file pointer and 
the second kind ones have value fields (sval and dval).
This is basically why the conversion operators deal with those 2 cases.
It would have been possible but somehow clumpsy to use an abstract
FitsKey class with 2 derived classes. */

/* we read numbers as strings because sometimes the same key (or equivalent keys) 
appears as a string or a number on different telescopes. but reading a number
as a string works with cfitsio routines. */

FitsKey::operator int() /* const */
{
if (fptr) /* genuine key to be read on a file */
  {
    int status=0;
    char value[256];
     fits_read_key(fptr, TSTRING, keyName, &value, NULL, &status);
    if (status == KEY_NO_EXIST)
      {
      if (warn) 
	cout << "trying to read key " << keyName 
	     << " which is not there " << endl;
      return 0;
      }
    return atoi(value);
  }
else return int(dval); /* "virtual" key */
}

FitsKey::operator float() /* const */  /* may be useless, because double would do the work?? */
{
if (fptr)
  {
    int status=0;
    char value[256];
    fits_read_key(fptr, TSTRING, keyName, &value, NULL, &status);
    if (status == KEY_NO_EXIST)
      {
      if (warn) 
	cout << "trying to read key " << keyName 
	     << " which is not there " << endl;
      return 0;
      }
    return float(atof(value));
  }
else return float(dval);
}

FitsKey::operator double() /* const */
{
if (fptr)
  {
    int status=0;
    char value[256];
    fits_read_key(fptr, TSTRING, keyName, &value, NULL, &status);
    if (status == KEY_NO_EXIST)
      {
      if (warn) 
	cout << "trying to read key " << keyName 
	     << " which is not there " << endl;
      return 0;
      }
    return atof(value);
  }
else return dval;
}

FitsKey::operator bool() /* const */
{
  if (fptr)
    {
      int status=0;
      // have to go through an int because cfitsio uses int for "logicals"
      int value = 0; 
      fits_read_key(fptr, TLOGICAL, keyName, &value, NULL, &status);
      if (status == KEY_NO_EXIST)
	{
	  if (warn) 
	    cout << "trying to read key " << keyName 
		 << " which is not there " << endl;
	  return false;
      }
      return bool(value);
    }
  else return false;
}



FitsKey::operator string() /* const */
{
if (fptr)
  {
  int status=0;
  char a_C_string[80];
  a_C_string[0] = '\0'; // default = empty string
  fits_read_key(fptr, TSTRING, keyName, a_C_string, NULL, &status);
  if (status == KEY_NO_EXIST)
    {
      if (warn) 
	cout << "trying to read key " << keyName 
	     << " which is not there " << endl;
    }
  return string(a_C_string);
  }
else return sval;
}


ostream& operator <<(ostream& stream, const FitsKey &This)
{
if (This.fptr)
  {
  char a_string[80];
  char a_comment[80];
  int status = 0;
  char s_key_name[80];
  strcpy(s_key_name, This.keyName);
  /* read the data card in the FITS header as a string, and print out Key and key value */
  fits_read_key(This.fptr, TSTRING, s_key_name, a_string, a_comment, &status);
  //stream << This.keyName << ' ' << a_string;
  // If we want to print only the value to use it in a shell for exemple
  stream << a_string;
  }
else
  {
  stream << This.sval;
  }
return stream;
}



/* error reporting and exception throwing:
   the very serious errors now throw a FitsException
   that calling routine may or may not catch.
   Since this is new (by oct 2006), most of the
   current code does not catch anything, including mains.
   This will change as time goes on...
*/

#define CHECK_STATUS(status, message, return_statement) \
if (status)\
  {\
  cfitsio_report_error(status, message);\
  cerr << " for file '" << this->fileName << "'" << endl;\
  return_statement;\
  }

static void cfitsio_report_error(int status, const string &Message)
{
cerr << Message ;
fits_report_error(stderr, status);
}

static void print_cfitsio_long_error_message(ostream &s, int status, const FitsHeader &Head, const char *Caller)
{
  char errmess[80];
  fits_get_errstatus(status, errmess);  /* get the error description */
  s << "ERROR: from " << Caller << ", for file " <<  Head.FileName() << endl
    << " FITSIO status = " << status 
    << ": " << errmess << std::endl;
  while ( fits_read_errmsg(errmess) )  /* get error stack messages */
    s << errmess << endl;;
}

static string one_line_error_message(int status, const FitsHeader &Head)
{
  char errmess[80];
  fits_get_errstatus(status, errmess);
  char temp[256];
  sprintf(temp,"(IO) %s / file=%s",errmess, Head.FileName().c_str());
  return string(temp);
}

#define CHECK_STATUS_AND_THROW_EXCEPTION(status, Head, message)	\
if (status)\
  {\
    print_cfitsio_long_error_message(std::cerr,status, Head, message);\
    throw(FitsException(one_line_error_message(status, Head)));\
  }\




static const char *FileModeName(const FitsFileMode Mode)
{
if (Mode == RO) return "ReadOnly";
if (Mode == RW) return "ReadWrite";
return " A bad Mode was provided";
}


/*************** FitsHeader ***********/

/* implementation choice : a fits image, which needs a fits header, always have
an associated file. if the file does not exist, it is created. */


FitsHeader::FitsHeader()
{
  cout << " gone through the default FitsHeader constructor: this should not happen " << endl;
}


FitsHeader::FitsHeader(const string &FileName, const FitsFileMode Mode)
{
#ifdef DEBUG
  cout << " > FitsHeader::FitsHeader(const string &FileName, const FitsFileMode Mode) " 
       << "FileName = " << FileName
       << endl;
#endif

  int status = 0;
  fileName = FileName;
  fileModeAtOpen = Mode;
  compressedImg = false; 
  fptr = 0;
  
  telInst = NULL;
  writeEnabled = true; // only relevant if mode != RO
 
  /* handle here ascii header files (which differ from FITS by the
     fact that they are split in lines and lines are not padded with
     blanks. The strategy is to transform them into genuine fits
     files, but in memory (since cfitsio enables that).

     So the strategy is simple: detect wether the file is an ascii file,
     and if so copy it into a memory fits file. Then, everything works as usual,
     except that:
     - it is a readonly mechanism
     - it does not handle compressed (gzipped) files
     - do not expect pixels in these files....
  */

  bool is_ascii=false;
  if (!strstr(fileName.c_str(),".gz")) // is NOT a compressed file 
    {
      FILE * f = fopen(fileName.c_str(), "r");
      if (f)
	{
	  char line[96];
	  // try to read (at most) 82 chars, to see if the header is ascii or fits
	  is_ascii = ((fgets(line,82,f)) && strlen(line)<=80);
	  fclose(f);
	}
    }


  if (is_ascii)
    {
      if (Mode == RW) 
	cerr << " The file " << fileName 
	     << " appears to be an ascii file and is requested to be opened RW " << endl
	     << " The code does not handle that: the file will not be modified " << endl;      
      FILE * f = fopen(fileName.c_str(), "r");
      string memName="mem://"+BaseName(fileName); // hopefully unique name
      fits_create_file(&fptr,memName.c_str(), &status);
      CHECK_STATUS_AND_THROW_EXCEPTION(status,*this,"fits_create_file(mem)");
      char line[96];
      while (fgets(line,82,f))
	{
	  char *last_char = line+strlen(line)-1;
	  if (last_char<line) /* empty_line */ continue;
	  if (strcmp(line,"END")==0) break;
	  if (*last_char=='\n') *last_char='\0'; // remove '\n' (mandatory)
	  fits_write_record(fptr,line,&status); // copy
	  CHECK_STATUS(status,line,);
	}
      fclose(f);
      return; // if it is an ascii input file, we are done.
    }

  // genuine fits file.
  else fits_open_file(&fptr, fileName.c_str(), int(Mode), &status);
 /* there seems to be a bug in cfitsio, when one requests some
    preprocessing before actually accessing the data, e.g. adressing
    subimages through actual_filename[imin:imax,jmin:jmax]. Then the
    returned pointer has RW mode although RO was requested. this is
    handled using fileModeAtOpen, rather than trying to hack cfitsio
    internals
 */
 

  if (status == 0)
    {
      int hdutype;
      fits_get_hdu_type(fptr, &hdutype, &status);
      if ( (hdutype == IMAGE_HDU) &&(int(KeyVal("NAXIS")) == 0 ))
        {// this is a NULL IMAGE_EXTENSION, we move forward, because we expect
         // to get an actual image.
	  fits_movrel_hdu(fptr, 1, NULL, &status);
	  //CHECK_STATUS(status," movrel_hdu in FitsHeader::FitsHeader",);
	  //	  cout << " moving by 1 HDU since the image is a compressed table " 
	  //   << endl;
        }
      // compressed images have a ZIMAGE key, so:
      compressedImg = (HasKey("ZIMAGE") && bool(KeyVal("ZIMAGE")));
    }

  /* when reopening RW a file image which was just opened RW
     cfitsio tries to uncompress the image, and fails. this is 
     the error 414, which we just ignore
  */
  if (status == 414) status = 0;

 if (status == 0) 
   {
     existingFile = true;
     return; // the file exists and was opened
   }
 existingFile = false;
 if (Mode == RO) // no need to try to create the file
   CHECK_STATUS_AND_THROW_EXCEPTION(status,*this,"FitsHeader in mode RO");
 
 // from here on, we are in RW mode and could not open an existing file

 int first_status = status; // for eventual error reporting
 int second_status = create_file();

 if (second_status)
   {
     cerr << " when opening file : " << fileName;
     cerr << ", with mode " << FileModeName(Mode) <<  ',' << endl;
     cerr << " we got these successive errors:" << endl;
     fits_report_error(stderr, first_status);
     // now throw instaed of going on.
     throw(FitsException(one_line_error_message(second_status, *this)));
     // the alternative : report error (and possibly abort).
     // fits_report_error(stderr, second_status);
     // abort();
   }
}


int FitsHeader::create_file()
{
#ifdef DEBUG
  cout << " > FitsHeader::create_file() " << endl;
#endif

 int status = 0;
 fileModeAtOpen = RW;

 if (FileExtension(fileName) == "fz") 
   /* this is a Toads convention: requesting a (non-existent) file
      named fz means you want a compressed image. In this source file,
      compression refers to internal compression of the pixel buffer,
      not external compression of files using gzip (which cfitsio
      handles totally transparently).

      Notes about the way we use internal compression (RICE
      compression by default) is decribed in a long comment at the end
      of this source file.
   */

   {
     fits_create_file(&fptr, (fileName+"[compress]").c_str(), &status);
     compressedImg = false;
     minimum_header();// minimum_header refuses compressed images.
     compressedImg = true;

     // create a dummy compressed image 
     // (with 0 pixels, it is not created)
     long naxis[2];
     naxis[0] = naxis[1] = 1;
     fits_create_img(fptr,16,2,naxis, &status);
     /* now set NAXIS1 and NAXIS2 to 0 so that nobody believes that 
	there are pixels */
     ModKey("NAXIS1",0);
     ModKey("NAXIS2",0);
   }
 else
   {
     compressedImg = false;
     fits_create_file(&fptr, fileName.c_str(), &status);
     minimum_header();
   }
 return status;
}
 



void FitsHeader::minimum_header()
{
#ifdef DEBUG
  cout << " > FitsHeader::minimum_header() " << endl;
#endif

  if (compressedImg)
    {
      cerr << " should never come here with a compressed image " << std::endl;
      abort();
    }
AddKey("SIMPLE", true);
AddKey("BITPIX", 16);
AddKey("NAXIS", 2);
AddKey("NAXIS1",0);
AddKey("NAXIS2",0);
AddKey("ORIGIN","This horrible toads software");
}

FitsFileMode FitsHeader::FileMode() const
{
  return fileModeAtOpen; // rather than file_mode
}

#ifdef STORAGE
FitsHeader::FitsHeader(const FitsHeader &Header)
{
  cerr <<" Using a routine that badly supports compressed images ! " << endl;
  int status =0;
  telInst = 0;
  writeEnabled = Header.writeEnabled; // probably irrelevant
  fileName = Header.fileName;
  compressedImg = Header.compressedImg;
  fits_open_file(&fptr, Header.fileName.c_str(), RO, &status);
  fileModeAtOpen = RO;
  CHECK_STATUS(status," FitsHeader(fitsHeader &) ", return );
}
#endif

FitsHeader::FitsHeader(const FitsHeader &Template, 
		       const string &NewFileName)
{
#ifdef DEBUG
  cout << " > FitsHeader::FitsHeader(const FitsHeader &Template,const string &NewFileName) "
       << " NewFileName = " << NewFileName
       << endl;
#endif
  fileName = NewFileName;
  remove(fileName.c_str());
  compressedImg = false; // will be eventually set by create_file()
  create_file();
  if (Template.HasKey("ORIGIN")) if (HasKey("ORIGIN")) RmKey("ORIGIN");

  // copy by hand the keys that may require a translation
  /* if one or both of the images are compressed, the following
     keys need to be "translated" (BITPIX -> ZBITPIX), and our routines
     do that
  */

  int bitpix = Template.KeyVal("BITPIX"); AddOrModKey("BITPIX",bitpix);
  int naxis = Template.KeyVal("NAXIS"); AddOrModKey("NAXIS",naxis);
  int naxis1 = Template.KeyVal("NAXIS1"); AddOrModKey("NAXIS1",naxis1);
  int naxis2 = Template.KeyVal("NAXIS2"); AddOrModKey("NAXIS2",naxis2);

/* copy all the user keywords (not the structural keywords) */
  int nkeys = 0;
  int status = 0;

  fits_get_hdrspace(Template.fptr, &nkeys, NULL, &status); 

  for (int ii = 1; ii <= nkeys; ii++) 
    {
      char card[81];
      fits_read_record(Template.fptr, ii, card, &status);
      if (fits_get_keyclass(card) > TYP_CMPRS_KEY)
	fits_write_record(fptr, card, &status);
   }


/* the following Flush has the effect that the produced image file
cannot be smaller than the one associated with Template. this is disastrous.
calling fits_resize_img does not help. */
/* Flush(); */
 telInst = 0;
 writeEnabled = true;
 existingFile = false;
 CHECK_STATUS(status, " CopyHeader ", )
}


/*const ? */ VirtualInstrument* FitsHeader::TelInst() const
{
  // this is a lie, it is NOT const...
  if (telInst) return telInst;
  FitsHeader &pseudo_copy = (FitsHeader &) *this;
  pseudo_copy.telInst = SniffTelInst(*this);
  return telInst;
}


void FitsHeader::DeleteFile()
{
  int status = 0;
  fits_delete_file(fptr, &status);
  // we get a bad status when the file is corrupted, but we just don't care
  //   CHECK_STATUS(status, "DeleteFile",);
  fptr = NULL;
}

FitsHeader::~FitsHeader()
{int status = 0;
 if (fptr)
   {
     if (FileMode() == RW && !writeEnabled && !existingFile)
       {
	 cout << "deleting file " << fileName << endl;
	 fits_delete_file(fptr,&status);
       }
     else fits_close_file(fptr, &status);
     /* When closing a file that contains a compressed image, 
	cfitsio tries to uncompress it, perhaps with good reasons.
	This may fail, especially when there is no actual image 
	on the file yet. We then get error 414, which 
	seems to be harmless.
	Other problem here : when we provide cfitsio with ascii header files,
	we map them on memory files. Often, these "ascii header files" miss
	structural keys (SIMPLE, BITPIX, ...). Cfitsio complains when we 
	delete the memory file. We just don't care. 

     */

     if (status != 414 &&  
	 (status != 252 // 1st key not SIMPLE or XTENSION
	  || fileName.find("mem://") == 0) // this is a "memory file" which has not disk counterpart.
	 )
       CHECK_STATUS(status, " ~FitsHeader "+fileName, )
   }
 if (telInst) VirtualInstrumentDestructor(telInst);
}

void FitsHeader::ImageSizes(int &nx, int &ny) const
{
nx = KeyVal("NAXIS1"); 
ny = KeyVal("NAXIS2"); 
}

bool FitsHeader::SameImageSizes(const FitsHeader &Other) const
{
  int nx = KeyVal("NAXIS1");  
  int ny = KeyVal("NAXIS2");
  int onx = Other.KeyVal("NAXIS1");
  int ony = Other.KeyVal("NAXIS2");
  return ((nx == onx) && (ny == ony));
}


Point FitsHeader::ImageCenter() const
{
  return Point(0.5*double(KeyVal("NAXIS1")), 0.5*double(KeyVal("NAXIS2")));
}


// if you add keys in this "translator", you should also
// probably add them into FitsHeader::FitsHeader(OldHeader, NewName);
static std::string common_to_file_name(const std::string &KeyName)
{
      if (KeyName == "BITPIX" || KeyName == "NAXIS" || 
	  KeyName == "NAXIS1" || KeyName == "NAXIS2")
	return "Z"+KeyName; 
      return KeyName;
}

#ifdef STORAGE
static std::string file_to_common_name(const std::string &KeyName)
{
      if (KeyName == "ZBITPIX" || KeyName == "ZNAXIS" || 
	  KeyName == "ZNAXIS1" || KeyName == "ZNAXIS2")
	return KeyName.substr(2,KeyName.size()); 
      return KeyName;
}
#endif

FitsKey 
FitsHeader::KeyVal(const string &KeyName, const bool Warn) const
{
  if (KeyName.find("TOAD") != KeyName.npos && !HasActualKey(KeyName)) 
    {return ToadsKeyVal(*this,KeyName,Warn);}
  else
    {
      // trap calls to NAXIS* and BITPIX
      string keyName = KeyName;
      if (compressedImg) keyName = common_to_file_name(KeyName);
      return FitsKey(keyName, fptr, Warn);
    }
}


int FitsHeader::mod_key(int type, const string &ThisKeyName, void *Value, 
			const string &Comment) const
{
if (fptr)
  {
  int status = 0;

  std::string KeyName = ThisKeyName;
  if (compressedImg) KeyName = common_to_file_name(ThisKeyName);


  // have to copy because of missing const's in cfitsio
  char s_keyname[80];
  strcpy(s_keyname,KeyName.c_str());
  char comment[256];
  char *the_comment = NULL;
  if (Comment!= "" )
    {
      strcpy(comment,Comment.c_str());
      the_comment = comment;
    }
  fits_update_key(fptr, type, s_keyname, Value, the_comment, &status);
  CHECK_STATUS(status, " ModKey : " + string(KeyName) ,);
  return (!status);
  }
else return 0;
}


int FitsHeader::ModKey(const string &KeyName, const int Value, const string Comment) const
{
  int value = Value;
  return mod_key(TINT, KeyName, &value, Comment);
}

int FitsHeader::ModKey(const string &KeyName, const double Value, const string Comment) const
{
  double value = Value;
  return mod_key(TDOUBLE, KeyName, &value, Comment);
}

int FitsHeader::ModKey(const string &KeyName, const char *Value, const string Comment) const
{
  char value[256];
  strcpy(value,Value);
  return mod_key( TSTRING, KeyName, value, Comment);
}

int FitsHeader::ModKey(const string &KeyName, const bool Value, const string Comment) const
{
  // have to cast to an int because cfitsio goes through int type for logicals
  int value = (int) Value;
  return mod_key(  TLOGICAL, KeyName, &value, Comment);
}



int FitsHeader::write_key(const string &ThisKeyName, void* KeyVal, const string &Comment, const int type) 
{
  std::string KeyName = ThisKeyName;
  if (compressedImg) KeyName = common_to_file_name(ThisKeyName);

  if (!fptr || KeyName.length()==0 || !KeyVal) return 0;
  if (HasKey(KeyName)) {cerr << " trying to write a second " << KeyName << "keyword" << endl; return 0;}
  int status = 0;
  // have to copy because of missing const's in cfitsio prototypes.
  // very innefficient !
  char s_key_name[80];
  strcpy(s_key_name,KeyName.c_str());
  char comment[256];
  char *the_comment = NULL;
  if (Comment != "")
    {
      strcpy(comment, Comment.c_str());
      the_comment = comment;
    }
  fits_write_key(fptr, type, s_key_name, KeyVal, the_comment, &status);
  CHECK_STATUS(status,"AddKey : "+ string(KeyName),);
  if (FileMode() != RW)
    {
      cerr << " trying to write key " << KeyName << " in file " 
	   << fileName << " opened RO " << endl;
    }
  return (!status);
}


int FitsHeader::AddKey(const string &KeyName, const char* KeyVal, const string Comment)
{
  char keyval[256];
  strcpy(keyval,KeyVal);
  return write_key(KeyName, keyval, Comment, TSTRING);
}



int FitsHeader::AddKey(const string &KeyName, const double KeyVal, const string Comment)
{
  double keyval = KeyVal;
  return write_key(KeyName, &keyval, Comment,TDOUBLE);
}

int FitsHeader::AddKey(const string &KeyName, const int KeyVal, const string Comment)
{
  int keyval = KeyVal;
  return write_key(KeyName, &keyval, Comment, TINT);
}


int FitsHeader::AddKey(const string &KeyName, const bool KeyVal, const string Comment)
{
  // have to cast to an int because cfitsio goes through int type for logicals
  int keyval = (int) KeyVal;
  return write_key(KeyName, &keyval, Comment, TLOGICAL);
}


int FitsHeader::AddOrModKey(const string &KeyName, const char *Value, const string Comment)
{
if (fptr)
  {
  if (HasKey(KeyName)) return ModKey(KeyName,Value,Comment);
  else return AddKey(KeyName,Value,Comment);
  }
else return 0;
}



int FitsHeader::AddOrModKey(const string &KeyName, const int Value, const string Comment)
{
if (fptr)
  {
  if (HasKey(KeyName)) return ModKey(KeyName,Value,Comment);
  else return AddKey(KeyName,Value,Comment);
  }
else return 0;
}

int FitsHeader::AddOrModKey(const string &KeyName, const double Value, const string Comment)
{
if (fptr)
  {
  if (HasKey(KeyName)) return ModKey(KeyName,Value,Comment);
  else return AddKey(KeyName,Value,Comment);
  }
else return 0;
}

int FitsHeader::AddOrModKey(const string &KeyName, const bool Value, const string Comment)
{
if (fptr)
  {
  if (HasKey(KeyName)) return ModKey(KeyName,Value,Comment);
  else return AddKey(KeyName,Value,Comment);
  }
else return 0;
}

int FitsHeader::KeyMatch(const string &KeyPattern, FitsKeyArray &Array) const
{
  Array.clear();
  char accept[32];
  strcpy(accept,KeyPattern.c_str());
  // rewind the file
  int status = 0;
  char card[256];
  char value[256];
  fits_read_record(fptr,0, card, &status);
  do
    {
      char *accept_list = accept;
      fits_find_nextkey(fptr, &accept_list, 1, NULL, 0, card, &status);
      if (status != 0) break;
      // key the actual keyname
      char keyname[80];
      int namelen;
      fits_get_keyname(card, keyname, &namelen, &status);
      fits_parse_value(card, value, NULL, &status);
      Array.push_back(FitsKey(keyname, value));
    } while (true);
  return int(Array.size());
}


int FitsHeader::ModKeyName(const string &OldKeyName, const string &NewKeyName) const
{
if (fptr)
  {
    int status = 0;
    fits_modify_name(fptr, const_cast<char*>(OldKeyName.c_str()), 
		     const_cast<char*>(NewKeyName.c_str()), &status);
    CHECK_STATUS(status, "ModKeyName ",);
    return (!status);
  }
else return 0;
}

int FitsHeader::ModKeyComment(const string &KeyName, const string &NewComment) const
{
if (fptr)
  {
    int status = 0;
    fits_modify_comment(fptr, const_cast<char*>(KeyName.c_str()), 
			const_cast<char*>(NewComment.c_str()), &status);
    CHECK_STATUS(status, "ModKeyComment ",);
    return (!status);
  }
else return 0;
}

bool FitsHeader::HasActualKey(const string &KeyName, const bool Warn) const
{
/* did not find an easy way in cfitsio of only checking if a key is there or not */
int status = 0;
char a_C_string[80];
fits_read_key(fptr, TSTRING, const_cast<char*>(KeyName.c_str()), a_C_string, NULL, &status);
if (Warn && status == KEY_NO_EXIST)
    cout << " file : " << fileName << " has no " << KeyName << endl;   
return (status != KEY_NO_EXIST);
}

bool FitsHeader::HasKey(const string &ThisKeyName, const bool Warn) const
{
  std::string KeyName = ThisKeyName;
  if (compressedImg) KeyName = common_to_file_name(ThisKeyName);

if (HasActualKey(KeyName, Warn)) return true;
if (KeyName.find("TOAD") != KeyName.npos) 
  {
  string value = ToadsKeyVal(*this, KeyName, Warn);
  if (value != NOVAL) return true;
  }
return false;
}


int FitsHeader::RmKey(const string &ThisKeyName) const
{
  std::string KeyName = ThisKeyName;
  if (compressedImg) KeyName = common_to_file_name(ThisKeyName);

  int status = 0;
  char keyName[256];
  strncpy(keyName,KeyName.c_str(),256); // copy for constness
  fits_delete_key(fptr, keyName, &status);
  CHECK_STATUS(status, " RmKey " + KeyName,);
  return (!status);
}

int FitsHeader::NKeys() const
{
int nkeys = 0; int morekeys; int status = 0;
fits_get_hdrspace(fptr, &nkeys, &morekeys, &status);
return nkeys;
}

int  FitsHeader::AddCommentLine(const string &Comment)
{
  int status = 0;
  if (fptr==0) return 0;
  fits_write_comment(fptr, Comment.c_str(), &status);
  CHECK_STATUS(status,"AddCommentLine",);
  return (!status);
}

int  FitsHeader::AddHistoryLine(const string &History)
{
int status = 0;
if (fptr==0) return 0;

fits_write_history(fptr,const_cast<char*>(History.c_str()), &status);
CHECK_STATUS(status,"AddHistoryLine",);
return (!status);
}

bool FitsHeader::ReadCard(const std::string &KeyName, string &Card) const
{
  int status = 0;
  char keyName[80];
  strcpy(keyName,KeyName.c_str());
  char card[256];
  if (fits_read_card(fptr, keyName, card, &status) == 0)
    {
      Card = string(card);
      return true;
    }
  return false;
}

int FitsHeader::AddOrModCard(const string &KeyName, const string &Card)
{
  char keyName[80];
  strcpy(keyName,KeyName.c_str());
  char card[256];
  strncpy(card,Card.c_str(),80);
  int status = 0;
  // in fact, fits_update_card adds or updates
  fits_update_card(fptr, keyName, card, &status);
  CHECK_STATUS(status, "fits_update_card",);
  return (!status);
}

#include <fstream> // for ofstream

bool FitsHeader::AsciiDump(const string &AsciiFileName, 
			   const StringList &WhichKeys) const
{
  string card;
  ofstream s(AsciiFileName.c_str());
  for (StringCIterator i = WhichKeys.begin(); i != WhichKeys.end(); ++i)
    if (ReadCard(*i, card))
      {
	s << card << std::endl;
      }
  s << "END     " << std::endl;
  return true;
}



int FitsHeader::Flush()
{
  int status = 0;
  fits_flush_file(fptr,&status); /* to make sure that the new key
				      values will be used by fits_write_img */
  CHECK_STATUS(status,"Flush",);
  return (!status);
}


ostream & operator << (ostream &stream, const FitsHeader &Header)
{
int nkeys = Header.NKeys();
char card[256];
int status = 0;
for (int i=1; i<= nkeys; i++ ) /* fortran numbering of header lines  */
  {
  card[0] = '\0';
  fits_read_record(Header.fptr, i, card, &status);
  if (card[0]) stream << card << endl;
  }
return stream;
}


bool FitsHeader::CopyKey(const std::string &KeyName, FitsHeader &To) const
{
  // this routine does not handle "key translation", because I did
  // not find a type independant way to do it.
  char card[512];
  int status = 0;
  char keyName[80];
  strncpy(keyName, KeyName.c_str(), 80);
  // copy everything verbatim (name ,value, comment)
  fits_read_card(fptr, keyName, card, &status);

  To.AddOrModCard(KeyName, card);
  CHECK_STATUS(status,"CopyKey", );
  return (status == 0);
}


/* sums 2 headers by simply putting all keys together. If a key appears in both,
keeps the value in 'this'. */

void FitsHeader::Append_LowPriority(const FitsHeader& ToAppend)
{
int nkeys = ToAppend.NKeys();

int status = 0;
for (int i=1; i<= nkeys; i++ ) /* fortran numbering of header lines  */
  {
  char key_name[256], key_val[256], key_comment[256]; 
  char card[256];
  // get the header line split into pieces.
  fits_read_keyn(ToAppend.fptr, i, key_name, key_val,key_comment, &status);
/* copy COMMENTS and HISTORY even if there is some already */
  if (string(key_name) != "COMMENT" && string(key_name) != "HISTORY" 
       && HasKey(key_name)) continue; /* because of Low Priority */
  /* reget the whole line */
  fits_read_record(ToAppend.fptr, i, card, &status);
  fits_write_record(fptr, card, &status);
  CHECK_STATUS(status,"fits_write_record", );
  }
}
bool FitsHeader::SameChipFilterInst(const FitsHeader &Other, const bool Warn) const
{
  bool value = (string(this->KeyVal("TOADBAND")) == string(Other.KeyVal("TOADBAND"))
		&& int(this->KeyVal("TOADCHIP"))  == int(Other.KeyVal("TOADCHIP"))
		&& string(this->KeyVal("TOADINST"))  == string(Other.KeyVal("TOADINST")));
  if (Warn && !value)
    {
      cerr << this->fileName << " and " << Other.fileName << " do not refer to the same chip/filter/instrument " << endl;
    }
  return value;
}

bool FitsHeader::SameChipFilter(const FitsHeader &Other, const bool Warn) const
{
bool value = (string(this->KeyVal("TOADBAND")) == string(Other.KeyVal("TOADBAND"))
       && int(this->KeyVal("TOADCHIP"))  == int(Other.KeyVal("TOADCHIP")));
if (Warn && !value)
  {
    cerr << this->fileName << " and " << Other.fileName << " do not refer to the same chip/filter " << endl;
  }
return value;
}

bool FitsHeader::SameChipFilter(const string &OtherFitsName, const bool Warn) const
{
FitsHeader other(OtherFitsName);
return (other.IsValid() && SameChipFilter(other, Warn));
}

bool FitsHeader::SameFilter(const FitsHeader &Other, const bool Warn) const
{
bool value = (string(this->KeyVal("TOADBAND")) == string(Other.KeyVal("TOADBAND")));
if (Warn && !value)
  {
    cerr << this->fileName << " and " << Other.fileName << " do not refer to the same filter " << endl;
  }
return value;
}

bool FitsHeader::SameFilter(const string &OtherFitsName, const bool Warn) const
{
FitsHeader other(OtherFitsName);
return (other.IsValid() && SameFilter(other, Warn));
}

bool FitsHeader::SameChip(const FitsHeader &Other, const bool Warn) const
{
bool value = (int(this->KeyVal("TOADCHIP"))  == int(Other.KeyVal("TOADCHIP")));
if (Warn && !value)
  {
    cerr << this->fileName << " and " << Other.fileName << " do not refer to the same chip" << endl;
  }
return value;
}

bool FitsHeader::SameChip(const string &OtherFitsName, const bool Warn) const
{
FitsHeader other(OtherFitsName);
return (other.IsValid() && SameChip(other, Warn));
}

/********** 4 routines used for splitting fits files with extensions, or accessing raw data */

int FitsHeader::MoveHDU(int HowMany)
{
  int status = 0;
  fits_movrel_hdu(fptr, HowMany, NULL, &status);
  CHECK_STATUS(status,"MoveHDU",);
  return (!status);
}

int FitsHeader:: CopyCHDUTo(FitsHeader &OutHeader)
{
  int status = 0;
  fits_copy_hdu(fptr, OutHeader.fptr, 0, &status);
  CHECK_STATUS(status,"CopyHDUTo",);
  return (!status);
}

int FitsHeader::CopyDataTo(FitsHeader &OutHeader)
{
  int status = 0;
  fits_copy_data(fptr, OutHeader.fptr, &status);
  CHECK_STATUS(status,"CopyDataTo",);
  return (!status);
}  

int FitsHeader::NHDU() const
{
  int nhdu = 0;
  int status =0;
  fits_get_num_hdus(fptr,&nhdu, &status);
  CHECK_STATUS(status,"NHDU", return 0);
  return nhdu;
}
  

void FitsHeader::EnableWrite(const bool YesOrNo)
{
  writeEnabled = YesOrNo;
}


/* there is in principal no need to write such a routine, since
cfitsio handles subimages. But since rice compression and decompression
does not work with data translation (e.g. BITPIX = 16 read in float's)
we provide this routine which does the right stuff, for reading subimages
of compressed images. It then also works for sub images of uncompressed data...
*/

int FitsHeader::read_image(const int xmin, const int ymin, 
			  const int xmax, const int ymax,
			  float *data)
{
#ifdef DEBUG
  cout << " > FitsHeader::read_image(...) " << endl;
#endif
  float nullval = 0;
  int anynull;
  int status = 0;
  int bitpix = KeyVal("BITPIX"); 
  int nx = KeyVal("NAXIS1");
  //  int ny = KeyVal("NAXIS2"); // don't need?
  int nyread = ymax - ymin;
  int nxread = xmax - xmin;
  double bscale = HasKey("BSCALE")? KeyVal("BSCALE") : 1.0;
  double bzero = HasKey("BZERO")? KeyVal("BZERO") : 0.0;
  if (!compressedImg || bitpix == -32 || (bscale == 1.0 && bzero == 0.0) )
    {
      if (nxread == nx) // pixels to read are all in a big chunk
	{
	  fits_read_img(fptr, TFLOAT, nx*ymin+1, nx*nyread, &nullval,  data, 
			&anynull, &status);
	  CHECK_STATUS_AND_THROW_EXCEPTION(status,*this," FitsHeader::read_image(fits_read_img) ");
	}
      else
	{
	  // read line per line
	  float *pdata = data;
	  for (int j = ymin; j < ymax; ++j)
	    {
	      // to be tested....
	      fits_read_img(fptr, TFLOAT, nx*j+xmin+1, nxread, &nullval, 
			      pdata, &anynull, &status);
	      CHECK_STATUS_AND_THROW_EXCEPTION(status,*this," FitsHeader::read_image(fits_read_img/split_read) ");
	      pdata += nxread;
	    }
	}
    }
  else // compressed image
    {/* there is a bug in cfitsio for compressed images : one can only
        read them back using the encoded type.  by default cfitsio
        groups pixels by rows before compressing so reading whole rows
        that contain pixels we are interested in should no have a
        significant CPU hit
     */
      if (bitpix != 16) 
	{
	  std:: cerr << " cannot read back compressed images with bitpix != {16,-32}, contact developpers " << std::endl;
	  std::cerr << "image name " << FileName() << std::endl;
	  char mess[256];
	  sprintf(mess," (IO) cannot read back compressed images with bitpix=%d for file=%s",bitpix,FileName().c_str());
	  throw(FitsException(mess));
	  abort();
	}
      float *pdata = data;
      int nychunk = 100;
      int chunk = nx*nychunk; // far faster than chunk = nx;
      short *temp = new short[chunk];
      // annihilates data rescaling. Perhaps we have to set them back
      fits_set_bscale(fptr, 1., 0., &status);
      CHECK_STATUS(status," fits_set_bscale in FitsImage::FitsImage",);
      for (int j=ymin ; j < ymax; j+= nychunk)
	{
	  int ny2read = std::min(nychunk, ymax-j);
	  fits_read_img(fptr, TSHORT, j*nx+1, ny2read*nx, &nullval, temp, 
			&anynull, &status);
	  int nread = ny2read*nx;
	  if (nxread == nx) // we use all pixels we read
	    for (int i=0; i < nread; ++i) 
	      {
		*pdata = temp[i]*bscale+bzero; pdata++;
	      }
	  else
	    { // hope we are not in a hurry
	      for (int j = 0; j<ny2read; ++j)
		{
		  for (int i=xmin; i < xmax; ++i)
		    {
		    *pdata = temp[j*nx+i]*bscale+bzero; pdata++;
		    }
		}
	    }
	}// end loop on pixel chunks
      CHECK_STATUS_AND_THROW_EXCEPTION(status,*this," fits_read_img in FitsImage::FitsImage");
      delete [] temp;
      if (pdata-data != nxread*nyread)
	{
	  cout << " aie aie aie, on n'a pas pousse le pointeur d'ecriture" 
	       << std::endl
	       << " le bon nombre de fois " << std::endl
	       << " pdata(end) -data = " << pdata-data << std::endl
	       << " nxread * nyread" << nxread*nyread << std::endl;
	}
    }
  return status;
}



/* this routine belongs to FitsHeader to enable classes derived
   from FitsHeader to write images */
int FitsHeader::write_image(int FirstPix, int NPix, Pixel *Data, 
			    int nx, int bitpix, double bscale, double bzero)
{
#ifdef DEBUG
  cout << " > FitsHeader::write_image(FirstPix = " << FirstPix 
       << ", NPix=" << NPix << ", nx="<<nx<<", bitpix="<<bitpix 
       << ", bscale=" << bscale << ", bzero="<<bzero << endl;
#endif


  int status = 0;
  if (!compressedImg || bitpix == -32)
    {
      /* write the pixels straight ahead : either the image uses
      compression without format change (float->float), or no
      compression, which both work for fits_write_img. fits_write_img
      does not do compression and format change together properly. */
      fits_write_img(fptr, TFLOAT, FirstPix+1, NPix, Data, &status);
      CHECK_STATUS_AND_THROW_EXCEPTION(status,*this," fits_write_img in FitsHeader::write_image ");
    }
  else
    {
      /* there is a serious bug in cfitsio : when writing a 
	 compressed image, the requested bitpix (which becomes ZBITPIX
	 is ignored. so we do our conversion ourselves into the target
	 BITPIX
      */
      if (bitpix != 16 && bitpix != 8) 
	{
	  std::cerr << " unhandled condition in FitsHeader::write_image(), contact developpers" << std::endl;
	  abort();
	}

      fits_set_bscale(fptr, 1., 0., &status); // turn off automatic scaling 
      // write slice chunks (it is somehow faster)
      long ychunk = 100;
      long npix_chunk = nx*ychunk;
      Pixel *pixel = Data;
      const Pixel *pend = Data+NPix;
      int firstPix = FirstPix;
      if (bitpix == 16)
	{
	  short *temp = new short[npix_chunk];
	  while (pixel < pend)
	    {
	      long npix = std::min(npix_chunk, (long int)(pend-pixel));
	      ffr4fi2(pixel, npix, bscale, bzero, temp, &status);
	      fits_write_img(fptr, TSHORT, firstPix+1, npix, temp, &status);
	      firstPix+= npix;
	      pixel += npix;
	    }
	  delete [] temp;
	  CHECK_STATUS(status,
		   " fits_write_img (TSHORT) in FitsHeader::write_image ",);
	}

      if (bitpix == 8)
	{
	  unsigned char* temp = new unsigned char[npix_chunk];
	  while (pixel < pend)
	    {
	      int npix = std::min(npix_chunk, (long int)(pend-pixel));
	      ffr4fi1(pixel, npix, 1., 0., temp, &status);
	      fits_write_img(fptr, TBYTE, firstPix+1, npix, temp, &status);
	      firstPix+= npix;
	      pixel += npix;
	    }
	  delete [] temp;
	  CHECK_STATUS(status,
		   " fits_write_img (TBYTE) in FitsHeader::write_image ",);
	}
    }
  return status;
}

int FitsHeader::create_image(int bitpix, int nx, int ny)
{
  int status = 0;
  /* if bitpix < 0 (floating point), then the file might have been created 
     for a compressed image (i.e header + bin table extension), and we
     might be giving up the compression just now. The easy way is to delete
     the file and recreate it from scratch. This is not very expensive
     since we only get here when we write the full pixel stack. */
  if (compressedImg || bitpix < 0)  
    /* 
       what we do here almost reproduces what imcopy (from cfitsio
       sources) does. imcopy only works if the output file does not
       exist. To get compressed images, we did not find any other way
       than reconstructing the file from scratch. Notes about cfitsio
       and compressed images can be found at the end of this file
    */
   {

     /* store all the user keywords (not the structural keywords) */
     int nkeys = 0;
     fits_get_hdrspace(fptr, &nkeys, NULL, &status); 
     vector<string> cards;
     cards.reserve(nkeys); // allocate more than needed : fast push_back

     for (int ii = 1; ii <= nkeys; ii++) 
       {
	 char card[81];
	 fits_read_record(fptr, ii, card, &status);
	 if (fits_get_keyclass(card) > TYP_CMPRS_KEY)
	   cards.push_back(card);
       }

     // delete file
     fits_delete_file(fptr, &status);
     /* structural header keys may be inconsistent, and fits_delete_file
	checks all that .... but we just don't care : no CHECK_STATUS and reset
     */
     //CHECK_STATUS(status," fits_delete_file(compress) in FitsImage::Write",);
     status = 0;
     if (bitpix>0)
       fits_create_file(&fptr, (fileName+"[compress]").c_str(), &status);
     else
       fits_create_file(&fptr, (fileName).c_str(), &status);
     CHECK_STATUS(status," fits_create_file (compress or not) in FitsImage::Write",);
     /* noise bits only matter for floting point images (bitpix = -32 or -64)
	the default value of cfitsio is 4 which really alters images */
     fits_set_noise_bits(fptr,16, &status);
     CHECK_STATUS(status," fits_set noise_bits in FitsHeader::create_image",);

     long int npix[2]; 
     npix[0] = nx;
     npix[1] = ny;
     fits_create_img(fptr, bitpix, 2, npix, &status);
     CHECK_STATUS(status," fits_create_img (compress) in FitsImage::Write",);
     for (unsigned int ii = 0; ii < cards.size(); ++ii)
       fits_write_record(fptr, cards[ii].c_str(), &status);
   }
  return status;
}
/*****************  FitsImage ****************/

FitsImage::FitsImage(const string &FileName, const FitsFileMode Mode) : FitsHeader(FileName, Mode) , Image()
{
#ifdef DEBUG
  cout << " > FitsImage::FitsImage(const string &FileName, const FitsFileMode Mode) " 
       << "FileName = " << FileName 
       << " Mode = " << Mode 
       << endl;
#endif
  written = 0 ;
  if (!IsValid()) return; 
  nx = KeyVal("NAXIS1");
  ny = KeyVal("NAXIS2");
  Image::allocate(nx,ny);
  if (nx*ny) read_image(0,0,nx,ny,data);

  /* My present (06/04) guess is that all the mess that follows in
     this routine could be removed. Check with Kyan.
  */

#ifdef STORAGE
  // there is a problem in the INT data with UltraDas : saturated pixels read
  // out 89 : we fix it up when reading unprocessed raw data.
  // Since we write WRITEDAT in all output files, the correction has already been 
  // applied on files that have this key.
if ( !HasKey("WRITEDAT") && (IsOfKind<IntWfcNewDaq>(*this) || IsOfKind<Cfht12K>(*this))
       && ( fabs(double (KeyVal("BSCALE")) - 1.0) < 1.e-10 ) 
     && !strstr(fileName.c_str(),"dead")  && !strstr(fileName.c_str(),"satur")
          // correction should not be applied to CFH dead/satur maps
    )
  {
      Pixel *end = Image::end();
      int n_pix_changed = 0;
      Pixel max_value = Image::MaxValue();
      Pixel badValue;
      int date = KeyVal("MJD-OBS");
      
      //      Pixel badValue=-1;


      if (IsOfKind<IntWfcNewDaq>(*this)) badValue = 89;
      if (IsOfKind<Cfht12K>(*this) && date <  51460 && date > 51430) {badValue = 0;max_value=65535;}

      for (Pixel *p = Image::begin() ; p<end; ++p) 
	  {
	      if (*p == badValue)
		  {
		      n_pix_changed++;
		      *p = max_value;
		  }
	  }


      if (IsOfKind<IntWfcNewDaq>(*this))
	cout << " Applied the INT UltraDas correction for 89-pixels to " 
	   << n_pix_changed << " pixels. new value : " << max_value << endl ;
      if (IsOfKind<Cfht12K>(*this) && date <  51460 && date > 51430 ) //HD problem only for the cfht99b campagn
	cout << " Applied the CFHT 12K  for 0-pixels to " 
	   << n_pix_changed << " pixels. new value :  " << max_value << endl ;


      //      int date = KeyVal("MJD-OBS");
      if (IsOfKind<Cfht12K>(*this) && date > 51430 && date <  51460) //HD problem only for the cfht99b campagn
	{
	  badValue = 0;
	  max_value=65535;
	}


      if (badValue != -1)
	{
	  for (Pixel *p = Image::begin() ; p<end; ++p) 
	    {
	      if (*p == badValue)
		{
		  n_pix_changed++;
		  *p = max_value;
		}
	    }
	  cout << " Applied corection for "<< TelInstName(*this) <<" for 0-pixels to " 
	       << n_pix_changed << " pixels. new value :  " << max_value << endl ;
	}

  }
#endif
}


FitsImage::FitsImage(const string &FileName, const FitsHeader &a_fits_header, const Image & an_image) :
   FitsHeader(a_fits_header, FileName),
   Image(an_image)
{ 
#ifdef DEBUG
  cout << " > FitsImage::FitsImage(string FileName, const FitsHeader &a_fits_header, const Image & an_image) " 
       << "FileName = " << FileName
       << endl;
#endif
ModKey("NAXIS1", nx);
ModKey("NAXIS2", ny);
 written = 0 ;
}

FitsImage::FitsImage(const string &FileName, const FitsHeader &a_fits_header) :
   FitsHeader(a_fits_header, FileName)
{ 
#ifdef DEBUG
  cout << " > FitsImage::FitsImage(string FileName, const FitsHeader &a_fits_header) " 
       << "FileName = " << FileName
       << endl;
#endif
allocate(int(KeyVal("NAXIS1")), int (KeyVal("NAXIS2")));
 written = 0 ;
}

FitsImage::FitsImage(const string &FileName, const Image& an_image) :
   FitsHeader(FileName,RW), Image(an_image)
{
#ifdef DEBUG
  cout << " > FitsImage::FitsImage(const string &FileName, const Image& an_image) " << endl;
#endif
ModKey("NAXIS1", nx);
ModKey("NAXIS2", ny);
 written = 0 ;
}


FitsImage::FitsImage(const string &FileName, const int Nx, const int Ny) :
  FitsHeader(FileName,RW), Image(Nx,Ny)
{
  ModKey("NAXIS1", Nx);
  ModKey("NAXIS2", Ny);
  written = 0 ;
}
  

#include <time.h>
static string  local_time()
{
time_t now = time(NULL);
//struct tm *date = localtime(&now);
char *ptime = ctime(&now);
// eventually remove \n at the end...
char *plast = ptime + strlen(ptime)-1;
if (*plast == '\n') *plast = '\0';
return string(ptime);
}


bool FitsImage::SetWriteAsFloat()
{
if (fileModeAtOpen !=  RW)
  {
    cerr << "ERROR: SetSaveAsFloat requested for RO file : " << fileName << endl; 
    return false;
  }
else ModKey("BITPIX",-32);
return true;
}

void FitsImage::PreserveZeros()
{
   AddOrModKey("KEEPZERO",true," actual 0's should be reloaded as 0 " );
}


#include <limits.h> /* for SHRT_MIN and SHRT_MAX */ 



int FitsImage::Write(bool force_bscale) 
{
  if (!writeEnabled) return 0;
  if (FileMode() != RW) return 0;
  int bitpix = KeyVal("BITPIX");
  if (bitpix == 0) bitpix = 16;

  cout << " Writing " << nx << "x" << ny << " " << fileName << endl;
  if (compressedImg && bitpix <0)
    {
      cout << " Rice compression for floating point images is messy : " << fileName << " will not be compressed " << endl;
      compressedImg=false;
    }

  if (compressedImg && !(bitpix == 16 || bitpix == 8 ))
    {
      std::cerr << " for image " << FileName() << " bitpix = " << bitpix 
		<< std::endl;
      std::cerr <<" only know how to write compressed image with BITPIX==[8,16] for now, no compression ... " << std::endl;
      compressedImg=false;
    }

  create_image(bitpix,nx,ny);
  double bscale = 1. ;
  double bzero = 0. ;
  Pixel mini,maxi;
  MinMaxValue(&mini,&maxi);
  bool int_data = HasKey("INT_DATA") && bool(KeyVal("INT_DATA"));

  if (bitpix == 16 && !int_data)
    {
      bscale = KeyVal("BSCALE");
      bzero = KeyVal("BZERO");
      if (HasKey("KEEPZERO") && bool(KeyVal("KEEPZERO")))
	{// compute Bscale and BZero so that 0 remain 0 (used for weight maps)
	  if (maxi != 0)
	    {
	      double short_max = SHRT_MAX;
	      double span = short_max;
	      bscale = maxi/span;
	      bzero = 0;
	    }
	}	
      else if (!force_bscale)
	{
	  /* CFITSIO DOES NOT COMPUTE automatically BSCALE and BZERO 
	     according to BITPIX in most of the cases */
	  double short_min =  SHRT_MIN +2 ;
	  double short_max =  SHRT_MAX -2 ;
	  double span = short_max -  short_min ;
	  bscale = (maxi - mini)/span;
	  bzero = mini - bscale * (short_min);
	  // preserve the value 0 in quantization:
	  //	  if (bscale !=0)
	  //	    bzero -= (bzero/bscale-floor(bzero/bscale))*bscale;
	}
    }

  if (bitpix == 8)
    { 
      if (mini<0 && maxi > 255)
	{ 
	  // untested code
	  cout << " We have to compress data to accomodate bitpix = 8 (you should consider bitpix=16)" << endl;
	  double span = 256;
	  bscale = (maxi-mini)/span;
	  bzero =  mini;
	}
      else if (mini<0 || maxi > 255)
	{// we assume we only have integer values : only use bzero
	  if (mini<0) bzero = mini;
	  if (maxi>255) bzero = maxi-255;
	}
    }
	  
  if (bscale == 0) bscale = 1;
  cout << "Min: " << mini << "  Max: " << maxi << endl ;
  cout << " with BITPIX=" << bitpix << " BSCALE=" << bscale << " BZERO=" << bzero << endl;
  int status = 0;

  string stime = local_time();
  char * time =  (char *) stime.c_str();
  AddOrModKey("WRITEDAT", time, " when this file was written"); 
  /* WRITEDAT used as a tag to find if an image was ever read by this software 
     (which corrects when reading senseless ADC values) */
  
  
  fits_set_bscale(fptr, bscale, bzero, &status);
  CHECK_STATUS(status, " set_bscale ", );
  AddOrModKey("BSCALE", bscale);
  AddOrModKey("BZERO" ,bzero);  
  ModKey("NAXIS1", nx);
  ModKey("NAXIS2", ny);
  

/* This "Flush()" is because new values of BSCALE and BZERO are not
  read (and thus not used for writing) by fits_write_img. The
  diagnostic is that when reloaded, the image has a different mean and
  sigma. Still true with v2r440.
 */

  Flush();

  status = write_image(0,nx*ny,data, nx, bitpix, bscale, bzero);

  written = 1 ;
  return (!status);
} 

FitsImage::~FitsImage()
{
  // cout << " FitsImage destructor for " << FileName() << " mode : " << FileModeName(file_mode(fptr)) << endl;


if (IsValid() && FileMode() == RW)
  {
    if (written == 0)
      Write(false);
  }
// fits_close_file(fptr, &status); in FitsHeader destructor;
}

bool FitsImage::Trim(const Frame &Region)
{
  cout << "Trimming the overscan region of " << FileName() << endl;

 int Nx_Image = int(Region.Nx());
 int Ny_Image = int(Region.Ny());
 int X_0 = int(Region.xMin);
 int Y_0 = int(Region.yMin);
  
 Frame wholeFrame(*this,WholeSizeFrame);
 // to be able to extract a subimage, Region should be inside the Image:
 if (Nx() == Nx_Image && Ny() == Ny_Image || (Region*wholeFrame != Region))
   {
     cout << "Image " << FileName () << " already trimed : nothing done." << endl ;
     return false;
   }
 *((Image *) this) = Subimage(X_0,Y_0,Nx_Image,Ny_Image);
 // Correct the WCS if any:
 // may be we should check for CRPIX1_<n> (following WCS recommandations)
 if (HasKey("CRPIX1") && HasKey("CRPIX2"))
   {
     if (FileMode() == RW)
       {
	 cout << "Updating WCS of " << FileName() << " to account for trim " << endl;
	 double crpix1 = double(KeyVal("CRPIX1")) - X_0;
	 ModKey("CRPIX1", crpix1);
	 double crpix2 = double(KeyVal("CRPIX2")) - Y_0;
	 ModKey("CRPIX2", crpix2);
       }
     else
       cout << " FitsImage::Trim : trimming a FITS image opened RO " << FileName() << endl;
   }


 if (HasKey("DATASEC")) {
   if (FileMode() == RW) {
     char datasec[100];
     sprintf(datasec,"[1:%d,1:%d]",Nx(),Ny());
     ModKey("DATASEC",datasec,"Modified by Toads::FitsImage::Trim");
   }else{
     cout << " FitsImage::Trim : trimming a FITS image opened RO " << FileName() << endl;
   }

 }
 return true;
}






//! copy an image, with possible (de)compression. No bitpix change enabled, because it requires some thinking.

/* This routine is the main of imcopy from cfitsio sources, however with 
   a key modification : no compression of fp images: Rice can be a disaster there */


int ImageCopy(const std::string &InFileName, std::string OutFileName,
	      const bool NoCompressionOnOutput)
{

  if (FileExists(OutFileName)) remove(OutFileName.c_str());
    fitsfile *infptr, *outfptr;   /* FITS file pointers defined in fitsio.h */
    int status = 0, ii = 1, iteration = 0, single = 0, hdupos;
    int hdutype, bitpix, bytepix, naxis = 0, nkeys, datatype = 0, anynul;
    long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    long first, totpix = 0, npix;
    double *array, bscale = 1.0, bzero = 0.0, nulval = 0.;
    char card[81];
    // addition to imcopy
    {
      //capture the bitpix for first image
      FitsHeader head(InFileName);
      bitpix = head.KeyVal("BITPIX");
    }

    /* Open the input file and create output file */
    fits_open_file(&infptr, InFileName.c_str(), READONLY, &status);
    if (FileExtension(OutFileName) == "fz" && !NoCompressionOnOutput)     // addition to imcopy
      {
	if (bitpix>=0) // means integer  encoding
	  OutFileName = OutFileName+"[compress]";
	else
	  cout << " ImageCopy : the file " << OutFileName  << " will not be compressed because input file is not integer" << endl;
      }
    
    fits_create_file(&outfptr, OutFileName.c_str(), &status);

    if (status != 0) {    
        fits_report_error(stderr, status);
        return(status);
    }

    fits_get_hdu_num(infptr, &hdupos);  /* Get the current HDU position */

    /* Copy only a single HDU if a specific extension was given */ 
    if (hdupos != 1 || strchr(InFileName.c_str(), '[')) single = 1;

    for (; !status; hdupos++)  /* Main loop through each extension */
    {

      fits_get_hdu_type(infptr, &hdutype, &status);

      if (hdutype == IMAGE_HDU) {

          /* get image dimensions and total number of pixels in image */
          for (ii = 0; ii < 9; ii++)
              naxes[ii] = 1;

          fits_get_img_param(infptr, 9, &bitpix, &naxis, naxes, &status);

          totpix = naxes[0] * naxes[1] * naxes[2] * naxes[3] * naxes[4]
             * naxes[5] * naxes[6] * naxes[7] * naxes[8];
      }

      if (hdutype != IMAGE_HDU || naxis == 0 || totpix == 0) { 

          /* just copy tables and null images */
	if (! (hdutype == IMAGE_HDU && naxis == 0)) // no .. not NULL Images (addition to imcopy)
	      fits_copy_hdu(infptr, outfptr, 0, &status);

      } else {

          /* Explicitly create new image, to support compression */
          fits_create_img(outfptr, bitpix, naxis, naxes, &status);

          /* copy all the user keywords (not the structural keywords) */
          fits_get_hdrspace(infptr, &nkeys, NULL, &status); 

          for (ii = 1; ii <= nkeys; ii++) {
              fits_read_record(infptr, ii, card, &status);
              if (fits_get_keyclass(card) > TYP_CMPRS_KEY)
                  fits_write_record(outfptr, card, &status);
          }

          switch(bitpix) {
              case BYTE_IMG:
                  datatype = TBYTE;
                  break;
              case SHORT_IMG:
                  datatype = TSHORT;
                  break;
              case LONG_IMG:
                  datatype = TLONG;
                  break;
              case FLOAT_IMG:
                  datatype = TFLOAT;
                  break;
              case DOUBLE_IMG:
                  datatype = TDOUBLE;
                  break;
          }

          bytepix = abs(bitpix) / 8;

          npix = totpix;
          iteration = 0;

          /* try to allocate memory for the entire image */
          /* use double type to force memory alignment */
          array = (double *) calloc(npix, bytepix);

          /* if allocation failed, divide size by 2 and try again */
          while (!array && iteration < 10)  {
              iteration++;
              npix = npix / 2;
              array = (double *) calloc(npix, bytepix);
          }

          if (!array)  {
              printf("Memory allocation error\n");
              return(0);
          }

          /* turn off any scaling so that we copy the raw pixel values */
          fits_set_bscale(infptr,  bscale, bzero, &status);
          fits_set_bscale(outfptr, bscale, bzero, &status);

          first = 1;
          while (totpix > 0 && !status)
          {
             /* read all or part of image then write it back to the output file */
             fits_read_img(infptr, datatype, first, npix, 
                     &nulval, array, &anynul, &status);

             fits_write_img(outfptr, datatype, first, npix, array, &status);
             totpix = totpix - npix;
             first  = first  + npix;
          }
          free(array);
      }

      if (single) break;  /* quit if only copying a single HDU */
      fits_movrel_hdu(infptr, 1, NULL, &status);  /* try to move to next HDU */
    }

    if (status == END_OF_FILE)  status = 0; /* Reset after normal error */

    fits_close_file(outfptr,  &status);
    fits_close_file(infptr, &status);

    /* if error occurred, print out error message */
    if (status)
       fits_report_error(stderr, status);
    return(status);
}


/*
std::string DecompressImageIfNeeded(const std::string &InFileName,
				    const std::string &UniqueString,
				    const std::string &TempDir,
				    std::string &ToRemove)
*/
std::string DecompressImageIfNeeded(const std::string &InFileName,
				    const std::string &OutFileName,
				    std::string &ToRemove)
{// uncompress images for sextractor

  // zipped?
  if (IsZipped(InFileName))
    {
      remove(OutFileName.c_str());
      std::string command = "gunzip -c "+InFileName+" > "+OutFileName;
      ToRemove += " "+OutFileName;
      //DEBUG
      cout << "executing " << command << endl;
      if (0 != system(command.c_str()))
	{
	  std::cerr << " could not execute " << std::endl
		    << command << std::endl;
	  return "";
	}
      else return OutFileName;
    }
  else
    {
      FitsHeader head(InFileName);
      if (head.CompressedImage())
	{
	  ToRemove += " "+OutFileName;
	  //DEBUG
	  cout << "executing ImageCopy(" << InFileName <<','<< OutFileName 
	       << ')' << std::endl;
	  if (0 == ImageCopy(InFileName, OutFileName, true))
	    return OutFileName;
	  else return "";
	}
      else return InFileName;
    }
}


bool DecompressOrLinkImage(const std::string &InFileName,
			   const std::string &OutFileName)
{
  //if the file exists, nothing is gone work.
  remove(OutFileName.c_str());
  // do we have do decompress?
  std::string dummy;
  std::string decompImage = DecompressImageIfNeeded(InFileName, OutFileName, 
						    dummy);
  if (decompImage == "") return false; // it did not work;
  if (decompImage == InFileName) // no decompression
    {
      return (MakeRelativeLink(InFileName, OutFileName));
    }
  return true;
}






/*   Notes on using compressed images with cfitsio,
and the choosen call sequnce to handle them transparently.


cfitsio proposes to compress the pixel part using elaborate
algorithms. One important capability of the approach is
to tile images prior to compression, in order to be able
to retrieve subimages withour uncompressing the whole image.

From the very beginning, the FitsImage class is intended for
files that contain a single image (a single pixel array).

One should be aware that compressed images are stored as a binary
table, which cannot be in the primary HDU, by FITS regulations. So, a
single compressed image is a fits file with an extension, the primary
header being mainly dummy. WE decided to hide this difference between
uncompressed and compressed images by "skipping" this dummy header
for compressed images (i.e., when opening an existing file, we find
a main header describing a zero pixel image, we go to next HDU).
When we create a file for a compressed image. Using a binary
table has other drawbacks : NAXIS, NAXIS1 and NAXIS2 have other meanings
than the usual one, and the image sizes, and bitpix go into 
ZNAXIS, ZNAXIS1, ZNAXIS2, ZBITPIX. We decided to operate a "key translation",
i.e, reading or writing NAXIS1 in a compressed image adresses ZNAXIS1.


Cfitsio documentation proposes the compressing scheme and describes it
as "transparent". I do not fully agree. In fact very few
functionnalities of the library work with compressed images.

In cfitsio, the user is supposed to have its data somewhere,
to process it and eventually write the results. Although
it is what most of our programs do, the FitsImage class
does not enforce any order in which the data should be accessed,
and altered.

  To be more precise, when we create a file to write a fitsimage in,
we do not necessarily know its size when we go through the
constructor. The size of the written image has to be the actual size at the
time of actually writing (triggered by the destructor, most of the
 time).

  In the implementation choosen long ago, for sake of simplicity,
the idea was to let cfitsio handle all things that concern the fits
headers : when you read or write a fits key, you go through cfitsio 
routines. There is no image of the header on toads side. For the image data,
it is loaded in memory, and eventually written to disk. as a side remark,
if FitsImage reads the whole image in a single chunk, we have other classes
which enable to traverse monster images by slices.

  So the standard cfitsio interface for a FitsImage (opened RW, the ones 
  opened RO don't cause any stress )  is
1-  open file (possibly create it)
2-  read image if any
3- alter header keywords, (or create it if empty)
4- alter pixels, possibly change size of image (in memory)
5- write pixels
6- close file

This kind of sequence works for uncompressed images, in a rather
straighforward way:
1- fits_open_file(existing file) or fits_create_file
2- fits_read_img
3- fits_{read,write,delete}_key, and 10 dozens of other routines
4- our own Image code stuff.
5- 1 tell sizes to cfitsio 
       - write BITPIX, NAXIS1, MAXIS2 in header
       - or call fits_resize_img
       - or call fits_create_img, if there is no image yet.

   2 then fits_write_img does the job


For compressed images, this scheme does not work, for various reasons,
not all understood:
  fits_resize_img does not accept compressed images.
  fits_create_img does the whole job of creating all the needed data
for compression. It actually creates a new image extension,
and header keys will go into this extension. But to be able to
alter a given image, we have to be able to write a compressed image
in a file that was created by an other program. This is where
fits_resize_img would be very useful. I even tried to go
into undocumented routines, and only got crashes. I though
I had figured out what are the calls in fits_create_img that actually
load the needed structures with the actual sizes, but it did not work.

 I then tried to delete the current HDU, and call fits_create_img.
for an obscure reason, it does not work either.

  So, how does it work? very simple: remark that the main header is dummy.
The actual image header (the one that contains our "USER DATA",
in the first extension) is very small compared to pixels.
So when we want to write, we just 
 - collect all header keys
 - delete the file (fits_delete_file)
 - create it (fits_open_file, and tell cfitsio that we want compression)
 - create the image extension (fits_create_img).
 - copy header keys back.
 - write pixels.

Isn't it nice? In fact this is exactly what is done in imcopy
(the utility/example provided with cfitsion library), which only
works if the target file does not exist. If you happen to find
a nicer way to change existing compressed images, please drop a mail to
pierre.astier@in2p3.fr
*/

/*  SCALING and TYPE CONVERSION
    
  Second difficulty with compressed images, is the conversion between
the format used on disk an the one in memory. Typically, we use short
int (BITPIX = 16) on disk and float in memory. It seems that when
writing or reading compressed images, there is a confusion between the
data type the user wishes and the one that is on disk.

   To get back to imcopy, one should remark that in this code, the user
format coincides with the disk format. Furthermore, bscale and bzero
are ignored, in order to make sure that raw pixel values are copied.
as a trial, we changed imcopy so that requested data in float (for a
bitpix = 16 image, there is no information loss), and let the
automatic scaling do its job: Everything works for non compressed
images, and fails for compressed ones (i.e. input and ouput differ).
This means that the cfitsio code for compressed images only works for
user type identical to disk type, and data rescaling is not properly
done.

This is why in FitsHeader::read_image and FitsImage::Write, we have
those 2 operations (convertion + rescaling).  We wrote
FitsHeader::read_image to be able to read subimages, because there are
many places not in this source file where fits_read_img is called, to
read a subimage (e.g. in FitsSlice). In order to accomodate compressed
images, we have to use our own wrapping of fits_read_img.


P.A 24/06/04

Note from  March 2010:

I tried again with cfitsio v3r006, and the fact that type conversion, rescaling, and 
compression do not work all at once when writing is still true. So 
FitsHeader::read_image and FitsHeader::write_image are still needed.

I tried with version 3.24 (from jan 2010) and the bug seems to be gone.
Might be worth to upgrade and simplify our code here.

*/


