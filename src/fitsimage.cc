#include <iostream>
#include <fitsio.h>
#include "fitsimage.h"
#include "fitstoad.h" 
#include "fileutils.h"
#include "frame.h"
#include "stringlist.h"

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



#define CHECK_STATUS(status, message, return_statement) \
if (status)\
  {\
  cfitsio_report_error(status, message);\
  cerr << " for file " << this->fileName << endl;\
  return_statement;\
  }

static void cfitsio_report_error(int status, const string &Message)
{
cerr << Message ;
fits_report_error(stderr, status);
}

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


FitsHeader::FitsHeader(const string &FileName, const FitsFileMode Mode, 
		       bool EmptyFile)
{
int status = 0;
fileName = FileName;
fileModeAtOpen = Mode;
fptr = 0;
  /* if one wants an empty file, delete any existing file unless
     cfitsio will append the data to it. This forbids to create with
     this software files with several images. */
if (EmptyFile) remove(FileName.c_str()); 

#ifdef DEBUG
 cout <<  " FitsHeader::FitsHeader(string FileName =" << FileName << " , int Mode =" << Mode << " ) " << endl;
#endif /* DEBUG */ 
 telInst = NULL;
 writeEnabled = true; // only relevant if mode != RO
 
 string temp_file_name=FileName;
 fits_open_file(&fptr, temp_file_name.c_str(), int(Mode), &status);
 /* there seems to be a bug in cfitsio, when one requests some
    preprocessing before actually accessing the data, e.g. adressing
    subimages through actual_filename[imin:imax,jmin:jmax]. Then the
    returned pointer has RW mode although RO was requested. this is
    handled using fileModeAtOpen, rather than trying to hack cfitsio
    internals
 */

 if (status == 0) return;
 if (Mode == RO) // no need to try to create the file
   CHECK_STATUS(status,"FitsHeader", return);
 
 // form here on, we are in RW mode and could not open an existing file
 int first_status = status; // for eventual error reporting
 status = 0;
 fits_create_file(&fptr, temp_file_name.c_str(), &status);   
 if (status)
   {
     cerr << " when opening file : " << fileName;
     cerr << ", with mode " << FileModeName(Mode) <<  ',' << endl;
     cerr << " we got these successive errors:" << endl;
     fits_report_error(stderr, first_status);
     fits_report_error(stderr, status);
   }
 else
   {
     if (!EmptyFile) minimum_header();
   }
}

void FitsHeader::minimum_header()
{
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

FitsHeader::FitsHeader(const FitsHeader &Header)
{int status =0;
  telInst = 0;
  writeEnabled = Header.writeEnabled; // probably irrelevant
  fileName = Header.fileName;
  fits_open_file(&fptr, Header.fileName.c_str(), RO, &status);
  fileModeAtOpen = RO;
  CHECK_STATUS(status," FitsHeader ", return );
}


FitsHeader::FitsHeader(const FitsHeader &Template, 
		       const string &NewFileName)
{
int status = 0;
#ifdef DEBUG
 cout << " FitsHeader::FitsHeader(const FitsHeader &Template= " << Template.fileName << 
   ", string &NewFileName = " << NewFileName << " ) " << endl;
#endif
fileName = NewFileName;
remove(NewFileName.c_str());
fits_create_file(&fptr, NewFileName.c_str(), &status);
fileModeAtOpen = RW;
fits_copy_header(Template.fptr, fptr, &status);
/* the following Flush has the effect that the produced image file
cannot be smaller than the one associated with Template. this is disastrous.
calling fits_resize_img does not help. */
/* Flush(); */
 telInst = 0;
 writeEnabled = true;
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
     if (FileMode() == RW && !writeEnabled)
       {
	 cout << "deleting file " << fileName << endl;
	 fits_delete_file(fptr,&status);
       }
     else fits_close_file(fptr, &status);
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

FitsKey 
FitsHeader::KeyVal(const string &KeyName, const bool Warn) const
{
  if (KeyName.find("TOAD") != KeyName.npos && !HasActualKey(KeyName)) 
    {return ToadsKeyVal(*this,KeyName,Warn);}
  else
    return FitsKey(KeyName, fptr, Warn);
}

FitsKey 
FitsHeader::KeyVal_Secure(const string &KeyName) const
{
  if (HasKey(KeyName)) 
    return KeyVal(KeyName);
  else
    {
      cerr << "No key " << KeyName << " in header " << fileName <<  endl ;
      exit(0);
    }
}




// has to be a member function because CHECK_STATUS uses this->fileName
int FitsHeader::mod_key(int type, const string &KeyName, void *Value, const string &Comment) const
{
if (fptr)
  {
  int status = 0;

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



int FitsHeader::write_key(const string &KeyName, void* KeyVal, const string &Comment, const int type) 
{
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

bool FitsHeader::HasKey(const string &KeyName, const bool Warn) const
{
  
if (HasActualKey(KeyName, Warn)) return true;
if (KeyName.find("TOAD") != KeyName.npos) 
  {
  string value = ToadsKeyVal(*this, KeyName, Warn);
  if (value != NOVAL) return true;
  }
return false;
}


int FitsHeader::RmKey(const string &KeyName) const
{
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
  char card[512];
  int status = 0;
  char keyName[80];
  sprintf(keyName, "%s", KeyName.c_str());
  // copy everything verbatim (name ,value, comment)
  fits_read_card(fptr, keyName, card, &status);
  fits_write_record(To.fptr, card, &status);
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

/*****************  FitsImage ****************/

FitsImage::FitsImage(const string &FileName, const FitsFileMode Mode) : FitsHeader(FileName, Mode) , Image()
{
#ifdef DEBUG
  cout << " FitsImage::FitsImage(string FileName=" << FileName << ", int Mode=" << Mode << ") " << endl;
#endif
  if (!IsValid()) return; 
nx = KeyVal("NAXIS1");
ny = KeyVal("NAXIS2");
Image::allocate(nx,ny);
float nullval = 0;
int anynull;
int status = 0;
written = 0 ;
fits_read_img(fptr, TFLOAT, 1, nx*ny, &nullval,  data, &anynull, &status);

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
CHECK_STATUS(status," FitsImage ", );
}

FitsImage::FitsImage(const FitsHeader &Header) : FitsHeader(Header), Image()
{
#ifdef DEBUG
  cout << " FitsImage::FitsImage(string FileName=" << FileName << ", int Mode=" << Mode << ") " << endl;
#endif 
nx = KeyVal("NAXIS1");
ny = KeyVal("NAXIS2");
Image::allocate(nx,ny);
float nullval = 0;
int anynull;
int status = 0;
 written = 0 ;
fits_read_img(fptr, TFLOAT, 1, nx*ny, &nullval,  data, &anynull, &status);
CHECK_STATUS(status," FitsImage ", );
}

FitsImage::FitsImage(const string &FileName, const FitsHeader &a_fits_header, const Image & an_image) :
   FitsHeader(a_fits_header, FileName),
   Image(an_image)
{ 
#ifdef DEBUG
  cout << "  FitsImage::FitsImage(string FileName, const FitsHeader &a_fits_header, const Image & an_image" << endl;
#endif
ModKey("NAXIS1", nx);
ModKey("NAXIS2", ny);
 written = 0 ;
}

FitsImage::FitsImage(const string &FileName, const FitsHeader &a_fits_header) :
   FitsHeader(a_fits_header, FileName)
{ 
#ifdef DEBUG
  cout << "  FitsImage::FitsImage(string FileName, const FitsHeader &a_fits_header, const Image & an_image" << endl;
#endif
allocate(int(KeyVal("NAXIS1")), int (KeyVal("NAXIS2")));
 written = 0 ;
}

FitsImage::FitsImage(const string &FileName, const Image& an_image) :
   FitsHeader(FileName,RW), Image(an_image)
{
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
  int status = 0;
  cout << " Writing " << nx << "x" << ny << " " << fileName << endl;
  if (FileMode() != RW) return 0;
  int bitpix = KeyVal("BITPIX");
  double bscale = 1. ;
  double bzero = 0. ;

  if (bitpix == 16)
    {
      if (HasKey("KEEPZERO") && bool(KeyVal("KEEPZERO")))
	{// compute Bscale and BZero so that 0 remain 0 (used for weight maps)
	  Pixel min,max;
	  MinMaxValue(&min,&max);
	  if (max != 0)
	    {
	      double short_max = SHRT_MAX;
	      double span = short_max;
	      bscale = max/span;
	      // this in fact assumes that the content is >=0
	      bzero = 0;
	      //	    cout << bzero/bscale << endl;
	    }
	}	
      else if (!force_bscale)
	{
	  /* CFITSIO DOES NOT COMPUTE automatically BSCALE and BZERO 
	     according to BITPIX in most of the cases */
	  Pixel min,max;
	  MinMaxValue(&min,&max);
	  double short_min =  SHRT_MIN +2 ;
	  double short_max =  SHRT_MAX -2 ;
	  double span = short_max -  short_min ;
	  bscale = (max - min)/span;
	  bzero = min - bscale * (short_min);

	}
    }

  cout << " with BITPIX=" << bitpix << " BSCALE=" << bscale << " BZERO=" << bzero << endl;

  AddOrModKey("BSCALE", bscale);
  AddOrModKey("BZERO" ,bzero);  
  ModKey("NAXIS1", nx);
  ModKey("NAXIS2", ny); 

  string stime = local_time();
  char * time =  (char *) stime.c_str();
  AddOrModKey("WRITEDAT", time, " when this file was written"); 
  /* WRITEDAT used as a tag to find if an image was ever read by this software 
     (which corrects when reading senseless ADC values) */

/* This "Flush()" is because new values of BSCALE and BZERO are not
  read (and thus not used for writing) by fits_write_img. The
  diagnostic is that when reloaded, the image has a different mean and
  sigma. Still true with v2r440.
 */

 Flush();

 /* actually write it */
  fits_write_img(fptr, TFLOAT, 1, nx*ny, data, &status);
  CHECK_STATUS(status," WriteImage ", );
  written = 1 ;
  return (!status);
} 

int FitsImage::Write(const double &Bscale, const double &Bzero) 
{
  if (!writeEnabled) return 0;
  int status = 0;
  cout << " Writing " << nx << "x" << ny << " " << fileName << endl;
  if (FileMode() != RW) return 0;
  int bitpix = KeyVal("BITPIX");
  cout << " with BITPIX=" << bitpix << " BSCALE=" << Bscale << " BZERO=" << Bzero << endl;
  AddOrModKey("BSCALE", Bscale);
  AddOrModKey("BZERO" ,Bzero);  
  ModKey("NAXIS1", nx);
  ModKey("NAXIS2", ny); 

  string stime = local_time();
  char * time =  (char *) stime.c_str();
  AddOrModKey("WRITEDAT", time, " when this file was written"); 
  /* WRITEDAT used as a tag to find if an image was ever read by this software 
     (which corrects when reading senseless ADC values) */

/* This "Flush()" is because new values of BSCALE and BZERO are not
  read (and thus not used for writing) by fits_write_img. The
  diagnostic is that when reloaded, the image has a different mean and
  sigma. Still true with v2r440.
 */

 Flush();

 /* actually write it */
  fits_write_img(fptr, TFLOAT, 1, nx*ny, data, &status);
  CHECK_STATUS(status," WriteImage ", );
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

void FitsImage::Trim(const Frame &Region)
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
     return;
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
}
