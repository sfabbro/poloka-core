#include <iostream>
#include <cmath>
#include <algorithm>
#include <map>

#include "fileutils.h"  // StringTo{Lower,Upper}
#include "fitsimage.h"
#include "polokaexception.h"
#include "astroutils.h" /* for conversion routines */
//#include "../src/gtransfo.h"

#include "frame.h"

#ifndef M_PI
#define     M_PI            3.14159265358979323846  /* pi */
#endif

/*! \file 
  \brief TOAD keys and fits header "standardization". */

/*! \page fitstoad_long Header Standardization

We explain here how headers are presented with a uniform interface to the
other parts of the code.

\subsection fitstoad_implementation  Implementation description and usage of fitstoad

  fitstoad.cc is the main file for presenting Fits headers in a uniform fashion
to the rest of the code.

  Seb implemented a first version of these functionnalities. Although it
did the job, it was hard to scale it up. It also missed a default
behavior, functionnality per functionnality.

  The (2nd generation) choosen implementation is the following:
We have a VirtualInstrument class that interfaces the code
to actual instrument classes. This VirtualInstrument class
has about 20 (virtual) routines, most of them being the translation
of the famous TOAD keys. I added a few missing routines (which you
could locate as big switches in the "old" code, in usnoutils.cc and
superflat.cc).
  fitstoad.cc contains the default implementation (VirtualInstrument)
and the interface to external world. To make the things really tight,
I even did not put the VirtualInstrument class in fitstoad.h . I
realize that this is arguable, but enforces a rigorous design.
  Then, the actual instruments are in a separate directory,
and if you remove a file from this directory, your code still
compiles and an instrument has disappeared. The consequence is that adding
an instrument involves the coding of one file, and NO MODIFICATION elsewhere.
This was the baseline of the new design. We also allow to use all C++ resources
about class derivation you may want to use. For example, if some changes
in the DAQ of a telescope involve minor changes, you should
implement them as derived class of something that already worked.
Implement the derived class in the same file as the base class
you want to derive it from, in order to control the ordering
in which the main "type" switch tries to match a given header
to the provided classes.
  For most of your actual headers, the default implementation 
(The VirtualInstrument routines) will provide the right behavior.
I coded a test routine called if your require the pseudo key TOADALL.
(header -k TOADALL your_fitsfile invokes it)


    The provided instruments are in the directory ../telinst/ 
*.cc in ../telinst are the files included in fitstoad.cc.
Adding files there does not require modifying the Makefile
(or mgr/requirements file). Some code is automatically
generated from parsing these files (using grep and sed)
and generates alltelinst.cc and alltelinst.h
The "code generator" is a makefile section located at the end of
fitstoad.cc, and extracted by the main Makefile.

This code generator collects the Acceptor routines,
and put all them in an array (AcceptorsArray). SniffTelInst tries all
components of the array and keeps the first match.
An Unknown class (placed at the end of the array)
ensures that there will be a successful match. 
You have a single way to control the ordering of
trials for related instruments. For classes in the
same file, "Acceptor" routines will be called
from last to first. If you derive an existing class,
you should put the derived class after the base class,
in the SAME file.

\subsection newtelinst Coding a new "tel/inst"
Here are the guidelines for the actual implementation of a new tel/inst.
You may use already coded stuff as an example.
There are several related instruments coded using derived classes 
in intwfc.cc.
you should first try a "header -k TOADALL <samplefile.fits>"
to see what does not work for your case(s) in the default implementation.

To code <yourtelinst>.cc:
- the best thing to do is to copy/paste one file that exists.
- in general you do not need include files: your code will be included
in fitstoad.cc after the utility routines you may need. If you refer
to other toads parts (such as wcsutils), include the correct include file.

- your class should derive (directly or not) from VirtualInstrument.
- In most of the cases, your class will contain an Acceptor routine,
i.e. the routine that indicates that
a given header belongs to the class you are coding. In this
likely case, put a comment containing 
    TYPE_SELECTOR on the same line as the class declaration.
   See already coded routines for the signature of Acceptor.
   DO NOT USE "TOAD*" keys IN THE CODE OF THE Acceptor ROUTINE.
  You may easily figure out why!

- you have to define the 3 name routines (TelInstName, TelName, InstName)
-  you have to define the toadkey translators for the keys which do not
  work by default. Use extensively the SIMPLE_TRANSLATOR and RETURN_A_VALUE
macros. They protect you from trivial bugs, and will ease any future changes,
if needed.

Not all TOAD keys are necessary for every task.
TOADPIXS is important: return something correct.
TOADEXPO is easy to map, do it if default does not apply.
There are routines provided to decode dates.


\subsection functions_for_flatfielding   Routines used  for flatfielding

  There are 4 routines related to flafielding:
\code
  Frame OverscanRegion(const FitsHeader &Head, const int Iamp) const;
  Frame TotalIlluRegion(const FitsHeader &Head) const;
  Frame IlluRegion(const FitsHeader &Head, const int Iamp) const;
  Frame AmpRegion(const FitsHeader &Head, const int Iamp) const;
\endcode

Do not bother with those if you do not plan to flatfield with TOADS.

\subsection guess_wcs  Guessing a rough WCS.

  The GuessLinWCS is used by matchusno, to figure which patch of the 
usno catalog to use. If you have DbImages.
try "matchusno -n <your_db_image>" to see if it works by default.
The default WCS just assumes north-up east-left and that
RA and DEC refer to the center of the image (they are read
via RaDec2000, which cooks up values from TOADRASC, TOADDECL and TOADEQUI)

\subsection test_fitstoad Testing how it works
To test your implementation use :
header -k TOADALL <your_fitsfile>

\subsection type_tester Type testers:

If you want to generate a  "type tester" for your class,
because there are places in the code where you would
like to put some instrument dependent code,
you have to place a comment containing TYPE_TESTER
on the line that declares your class (as for 
TYPE_SELECTOR)

  You then have a routine
\code
bool IsOfKind<YourBeautifulInstrument>(const FitsHeader &) 
\endcode
which is generated. (instanciated in the C++ jargon)

  you may now use it (with #include "fitstoad.h"):
\code
  if (IsOfKind<YourBeautifulInstrument>(head))
    {
       .... // specific processing
    }
\endcode

You should however try not to resort to such a scheme,
because it is breaking the instrument abstraction.
This may just mean that Virtual Instrument misses 
one functionnality, that you may consider addding.

Pierre.

*/
  





// handy routines:

static bool CheckKeyToUpper(const FitsHeader &Header, const string &KeyName, const string &ExpectedVal)
{
  return (Header.HasKey(KeyName.c_str()) && (StringToUpper(Header.KeyVal(KeyName.c_str())) == ExpectedVal));
}


static bool CheckKey(const FitsHeader &Header, const string &KeyName, const string &ExpectedVal)
{
  return (Header.HasKey(KeyName.c_str()) && (string(Header.KeyVal(KeyName.c_str())) == ExpectedVal));
}


static string end_of_string( const string &S)
{
  return S.substr(S.length()-1,1);
}

static bool fits_imregion_to_frame(const string &FitsStuff, Frame &ToadsStuff)
{
  int x0,x1,y0,y1;
  if (sscanf(FitsStuff.c_str(),"[%d:%d,%d:%d]",&x0,&x1,&y0,&y1) == 4)
    {
      if (x0>x1) swap(x0,x1);
      if (y0>y1) swap(y0,y1);
      ToadsStuff = Frame(Point(x0-1,y0-1), Point(x1-1, y1-1));
      return true;
    }
  else
    {
      cout << " ERROR : cannot decode 4 integers in " << FitsStuff << endl;
      ToadsStuff = Frame();
      return false;
    }
}




// Date hacking routines

#ifdef FUTURE
// code this hopping anyway that it will remain useless
// Call via e.g   decode_date(my_string, YY,MM,DD) or decode_date(my_string,DD,MM,YY)

typedef enum {YY=0,MM=1,DD=2} DateEnum;

bool string decode_date(const string &FitsDate,const DateEnum n1, const DateEnum n2, const DateEnum n3)
{
  int n[3];
  if ((sscanf(FitsDate.c_str(),"%d/%d/%d",&n[n1],&n[n2],&n[n3]) != 3) &&
      (sscanf(FitsDate.c_str(),"%d-%d-%d",&n[n1],&n[n2],&n[n3]) != 3)
      )
    return string;
  if (n[YY] < 1949) n[YY] += (n[YY] < 50) ?  2000 : 1900;
  char date_string[64];
  sprintf(date_string,"%d/%d/%d",n[DD],n[MM],n[YY]);
  return string(date_string);
}

#endif /*FUTURE */ 
  



/* we could write a routine that uses the fits comment to decode dates with 2 digits years.
   hope it is not needed. Seems that FITS standard forbids it for yy>=2000 */







#include "virtualinstrument.h"
#include "virtualinstrument.cc"
#include "alltelinst.h"


#ifdef INSTRUMENTWCS
//used to compute the rotation/flip of coordinates given North and East directions
// on the fits image. This provides if needed the 3rd argument of ComputeLinWCS.
typedef enum {Up,Down,Right,Left} AxisDir;

static GtransfoLin RotationFlip(const AxisDir NorthDir, const AxisDir EastDir)
{
  double a11 = 0;
  double a12 = 0;
  double a21 = 0;
  double a22 = 0;
  // tested with NorthDir=Down,EastDir=Left, and NorthDir=Left,EastDir=Down.
  // It should then work for other cases.
  switch (NorthDir)
    {
    case Up    : a22 =  1; break;
    case Down  : a22 = -1; break;
    case Right : a21 =  1; break;
    case Left  : a21 = -1; break;
    }
  switch (EastDir)
    {
    case Up    : a12 =  1; break;
    case Down  : a12 = -1; break;
    case Right : a11 =  1; break;
    case Left  : a11 = -1; break;
    }
  GtransfoLin rotFlip(0,0,a11,a12,a21,a22);
  if (fabs(rotFlip.Determinant()) != 1.)
    {
      cerr << " RotationFlip computes a non unitary transfo :" << endl 
	   << rotFlip << endl;
    }
  return rotFlip;
}

//! a handy routine to compute aWCS given a RaDec reference, and a possible rotation and flip.     
static bool ComputeLinWCS(const FitsHeader &Head, 
		     const Point &CrPix, 
		     const GtransfoLin &RotFlip, 
		     TanPix2RaDec &WCS)
{
  double pixscale = Head.KeyVal("TOADPIXS");
  if (pixscale == 0)
    {
      cerr << " NO TOADPIXS in file " << Head.FileName() << " : cannot guess a WCS " << endl;
      return false;
    }
  double ra,dec;
  RaDec2000(Head, ra, dec);
  
  GtransfoLin cd = GtransfoLinScale(pixscale/(3600), pixscale/3600.)
    *RotFlip
    *GtransfoLinShift(-CrPix.x, -CrPix.y);
  WCS = TanPix2RaDec(cd, Point(ra,dec));
  return true;
}


bool GuessLinWCS(const FitsHeader &Header, TanPix2RaDec &Guess)
{
  cout  << " trying default GuessLinWCS" << endl;
  return ComputeLinWCS (Head, Head.ImageCenter(), GtransfoIdentity(), Guess);

  //VirtualInstrument *p = Header.TelInst();   
 
  //if (p->GuessLinWCS(Header, Guess)) return true;
  
  // the tel/inst specific procedure failed. try the default one ...
  //.... but only if it is different 
  // the following test is refused by the compiler
  //  if (&(p->GuessLinWCS) == &(p->VirtualInstrument::GuessLinWCS)) return false;
  
  //return (p->VirtualInstrument::GuessLinWCS(Header, Guess));
}
#endif

// wrapper routines. They enable not to propagate the structure of TelInst
// and provide a simple user interface.

string TelescopeName(const FitsHeader &Head)
{
  return Head.TelInst()->TelName();
}  


string InstrumentName(const FitsHeader &Head)
{
  return Head.TelInst()->InstName();
}

string TelInstName(const FitsHeader &Head)
{
  return Head.TelInst()->TelInstName();
}


template<class Inst> bool IsOfKind(const FitsHeader &Head)
{
  return (dynamic_cast<Inst *>(Head.TelInst()) != NULL);
}

template<class Inst> bool IsOfKind(const string &FitsFileName)
{
  FitsHeader Head(FitsFileName);
  return (dynamic_cast<Inst *>(Head.TelInst()) != NULL);
}


Frame OverscanRegion(const FitsHeader &Head, const int iAmp)
{
  return Head.TelInst()->OverscanRegion(Head, iAmp);
}


Frame IlluRegion(const FitsHeader &Head, const int Iamp)
{
  return Head.TelInst()->IlluRegion(Head, Iamp);
}

Frame TotalIlluRegion(const FitsHeader &Head)
{
  return Head.TelInst()->TotalIlluRegion(Head);
}

Frame AmpRegion(const FitsHeader &Head, const int Iamp)
{
  return Head.TelInst()->AmpRegion(Head, Iamp);
}

double AmpGain(const FitsHeader &Head, const int iAmp)
{
  return Head.TelInst()->AmpGain(Head, iAmp);
}

typedef VirtualInstrument *(*AcceptorType)(const FitsHeader &);

/* 
   alltellinst .cc contains the implementation of all other (real!)
   instruments implementation of other instruments : 
   this file is generated by a makefile located at the end of this fitstoad.cc 
   Guidelines for the implementation of a new instrument : see 
   the (long) comment at the beginning of this file.
*/
#define VIRTUAL_INSTRUMENTS
#include "alltelinst.cc"

/* AcceptorsArray is declared in alltelinst.cc , 
generated by a makefile at the end of this source file */
static int NumberOfAcceptors = sizeof(AcceptorsArray)/sizeof(AcceptorsArray[0]);

#ifdef NEWTELINST
ostream& operator << (ostream &stream, const TelInstEnum &telInst)
{
  stream << " " << telescope_name(telInst) 
	 << " with " << instrument_name(telInst) << endl;
  return stream; 
}
#endif


#ifdef NEW_TELINST
// no trivial way to do that.... with the implementation we have
static void toadinstlist()
{
  cout << endl 
       << "List of Available instruments :   " << endl
       << "----------------------------------" << endl;
  for (int i=0; i<NumberOfTelInst; ++i)
    {
      cout << AvailableInstruments[i]->TelInstName << endl;
    }
}
#endif

VirtualInstrument  *SniffTelInst(const FitsHeader &Head)
{
  VirtualInstrument *p;
  /* scan reversly so that several acceptors in the same source file 
     are tried from last to first (derived classes before base classes).
     The alternative would have beeen to use "tac" in the code generator,
     but it is Linux specific... */
  for (int i=NumberOfAcceptors-2; i>=0; --i)
    if ((p=AcceptorsArray[i](Head)))
      {return p;}
  p = AcceptorsArray[NumberOfAcceptors-1](Head);
  if (p) return p;
  // should never happen
  cerr << " we miss a default class in SniffTelInst. Something went wrong " << endl;
  return new Unknown;
}

			     

typedef struct ToadsKeyRec {
  string Name;
  string Comment;
  string DefaultKey;
  FitsKey (VirtualInstrument::*key_translator)(const FitsHeader&, const bool) const;
};


// Be careful : I did not find any way to ensure that the documented default is the implemented one
static ToadsKeyRec ToadsKeys[] = {
  {"TOADPIXS", " PIXel Size in arcseconds. (double)", "PIXSCALE", &VirtualInstrument::TOADPIXS},
  {"TOADINST",  "INSTrument that recorded the image ", "Instrument Name", &VirtualInstrument::TOADINST},
  {"TOADTELE",  "TELEscope  ", "Telescope Name", &VirtualInstrument::TOADTELE},
  {"TOADFILT",  "FILTer description (string)", "FILTER", &VirtualInstrument::TOADFILT},
  {"TOADEXPO", "EXPOsure time in seconds (double)", "EXPTIME", &VirtualInstrument::TOADEXPO},
  {"TOADGAIN", " GAIN of the CCD in e-/ADU (double)", "GAIN", &VirtualInstrument::TOADGAIN},
  {"TOADRDON", " ReaDOut Noise of the electronics in e- (double)", "RDNOISE", &VirtualInstrument::TOADRDON},
  {"TOADRASC", " Right ASCension of the target in h:mn:sec(string)", "RA", &VirtualInstrument::TOADRASC},
  {"TOADDECL", " DECLination of the target in deg:mn:sec (string)", "DEC", &VirtualInstrument::TOADDECL},
  {"TOADEQUI", " EQUInox of the coordinates for the target  (double)", "EQUINOX", &VirtualInstrument::TOADEQUI},
  {"TOADAIRM", " AIRMass number for this exposure(double)" , "AIRMASS", &VirtualInstrument::TOADAIRM},
  {"TOADUTIM", " Universal TIMe at start of observation (double)", "UT", &VirtualInstrument::TOADUTIM},
  {"TOADDATE", " DATE of observation in day/month/year (string)", "DATE-OBS & DATE", &VirtualInstrument::TOADDATE},
  {"TOADMJD", " Mod julian date of obs (double)", "TOADMJD", &VirtualInstrument::TOADMJD},
  {"TOADMMJD", " DATE of observation in days since January first 2003 (double)", "TOADDATE", &VirtualInstrument::TOADMMJD},
  //{"TOADSCAN", " over/under SCAN region of the ccd [Xstart,Nx;Ystart,Ny](string)","OVERSCAN", &VirtualInstrument::TOADSCAN},
  //{"TOADILLU", " ILLUminated region of the ccd  [Xstart,Nx;Ystart,Ny](string)","NOVAL", &VirtualInstrument::TOADILLU},
  {"TOADBAND", "Short filter BAND (UBVRIJKLGZ) (string)","ToadBand(KeyVal(\"TOADFILT\"))", &VirtualInstrument::TOADBAND},
  {"TOADTYPE", " TYPE of data (OBJECT, DOMEFLAT, SKYFLAT, DARK, BIAS, FOCUS) (string)", "IMAGETYP", &VirtualInstrument::TOADTYPE},
  {"TOADOBJE", " OBJEct name-usually the field or the target name (string)", "OBJECT", &VirtualInstrument::TOADOBJE},
  {"TOADNAMP", " Number of AMPlifiers with one CCD (double)", "1", &VirtualInstrument::TOADNAMP},
  {"TOADCHIP", " CHIP number in case of multiple chips CCD (string)", "1", &VirtualInstrument::TOADCHIP},
  {"TOADPZPT", " CHIP number in case of multiple chips CCD (double)", "ZEROUSNO", &VirtualInstrument::TOADPZPT}
};


static int ToadsKeyCount  = sizeof(ToadsKeys)/sizeof(ToadsKeys[0]);

static void toadkeylist()
{
  cout  << "List of available TOADKEYS (8 characters):" << endl;
  for (int i=0; i<ToadsKeyCount; ++i)
    {
      ToadsKeyRec &key = ToadsKeys[i];
      cout << key.Name << " " << key.Comment << " " << key.DefaultKey << endl;
    }
  cout  << "-------------------------------------------------------------" << endl;
}

static void (*TestWCS)(const FitsHeader &Head) = NULL;

int SetTestWCS(void (*toto)(const FitsHeader &Head))
{
  TestWCS=toto;
  return 0;
}

static void TestTelInst(const FitsHeader &Head)
{
  cout << endl << " test of fitstoad  for header " << Head.FileName() << endl;
  VirtualInstrument *p = Head.TelInst();
  cout << "Telescope " << p->TelName() << " | Instrument " << p->InstName() << " | TelInst " << p->TelInstName() << endl;
  cout << " ------------------------------------------------" << endl;
  cout << " test of TOADKEYS " << endl;
  for (int i=0; i<ToadsKeyCount; ++i)
    {
      ToadsKeyRec &key = ToadsKeys[i];
      cout << key.Name << ": " << Head.KeyVal(key.Name.c_str(), true) << endl;
    }
  cout << " ------------------------------------------------" << endl;
  cout << " illu/bias regions " << endl;
  int namp = Head.KeyVal("TOADNAMP");
  for (int iamp=1; iamp<=namp; ++iamp)
    {
      Frame illu = IlluRegion(Head, iamp);
      Frame bias = OverscanRegion(Head,iamp);
      cout << " iamp " << iamp << " illu " << illu << " | overscan " << bias << endl;
      cout << " amp region once trimmed " << AmpRegion(Head,iamp);
    }
  cout << " total illu region " << TotalIlluRegion(Head) << endl;
  cout << " ------------------------------------------------" << endl;
  if (TestWCS) {
    TestWCS(Head);
  }
  cout << " ------------------------------------------------" << endl;
}


static int key_rank(const string &KeyName)
{
  for (int i=0; i < ToadsKeyCount; ++i)
    if (ToadsKeys[i].Name == KeyName) return i;
  return ToadsKeyCount;
}



FitsKey ToadsKeyVal(const FitsHeader &Head, const string &KeyName, const bool Warn)
{
  if (KeyName=="TOADHELP") {toadkeylist(); return FitsKey(KeyName,NOVAL);}
  if (KeyName=="TOADALL") {TestTelInst(Head); return FitsKey(KeyName,NOVAL);}
  

  VirtualInstrument *p = Head.TelInst();
  int r = key_rank(KeyName);
  if (r==ToadsKeyCount) // undefined ToadKey
    {
      cerr << " NO ToadsKey called " << KeyName << " for file " << Head.FileName() << endl;
      return FitsKey("KeyName",NOVAL);
    }
  return (p->*ToadsKeys[r].key_translator)(Head,Warn);
}



// what follows is extracted by make and executed to generate alltelinst.cc and alltelinst.h
// you cannot put any make comments because this triggers a C preprocessor error....

#ifdef FOR_MAKE

repo = ../telinst/
telinst_sources = $(repo)*.cc

all : alltelinst.cc alltelinst.h fitstoad.cc

	
classname = \(.*class *\)\([^ :]*\)\(.*\)
comment = classname is \2 in sed second field of 's' statement

alltelinst.cc : $(telinst_sources) fitstoad.cc
	\rm -f $@.new
	echo " /* Automatically generated by make from" $(telinst_sources) "*/" >> $@.new
	echo " /* makefile at the end of fitstoad.cc */" >> $@.new
	for telinst in ${telinst_sources};  do 	echo "#include \"$$telinst\"" >> $@.new; done
	echo "#ifndef USE_WCS"  >> $@.new
	echo "static AcceptorType AcceptorsArray [] = {" >> $@.new
	grep TYPE_SELECTOR ${telinst_sources}  | sed 's/${classname}/\&\2::Acceptor,/' >> $@.new
	echo "&Unknown::Acceptor};" >> $@.new
	echo '/* instanciation of templates IsOfKind */' >> $@.new
	grep TYPE_TESTER ${telinst_sources}  | sed 's/${classname}/template bool IsOfKind<\2>(const FitsHeader \&);/' >> $@.new
	grep TYPE_TESTER ${telinst_sources}  | sed 's/${classname}/template bool IsOfKind<\2>(const string \&);/' >> $@.new
	echo "#endif"  >> $@.new
	-diff $@.new $@ || mv -f $@.new $@

alltelinst.h : $(telinst_sources) fitstoad.cc
	\rm -f  $@.new
	echo " /* Automatically generated by make from" $(telinst_sources) "*/" >> $@.new
	echo " /* makefile at the end of fitstoad.cc */" >> $@.new
	echo '#ifndef ALLTELINST__H' >> $@.new
	echo '#define ALLTELINST__H' >> $@.new
	grep TYPE_TESTER ${telinst_sources} | sed 's/${classname}/class \2;/'   >> $@.new
	echo '#endif' >> $@.new
	diff $@.new $@ || mv -f $@.new $@

fitstoad.cc : $(telinst_sources)
	  touch $@
#endif




