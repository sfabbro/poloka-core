#include <iostream>
#include <cmath>
#include <algorithm>
#include <map>

#include "fileutils.h"  // StringTo{Lower,Upper}
#include "fitsimage.h"
#include "gtransfo.h" /* needed for  GuessLinWCS in TestTelInst */
#include "astroutils.h" /* for conversion routines */
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
  






#include "alltelinst.h"

typedef map<string,char> StringCharMap;


static string ToadBand(const string &keyval) 
{
  // This routine creates (only once) a map of filters to simple bands.
  static bool called = false;
  static StringCharMap ToadBandMap;
  if (!called)
    {
      called = true;
      string X = "UBVRIJHKLGZ";
      char F;
      string sF;
      for (unsigned int i=0;i<X.length();i++)
	{
	  F = X[i];
	  sF = string(&F,1);
	  ToadBandMap["BESS_"+sF] = F;
	  ToadBandMap["BESS"+sF] = F;
	  ToadBandMap[sF+"_BESS"] = F;
	  ToadBandMap[sF+"BESS"] = F;
          ToadBandMap["guessed_"+sF] = F;
	  ToadBandMap[sF+"_guessed"] = F;
	  ToadBandMap[sF+"guessed"] = F;
	  ToadBandMap["HARRIS_"+sF] = F;
	  ToadBandMap["HARRIS"+sF] = F;
	  ToadBandMap[sF+"_HARRIS"] = F;
	  ToadBandMap[sF+"HARRIS"] = F;
	  ToadBandMap[sF+"_SLOAN"] = F;
	  ToadBandMap[sF+"SLOAN"] = F;
	  ToadBandMap["SLOAN"+sF] = F;
	  ToadBandMap["SLOAN_"+sF] = F;
	  ToadBandMap[sF+"_STROMGREN"] = F;
	  ToadBandMap[sF+"STROMGREN"] = F;
	  ToadBandMap["STROMGREN"+sF] = F;
	  ToadBandMap["STROMGREN_"+sF] = F;
	  ToadBandMap[sF+"CTIO"] = F;
	  ToadBandMap[sF+"_GUNN"] = F;
	  ToadBandMap[sF+"GUNN"] = F;
	  ToadBandMap["GUNN_"+sF] = F;
	  ToadBandMap["GUNN"+sF] = F;
	  ToadBandMap[sF+"_WIDE"] = F;
	  ToadBandMap[sF+"WIDE"] = F;
	  ToadBandMap["WIDE_"+sF] = F;
	  ToadBandMap["WIDE"+sF] = F;
	  ToadBandMap[sF] = F;
	  ToadBandMap["F850LP"] = 'Z';
	  ToadBandMap["F814W"] = 'I';
	  ToadBandMap["F625W"] = 'R';
	  ToadBandMap["F675W"] = 'R';
	  ToadBandMap["F555W"] = 'V';
	  ToadBandMap["F110W"] = 'J';
	  ToadBandMap[sF+"#810"] = F;
	  ToadBandMap[sF+"#811"] = F;
	  ToadBandMap[sF+"#812"] = F;
	  ToadBandMap[sF+"#813"] = F;
	  ToadBandMap[sF+"#814"] = F;
	  ToadBandMap[sF+"#815"] = F;
	  ToadBandMap["SLOGUN"+sF] = F;
	  ToadBandMap["6 "+sF] = F;
	  ToadBandMap["4 "+sF] = F;
	  ToadBandMap[sF+"s"] = F;
	  ToadBandMap["1"+sF] = F;
	  ToadBandMap[sF+"1"] = F;
	  ToadBandMap[sF+"3"] = F;
	  ToadBandMap["1"+sF+"#7"] = F;
	  ToadBandMap["2"+sF+"#74"] = F;
	  ToadBandMap["3"+sF+"#75"] = F;
	  ToadBandMap["4"+sF+"#76"] = F;
	  ToadBandMap["5"+sF+"#12"] = F;
	  ToadBandMap["Sloang`"] = 'G';
	  ToadBandMap["Sloanr`"] = 'R';
	  ToadBandMap["Sloani`"] = 'I';
	  ToadBandMap["Sloanu`"] = 'U';
	  ToadBandMap["Sloanz`"] = 'Z';
	  ToadBandMap[sF+"Johnson"] = F;
	}  
    }
  StringCharMap::iterator match = ToadBandMap.find(keyval);
  if (match == ToadBandMap.end()) match = ToadBandMap.find(StringToUpper(keyval));
  if (match == ToadBandMap.end())
    {
      cerr << " Did not find the associate band with filter " << keyval << endl;
      return keyval;
    }    
  char result = (*match).second; /*ToadBandMap[keyval];*/
  return string(&result,1);
}


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
      cerr << " cannot decode 4 integers in " << FitsStuff << endl;
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
  

static string toads_date(const string &FitsDate)
{
  int yy,mm,dd;
  string result="";
  if ((sscanf(FitsDate.c_str(),"%d/%d/%d",&yy,&mm,&dd) != 3) &&
      (sscanf(FitsDate.c_str(),"%d-%d-%d",&yy,&mm,&dd) != 3)
      )
    return result;
  bool ok = false;
  if (dd > 31) {swap (yy,dd); ok = true;}// done
  else if (yy > 31) ok = true;
  if (!ok) return result;
  // convert 2 digits dates
  if (yy < 1949) yy += (yy>50)? 1900 : 2000;
  char date_string[16];
  sprintf(date_string,"%02d/%02d/%04d",dd,mm,yy);
  return string(date_string);
}


/* we could write a routine that uses the fits comment to decode dates with 2 digits years.
   hope it is not needed. Seems that FITS standard forbids it for yy>=2000 */




// go through these simple "routines" in case we want to change

#define TRANSLATOR_DEC(RoutineName) \
      FitsKey RoutineName( const FitsHeader &Head, const bool Warn) const

#define SIMPLE_TRANSLATOR(RoutineName,KeyTag) \
       TRANSLATOR_DEC(RoutineName) { return Head.KeyVal(KeyTag,Warn);}

#define RETURN_A_VALUE(RoutineName,Value) \
        TRANSLATOR_DEC(RoutineName) { return FitsKey(#RoutineName,Value);}



// This is the class of which to derive the actual instruments.

class VirtualInstrument
{
public :
  // the routine that says if a given header belongs to this class
  // not in the base class
  //  static VirtualInstrument *Acceptor(const FitsHeader &);

  //! -
  virtual string TelInstName() const = 0;
  virtual string TelName() const = 0;
  virtual string InstName() const = 0;

  //! try to guess a (linear) WCS.
  virtual bool GuessLinWCS(const FitsHeader &Head, TanPix2RaDec &Guess) const; 

  //! Which region on the sky
  virtual Frame SkyRegion(const FitsHeader &Head) const;

  //! default uses BIASSEC
  virtual Frame OverscanRegion(const FitsHeader &Head, const int Iamp) const;

  //! default : uses DATASEC. To be overwritten if not applicable
  virtual Frame TotalIlluRegion(const FitsHeader &Head) const;

  //! region corresponding to amplifier Iamp (only one amp -> iamp = 1)
  virtual Frame IlluRegion(const FitsHeader &Head, const int Iamp) const;

  //! region corresponding to amplifier Iamp once image trimed (only one amp -> iamp = 1)
  virtual Frame AmpRegion(const FitsHeader &Head, const int Iamp) const;

  //! gain corresponding to amplifier Iamp (only ona amp -> iamp = 1)
  virtual double AmpGain(const FitsHeader &Head, const int Iamp) const;

  //  virtual Frame GuessSkyArea(const FitsHeader &, const int IAmp) const = 0;

  // the default translators
  // if you change here, change also the documented default in ToadsKeys array...
  // I would advise NOT to change these defaults: derived classes rely on it
  virtual SIMPLE_TRANSLATOR(TOADPIXS,"PIXSCALE")
  virtual TRANSLATOR_DEC(TOADINST) { return FitsKey("TOADINST",TelInstName());}
  virtual TRANSLATOR_DEC(TOADTELE) { return FitsKey("TOADTELE",TelName());}
  virtual SIMPLE_TRANSLATOR(TOADFILT,"FILTER");
  virtual SIMPLE_TRANSLATOR(TOADEXPO,"EXPTIME");
  virtual SIMPLE_TRANSLATOR(TOADGAIN,"GAIN");
  virtual SIMPLE_TRANSLATOR(TOADRDON,"RDNOISE");
  virtual SIMPLE_TRANSLATOR(TOADRASC,"RA");
  virtual SIMPLE_TRANSLATOR(TOADDECL,"DEC");
  virtual SIMPLE_TRANSLATOR(TOADEQUI,"EQUINOX");
  virtual SIMPLE_TRANSLATOR(TOADAIRM,"AIRMASS");
  virtual SIMPLE_TRANSLATOR(TOADUTIM,"UT")
  virtual FitsKey   TOADDATE(const FitsHeader &Head, const bool Warn) const
  { string result = toads_date(string(Head.KeyVal("DATE-OBS")));
    if (result.length() == 0) result = toads_date(string(Head.KeyVal("DATE")));
    if (result.length() == 0) result = toads_date(string(Head.KeyVal("DATEOBS")));
    if (Warn && result.length() == 0) cout << " no DATE-OBS nor DATE nor DATEOBS " << endl;
    return FitsKey("TOADDATE", result);
  }
  virtual FitsKey   TOADMMJD(const FitsHeader &Head, const bool Warn) const
  { 
    double val = ModifiedModifiedJulianDay(Head);
    return FitsKey("TOADMMJD", val);
  }
  
  virtual FitsKey   TOADBAND(const FitsHeader &Head, const bool Warn) const 
     { 
	 string filt_str = Head.KeyVal("TOADFILT");
	 RemovePattern(filt_str," ");	 
	 return FitsKey("TOADBAND",ToadBand(filt_str));
     }
  virtual SIMPLE_TRANSLATOR(TOADTYPE,"IMAGETYP");
  virtual SIMPLE_TRANSLATOR(TOADOBJE,"OBJECT"); 
  virtual RETURN_A_VALUE(TOADNAMP,1);
  virtual RETURN_A_VALUE(TOADCHIP,1);
  virtual SIMPLE_TRANSLATOR(TOADPZPT,"ZEROUSNO");
  virtual ~VirtualInstrument() {};

  /* if you add things in the functionnalities, update the TestTelInst function
  accordingly */

};



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



// default Implementations!
bool VirtualInstrument::GuessLinWCS(const FitsHeader &Head, 
				    TanPix2RaDec &Guess) const
{
  cout  << " trying default GuessLinWCS" << endl;
  return ComputeLinWCS (Head, Head.ImageCenter(), GtransfoIdentity(), Guess);
}

Frame  VirtualInstrument::SkyRegion(const FitsHeader &Head) const
{
  int nx,ny;
  Head.ImageSizes(nx,ny);
  TanPix2RaDec Pix2RaDec;
  if (!GuessLinWCS(Head, Pix2RaDec)) return Frame();
  cout << " Lin WCS Guess " << Pix2RaDec << endl;
  Point p00 = Pix2RaDec.apply(Point(0,0));
  Point p11 = Pix2RaDec.apply(Point(nx,ny));
  return Frame(p00,p11);
  
}  


Frame VirtualInstrument::AmpRegion(const FitsHeader &Head, const int Iamp) const
{
  Frame illuRegion = IlluRegion(Head, Iamp);
  Frame totIlluRegion = TotalIlluRegion(Head);
  return illuRegion.ApplyTransfo(GtransfoLinShift(-totIlluRegion.xMin, - totIlluRegion.yMin));
}

Frame VirtualInstrument::IlluRegion(const FitsHeader &Head, const int Iamp) const
{
  if (Iamp != 1)
    {
      cerr << "using default IlluRegion for iamp != 1. expect serious troubles " << endl;
      cerr << Head.FileName() << ' ' << "tel/inst" << TelInstName() << endl;
    }
  return TotalIlluRegion(Head);
}

Frame VirtualInstrument::TotalIlluRegion(const FitsHeader &Head) const
{
  Frame result;
  if (Head.HasKey("DATASEC"))
    fits_imregion_to_frame(string(Head.KeyVal("DATASEC")), result);
  else result = Frame(Head,WholeSizeFrame);
  return result;
}

Frame VirtualInstrument::OverscanRegion(const FitsHeader &Head, const int Iamp) const
{
  if (Head.HasKey("BIASSEC"))
    {
      string keyval = Head.KeyVal("BIASSEC");
      Frame result;
      fits_imregion_to_frame(keyval, result);
      return result;
    }
  else return Frame();
}

double VirtualInstrument::AmpGain(const FitsHeader &Head, const int Iamp) const
{
  double gain = Head.KeyVal("TOADGAIN");
  return gain;
}

void VirtualInstrumentDestructor(VirtualInstrument *p)
{
  delete p;
}

/******************** implementation of Unknown intrument *************/

class Unknown : public VirtualInstrument {

public :
  static VirtualInstrument *Acceptor(const FitsHeader &Head) {return new Unknown();};
  string TelInstName() const {return "Unknown";};
  string InstName() const {return "Unknown";};
  string TelName() const {return "Unknown";};

};
    


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


bool GuessLinWCS(const FitsHeader &Header, TanPix2RaDec &Guess)
{
  VirtualInstrument *p = Header.TelInst();  
  if (p->GuessLinWCS(Header, Guess)) return true;
  // the tel/inst specific procedure failed. try the default one ...
  //.... but only if it is different 
  // the following test is refused by the compiler
  //  if (&(p->GuessLinWCS) == &(p->VirtualInstrument::GuessLinWCS)) return false;
  return (p->VirtualInstrument::GuessLinWCS(Header, Guess));
}

Frame SkyRegion(const FitsHeader &Header)
{
  VirtualInstrument *p = Header.TelInst();
  return p->SkyRegion(Header);
}
  

typedef VirtualInstrument *(*AcceptorType)(const FitsHeader &);

/* 
   alltellinst .cc contains the implementation of all other (real!)
   instruments implementation of other instruments : 
   this file is generated by a makefile located at the end of this fitstoad.cc 
   Guidelines for the implementation of a new instrument : see 
   the (long) comment at the beginning of this file.
*/
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
  TanPix2RaDec wcs;
  cout << " test GuessLinWCS " << endl;
  if (GuessLinWCS(Head,wcs))
    {
      cout << " considered successful " << endl;
      cout << wcs << endl;
    }
  else cout << " return false " << endl;
  cout << " test SkyRegion " << endl << SkyRegion(Head) << endl;
  cout << " ------------------------------------------------" << endl;
}

#ifdef STORAGE
void FitsHeader::SniffTelInst() 
{
  //instruments with the INSTRUME key
  string KEYINST = KeyVal("INSTRUME");
  KEYINST = StringToUpper(KEYINST);  
  if (KEYINST == "BTC") {TelInst = CtioBtc;return;}
  if (KEYINST == "PFCCD") {TelInst = CtioPfc;return;}
  if (KEYINST == "EFOSC"){TelInst = EsoEfosc;return;}
  if (KEYINST == "CASS. AUX. PORT") {TelInst = WhtCat;return;}
  if (KEYINST == "JAG-CCD") {TelInst = JktJag;return;}
  if (KEYINST == "FORS1") {TelInst = VltAntuFors1;return;}
  if (KEYINST == "WFPC2") {TelInst = HstWfpc2;return;}
  if (KEYINST == "NICMOS") {TelInst = HstNicmos;return;}
  if (KEYINST == "PFCU") {TelInst = IntTek3;return;} 
  if (KEYINST == "ALFOSC-FASU") {TelInst = NotAlfosc;return;}
  if (KEYINST == "WFC")  
    {
      double jd = KeyVal("JD");
      if (HasKey("JD") &&  jd > 2451410.5822505) /* roughly Aug 20th 1999 when the Daq was changed */
          { TelInst = IntWfcNewDaq; return;}
      if (HasKey("CCDTYPE"))
        {
	  string keyccd  = KeyVal("CCDTYPE");
	  keyccd = StringToLower(keyccd);
	  if (keyccd=="loral") {TelInst = Int1Wfc;return;}
	  if (keyccd=="eev42") {TelInst = Int4Wfc;return;}
        }
    }

  //instruments with the DETECTOR key
  KEYINST = KeyVal("DETECTOR");
  if (KEYINST == "WIYNMiniMos") {TelInst = WIYNMiniMos;return;}  
  if (KEYINST == "Mosaic2") {TelInst = CtioMosaic2;return;} // $CHECK$ nrl
  KEYINST = StringToUpper(KEYINST);
  if (KEYINST == "S2KB") {TelInst = WiynS2kb;return;}  
  if (KEYINST == "CFH12K") {TelInst = Cfht12K;return;}  


  //instruments with another key
  KEYINST = KeyVal("TELESCOP");
  KEYINST = StringToUpper(KEYINST);
  RemovePattern(KEYINST," ");
  if (KEYINST == "KECKI") {TelInst = KeckIILris;return;}  
  if (KEYINST == "KECKII") {TelInst = KeckIIesi;return;}  
  KEYINST = KeyVal("DETSTAT");
  if (KEYINST=="VME-HP1000 V2.1") {TelInst=EsoEfosc;return;}//sounds quite specific

  //Hack.
  TelInst=Unknown;
}
#endif



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

// OLd FitsToad stuff. When you reimplement a tel/inst
// in the new scheme, put the old code into ifdefs at the end of the your new file,
// and remove it from here....

#ifdef NEW_TELINST
//! a substitute for TOADSCAN and TOADILLU when ccd is read out via several amps. WhichInfo is ILLU or SCAN.
Frame FitsHeader::OverIlluRegions(const IlluOrScan WhichInfo, const int IAmp) const
{
  if (AvailableInstruments[TelInst]->OverIlluRegions)
    return AvailableInstruments[TelInst]->OverIlluRegions(*this,WhichInfo, IAmp);
  else
    switch (WhichInfo)
      {
	case (SCAN)
	  {
	    int x0, y0, nx, ny;
	    string keyval = KeyVal("TOADSCAN");
	    if (sscanf(keyval.c_str(),"[%d,%d;%d,%d]", &x0, &nx, &y0, &ny) == 4) 
	      return Frame(x0-1, y0-1, x0+nx-1, y0+ny-1);	
	    cerr << " cannot convert " << keyval << " to overscan region for file " << FileName() << endl
		 << " tel/inst " << TelInstName() << endl;
	    return Frame();
	    break;
	  }
	case (ILLU)
	  {  
	    int x0, y0, nx, ny;
	    string keyval = KeyVal("TOADILLU");
	    if (sscanf(keyval.c_str(),"[%d,%d;%d,%d]", &x0, &nx, &y0, &ny) == 4) 
	      return Frame(x0-1, y0-1, x0+nx-1, y0+ny-1);
	    cerr << " cannot convert " << keyval << " to illuminated region for file " << FileName() << endl
	     << " tel/inst " << TelInstName() << endl;
	    return Frame();
	    break;
	  }
      }
}



//***************************************************************
//       Formatting routines for toadskeys
//****************************************************************

//***************************************************************
//       Formattting routines for toadskeys
//****************************************************************
FitsKey FitsHeader::CtioMosaic2Format(const string &KeyName) const
{
  if (KeyName == "TOADPIXS") {
    double xpxsz = KeyVal("PIXSIZE1") ;
    double ypxsz = KeyVal("PIXSIZE2") ;
    return FitsKey(KeyName,0.5*(xpxsz+ypxsz));
  }
  if (KeyName == "TOADINST") return KeyVal("DETECTOR");
  if (KeyName == "TOADFILT") return KeyVal("FILTER");
  if (KeyName == "TOADEXPO") return KeyVal("EXPTIME");
  if (KeyName == "TOADGAIN") return KeyVal("GAIN");
  if (KeyName == "TOADRDON") return KeyVal("RDNOISE");
  if (KeyName == "TOADRASC") return KeyVal("RA");
  if (KeyName == "TOADDECL") return KeyVal("DEC");
  
  if (KeyName == "TOADEQUI") return KeyVal("EQUINOX");
  if (KeyName == "TOADAIRM") return KeyVal("AIRMASS");
  if (KeyName == "TOADUTIM") return KeyVal("TIME-OBS");
  if (KeyName == "TOADDATE") 
    {
      int dd,mm,yy;
      string keyval = KeyVal("DATEOBS");
      if (sscanf(keyval.c_str(),"%d-%d-%d",& yy,&mm,&dd) == 3)
 	{
          char date_string[64];
	  sprintf(date_string,"%d/%d/%d",dd,mm,yy);
 	  return FitsKey(KeyName,string(date_string));
 	}
    };
  if (KeyName == "TOADSCAN")
    { 
      string keyval = KeyVal("BIASSEC");
      int x0,x1,y0,y1,nx,ny;
      if (sscanf(keyval.c_str(),"[%d:%d,%d:%d]",&x0,&x1,&y0,&y1) == 4)
	{
	  nx = x1-x0+1;
	  ny = y1-y0+1;	
          char sec_charstar[64];
          sprintf(sec_charstar,"[%d,%d;%d,%d]",x0,nx,y0,ny);
	  return FitsKey(KeyName,string(sec_charstar));
	}
    }
  if (KeyName == "TOADILLU")
    {
      string keyval = KeyVal("DATASEC");
      int x0,x1,y0,y1,nx,ny;
      if (sscanf(keyval.c_str(),"[%d:%d,%d:%d]",&x0,&x1,&y0,&y1) == 4)
	{
	  nx = x1-x0+1;
	  ny = y1-y0+1;
          char sec_charstar[64];
	  sprintf(sec_charstar,"[%d,%d;%d,%d]",x0,nx,y0,ny);
	  return FitsKey(KeyName, string(sec_charstar));
	}
    }
  if (KeyName == "TOADTYPE") return KeyVal("OBSTYPE");
  if (KeyName == "TOADOBJE") return KeyVal("OBJECT") ;
  if (KeyName == "TOADCHIP") {return KeyVal("CHIPID");} // 
  if (KeyName == "TOADBAND")
    {
      string Filter = KeyVal("FILTER");
      return FitsKey(KeyName,Filter.substr(0,1)) ;
    }
  if (KeyName == "TOADNAMP") return FitsKey(KeyName,1);
  
  return FitsKey(KeyName,NOVAL);
  
}


FitsKey FitsHeader::CtioPfcFormat(const string &KeyName) const
{
  if (KeyName == "TOADPIXS") 
    {
      double xpix = KeyVal("XPIXSIZE");
      double ypix = KeyVal("YPIXSIZE");
      return FitsKey(KeyName,(xpix+ypix)/2);
    }
  if (KeyName == "TOADINST") return KeyVal("INSTRUME");
  if (KeyName == "TOADFILT") return KeyVal("FNAME2");
  if (KeyName == "TOADEXPO") return KeyVal("EXPTIME");
  if (KeyName == "TOADGAIN")  
    {
      if (HasKey("GAIN")) return KeyVal("GAIN");
      else
	 {
	   double gain1 = KeyVal("GTGAIN11");
	   double gain2 = KeyVal("GTGAIN12");
	   double gain3 = KeyVal("GTGAIN21");
	   double gain4 = KeyVal("GTGAIN22");
	   return FitsKey(KeyName,(gain1+gain2+gain3+gain4)/4.0);
	 }
    }
  if (KeyName == "TOADRDON")
    {
      double ron1 = KeyVal("GTRON11");
      double ron2 = KeyVal("GTRON12");
      double ron3 = KeyVal("GTRON21");
      double ron4 = KeyVal("GTRON22");
      return FitsKey(KeyName,(ron1+ron2+ron3+ron4)/4.0);
    }
  if (KeyName == "TOADRASC") return KeyVal("RA");
  if (KeyName == "TOADDECL") return KeyVal("DEC");
  if (KeyName == "TOADEQUI") return KeyVal("EPOCH");
  if (KeyName == "TOADAIRM") return KeyVal("AIRMASS");
  if (KeyName == "TOADUTIM") return KeyVal("UT");
  if (KeyName == "TOADDATE") 
    {
      int dd,mm,yy;
      string keyval = KeyVal("DATE-OBS");
      if (sscanf(keyval.c_str(),"%d/%d/%d",&dd,&mm,&yy) == 3)
 	{
          char date_string[64];
	  sprintf(date_string,"%d/%d/%d",dd,mm,yy+1900);
 	  return FitsKey(KeyName,string(date_string));
 	}
    }
  if (KeyName == "TOADSCAN") return FitsKey(KeyName,"Has 4 scan regions");
  if (KeyName == "TOADILLU") return FitsKey(KeyName,"Has 4 illu regions") ;
  if (KeyName == "TOADTYPE") return KeyVal("IMAGETYP");
  if (KeyName == "TOADOBJE") return KeyVal("OBJECT");
  if (KeyName == "TOADCHIP") return FitsKey(KeyName,1);
  if (KeyName == "TOADBAND")
    {
      string filter = KeyVal("FILTER2");
      RemovePattern(filter," " );
      return FitsKey(KeyName,ToadBand(StringToUpper(filter)));	  
    }
  if (KeyName == "TOADNAMP") return FitsKey(KeyName,4);
  
  return FitsKey(KeyName,NOVAL);
  
}


FitsKey FitsHeader::WIYNMiniMosFormat(const string &KeyName) const
{
  if (KeyName == "TOADPIXS") return KeyVal("PIXSCAL1");
  if (KeyName == "TOADINST") return KeyVal("DETECTOR");
  if (KeyName == "TOADFILT") return KeyVal("FILTER");
  if (KeyName == "TOADEXPO") return KeyVal("EXPTIME");
  if (KeyName == "TOADGAIN") return KeyVal("GAIN");
  if (KeyName == "TOADRDON") return KeyVal("RDNOISE");
  /*
    Below is an example of problem associated with the wiyn: watch coordinates and equinox, they don't match.
    We then choose the RA and DEC keys, they seem more exact.
    RA      = '28:55:57.70'        / Telescope RA
    DEC     = '-3:36:57.68'        / Telescope DEC
    RAOFFST = '00:00:0.00'         / Telescope RA offset
    DECOFFST= '00:00:0.00'         / Telescope DEC offset
    RADECSYS= 'FK5'                / coordinate system
    TARGRA  = '04:56:11.60'        / right ascension
    TARGDEC = '-3:41:26.00'        / declination
    EQUINOX =               1998.1 / equinox of position
    EPOCH   =               1998.1 / same as EQUINOX (for back compat.)
  */
  if (KeyName == "TOADRASC") 
    {
      string rastr = KeyVal("RA");
      double radeg = RaStringToDeg(rastr);
      while (radeg > 360) radeg -= 360;
      return FitsKey(KeyName,RaDegToString(radeg));
    }
  if (KeyName == "TOADDECL") return KeyVal("DEC");
  if (KeyName == "TOADEQUI") /* return KeyVal("EQUINOX"); */
    {
      string radecsys = KeyVal("RADECSYS");
      //should be 2000 but works only with 1950...should really learn the coordinate systems
      if (strstr(radecsys.c_str(),"FK5")) return FitsKey(KeyName,2000.0);
      if (strstr(radecsys.c_str(),"FK4")) return FitsKey(KeyName,2000.0); // never seen that, according to RA and DEC keys
      return FitsKey(KeyName,2000.0);
    } 
  if (KeyName == "TOADAIRM") return KeyVal("AIRMASS");
  if (KeyName == "TOADUTIM") 
     {
      string keyval = KeyVal("DATE-OBS");
      return KeyVal(strstr(keyval.c_str(),"T"));
     }
  if (KeyName == "TOADDATE") 
    {  
      int dd,mm,yy;
      string keyval = KeyVal("DATE-OBS");
      if (sscanf(keyval.c_str(),"%d-%d-%dT",&dd,&mm,&yy) == 3)
 	{
 	  yy = yy + 1900;
          char date_string[64];
	  sprintf(date_string,"%d/%d/%d",dd,mm,yy);
 	  return FitsKey(KeyName,string(date_string));
 	}
    }
  if (KeyName == "TOADSCAN") 
    {
      string keyval = KeyVal("BIASSEC");
      int x0,x1,y0,y1,nx,ny;
      if (sscanf(keyval.c_str(),"[%d:%d,%d:%d]",&x0,&x1,&y0,&y1) == 4)
	{
	  nx = x1-x0+1;
	  ny = y1-y0+1;
	  char sec_charstar[64];
	  sprintf(sec_charstar,"[%d,%d;%d,%d]",x0,nx,y0,ny);
	  return FitsKey(KeyName, string(sec_charstar));
	}
    }
  if (KeyName == "TOADILLU") 
    {
      string keyval = KeyVal("DATASEC");
      int x0,x1,y0,y1,nx,ny;
      if (sscanf(keyval.c_str(),"[%d:%d,%d:%d]",&x0,&x1,&y0,&y1) == 4)
      {
	nx = x1-x0+1;
	ny = y1-y0+1;
	char sec_charstar[64];
	sprintf(sec_charstar,"[%d,%d;%d,%d]",x0,nx,y0,ny);
	return FitsKey(KeyName, string(sec_charstar));
      }
    }
  if (KeyName == "TOADTYPE") return KeyVal("OBSTYP");
  if (KeyName == "TOADOBJE") return KeyVal("OBJECT");
  if (KeyName == "TOADCHIP") return KeyVal("IMAGEID");
  if (KeyName == "TOADBAND")
    {
      string filter = KeyVal("FILTER");
      return FitsKey(KeyName,ToadBand(StringToUpper(filter)));	  
    }
// although there are 2 per chip but when splitted, 1 image (1x4k) corresponds to 1 amp
  if (KeyName == "TOADNAMP") return FitsKey(KeyName,1); 


  return FitsKey(KeyName,NOVAL);      
}

//*********************************************************************
FitsKey FitsHeader::IntTek3Format(const string &KeyName) const
{
  if (KeyName == "TOADPIXS") return FitsKey(KeyName,0.587);
  if (KeyName == "TOADINST") return KeyVal("INSTRUME");
  if (KeyName == "TOADFILT") return FitsKey(KeyName, string(KeyVal("AGCOLFLP")));
  if (KeyName == "TOADEXPO") return KeyVal("EXPTIME");
  if (KeyName == "TOADGAIN") return KeyVal("GAIN");
  if (KeyName == "TOADRDON") return KeyVal("READNOIS");
  if (KeyName == "TOADRASC") return KeyVal("RA");
  if (KeyName == "TOADDECL") return KeyVal("DEC");
  if (KeyName == "TOADEQUI") 
    {
      double epoch = 0;
      char dumb;
      string keyval = KeyVal("EQUINOX"); 
      if (sscanf(keyval.c_str(),"%c%lf",&dumb,&epoch) != 2)
	{cerr << "Needs epoch "<< epoch 
	      <<" in correct INT-Tek3 format : " << KeyVal("EQUINOX") <<endl;}
      return FitsKey(KeyName,epoch);
    }
  if (KeyName == "TOADAIRM") return KeyVal("AIRMASS");
  if (KeyName == "TOADUTIM") return KeyVal("UTSTART");
  if (KeyName == "TOADDATE")
    {
      int dd,mm,yy;
      string keyval = KeyVal("DATE-OBS");
      if (sscanf(keyval.c_str(),"%d/%d/%d",&dd,&mm,&yy) == 3)
 	{ 
	  yy +=1900;  
	  char date_charstar[64];
	  sprintf(date_charstar,"%d/%d/%d",dd,mm,yy);
	  return FitsKey(KeyName,string(date_charstar));
	}
    } 

  if (KeyName == "TOADSCAN") return FitsKey(KeyName,"[1079,42;1,1024]");
  if (KeyName == "TOADILLU") return FitsKey(KeyName,"[51,1024;1,1024]");
  if (KeyName == "TOADTYPE") return KeyVal("IMAGETYP");
  if (KeyName == "TOADOBJE") return KeyVal("OBJECT");
  if (KeyName == "TOADCHIP") return FitsKey(KeyName,1);
  if (KeyName == "TOADBAND") return FitsKey(KeyName,"R");
    /*
    {
      string filter = KeyVal("AGCOLFLP");
      int which;
      char band;
      sscanf(filter.c_str(),"%d (%c)   ",&which,&band);
      return FitsKey(KeyName,band);
      }*/
  if (KeyName == "TOADNAMP") return FitsKey(KeyName,1);
  return FitsKey(KeyName,NOVAL);      
}




//*********************************************************************
FitsKey FitsHeader::WhtCatFormat(const string &KeyName) const
{
  if (KeyName == "TOADINST") return KeyVal("INSTRUME");
  if (KeyName == "TOADPIXS") return FitsKey(KeyName,0.108);
  if (KeyName == "TOADFILT") return KeyVal("CAGAXFLT");
  if (KeyName == "TOADEXPO") return KeyVal("EXPTIME");
  if (KeyName == "TOADGAIN") return KeyVal("GAIN");
  if (KeyName == "TOADRDON") return KeyVal("READNOIS");
  if (KeyName == "TOADRASC") return KeyVal("RA");
  if (KeyName == "TOADDECL") return KeyVal("DEC");
  if (KeyName == "TOADEQUI") 
    {
      double epoch = 0;
      char dumb;
      string keyval = KeyVal("EQUINOX"); 
      if (sscanf(keyval.c_str(),"%c%lf",&dumb,&epoch) != 2)
	{cerr << "Needs epoch "<< epoch <<" in correct WHT-CAT format : " << KeyVal("EQUINOX") <<endl;}
      return FitsKey(KeyName,epoch);
    }
  if (KeyName == "TOADAIRM") return KeyVal("AIRMASS");
  if (KeyName == "TOADUTIM") return KeyVal("UTSTART");
  if (KeyName == "TOADDATE")
    {
      int dd,mm,yy;
      string keyval = KeyVal("DATE-OBS");
      char date_charstar[64];
      if (sscanf(keyval.c_str(),"%d/%d/%d",&yy,&mm,&dd) == 3)
 	{ 
	  sprintf(date_charstar,"%d/%d/%d",dd,mm,yy);
	  return FitsKey(KeyName,string(date_charstar));
	}
    } 
  if (KeyName == "TOADSCAN") 
    {
      char sec_charstar[64];
      sprintf(sec_charstar,"[%d,%d;%d,%d]",1,47,1,1024);
      return FitsKey(KeyName,string(sec_charstar));
    }
  if (KeyName == "TOADILLU") 
    {
      char sec_charstar[64];
      sprintf(sec_charstar,"[%d,%d;%d,%d]",51,1024,1,1024);
      return FitsKey(KeyName,string(sec_charstar));
    }
  if (KeyName == "TOADTYPE") return KeyVal("OBSTYPE");
  if (KeyName == "TOADOBJE") return KeyVal("OBJECT");
  if (KeyName == "TOADCHIP") return FitsKey(string("TOADCHIP"),"1");
  if (KeyName == "TOADBAND")
    {
      string filter = KeyVal("CAGAXFLT");
      return FitsKey(KeyName,ToadBand(StringToUpper(filter)));	  
    }
  if (KeyName == "TOADNAMP") return FitsKey(KeyName,1);
  return FitsKey(KeyName,NOVAL);      
}


//*********************************************************************
// Let's face it EFOSC is a piece of shit
FitsKey FitsHeader::EsoEfoscFormat(const string &KeyName) const
{ 
  if (KeyName == "TOADINST") return FitsKey(KeyName,"EFOSC");
  if (KeyName == "TOADPIXS") return FitsKey(KeyName,0.604);
  if (KeyName == "TOADFILT") 
    {
      int filter = KeyVal("FILTRNR");
      if (filter==4) return FitsKey(KeyName,"R");	  
      if (filter==10) return FitsKey(KeyName,"I");  
    }
  if (KeyName == "TOADEXPO") 
    {
      double t0 =  KeyVal("TM-START");
      double t1 =  KeyVal("TM-END");
      return FitsKey(KeyName,t1-t0);
    }
  if (KeyName == "TOADGAIN") return FitsKey(KeyName, 3.5);
  if (KeyName == "TOADRDON") return FitsKey(KeyName, 8.2);
  if (KeyName == "TOADRASC") return FitsKey(KeyName, RaDegToString(double(KeyVal("POSTN-RA"))));
  if (KeyName == "TOADDECL") return FitsKey(KeyName, DecDegToString(double(KeyVal("POSTN-DE"))));
  if (KeyName == "TOADEQUI") return FitsKey(KeyName,2000.00);
  if (KeyName == "TOADAIRM") return KeyVal("AIRMASS");
  if (KeyName == "TOADUTIM") 
    {
      int aftermidnight = int(KeyVal("TM-START"));
      int hh = aftermidnight/3600;
      int mm = (aftermidnight % 3600)/60;
      int ss = (aftermidnight % 3600) % 60;
      char chartime[12];
      sprintf(chartime,"%d:%d:%d",hh,mm,ss);
      return FitsKey(KeyName,string(chartime));
    }
  if (KeyName == "TOADDATE") 
    {
      int dd,mm,yy;
      string keyval = KeyVal("DATE-OBS");
      char date_charstar[64];
      if (sscanf(keyval.c_str(),"%d/%d/%d",&dd,&mm,&yy) == 3)
 	{
 	  yy = yy + 1900;
	  sprintf(date_charstar,"%d/%d/%d",dd,mm,yy);
 	  return FitsKey(KeyName,string(date_charstar));
 	}; 
    }
  if (KeyName == "TOADSCAN") return FitsKey(KeyName,string("[1,49;1,512]"));
  if (KeyName == "TOADILLU") return FitsKey(KeyName,string("[51,512;1,512]"));
  if (KeyName == "TOADTYPE") return KeyVal("EXPTYP");
  if (KeyName == "TOADOBJE") return KeyVal("OBJECT");
  if (KeyName == "TOADBAND")
    {
      int filter = KeyVal("FILTRNR");
      if (filter==4) return FitsKey(KeyName,"R");	  
      if (filter==10) return FitsKey(KeyName,"I");  
    }
  if (KeyName == "TOADCHIP") return FitsKey(KeyName,"1");
  if (KeyName == "TOADNAMP") return FitsKey(KeyName,1);
  return FitsKey(KeyName,NOVAL);      
}


//*********************************************************************
FitsKey FitsHeader::JktJagFormat(const string &KeyName) const
{
  if (KeyName == "TOADINST") return KeyVal("INSTRUME");
  if (KeyName == "TOADPIXS") return FitsKey(KeyName,0.33);
  if (KeyName == "TOADFILT") return FitsKey(KeyName, string(KeyVal("JAGFBAND"))+string(KeyVal("JAGFSYS")));
  if (KeyName == "TOADEXPO") return KeyVal("EXPTIME");
  if (KeyName == "TOADGAIN") return KeyVal("GAIN");
  if (KeyName == "TOADRDON") return KeyVal("READNOIS");
  if (KeyName == "TOADRASC") return KeyVal("RA");
  if (KeyName == "TOADDECL") return KeyVal("DEC");
  if (KeyName == "TOADEQUI") 
    {
      double epoch = 0;
      char dumb;
      string keyval = KeyVal("EQUINOX"); 
      if (sscanf(keyval.c_str(),"%c%lf",&dumb,&epoch) != 2)
	{cerr << "Needs epoch in correct JKT-JAG format" << KeyVal("EQUINOX")<<endl;}
      return FitsKey(KeyName,epoch);
    }
  if (KeyName == "TOADAIRM") return KeyVal("AIRMASS");
  if (KeyName == "TOADUTIM") return KeyVal("UTSTART");
  if (KeyName == "TOADDATE")
    {
      int dd,mm,yy;
      string keyval = KeyVal("DATE-OBS");
      if (sscanf(keyval.c_str(),"%d-%d-%d",&yy,&mm,&dd) == 3)
 	{ 
	  char date_charstar[64];
	  sprintf(date_charstar,"%d/%d/%d",dd,mm,yy);
	  return FitsKey(KeyName,string(date_charstar));
	}
    } 
  if (KeyName == "TOADSCAN") 
    {
      string keyval = KeyVal("BIASSEC");
      int x0,x1,y0,y1,nx,ny;
      if (sscanf(keyval.c_str(),"[%d:%d,%d:%d]",&x0,&x1,&y0,&y1) == 4)
	{
	  nx = x1-x0+1;
	  ny = y1-y0+1;
	  char sec_charstar[64];
	  sprintf(sec_charstar,"[%d,%d;%d,%d]",x0,nx,y0,ny);
	  return FitsKey(KeyName,string(sec_charstar));
	}
    }
  if (KeyName == "TOADILLU") 
    {

      return FitsKey(KeyName,string("[350,1200;550,1200]"));

     /* string keyval = KeyVal("TRIMSEC");
      int x0,x1,y0,y1,nx,ny;
      if(sscanf(keyval.c_str(),"[%d:%d,%d:%d]",&x0,&x1,&y0,&y1) == 4)
	{
	  nx = x1-x0+1;
	  ny = y1-y0+1;
	  char sec_charstar[64];
	  sprintf(sec_charstar,"[%d,%d;%d,%d]",x0,nx,y0,ny);
	  return FitsKey(KeyName,string(sec_charstar)) ;
	}*/ 
    }
  if (KeyName == "TOADTYPE") return KeyVal("IMAGETYP");
  if (KeyName == "TOADOBJE") return KeyVal("OBJECT");
  if (KeyName == "TOADCHIP") 
    if (TelInst == IntWfcNewDaq) return KeyVal("IMAGEID");
    else return KeyVal("CCDCHIP");
  if (KeyName == "TOADBAND")
    {
      string filter = KeyVal("JAGFBAND");
      return FitsKey(KeyName,ToadBand(StringToUpper(filter)));	  
    }
  if (KeyName == "TOADNAMP") return FitsKey(KeyName,1);
  return FitsKey(KeyName,NOVAL);      
}


//*********************************************************************
FitsKey FitsHeader::KeckIILrisFormat(const string &KeyName) const
{
  if (KeyName == "TOADNAMP") return FitsKey(KeyName,2);
  if (KeyName == "TOADINST") return KeyVal("PONAME");
  if (KeyName == "TOADPIXS") return FitsKey(KeyName,0.215);
  if (KeyName == "TOADFILT") return KeyVal("REDFILT");
  if (KeyName == "TOADEXPO") return KeyVal("EXPOSURE");
  if (KeyName == "TOADGAIN") 
     {
	if (HasKey("GAIN")) return KeyVal("GAIN");
	else return FitsKey(KeyName,(1.97+2.1)/2.0);
     }
  if (KeyName == "TOADRDON") return FitsKey(KeyName,(6.3+6.6)/2.0);
  if (KeyName == "TOADRASC") return KeyVal("RA");
  if (KeyName == "TOADDECL") return KeyVal("DEC");
  if (KeyName == "TOADEQUI") return KeyVal("EQUINOX");
  if (KeyName == "TOADAIRM") return KeyVal("AIRMASS");
  if (KeyName == "TOADUTIM") return KeyVal("UT");
  if (KeyName == "TOADDATE") 
    {
      string keyval;
      if (HasKey("DATE-OBS")) keyval = KeyVal("DATE-OBS");
      else keyval= KeyVal("DATE");
      int dd,mm,yy;
      if (sscanf(keyval.c_str(),"%d/%d/%d",&dd,&mm,&yy) == 3)
 	{ 
	  char date_charstar[64];
	  sprintf(date_charstar,"%d/%d/%d",dd,mm,yy+1900);
	  return FitsKey(KeyName,string(date_charstar));
	}
      if (sscanf(keyval.c_str(),"%d-%d-%d",&yy,&mm,&dd) == 3)
 	{ 
	  char date_charstar[64];
	  sprintf(date_charstar,"%d/%d/%d",dd,mm,yy);
	  return FitsKey(KeyName,string(date_charstar));
	}
    } 

  if (KeyName == "TOADSCAN") 
    {
      //Only considering the right side of overscan
      string sec  = "[2091,80;1,2048]";
      return FitsKey(KeyName,sec);
    }
  if (KeyName == "TOADILLU") 
    {
      string sec  = "[240,1640;1,2048]";
      return FitsKey(KeyName,sec);
    }
  if (KeyName == "TOADTYPE") return KeyVal("IMAGETYP");
  if (KeyName == "TOADOBJE") return KeyVal("OBJECT");
  if (KeyName == "TOADCHIP") return FitsKey(KeyName,1);
  if (KeyName == "TOADBAND")
    {
      string filter = KeyVal("REDFILT");
      return FitsKey(KeyName,ToadBand(StringToUpper(filter)));	  
    }
  return FitsKey(KeyName,NOVAL);      
}
//*********************************************************************
FitsKey FitsHeader::HstWfpc2Format(const string &KeyName) const
{
  if (KeyName == "TOADNAMP") return FitsKey(KeyName,1);
  if (KeyName == "TOADINST") return KeyVal("INSTRUME");
  if (KeyName == "TOADPIXS") return FitsKey(KeyName,0.1);
  if (KeyName == "TOADFILT") return KeyVal("FILTNAM1");
  if (KeyName == "TOADEXPO") return KeyVal("EXPTIME");
  if (KeyName == "TOADGAIN") return KeyVal("ATODGAIN");
  // check the crap with PHOTFLAM or PHOTZPT for a proper calib.
  if (KeyName == "TOADRDON") return FitsKey(KeyName,0.72);  
  if (KeyName == "TOADRASC") 
      {
	  double rastr=KeyVal("RA_TARG");
	  return FitsKey(KeyName,RaDegToString(rastr));
      }
  if (KeyName == "TOADDECL")       
      {
	  double decstr=KeyVal("DEC_TARG");
	  return FitsKey(KeyName,DecDegToString(decstr));
      }

  if (KeyName == "TOADEQUI") return KeyVal("EQUINOX");
  if (KeyName == "TOADAIRM") return FitsKey(KeyName,0.0);
  if (KeyName == "TOADUTIM") return KeyVal("EXPSTART");
  if (KeyName == "TOADDATE") return KeyVal("DATE-OBS");
  if (KeyName == "TOADSCAN") 
    {
      string sec  = "not yet implemented";
      return FitsKey(KeyName,sec);
    }
  if (KeyName == "TOADILLU") 
    {
      string sec  = "not yet implemented";
      return FitsKey(KeyName,sec);
    }
  if (KeyName == "TOADTYPE") return KeyVal("IMAGETYP");
  if (KeyName == "TOADOBJE") return KeyVal("TARGNAME");
  if (KeyName == "TOADCHIP") return KeyVal("DETECTOR");
  if (KeyName == "TOADBAND")
    {
      string filter = KeyVal("FILTNAM1");
      return FitsKey(KeyName,ToadBand(StringToUpper(filter)));	  
    }
  return FitsKey(KeyName,NOVAL);      
}

//*********************************************************************
FitsKey FitsHeader::HstNicmosFormat(const string &KeyName) const
{
  if (KeyName == "TOADNAMP") return FitsKey(KeyName,1);
  if (KeyName == "TOADINST") return KeyVal("INSTRUME");
  if (KeyName == "TOADPIXS") return FitsKey(KeyName,0.1);
  if (KeyName == "TOADFILT") return KeyVal("FILTER");
  if (KeyName == "TOADEXPO") return KeyVal("EXPTIME");
  if (KeyName == "TOADGAIN") return KeyVal("ADCGAIN");
  if (KeyName == "TOADRDON") return FitsKey(KeyName,0.72);
  if (KeyName == "TOADRASC") 
      {
	  double rastr=KeyVal("RA_TARG");
	  return FitsKey(KeyName,RaDegToString(rastr));
      }
  if (KeyName == "TOADDECL")
      {
	  double decstr=KeyVal("DEC_TARG");
	  return FitsKey(KeyName,DecDegToString(decstr));
      }
  if (KeyName == "TOADEQUI") return KeyVal("EQUINOX");
  if (KeyName == "TOADAIRM") return FitsKey(KeyName,0.0);
  if (KeyName == "TOADUTIM") return KeyVal("EXPSTART");
  if (KeyName == "TOADDATE") 
    {
      string keyval = KeyVal("DATE-OBS");
      int dd,mm,yy;
      if (sscanf(keyval.c_str(),"%d/%d/%d",&dd,&mm,&yy) == 3)
 	{ 
	  char date_charstar[64];
	  sprintf(date_charstar,"%d/%d/%d",dd,mm,yy+1900);
	  return FitsKey(KeyName,string(date_charstar));
	}      
    }
  if (KeyName == "TOADSCAN") 
    {
      string sec  = "not yet implemented";
      return FitsKey(KeyName,sec);
    }
  if (KeyName == "TOADILLU") 
    {
      string sec  = "not yet implemented";
      return FitsKey(KeyName,sec);
    }
  if (KeyName == "TOADTYPE") return KeyVal("IMAGETYP");
  if (KeyName == "TOADOBJE") return KeyVal("TARGNAME");
  if (KeyName == "TOADCHIP") return KeyVal("CAMERA");
  if (KeyName == "TOADBAND")
    {
      string filter = KeyVal("FILTER");
      return FitsKey(KeyName,ToadBand(StringToUpper(filter)));	  
    }
  return FitsKey(KeyName,NOVAL);      
}

//*********************************************************************
FitsKey FitsHeader::StandardFormat(const string &KeyName) const
{
  if (KeyName == "TOADINST") return KeyVal("INSTRUME");
  if (KeyName == "TOADPIXS") return KeyVal("PIXSCALE");
  if (KeyName == "TOADFILT") return KeyVal("FILTER");
  if (KeyName == "TOADEXPO") return KeyVal("EXPTIME");
  if (KeyName == "TOADGAIN") return KeyVal("GAIN");
  if (KeyName == "TOADRDON") return KeyVal("RDNOISE");
  if (KeyName == "TOADRASC") return KeyVal("RA");
  if (KeyName == "TOADDECL") return KeyVal("DEC");
  if (KeyName == "TOADEQUI") return KeyVal("EQUINOX");
  if (KeyName == "TOADAIRM") return KeyVal("AIRMASS");
  if (KeyName == "TOADUTIM") return KeyVal("UT");
  if (KeyName == "TOADDATE")
    {
      int dd,mm,yy;
      string keyval = KeyVal("DATE-OBS");
      if (sscanf(keyval.c_str(),"%d-%d-%d",&yy,&mm,&dd) == 3)
	{ 
	  char date_charstar[64];
	  sprintf(date_charstar,"%d/%d/%d",dd,mm,yy);
	  return FitsKey(KeyName,string(date_charstar));
	}
      } 
  if (KeyName == "TOADSCAN") 
    {
      string keyval = KeyVal("BIASSEC");
      int x0,x1,y0,y1,nx,ny;
      if (sscanf(keyval.c_str(),"[%d:%d,%d:%d]",&x0,&x1,&y0,&y1) == 4)
	{
	  nx = x1-x0+1;
	  ny = y1-y0+1;
	  char sec_charstar[64];
	  sprintf(sec_charstar,"[%d,%d;%d,%d]",x0,nx,y0,ny);
	  return FitsKey(KeyName,string(sec_charstar));
	}
    }
    if (KeyName == "TOADILLU") 
      {
	string keyval = KeyVal("DATASEC");
	int x0,x1,y0,y1,nx,ny;
	if (sscanf(keyval.c_str(),"[%d:%d,%d:%d]",&x0,&x1,&y0,&y1) == 4)
	  {
	    nx = x1-x0+1;
	    ny = y1-y0+1;
	    char sec_charstar[64];
	    sprintf(sec_charstar,"[%d,%d;%d,%d]",x0,nx,y0,ny);
	    return FitsKey(KeyName,string(sec_charstar));
	  }
      }
    if (KeyName == "TOADTYPE") return KeyVal("OBSTYPE");
    if (KeyName == "TOADOBJE") return KeyVal("OBJECT");
    if (KeyName == "TOADCHIP") return KeyVal("IMAGEID");
    if (KeyName == "TOADBAND")
      {
	string filter = KeyVal("FILTER");
	return FitsKey(KeyName,ToadBand(StringToUpper(filter)));	  
      }
    if (KeyName == "TOADNAMP") return FitsKey(KeyName,1);

    return FitsKey(KeyName,NOVAL);      
}

FitsKey FitsHeader::NotAlfoscFormat(const string &KeyName) const
{
  if (KeyName == "TOADINST") return KeyVal("INSTRUME");
  if (KeyName == "TOADPIXS") return FitsKey(KeyName,0.188);
  if (KeyName == "TOADFILT") return FitsKey(KeyName, string(KeyVal("FILTER")));
  if (KeyName == "TOADEXPO") return KeyVal("EXPTIME");
  if (KeyName == "TOADRASC")
    {
     
      double radeg = KeyVal("RA");
      while (radeg > 360) radeg -= 360;
      return FitsKey(KeyName,RaDegToString(radeg));
    }

  if (KeyName == "TOADDECL")
    {
      double decdeg = KeyVal("DEC");
      return FitsKey(KeyName,DecDegToString(decdeg));
    }

  if (KeyName == "TOADRDON") return FitsKey(KeyName,6);
  if (KeyName == "TOADGAIN") return FitsKey(KeyName,1.06);
  if (KeyName == "TOADEQUI") 
    {
      double epoch = 0;
      char dumb;
      string keyval = KeyVal("EQUINOX"); 
      if (sscanf(keyval.c_str(),"%c%lf",&dumb,&epoch) != 2)
	{cerr << "Needs epoch in correct NOT-ALFOSC format" << KeyVal("EQUINOX")<<endl;}
      return FitsKey(KeyName,epoch);
    }
  if (KeyName == "TOADAIRM") return KeyVal("AIRMASS");
  if (KeyName == "TOADUTIM")
    {
      double utdeci = KeyVal("UT");
      return FitsKey(KeyName,UtDeciToString(utdeci));
    }

  if (KeyName == "TOADDATE") 
    {
      int dd,mm,yy;
      string keyval = KeyVal("DATE-OBS");
      if (sscanf(keyval.c_str(),"%d-%d-%d",&yy,&mm,&dd) == 3)
 	{ 
	  char date_charstar[64];
	  sprintf(date_charstar,"%d/%d/%d",dd,mm,yy);
	  return FitsKey(KeyName,string(date_charstar));
	}
    } 
  if (KeyName == "TOADSCAN") return KeyVal("OVERSCAN");
/*    {
      string keyval = KeyVal("OVERSCAN");
      int x0,x1,y0,y1,nx,ny;
      if (sscanf(keyval.c_str(),"[%d:%d,%d:%d]",&x0,&x1,&y0,&y1) == 4)
 	{
 	  nx = x1-x0+1;
 	  ny = y1-y0+1;
 	  char sec_charstar[64];
 	  sprintf(sec_charstar,"[%d,%d;%d,%d]",x0,nx,y0,ny);
 	  return FitsKey(KeyName,string(sec_charstar));
 	}
    } */
  if (KeyName == "TOADILLU") 
    {

/*      return FitsKey(KeyName,string("[253,1284;536,1284]"));*/

      string keyval = KeyVal("TRIM");
      int x0,x1,y0,y1,nx,ny;
      if (sscanf(keyval.c_str(),"[%d:%d,%d:%d]",&x0,&x1,&y0,&y1) == 4)
        {
	  nx = x1-x0+1;
	  ny = y1-y0+1;
	  char sec_charstar[64];
	  sprintf(sec_charstar,"[%d,%d;%d,%d]",x0,nx,y0,ny);
	  return FitsKey(KeyName,string(sec_charstar)) ;
	} 
    }
  if (KeyName == "TOADTYPE") return KeyVal("IMAGETYP");
  if (KeyName == "TOADOBJE") return KeyVal("OBJECT");
  if (KeyName == "TOADCHIP") return FitsKey(KeyName,7);    
  if (KeyName == "TOADBAND")
    {
      string filter = KeyVal("FILTER");
      return FitsKey(KeyName,ToadBand(StringToUpper(filter)));	  
    }
  return FitsKey(KeyName,NOVAL);      
}
#endif

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
	echo "static AcceptorType AcceptorsArray [] = {" >> $@.new
	grep TYPE_SELECTOR ${telinst_sources}  | sed 's/${classname}/\&\2::Acceptor,/' >> $@.new
	echo "&Unknown::Acceptor};" >> $@.new
	echo '/* instanciation of templates IsOfKind */' >> $@.new
	grep TYPE_TESTER ${telinst_sources}  | sed 's/${classname}/template bool IsOfKind<\2>(const FitsHeader \&);/' >> $@.new
	$(CMTROOT)/mgr/cmt check_files $@.new $@

alltelinst.h : $(telinst_sources) fitstoad.cc
	\rm -f  $@.new
	echo " /* Automatically generated by make from" $(telinst_sources) "*/" >> $@.new
	echo " /* makefile at the end of fitstoad.cc */" >> $@.new
	echo '#ifndef ALLTELINST__H' >> $@.new
	echo '#define ALLTELINST__H' >> $@.new
	grep TYPE_TESTER ${telinst_sources} | sed 's/${classname}/class \2;/'   >> $@.new
	echo '#endif' >> $@.new
	$(CMTROOT)/mgr/cmt check_files $@.new $@

fitstoad.cc : $(telinst_sources)
	  touch $@
#endif




