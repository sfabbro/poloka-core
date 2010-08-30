#include <math.h>
#include <string.h>

#include "fitsimage.h"
#include "gtransfo.h"
#include "wcsutils.h"
#include "astroutils.h"
#include "frame.h"
#include "fitstoad.h"
#include "basestar.h" /* for MEMPIX2DISK */
#include "stringlist.h"
#include "imageutils.h"

#ifndef M_PI
#define     M_PI            3.14159265358979323846  /* pi */
#endif


/* some history about Toads and WCS:
   to begin with, I (Pierre astier) have to admit I thought
   that all the projection stuff was useless for the size of 
   images we usually handle, and that a linear approximation
   of the (de)projection was introducing less error than the 
   uncertainty of USNO astrometry. Roughly speaking, this
   is true, but finding the correct linear transformation
   from the CD's, CRPIX's, and CRVAL's is not as simple as
   taking the CD's and CRPIx's for the linear part, and
   taking the projection as just a contraction by cos(dec)
   along x. This model will be called LinearProjApproximation.
   The fact that it gives wrong results if interpreted by
   official WCS readers (such as WCStools) is not a serious
   burden for us, as long as all astrometry goes through
   this package. from now on (may 03), we will anyway encode
   normal WCS's.
   
   About optical distortions, our first implementation
   consisted in fitting a cubic correction on top of the
   linear stuff above, the cubic correction being applied 
   before:

      read the linear part (WCSTransfoFromHeader ), 
      read the cubic part (WCS3TransfoFromHeader)
      and convert from pixels to a,d by using
      Pix2RaDec = GtransfoCompose(&linPix2RaDec, &cubCorr);

   A lot of headers have been encoded with this convention, so we have
   to keep the code that reads and interprets them. These headers have
   WCS3 keys. Then there is a bug that lasted a few weeks in march 03,
   where we have headers with WCS encoded using official WCS tools,
   but also have WCS3keys (a bug in WCSCopy, where the WCS3 keys were
   not removed if present in the target and not in the source). By
   chance those headers all miss CD keys (they use CDELT instead),
   which enables to trigger the "official" behaviour:

   From now (may 03) on, we will encode WCS's the normal way, with corrections
   on top of the linear + deprojection standard steps. Following SWarp and 
   Atrometrix, we implement this correction between these 2 steps. The
   encoding follows the proposal of "paper IV" :
   http://www.atnf.csiro.au/people/Mark.Calabretta/WCS/dcs_20020312.ps.gz

*/




typedef enum { LinearProjApproximation, TanProjWithoutCorr, TanProjWithCorr, 
       Unknown, NoWCS} WCSKind;

static bool HasWCS3Correction(const FitsHeader &Head);

static void RemoveWCS3Correction(FitsHeader &Head);

/*************** reading of LinearProjApproximation **********************/


#ifdef STORAGE
int WCSTransfo2Header(const string &FitsImageName, const GtransfoLin &Pix2RaDec)
{
   FitsHeader header(FitsImageName, RW);
   if (!header.IsValid()) return 0;
   return WCSTransfo2Header(header,Pix2RaDec);
} 

// This put the right keywords in the header of the image, for TAN projection.
// one trick is that Pix2RaDec actually converts to Ra and Dec: when Dec != 0
// the matrix of the transformation is not an isometry (because of cos(dec)).
// The WCS conventions want the CD matrix to be an isometry
// if the pixels are squares. so we have to account for it
// when writing the CD stuff

int WCSTransfo2Header(FitsHeader &header, const GtransfoLin &Pix2RaDec)
{ 

   double ra,dec;
   GetRaDecInDeg(header, ra, dec);
   double cosdec = cos(dec*M_PI/180.);

   // type of WCS 
   header.AddOrModKey("CTYPE1", "RA---TAN");
   header.AddOrModKey("CTYPE2", "DEC--TAN");

   double ra_tan, dec_tan;
   if (IsOfKind<Int4Wfc>(header)) 
     {
       double xapoff = header.KeyVal("XAPOFF");
       double yapoff = header.KeyVal("YAPOFF");
       ra_tan = ra - yapoff/cosdec; 
       dec_tan = dec + xapoff;
     }
   else
     {
     ra_tan = ra;
     dec_tan = dec;
     }

   // The coefficients of the matrix. cosdec is not a bug !

   header.AddOrModKey("CD1_1", Pix2RaDec.A11()*cosdec, " Transformation matrix : a11 (deg/pix)");
   header.AddOrModKey("CD2_2", Pix2RaDec.A22(),        " Transformation matrix : a22 (deg/pix)");
   header.AddOrModKey("CD1_2", Pix2RaDec.A12()*cosdec, " Transformation matrix : a12 (deg/pix)");
   header.AddOrModKey("CD2_1", Pix2RaDec.A21(),        " Transformation matrix : a21 (deg/pix)");


   header.AddOrModKey("CRVAL1", ra_tan, "Rotator x coordinates (deg)");
   header.AddOrModKey("CRVAL2", dec_tan, "Rotator y coordinates (deg)");


   Point OptCentreDeg(ra_tan, dec_tan);
   GtransfoLin Deg2Pix;
   Deg2Pix = Pix2RaDec.invert();
   Point OptCentrePix = Deg2Pix(OptCentreDeg);
   
   header.AddOrModKey("CRPIX1", OptCentrePix.x + MEMPIX2DISK, "Rotator x coordinates (pix)");
   header.AddOrModKey("CRPIX2", OptCentrePix.y + MEMPIX2DISK, "Rotator y coordinates (pix)");
   header.AddOrModKey("FITSCONV",true,"CRPIX1 and CRPIX2 in standard FITS pixel coords (1,1)");
   return 1;
} 

#endif /* STORAGE */



// To read a transformation written in the header of the image "FitsImageName".
//    Return the transformation
bool WCSLinTransfoFromHeader(const string &FitsImageName, GtransfoLin &Pix2RaDec)
{
  FitsHeader header(FitsImageName);
  return  WCSLinTransfoFromHeader(header, Pix2RaDec);
}


bool WCSLinTransfoFromHeader(const FitsHeader& Header, GtransfoLin &Pix2RaDec)
{
   Gtransfo *wcs;
   if (!WCSFromHeader(Header,wcs)) return false;
   Frame frame(Header);
   // scale the derivation step with the size of the frame. Not critical hopefully
   Pix2RaDec = wcs->LinearApproximation(frame.Center(), sqrt(frame.Area())/3.);
   delete wcs;
   return true;
}

static bool OldLinWCSFromHeader(const FitsHeader& Header, GtransfoLin &Pix2RaDec)
{
   if (!HasLinWCS(Header)) return false;
   string ctype1 = Header.KeyVal("CTYPE1");
   if ( ctype1 != "RA---TAN")
     {
       cerr << " bad WCS type :" << ctype1 << "don\'t know yet how to read this " << endl;
       return false;
     }
 
   // ("CTYPE2", "DEC--TAN");
 
   double ra,dec;
   GetRaDecInDeg(Header, ra, dec);
   double cosdec = cos(dec*M_PI/180.);
   double ra_rot = Header.KeyVal("CRVAL1");
   double dec_rot = Header.KeyVal("CRVAL2");
   double x_rot = Header.KeyVal("CRPIX1");
   double y_rot = Header.KeyVal("CRPIX2");
   /* we want a WCS to be interpreted in TOADs ccordinates, i.e.
   starting at (0,0). So
    - if the WCS was NOT written by Toads, it misses cubic corrections
    - if it was wriiten by Toads, then we have to shift only if the WCS
    was written after the "normalization" */
   if (!HasWCS3Correction(Header) ||  Header.HasKey("FITSCONV"))
     {
       x_rot -= MEMPIX2DISK;
       y_rot -= MEMPIX2DISK;
     }
       
   double a11,a12,a21,a22;
   // aij are sometime absent (e.g. from swarp resampled images)
   if (Header.HasKey("CD1_1") && Header.HasKey("CD1_2") &&
       Header.HasKey("CD2_1") && Header.HasKey("CD2_2"))
     {
       a11 = double(Header.KeyVal("CD1_1"))/cosdec;
       a12 = double(Header.KeyVal("CD1_2"))/cosdec;
       a21 = Header.KeyVal("CD2_1");
       a22 = Header.KeyVal("CD2_2");
     }
   else if (Header.HasKey("CDELT1") && Header.HasKey("CDELT2"))
     {
       a11 = double(Header.KeyVal("CDELT1"))/cosdec;
       a22 = Header.KeyVal("CDELT2");
       a12 = 0;
       a21 = 0;
     }
   else return false;

   GtransfoLin rotateur(ra_rot, dec_rot, -a11, -a12, -a21, -a22);
   Point rotpix(x_rot, y_rot);
   Point delta = rotateur.apply(rotpix);
   
   Pix2RaDec = GtransfoLin(delta.x, delta.y, a11, a12, a21, a22);
   return true;
}

/********** handling of TanProjWithCorr, and TanProjWithoutCorr **********/

/* check that there are DJ or PV keys in the header ******/
static bool HasTanCorrection(const FitsHeader &Head, string &CorrTag, 
			     FitsKeyArray *CorrKeys = NULL)
{
  if (Head.HasKey("DJ1_0")) CorrTag = "DJ"; // only Toads
  else if (Head.HasKey("PV1_0")|| Head.HasKey("PV1_1"))   // swarp & co:
    CorrTag = "PV"; 
  else if (Head.HasKey("DV1_0") || Head.HasKey("DV1_1")) // future names
    CorrTag = "DV"; 
  else return false;
  FitsKeyArray array;
  string selector = CorrTag+"?_#"; // maches e.g. DJ2_10
  int count = Head.KeyMatch(selector.c_str(), (CorrKeys)? 
			    (*CorrKeys) : array );
  if (count >= 12) 
    /* min 6 for x and 6 for y for a correction beyond 
       first order, but some may be 0? */
    return true;
  else
    {
    if (CorrKeys) CorrKeys->clear();
    return false;
    }
}


bool TanLinWCSFromHeader(const FitsHeader &Head, TanPix2RaDec &TanWcs, 
			 const bool Warn)
{
  string ctype1 = Head.KeyVal("CTYPE1");
  string ctype2 = Head.KeyVal("CTYPE2");
  if ((ctype1 != "RA---TAN" || ctype2 != "DEC--TAN") && (ctype1 != "RA---TNX" || ctype2 != "DEC--TNX"))
    {
      if (Warn) cerr << " No TAN WCS in " << Head.FileName () << endl;
      return false;
    }
  if (!HasLinWCS(Head))
    {
      if (Warn) cerr << "No Lin WCS in " << Head.FileName() << endl;
      return false;
    }
  double crpix1, crpix2;
  crpix1 = Head.KeyVal("CRPIX1",true);
  crpix2 = Head.KeyVal("CRPIX2", true);
  /* this WCS is expressed with the bottom lefmost pixel = 1. 
     Internally, we use 0 */
  crpix1 -= MEMPIX2DISK;
  crpix2 -= MEMPIX2DISK;
  GtransfoLin pix2ThetaPhi = GtransfoLinShift(-crpix1,-crpix2);
  if (Head.HasKey("CD1_1"))
    {
      double cd11 = Head.KeyVal("CD1_1",true);
      double cd12 = Head.KeyVal("CD1_2",true);
      double cd21 = Head.KeyVal("CD2_1",true);
      double cd22 = Head.KeyVal("CD2_2",true);
      GtransfoLin rot(0,0, cd11,cd12,cd21,cd22);
      pix2ThetaPhi = rot*pix2ThetaPhi;
    }
  else if (Head.HasKey("CDELT1"))
    {
      double cdelt1 = Head.KeyVal("CDELT1",true);
      double cdelt2 = Head.KeyVal("CDELT2",true);
      GtransfoLinScale scale(cdelt1,cdelt2);
      pix2ThetaPhi = scale*pix2ThetaPhi;
    }
  // following WCS conventions, all angles are expressed in degrees
  double crval1 = Head.KeyVal("CRVAL1", true);
  double crval2 = Head.KeyVal("CRVAL2", true);
  TanWcs = TanPix2RaDec(pix2ThetaPhi, Point(crval1, crval2));
  return true;
}


typedef struct
{
  const char *fitsName;
  const char *coeffName;
}   DistortionMapping;


/* pretty obsolete : unused since 2002 */
static const DistortionMapping WCS3Names[] =
  {
    {"DX","dx"   },
    {"DY",  "dy"   },
    {"A11", "a11"  },
    {"A12", "a12"  },
    {"A21", "a21"  },
    {"A22", "a22"  },
    {"A1X2","a1x2" },
    {"A1XY","a1xy" },
    {"A1Y2","a1y2" },
    {"A2X2","a2x2" },
    {"A2XY","a2xy" },
    {"A2Y2","a2y2" },
    {"A1X3","a1x3" },
    {"A1X2Y", "a1x2y"},
    {"A1XY2", "a1xy2"},
    {"A1Y3", "a1y3" },
    {"A2X3","a2x3" },
    {"A2X2Y", "a2x2y"},
    {"A2XY2", "a2xy2"},
    {"A2Y3","a2y3" },
    {"end",NULL}    /* end marker */
  };


/* This exists in some old images , but is "bugged":
   for the y transfo, x and y "should" be swapped,
   according to the convention : the correct version 
   is the DV one below. 
   Still in the code for backward compatibility.
*/
static const DistortionMapping DJNames [] = 
  {
    {"1_0",  "dx"   },
    {"2_0",  "dy"   },
    {"1_1", "a11"  },
    {"1_2", "a12"  },
    {"2_1", "a21"  },
    {"2_2", "a22"  },
    {"1_4","a1x2" },  /* 1_3 and 2_3 are coeffs for r */
    {"1_5","a1xy" },
    {"1_6","a1y2" },
    {"2_4","a2x2" },
    {"2_5","a2xy" },
    {"2_6","a2y2" },
    {"1_7","a1x3" },
    {"1_8", "a1x2y"},
    {"1_9", "a1xy2"},
    {"1_10", "a1y3" },
    {"2_7","a2x3" },
    {"2_8", "a2x2y"},
    {"2_9", "a2xy2"},
    {"2_10","a2y3" },
    {"end",NULL} /* end marker */
  };

/* 
   This mapping of the polynomial coefficients was proposed in a draft
   paper about distorions handling in WCS's (Calbretta and Greisen),
   but it disappeared in a later version.  It is no longer in the
   current version (by May 03) of the paper (which is really a draft).
   The actual mapping adopted is the one from Astrometrix (described
   in the documentation of the package), see
   http://terapix.iap.fr/soft/.

   The difference with the "DJ" stuff is that for the second polynomial
   (the DV2_# keys), the role of x and y are swapped with respect to the 
   first polynomial.

   This mapping essentially exists with "PV" as the key head. although
   the Calabretta and Greisen paper suggests to use "DV". swarp 2.10 
   searches for PV but not DV.

*/


static const DistortionMapping DVNames[] =   
  {
    {"1_0",  "dx"},
    {"2_0",  "dy"},
    {"1_1", "a11"},
    {"1_2", "a12"},
    {"2_1", "a22"},
    {"2_2", "a21"},
    {"1_4", "a1x2"},  /* 1_3 and 2_3 are coeffs for r */
    {"1_5", "a1xy"},
    {"1_6", "a1y2"},
    {"2_4", "a2y2"},
    {"2_5", "a2xy"},
    {"2_6", "a2x2"},
    {"1_7", "a1x3"},
    {"1_8", "a1x2y"},
    {"1_9", "a1xy2"},
    {"1_10","a1y3"},
    {"2_7", "a2y3" },
    {"2_8", "a2xy2"},
    {"2_9", "a2x2y"},
    {"2_10","a2x3" },
    {"end",NULL} /* end marker */
  };



/* routine the decodes a correction polynomial , independent from the 
actual WCS type */
static GtransfoPoly ReadWCSCorrection(const FitsHeader &Head,
				      const string &corrTag,
				      FitsKeyArray &corrections)
{
  /* first loop on Fits keys to identify the corrections scheme 
     and evaluate correction degree */
  int maxj = 0;
  for (unsigned k = 0; k < corrections.size(); ++k)
    {
      FitsKey &key = corrections[k];
      string keyName = key.KeyName();
      // decoding routine
      // analyse the keyname to guess the degree.
      int i,j;
      char toto[8];
      if (sscanf(keyName.c_str(),"%2s%d_%d",toto,&i,&j) == 3)
	if (j>maxj) maxj = j;
    }
  // if the highest j is 6, then second degree is enough
  GtransfoPoly corr((maxj==6) ? 2 : 3);
  const DistortionMapping * assoc = 0;
  if (corrTag == "DJ") assoc = DJNames;
  else if (corrTag == "DV" || corrTag == "PV") assoc = DVNames;
  else if (corrTag == "WCS3") assoc = WCS3Names;
  unsigned tagLen = corrTag.length();
  for (unsigned k = 0; k < corrections.size(); ++k)
    {
      FitsKey &key = corrections[k];
      /* remove PV, DJ or DV at the beginning of key names (they are
	 present by construction, and absent from the "assoc" tables)
      */
      const char *fitsName = key.KeyName().c_str()+tagLen;
      double val = key; 
      size_t j;
      for (j=0; assoc[j].coeffName ; ++j)
	{
	  if (strcmp(assoc[j].fitsName,fitsName) == 0) break;
	}
      if (!assoc[j].coeffName) // reached the end
	{
	  cerr << " ReadWCSCorrection : no item named " 
	       << fitsName << '(' << key.KeyName() << ')' << endl;
	  continue;
	}
      corr.Coeff(assoc[j].coeffName) = val;
    }
  return corr;
}


bool TanWCSFromHeader(const FitsHeader &Head, TanPix2RaDec &TanWcs,
		      const bool Warn)
{
  // read the linear (CD's, CRPIX's) and (de)projection (CRVAL's) parts
  if (!TanLinWCSFromHeader(Head, TanWcs,Warn)) return false;
  // now go for non-linear distortions :
  FitsKeyArray corrections;
  /* there has been a whole history of those distortions. Their
     meanings depend on the "tag" value: */
  string corrTag; // PV, DJ, DV...
  if (HasTanCorrection(Head, corrTag, &corrections))
    {
      GtransfoPoly  corr = ReadWCSCorrection(Head,corrTag, corrections);
      TanWcs.SetCorrections(&corr);
    }
  return true;
}

// remove non linear corrections that apply to Tangent point WCS's
static void RemoveTanCorrection(FitsHeader &Head)
{
  static const string keyHead[] = {"DJ","PV","DV"};
  static const DistortionMapping* assoc[] = {DJNames, DVNames, DVNames};
  static int nHead = sizeof(keyHead)/sizeof(keyHead[0]);;
  for (int i=0; i<nHead; ++i)
    {      
      for (unsigned k =0; assoc[i][k].coeffName; ++k)
	{
	  string keyName = keyHead[i]+assoc[i][k].fitsName;
	  if (Head.HasKey(keyName)) Head.RmKey(keyName);
	}
    }
}
  
#define TOADS_WCS_VERSION "1.2"

int TanWCS2Header(const string &FitsImageName, const TanPix2RaDec &TanWcs)
{
   FitsHeader header(FitsImageName, RW);
   if (!header.IsValid()) return 0;
   TanWCS2Header(header,TanWcs);
   return 1;
} 


void TanWCS2Header(FitsHeader &Head, const TanPix2RaDec &TanWcs)
{
  // clear a few keys:
  // old toads WCS distortion keys:
  RemoveWCS3Correction(Head);
  // new Tan plane distortion keys:
  RemoveTanCorrection(Head);
  // wrong tags were introduced at some point. remove them to avoid confusion.
  // but this happened to be (almost) harmless
  if (Head.HasKey("CRTYPE1")) Head.RmKey("CRTYPE1");
  if (Head.HasKey("CRTYPE2")) Head.RmKey("CRTYPE2");

  // filling
  Head.AddOrModKey("CTYPE1","RA---TAN");
  Head.AddOrModKey("CTYPE2","DEC--TAN");
  Point crpix = TanWcs.CrPix();
  Head.AddOrModKey("CRPIX1",crpix.x+MEMPIX2DISK);
  Head.AddOrModKey("CRPIX2",crpix.y+MEMPIX2DISK);
  Point crval = TanWcs.TangentPoint();
  Head.AddOrModKey("CRVAL1",crval.x);
  Head.AddOrModKey("CRVAL2",crval.y);
  const GtransfoLin linPart = TanWcs.LinPart();
  Head.AddOrModKey("CD1_1",linPart.A11());
  Head.AddOrModKey("CD1_2",linPart.A12());
  Head.AddOrModKey("CD2_1",linPart.A21());
  Head.AddOrModKey("CD2_2",linPart.A22());
  // since we write the CD stuff down, we should remove the 
  // CDELT which may not coexist with CD's
  if (Head.HasKey("CDELT1")) Head.RmKey("CDELT1");
  if (Head.HasKey("CDELT2")) Head.RmKey("CDELT2");

  // corrections ?
  const GtransfoPoly *corr = TanWcs.Corr();
  if (corr)
    {
      /* We chose to follow the PV encoding, because 
	 I am sure (P.A) that swarp reads it (version 2.10)
      */
      const DistortionMapping *assoc = DVNames;
      for (unsigned k=0; assoc[k].coeffName; ++k)
	{
	  string keyName = "PV"+string(assoc[k].fitsName);
	  const char *coeffName = assoc[k].coeffName;
	  if (!corr->HasCoeff(coeffName)) continue;
	  double val = corr->Coeff(coeffName);
	  Head.AddOrModKey(keyName, val);
	}
    }
  Head.AddOrModKey("WCSVERS", TOADS_WCS_VERSION);
}

/************** returns the WCS (of the right kind) ***********/

static bool WCS3TransfoFromHeader(const FitsHeader& Header, 
				  GtransfoPoly &CubicCorr, 
				  const bool Warn);



/*! This routine combines a "Linear" WCS (if any) with a cubic correction
  (if any). The result is exact if the image was matched using matchusno.
  For WCS encoded by other tools (such as SWarp), the result is now also exact,
  using the TanPix2RaDec Gtransfo */
bool WCSFromHeader(const FitsHeader &Head, Gtransfo* &Pix2RaDec)
{
  Pix2RaDec = NULL;
  if (!Head.IsValid() || !HasLinWCS(Head))
    {
      cerr << " do not find the expected WCS in " << Head.FileName() << endl;
      return false;
    }
  /* there was a bug in CopyWCS : the cubic part of the WCS was not
     deleted, when copying from a header that did not have any cubic 
     correction. In practise, only WCS from images from SWarp were copied,
     letting WCS3 keys in the header. But since SWarp headers do not contain 
     the CD matrix, but CDELT keys instead, there is a possible correction:
     ignore the cubic terms if CD is absent. When this bug was corrected
     we introduced the key WCSVERS in the headers, to simplify the handling
     of possible  present or future bugs .
   */
     
  if ( (Head.HasKey("CDELT1") && !Head.HasKey("CD1_1")) 
       || Head.HasKey("WCSVERS") 
       ||!HasWCS3Correction(Head))
    {
       TanPix2RaDec wcs;
       bool ok = TanWCSFromHeader(Head, wcs);
       Pix2RaDec = (ok)? wcs.Clone() : NULL;
       return ok;
    }
  else
    {
      GtransfoLin linPix2RaDec;
      if (!OldLinWCSFromHeader(Head,linPix2RaDec)) return false;
      GtransfoPoly cubCorr(3);
      if (WCS3TransfoFromHeader(Head, cubCorr, /*warn = */ false))
	{
	  Pix2RaDec = GtransfoCompose(&linPix2RaDec, &cubCorr);
	}
      else Pix2RaDec = new GtransfoLin(linPix2RaDec);
      return true;
    }
  return true;
}

bool WCSFromHeader(const string &FitsName, Gtransfo* &Pix2RaDec)
{
  FitsHeader head(FitsName);
  return WCSFromHeader(head,Pix2RaDec);
}


static const char* WCSKeys[] = 
  {"CTYPE1","CTYPE2","CRVAL1","CRVAL2","CRPIX1","CRPIX2",
   "CD1_1","CD1_2","CD2_1","CD2_2",
   "CDELT1","CDELT2","CUNIT1","CUNIT2","EQUINOX","RADECSYS"};

static int NWCSKeys = sizeof(WCSKeys)/sizeof(WCSKeys[0]);

//! returns key names that describe the WCS in this header
void WCSKeyList(const FitsHeader &Head, StringList &KeyNames)
{
  KeyNames.clear();
  for (int i=0; i< NWCSKeys; ++i) 
    {
    const char *keyName = WCSKeys[i];
    if (Head.HasKey(keyName)) KeyNames.push_back(keyName);
    }
  FitsKeyArray keys;
  string dummy;
  if (!HasTanCorrection(Head, dummy,  &keys)) // try modern corr scheme
    Head.KeyMatch("WCS3*", keys); // or old one.
  for (unsigned k =0; k<keys.size(); ++k)
    {
      FitsKey &key = keys[k];
      KeyNames.push_back(key.KeyName());
    }
}


//#define DEBUG

bool CopyWCS(const FitsHeader &FromHeader, FitsHeader &ToHeader)
{
  if (!HasLinWCS(FromHeader)) return false;
  StringList fromKeys;
  WCSKeyList(FromHeader, fromKeys);
  StringList toKeys;
  WCSKeyList (ToHeader, toKeys);
  /* rather than removing all keys from To and inserting all keys 
     from From, we update and remove: this preserves original 
     locations in the header */
  for (StringCIterator i = toKeys.begin(); i != toKeys.end(); ++i)
    {
      const std::string &toKey = *i;
      if (fromKeys.Locate(toKey) == fromKeys.end()) 
	{
	  ToHeader.RmKey(toKey);
#ifdef DEBUG
	  cout << " removing " << toKey << " from " << ToHeader.FileName() 
	       << endl;
#endif
	}
    }
  for (StringCIterator i = fromKeys.begin(); i != fromKeys.end(); ++i)
    {
      const std::string &fromKey = *i;
      FromHeader.CopyKey(fromKey, ToHeader);
#ifdef DEBUG
      cout << " inserting " << fromKey << " into " << ToHeader.FileName() << endl;
#endif
    }
  return true;
}



bool WCSTransfoBetweenHeader(const FitsHeader &Header1, const FitsHeader &Header2, GtransfoLin &Transfo1to2)
{
  if ( HasLinWCS(Header1) && HasLinWCS(Header2))
    {
      GtransfoLin transfo1;
      GtransfoLin transfo2;
      if (WCSLinTransfoFromHeader(Header1,transfo1) && WCSLinTransfoFromHeader(Header2,transfo2))
	Transfo1to2 = transfo2.invert()*transfo1;
      return true;
    }
  return false;
}

double PixelSize(const FitsHeader& Head)
{
  double sizx,sizy;
  GetPixelSize(Head,sizx,sizy);
  return (sizx+sizy)/2.;
}

double sqr(const double x) {return x*x;}

void GetPixelSize(const FitsHeader& Head, double &SizeX, double &SizeY)
{
  if (!HasLinWCS(Head)) {SizeX=SizeY=Head.KeyVal("TOADPIXS");return;}
  if (Head.HasKey("CDELT1") && Head.HasKey("CDELT2"))
    {
      SizeX = double(Head.KeyVal("CDELT1")) * 3600;
      SizeY = double(Head.KeyVal("CDELT2")) * 3600;
      return;
    }
  Gtransfo *Pixtodeg;
  WCSFromHeader(Head, Pixtodeg);
  double nx = Head.KeyVal("NAXIS1");
  double ny = Head.KeyVal("NAXIS2");
  double ra0,dec0,ra1,dec1,ra2,dec2,ra3,dec3;
  Pixtodeg->apply(0,ny/2,ra0,dec0);
  Pixtodeg->apply(nx,ny/2,ra1,dec1);
  Pixtodeg->apply(nx/2,0,ra2,dec2);
  Pixtodeg->apply(nx/2,ny,ra3,dec3);
  double dec = DecStringToDeg(string(Head.KeyVal("TOADDECL")));
  double cosdec = cos(M_PI*dec/180);
  SizeX = sqrt( sqr((ra1-ra0)*cosdec) + sqr(dec1-dec0) )*3600./nx;
  SizeY = sqrt( sqr((ra3-ra2)*cosdec) + sqr(dec3-dec2) )*3600./ny;
  delete Pixtodeg;
}

double Arcmin2Area(const Frame &aFrame,const FitsHeader &Header)
{
  Gtransfo* pix2radec;
  if (WCSFromHeader(Header, pix2radec))
    {
      Frame degframe = ApplyTransfo(aFrame,*pix2radec);
      delete pix2radec;
      double ra,dec;
      RaDecFromWCS(Header,ra,dec);
      double cosdec = cos(M_PI*dec/180);
      return degframe.Area()*3600*cosdec;
    }
  return aFrame.Area()*sqr(PixelSize(Header))/3600.;  
}

void RaDecFromWCS(const FitsHeader &Header, double &Ra, double &Dec)
{
  Gtransfo *pix2RaDec = 0;
  if (WCSFromHeader(Header, pix2RaDec))
    {
      double xc = double(Header.KeyVal("NAXIS1"))/2.;
      double yc = double(Header.KeyVal("NAXIS2"))/2.;
      pix2RaDec->apply(xc,yc,Ra,Dec);
      delete pix2RaDec;
    }
   
}

double Arcmin2Area(const FitsHeader &Header)
{Frame full(Header); return Arcmin2Area(full,Header);}

double Arcmin2Overlap(const FitsHeader& Head1,const FitsHeader& Head2)
{
  //check overlap with WCS 
  GtransfoLin transfo1to2;
  if (WCSTransfoBetweenHeader(Head1,Head2,transfo1to2))
    {
      Frame frame1(Head1); // boundaries of the image
      Frame frame2(Head2); // boundaries of the image
      Frame frame1in2 = ApplyTransfo(frame1,transfo1to2); //assume simple rotation      
      frame1in2 *= frame2;      
      if ((frame1in2.xMin < frame1in2.xMax) &&
	  (frame1in2.yMin < frame1in2.yMax) &&
	  (frame1in2.Area() > 10) )// assume that bad WCS transfo gives you a too small area
	{
	  //cout << frame1in2 << endl;
	  return Arcmin2Area(frame1in2,Head2);
	} 
      return 0;
    }

  //if no WCS or bad one, check overlap brutally by assuming RA and DEC 
  // at center and a square area around it
  double siz1,siz2,ra1,dec1,ra2,dec2;
  RaDec2000(Head1,ra1,dec1); //in deg
  RaDec2000(Head2,ra2,dec2);
  double nx=Head1.KeyVal("NAXIS1");
  double ny=Head1.KeyVal("NAXIS2");
  double npix1 = min(nx,ny);
  nx=Head2.KeyVal("NAXIS1");
  ny=Head2.KeyVal("NAXIS2");
  double npix2 = min(nx,ny);
  siz1= PixelSize(Head1)*npix1/2; //in arcsec
  siz2= PixelSize(Head2)*npix2/2; //in arcsec
  double cosdec1=cos(M_PI*dec1/180.);
  double cosdec2=cos(M_PI*dec2/180.);
  ra1 *= 3600;
  dec1 *= 3600;
  ra2 *= 3600;
  dec2 *= 3600;
  Frame frame1((ra1-siz1)*cosdec1,dec1-siz1,(ra1+siz1)*cosdec1,dec1+siz1);
  Frame frame2((ra2-siz2)*cosdec2,dec2-siz2,(ra2+siz2)*cosdec2,dec2+siz2);
  Frame intersect = frame1 * frame2;
  double area = 0.;
  if ((intersect.xMin < intersect.xMax) && (intersect.yMin < intersect.yMax))
    area=intersect.Area()*3600;
  return area;
}
  
// We write here the 3rd order polynomial transformation from pix to Ra and Dec.
// This corresponds by no mean to an official WCS implementation.



// To read a (cubic) transformation written in the header of the image "FitsImageName".
//    Return the transformation
static bool read_old_cubic_corr(const FitsHeader& Header, 
				GtransfoPoly &CubicCorr,
				const bool Warn,
				const string &KeyHead)
{
  bool correct = true;
  for (unsigned k=0; WCS3Names[k].coeffName; ++k)
    {
      string keyName = KeyHead+string(WCS3Names[k].fitsName);
      //       const char *keyName = sKeyName.c_str();
      if (!Header.HasKey(keyName))
	{
	  if (Warn) 
	    cerr << " read_old_cubic_corr :: no key named " 
		 <<  keyName << endl;
	  correct = false;
	  continue;
	}
      double val = Header.KeyVal(keyName);
      CubicCorr.Coeff(WCS3Names[k].coeffName) = val;
    }
  return  correct;
}



bool WCS3TransfoFromHeader(const FitsHeader& Header, GtransfoPoly &CubicCorr, 
			   const bool Warn)
{
  /* due to a (temporary) bug, there exists a set of images with OLD WCS's
  (WCS3 keys) encoded with the new (i_j) scheme. So search for both of them. 
  */
  if (read_old_cubic_corr(Header, CubicCorr, Warn, "WCS3"))
    return true;
  return read_old_cubic_corr(Header, CubicCorr, Warn, "DJ"); 
}


static void RemoveWCS3Correction(FitsHeader &Head)
{  
  for (unsigned k=0; WCS3Names[k].coeffName; ++k)
    {
      string keyName("WCS3"+string(WCS3Names[k].fitsName));
      if (Head.HasKey(keyName)) Head.RmKey(keyName);
    }
  // due to a (temporary) bug, there are old WCS's (with WCS3 keys)
  // encoded with  the new (i_j) tagging scheme. So remove them here:
  for (unsigned k=0; DJNames[k].coeffName; ++k)
    {
      string keyName("DJ"+string(DJNames[k].fitsName));
      if (Head.HasKey(keyName)) Head.RmKey(keyName);
    }

}



bool WCS3TransfoFromHeader(const string &FitsImageName, GtransfoPoly &Pix2RaDec,
			   const bool Warn)
{
  FitsHeader header(FitsImageName);
  return  (header.IsValid() && WCS3TransfoFromHeader(header, Pix2RaDec, Warn));
}

static bool HasWCS3Correction(const FitsHeader &Head)
{
  GtransfoPoly cubCorr(3);
  return WCS3TransfoFromHeader(Head, cubCorr, /*warn = */ false);
}


/********** Routines that analyze WCS ***************/

// check the presence of basic keys
bool HasLinWCS(const FitsHeader &Header)
{
  return( (Header.HasKey("CTYPE1") 
	  && Header.HasKey("CTYPE2") 
	  && Header.HasKey("CRVAL1") 
	  && Header.HasKey("CRVAL2") 
	  && Header.HasKey("CRPIX1") 
	  && Header.HasKey("CRPIX2"))

	  && (( Header.HasKey("CD1_1") && Header.HasKey("CD1_2") 
		&& Header.HasKey("CD2_1") && Header.HasKey("CD2_2") )
	      ||
	      (Header.HasKey("CDELT1") &&  Header.HasKey("CDELT1")))
	  );
}


bool UpdateRaDec(FitsHeader &Header)
{
  Gtransfo *pix2RaDec;
  if (!WCSFromHeader(Header, pix2RaDec)) return false;
  Frame frame(Header);
  Point middleRaDec = pix2RaDec->apply(frame.Center());
  Header.AddOrModKey("TOADRASC",RaDegToString(middleRaDec.x).c_str(),"Updated RA from WCS and usable part");
  Header.AddOrModKey("TOADDECL",DecDegToString(middleRaDec.y).c_str(),"Updated DEC from WCS and usable part");
  delete pix2RaDec;
  return true;
}

//////////////////////////////////////////////////



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


//! a handy routine to compute a WCS given a RaDec reference, and a possible rotation and flip.     
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
  cout << " Found that one: " << endl;
  cout << WCS << endl;
  return true;
}



// include alltelinst here
typedef bool (*GuessLinWCS_Type)(const FitsHeader &Head, TanPix2RaDec &Guess);
static std::map<string,GuessLinWCS_Type> WCS_Functions;
static int Add2GuessWCSMap(const string &Name, GuessLinWCS_Type Function) {
  WCS_Functions[Name]=Function;
  return 0;
}


#define USE_WCS
#ifdef VIRTUAL_INSTRUMENTS
#undef VIRTUAL_INSTRUMENTS
#endif
#include "alltelinst.cc"


bool GuessLinWCS(const FitsHeader &Header, TanPix2RaDec &Guess)
{
  string name = TelInstName(Header);
  std::map<string,GuessLinWCS_Type>::const_iterator iter = WCS_Functions.find(name);
  if (iter!=WCS_Functions.end()) {
    cout << "yes, the instrument " << name << " has a specific WCS function" << endl;
    if ( iter->second(Header,Guess) )
      return true;
  }
  
  cout  << " trying default GuessLinWCS" << endl;
  
  // the tel/inst specific procedure failed. try the default one ...
  
  if (HasLinWCS(Header)) return TanLinWCSFromHeader(Header,Guess);
  cout << " failed. now trying simple shift\n";
  return ComputeLinWCS(Header, Header.ImageCenter(), GtransfoIdentity(), Guess);  
}

static void TestWCSimplementation(const FitsHeader &Head) {
  TanPix2RaDec wcs;
  cout << " test GuessLinWCS " << endl;
  if (GuessLinWCS(Head,wcs))
    {
      cout << " considered successful " << endl;
      cout << wcs << endl;
    }
  else cout << " return false " << endl;
  cout << " test SkyRegion " << endl << SkyRegion(Head) << endl;
}

// tell fitstoad to use this function to test WCS
static int tata = SetTestWCS(TestWCSimplementation);

Frame SkyRegion(const FitsHeader &Header)
{
  int nx,ny;
  Header.ImageSizes(nx,ny);
  TanPix2RaDec pix2RaDec;
  if (!GuessLinWCS(Header, pix2RaDec)) return Frame();
  //  cout << " Lin WCS Guess " << Pix2RaDec << endl;
  return ApplyTransfo(Frame(Header), pix2RaDec, LargeFrame);
}


bool  CheckWCSIsAccurate( const FitsHeader & header )
{
  return ( (header.HasKey("WCSVERS")) || 
	   ("Swarp" == string(header.KeyVal("TOADINST"))) );
}
