#include <math.h>
#include <iostream>
#include <fstream>
#include <stdio.h>

#include "fitsimage.h"
#include "gtransfo.h"
#include "basestar.h"
#include "sestar.h"
#include "fileutils.h"
#include "vutils.h"
#include "frame.h"
#include "dbimage.h"
#include "starmatch.h"
#include "astroutils.h" /* For RaDec2000 */
#include "usnoutils.h"
#include "wcsutils.h"
#include "fitstoad.h"
#include "listmatch.h"

#ifndef M_PI
#define     M_PI            3.14159265358979323846  /* pi */
#endif

#ifdef TO_BE_DONE
- in UsnoProcess, there is only one "collection".
  it would be better to do two of them, one with a large tolerance
  (say 5"),  then a fit and a second one with a smaller (1.5") tolerance.
- read the various cuts (those tolerances) in datacards.
- read from a "virtual catalog"
#endif



//! finds the zeropoint for magnitudes
/*! from a matchlist between USNO and an image catalog, given a filter. 
  Assumes that stars in the starMatchList contain fluxes. */
int GetUsnoZeroPoint(const StarMatchList *MatchingList, UsnoColor Color, double &zeropoint, double &errzero)
{
  int n = 0;
  cout << "entering GetUsnoZeroPoint with " << MatchingList->size() << " matches" << endl;

  double *zp = new double[MatchingList->size()];
  // loop over the list 
  // sacre bricolage :
  bool use_sat = (getenv("NOSAT") == 0);
  for(StarMatchCIterator smi = MatchingList->begin(); smi != MatchingList->end() ; smi++)
    {
      const SEStar *catstar = smi->s1;
      const BaseStar *usnostar = (smi->s2);
      //get rid of saturated and negative flux
      if ( (catstar->flux > 0) &&  (usnostar->flux > 0)) 
	{
	  if (!use_sat || !catstar->IsSaturated())
	    {
	      zp[n] = -2.5*log10(usnostar->flux) + 2.5*log10(catstar->flux);
	      n++;
	    }
	}
    }
  // zeropoint is the clipped-mean at 3sigmas with 4 iterations.
  cout << " GetUsnoZeroPoint : " << n << ' ' << " matches" << endl;
  if (n<3) {delete [] zp; return 0;} 
  zeropoint = clipmean(zp,n,errzero,3,4);
  cout << "  GetUsnoZeroPoint : " << n << " matches survived clipped mean " << endl;
  delete [] zp;

  // Hack for B band : SN zeropoint at z=0.4 is 0.7mag higher    
  if (Color==BColor) zeropoint += 0.7;
  cout << "Zeropoint = " << zeropoint << " +- " << errzero << endl; 
 
  return 1;
}



#ifdef __linux__  /* byte swap */
#define LOWFIRST
#endif

static inline void intswap(unsigned int &a_word)
{
  char tmp;
  char *p = reinterpret_cast<char *>(&a_word);
  tmp = p[0]; p[0] = p[3]; p[3] = tmp; tmp = p[2]; p[2] = p[1]; p[1] = tmp;
}

static inline void read_a_star(FILE *ifp, unsigned int &raword, unsigned int  &decword, unsigned int &magword)
{ 
  fread(&raword,4,1,ifp);
  fread(&decword,4,1,ifp);
  fread(&magword,4,1,ifp);
#ifdef LOWFIRST
  intswap(raword);
  intswap(decword);
  intswap(magword);
#endif
}

#ifdef COMMENT

/* found in read.use in the catalog directory */

12) The RA takes a full 32-bit integer as does the SPD.  The third 32-bit
integer has been packed according to the following format.

SQFFFBBBRRR   (decimal), where

  S = sign is - if there is a correlated GSC entry, + if not.
    Q = 1 if internal PMM flags indicate that the magnitude(s) might be
      in error, or is 0 if things looked OK.  As discussed in read.pht,
	the PMM gets confused on bright stars.  If more than 40% of the
	pixels in the image were saturated, our experience is that the
	image fitting process has failed, and that the listed magnitude
	can be off by 3 magnitudes or more.  The Q flag is set if either
	  the blue or red image failed this test.  In general, this is a
	  problem for bright (<12th mag) stars only.

	  FFF = field on which this object was detected.  In the north, we
	    adopted the MLP numbers for POSS-I.  These start at 1 at the
	    north pole (1 and 2 are degenerate) and end at 937 in the -30
	    degree zone.  Note that fields 723 and 724 are degenerate, and we
	    measured but omitted 723 in favor of 724 which corresponds to the
	    print in the paper POSS-I atlas.  In the south, the fields start
	    at 1 at the south pole and the -35 zone ends at 408.  To avoid
	    wasting space, the field numbers were not put on a common system.
	    Instead, you should use the following test.

            IF ((zone.le.600).and.(field.le.408)) THEN
	    south(field)
	      ELSE
	    north(field)
	      ENDIF

	    BBB = 10 times the blue magnitude.  The range 0 through 250 contains
	    reasonable magnitudes.  500 is reserved for a PMM flux estimator
	    that was exactly zero, and 501 through 750 are reserved for PMM
	    flux estimators that were negative.  Only the reasonable magnitudes
	    were calibrated: the weird ones are just as they came out of the PMM.
	    For northern fields, magnitudes are defined by the 103a-O emulsion
	    and filter, while southern fields are defined by the IIIa-J emulsion
	    and filter.

	    RRR = 10 times the red magnitude.  As above except that northern plates
	    are 103a-E emulsions and southern plates are IIIa-F emulsions.

#endif /* comment */


/**********************************************************************/
/* Does the actual work of reading a file                             */
/**********************************************************************/

//#define DEBUG_READ /* to actually debug the reading of a file */
// routine for the official catalogue distribution
static void readusno(const string &filebase,double minra,double maxra,
	   double mindec,double maxdec, UsnoColor Color,
	   BaseStarList &ApmList, const Gtransfo& T= GtransfoIdentity())
{

  char line[80];
  int idx[96],len[96],i,raidx0,pos;
  float junk;
  FILE *ifp=NULL;
  unsigned int  raword,decword,magword;
  unsigned int minraword,maxraword,mindecword,maxdecword;
  double mag;

string catname = string(filebase) + ".cat";
string accname = string(filebase) + ".acc";

#ifdef DEBUG_READ
cerr << "readusno: reading USNO file  " << filebase << endl;
#endif

/* Read the accelerator */
  
  if (!(ifp=fopen(accname.c_str(),"r"))) {
    cerr << "readusno error: Couldn't open " << accname << endl;
    return;
    }
  for (i=0 ; i<96 && !feof(ifp) ; ++i) {
    fgets(line,80,ifp);
    sscanf(line,"%f %d %d",&junk,&(idx[i]),&(len[i]));
    }
  fclose(ifp);
  if (i!=96) {
    cerr << "parsusno error: " << accname << " too short, only " << i << "  lines " << endl;
    return;
    }


/* Open the cat file and skip to the first index we'll need */

  if (!(ifp=fopen(catname.c_str(),"r"))) {
    cerr << "readusno error: Couldn't open " << catname << endl;
    return;
    }

if (!(ifp=fopen(catname.c_str(),"r"))) {
cerr << "readusno error: Couldn't open " << catname << endl;
return;
}

raidx0=(int)floor(minra/3.75);       /* Indexed every 3.75deg = 15min */
pos=(idx[raidx0]-1)*12;
#ifdef DEBUG_READ
fprintf(stderr,"radix0=%d\n",raidx0);
fprintf(stderr,"Seeking to position %d in file %s\n",idx[raidx0],catname.c_str());
#endif
  if (fseek(ifp,pos,SEEK_SET) < 0) {
    cerr << "readusno: Error setting position in " << catname << endl;
    fclose(ifp);
    return ;
  }

/* Figure out the min and max ra and dec words */

  minraword=(unsigned int)floor(minra*3600*100 + 0.5);
  maxraword=(unsigned int)floor(maxra*3600*100 + 0.5);
  mindecword=(unsigned int)floor(mindec*3600*100 + 0.5);
  maxdecword=(unsigned int)floor(maxdec*3600*100 + 0.5);


  int mag_fact; /* mag = (magword/mag_fact)%1000, where mag_fact = 1 for R and 1000 for B  */
  if (Color == RColor) mag_fact = 1; 
  else if (Color == BColor) mag_fact = 1000;
  else {cerr << " unknown color value requested in UsnoListGet " << endl; fclose(ifp); return;}

#ifdef DEBUG_READ
fprintf(stderr,"minraword=%d, maxraword=%d\n",minraword,maxraword);
fprintf(stderr,"mindecword=%d, maxdecword=%d\n",mindecword,maxdecword);
#endif

/* Read until we find get into the RA range */

raword=0;

while (raword<minraword && !feof(ifp)) read_a_star(ifp, raword, decword, magword);

#ifdef DEBUG_READ
  fprintf(stderr,"Should be at minraword; minraword=%d, raword=%d\n",
  minraword,raword);
#endif

/* Read as long as we are in the right RA range; if a star read
   is in the right dec range, then keep it */

  while (raword<=maxraword  && !feof(ifp)) 
  {
    read_a_star(ifp, raword, decword, magword);
    if (decword>=mindecword && decword<=maxdecword && raword<=maxraword ) 
    {
      mag=(double)((magword/mag_fact)%1000)/10.;
      double ra = double(raword)/3600./100.;
      double dec = double(decword)/3600./100.-90;
      if (mag<99.9)
        {
	  //cout << raword << ' ' << decword << ' ' << mag << endl;
          BaseStar *s = new BaseStar(ra,dec,mag);
          *s = T.apply(*s);
	  ApmList.push_back(s);
        }        
    }
}
  

 fclose(ifp);
}


// read an ascii file containing alpha,delta mag, 1 item per line.
static void read_ascii_usno(const string &FileName, double MinRa, double MaxRa,
	   double MinDec, double MaxDec, 
	   BaseStarList &ApmList, const Gtransfo& T= GtransfoIdentity())
{
  FILE *file = fopen(FileName.c_str(),"r");
  if (!file)
    {
      cerr << " cannot open (supposedly ascii USNO) file " << endl;
      return;
    }
  char line[512];
  while (fgets(line,512,file))
    {
      if (line[0] == '#') continue;
      char cra[32], cdec[32];
      double mag;
      if (sscanf(line,"%s %s %lf",cra,cdec,&mag) != 3)
	{
	  cerr << " cannot decode ra dec mag in :" << endl << line << endl;
	  break;
	}
      double ra = RaStringToDeg(cra);
      double dec = DecStringToDeg(cdec);
      if (( ra < MinRa) || (ra > MaxRa) || ( dec < MinDec) || ( dec > MaxDec))
	continue;
      BaseStar *s = new BaseStar(ra,dec,mag);
      *s = T.apply(*s);
      ApmList.push_back(s);
    }
  fclose(file);
}

void ConvertMagToFlux(BaseStarList *List)
{
  for (BaseStarIterator si = List->begin(); si != List->end(); ++si)
    {
      BaseStar *s = (*si);
      if (s->flux < 40)  s->flux = pow(10., -(s->flux)/2.5);
      else s->flux = 0;
    }
}

int UsnoRead(const Frame &W, UsnoColor Color, BaseStarList &ApmList)
{
  return UsnoRead(W.xMin, W.xMax, W.yMin, W.yMax, Color, ApmList);
}


//! ra's and dec's in degrees. Handles limits across alpha=0. (not with USNOFILE)
int UsnoRead(double minra, double maxra, 
             double mindec, double maxdec, UsnoColor Color, BaseStarList &ApmList)
{
  char *env_var = getenv("USNODIR");
  if (!env_var)
    {
      env_var = getenv("USNOFILE");
      if (env_var)
	{
	  read_ascii_usno(env_var,minra,maxra,mindec,maxdec, ApmList);
	  return ApmList.size();
	}
      else
	{
	  cerr << " you should define USNODIR or USNOFILE env var " 
	       << "to run this code " << endl;
	  return 0;
	}
    }
  string usnodir(env_var);
  int cat0,cat1;
  
  cout << " reading usno in window (" << minra << ',' << maxra << ") (" << mindec << ',' << maxdec << ")" << endl;

  /* Read parameters */
  
#ifdef DEBUG_READ
  fprintf(stderr,"Going to read USNO from directory %s\n",usnodir.c_str());
  fprintf(stderr,"minra=%lf, maxra=%lf\n",minra,maxra);
  fprintf(stderr,"mindec=%lf, maxdec=%lf\n",mindec,maxdec);
#endif

  /* Make sure the region is reasonable */

  if (minra == maxra || mindec == maxdec) return 0;
  if (minra > maxra) swap (minra,maxra);
  if (mindec > maxdec) swap(mindec,maxdec);
  if (maxra>=360.) return 0;
  if (mindec<-90. || maxdec>90.) return 0;

  /* Figure out which catalog files we need to read, and read them */

  mindec+=90.;
  maxdec+=90.;           /* Dec->spd */

  cat0=int(75*floor(mindec/7.5));
  cat1=int(75*floor(maxdec/7.5));
#ifdef DEBUG_READ
  fprintf(stderr,"cat0=%d,cat1=%d\n",cat0,cat1);
#endif
  for ( ; cat0<=cat1; cat0+=75) 
    {
      char cat_char[12];
      sprintf(cat_char,"zone%04d",cat0);
#ifdef DEBUG_READ
    fprintf(stderr,"Calling readusno on %s\n",cat_char);
#endif
    // some precaution if minra<0 and maxra >0 (have to read in 2 steps)
    string filename = usnodir+'/'+cat_char;
    if (minra<0 and maxra >0)
      {
	readusno(filename, minra+360.0,360.0,mindec,maxdec,Color, ApmList, 
		 GtransfoLinShift(-360, 0));
	readusno(filename, 0,maxra,mindec,maxdec,Color, ApmList);
      }
    else if (minra<0 and maxra <0)
      {
	readusno(filename, minra+360.0,maxra+360.0,mindec,maxdec,Color, ApmList, 
		 GtransfoLinShift(-360, 0));
      }
    else if (maxra>360. and minra<360.)
      {
	readusno(filename, minra, 360,mindec,maxdec,Color, ApmList);
	readusno(filename, 0,maxra-360,mindec,maxdec,Color, ApmList, 
		 GtransfoLinShift(360,0));
      }
    else if (maxra>360. and minra>360.)
      {
	readusno(filename, minra-360.,maxra-360,mindec,maxdec,Color, ApmList, 
		 GtransfoLinShift(360,0));
      }
    else readusno(filename, minra,maxra,mindec,maxdec,Color, ApmList);
    }

#ifdef DEBUG_READ
  fprintf(stderr,"Returning.\n");
#endif
  for (BaseStarIterator si = ApmList.begin(); si != ApmList.end(); ++si)
    {
     BaseStar &s = *(*si);
     if (s.x <minra || s.x > maxra) cout << " problem with usno star  : " << s ;
   }
  int count = ApmList.size();
  cout << " collected " << count  << " objects" << endl;  
  return count;
}



//! Guesses the transfo from sestarlist to catalog.  
/*! UsnoToPix is the transformation that was applied to 
  the USNO catalog to project it on the tangent plane, 
  with coordinates expressed in degrees.
    Hence, UsnoCat constains coordinates in degrees in the
  tangent plane, and fluxes (in dummy unit) instead of magnitudes.
  The  tangent point can be extracted from UsnoToPix.
  The found matches and the guessed transfo derived from 
  them are in the
  returned StarMatchList, to be deleted by the caller.*/
StarMatchList *FindMatchUsno(const string &FitsImageName, 
			       SEStarList &ImageList, 
			       Point &TangentPoint, 
			       BaseStarList &UsnoCat)
{
  
  FitsHeader header(FitsImageName);
  if (!header.IsValid()) return NULL;

  Frame usnoFrame = SkyRegion(header);
  if (usnoFrame.Area()==0) // something went obviously wrong
    {
      cerr << " cannot find where in the sky the image " <<  FitsImageName << " corresponds" << endl;
      return NULL;
    }
  usnoFrame = usnoFrame.Rescale(1.2); // add 20%  
  UsnoRead(usnoFrame, RColor, UsnoCat);

  ConvertMagToFlux( &UsnoCat);
  cout << " trying to match " << FitsImageName << " with usno R catalog " << endl;

  TanPix2RaDec Pix2RaDec;
  if (!GuessLinWCS(header, Pix2RaDec)) return NULL;  
  TanRaDec2Pix UsnoToPix = Pix2RaDec.invert();

  // convert the ( ra, dec) to pixel scale
  UsnoCat.ApplyTransfo(UsnoToPix);


  bool okShift = false;
  StarMatchList *match = NULL;
  if (!getenv("NO_SHIFT"))
    {
  // try to guess a simple shift, in case the guessed WCS is good:
      cout << " trying to guess a shift" << endl;
      GtransfoLin *shift = ListMatchupShift(*SE2Base(&ImageList), UsnoCat, 
					    GtransfoIdentity(), 200.,10.);
      match = ListMatchCollect(*SE2Base(&ImageList),
					      UsnoCat, shift, 10.);
      delete shift;
      unsigned shiftMatches = match->size();
      if (3*shiftMatches > UsnoCat.size() && shiftMatches > 25 ) // considered as good
	{
	  cout << " found " << match->size() << " matches with a simple shift" 
	   << endl;
	  okShift =true;
	}
      else
	{
	  cout << " only " << match->size() << " matches found "<< endl;
	  delete match; match = NULL;
	}
    }
  if (!okShift) // try the hard way
    {

      MatchConditions conditions;
      // areas : 
      double usnoWindowSizeInPix = usnoFrame.ApplyTransfo(UsnoToPix).Area();
      double imageSize = Frame(header).Area();
      
      conditions.NStarsL1 = 60;
      conditions.NStarsL2 = 60;
      conditions.NStarsL2 = int (conditions.NStarsL1 * (usnoWindowSizeInPix/imageSize));
      
      cout << " combinatorial search with lists of " << conditions.NStarsL1 << ' ' <<  conditions.NStarsL2 << " stars " << endl;
      
      match = MatchSearchRotShiftFlip(*SE2Base(&ImageList), UsnoCat, conditions);
      if (!match) return NULL;
    }

  // Here we check that our initial guess for the usnoWindow is accurate (i.e 
  // we took in the USNO catalog all the area covered by the image.
  GtransfoLin guessCorr = *dynamic_cast<const GtransfoLin *>(match->Transfo());
  cout << " correction to guessed WCS, found using " << match->size() 
       << " pairs" << endl<< guessCorr << endl;
  TanPix2RaDec newPix2RaDec = Pix2RaDec*guessCorr;
  cout << " updated lin WCS: " << endl << newPix2RaDec << endl;

  /* now we know where the image is on the sky. Check that 
     we collected catalog actually covers the image...
  */
  
  int nx,ny;
  header.ImageSizes(nx, ny);
  Frame newUsnoFrame(newPix2RaDec.apply(Point(0,0)), 
		     newPix2RaDec.apply(Point(nx,ny)));
  // intersection between rough window and more accurate one:
  Frame overlap = newUsnoFrame*usnoFrame; 
  // allow for 10 % loss
  if (match->Nused() > 4 && overlap.Area() < newUsnoFrame.Area()*0.9) 
    {
      // recollect USNO
      cout << " FindTransfoUsno: bad initial guess of window on USNO, recollecting" << endl;
      cout << " shifting by " 
	   << newUsnoFrame.Center()-usnoFrame.Center() << endl;
      newUsnoFrame = newUsnoFrame.Rescale(1.2); // add 20%
      UsnoCat.clear(); 
      UsnoRead(newUsnoFrame, RColor, UsnoCat);
      ConvertMagToFlux( &UsnoCat);
      UsnoToPix = newPix2RaDec.invert();

      // convert the ( ra, dec) to pixel units
      UsnoCat.ApplyTransfo(UsnoToPix); 
      // recollect , but no need for a combinatorial match:
      delete match;
      GtransfoIdentity id;
      match = ListMatchCollect(*SE2Base(&ImageList), UsnoCat, &id, 5);
      match->RefineTransfo(3.);
    }

  // express outputs in the tangent plane , in degrees
  // transfo from guEssed pixels to degrees
  GtransfoLin pixToDegrees(UsnoToPix.LinPart().invert());
  UsnoCat.ApplyTransfo(pixToDegrees);
  TangentPoint = Pix2RaDec.TangentPoint();
  // rebuild the match, since the UsnoList has changed (scale)
  delete match;
  GtransfoLin cdFound(newPix2RaDec.LinPart());
  match = ListMatchCollect(*SE2Base(&ImageList), UsnoCat, 
			   &cdFound, 2./3600.);
  match->RefineTransfo(5);
  return match;
}


static double object_size(const SEStar *object)
{
  return sqrt(object->Mxx()*object->Myy());
}

static void match_clean(StarMatchList &List)
{
  int count = List.size();
  double *sizes = new double[count];
  int i=0;
  for (StarMatchIterator si=List.begin(); si!= List.end(); ++si)
    {
      SEStar *object = (SEStar *) (BaseStar *) si->s1;
      sizes[i] = object_size(object);
      i++;
    }
  double cut = 2.0 * DArrayMedian(sizes,count);
  cout << "cut " << cut << endl;
  delete [] sizes;
  for (StarMatchIterator si=List.begin(); si!= List.end(); )
    {
      SEStar *object = (SEStar *) (BaseStar *) si->s1;
      if (object_size(object)>cut)
	{
	  si = List.erase(si);
	}
      else ++si;
    }
}

      
//! fills a file to be used when building finding charts, or to check the match.
void FillMatchFile(const DbImage &Image, const Gtransfo &Pix2RaDec)
{
  FillMatchFile(FitsHeader(Image.FitsImageName(Calibrated)), 
		SEStarList(Image.ImageCatalogName()),
		Pix2RaDec, Image.Dir()+"/match_usno.dat");
}

//! fills a file to be used when building finding charts, or to check the match.
void FillMatchFile(const FitsHeader &header, const SEStarList &imageList, const Gtransfo &Pix2RaDec, const string &MatchFileName)
{

  Frame pixFrame(header); // boundaries of the image
  Frame raDecFrame = pixFrame.ApplyTransfo(Pix2RaDec);
  Frame bigraDecFrame = raDecFrame.Rescale(1.1);
  BaseStarList catalog;
  UsnoRead( bigraDecFrame, RColor, catalog);

  // take a copy of the input list
  // match to the catalog
  StarMatchList *matches = ListMatchCollect(*SE2Base(&imageList), catalog, &Pix2RaDec, 2./3600.); /* not farther than 2", becomes wrong if abs(Dec) is too large.*/

  cout << " writing " << MatchFileName << endl;
  FILE  *matchstream = fopen(MatchFileName.c_str(), "w"); // use C IO's because formats are easy !

  // header : descriptor used by l2tup
  fprintf(matchstream,"# X_im : x position of the stars from the image after fitting the flux\n");
  fprintf(matchstream,"# Y_im : same for y \n");
  fprintf(matchstream,"# ra_im : same for ra \n");
  fprintf(matchstream,"# dec_im : same for dec \n");
  fprintf(matchstream,"# flux_im : flux on image \n");
  fprintf(matchstream,"# ra : ra from catalog (degrees)\n");
  fprintf(matchstream,"# dec : dec from catalog\n");
  fprintf(matchstream,"# MAG_cat : magnitude from the usno catalog\n");
  fprintf(matchstream,"# apflux : aperture flux\n");
  fprintf(matchstream,"# distsec : distance beween catalog and image in arcsec\n");
  fprintf(matchstream,"# end\n");

  for (StarMatchCIterator it = matches->begin(); it != matches->end(); ++it)
    {
      const SEStar *s1 = (const SEStar *) (const BaseStar *) it->s1; // since *s1 belongs to imageList ...
      const BaseStar *s2 = it->s2;
      Point ra_dec_im = Pix2RaDec(*s1);
      double distsec = s2->Distance(ra_dec_im)*3600.;

      fprintf(matchstream,"%8.2f %8.2f %12.6f %12.6f %12.1f %12.6f %12.6f %6.2f %12.1f %8.3f\n",
	      //             x     y     ra     dec   flux   ra      dec   mag  apflux
	      s1->x , s1->y , ra_dec_im.x, ra_dec_im.y, s1->flux, s2->x , s2->y, s2->flux, s1->Flux_aper(), distsec);
    }

  for (SEStarCIterator it = imageList.begin(); it != imageList.end();++it  )
    {
      const SEStar &star = **it;
      Point ra_dec_im = Pix2RaDec(star);
      fprintf(matchstream,"%8.2f %8.2f %12.6f %12.6f %12.1f %12.6f %12.6f %6.2f %12.1f %8.3f\n",
	      //             x     y     ra     dec   flux   ra      dec   mag  apflux
	      star.x, star.y, ra_dec_im.x, ra_dec_im.y, star.flux,ra_dec_im.x, ra_dec_im.y, -10., star.flux, 0.);
    } 

  fclose(matchstream);
  delete matches;
}

static void FillAllStarsFile(const DbImage &dbImage, const Gtransfo &Pix2RaDec)
{
  SEStarList imageList(dbImage.ImageCatalogName());
  string fileName = dbImage.Dir()+"/"+"all_coords.dat";
  cout << " writing " << fileName << endl;
  FILE  *matchstream = fopen(fileName.c_str(), "w"); // use C IO's because formats are easy !
  // header : descriptor used by l2tup
  fprintf(matchstream,"# X_im : x position of the stars from the image after fitting the flux\n");
  fprintf(matchstream,"# Y_im : same for y \n");
  fprintf(matchstream,"# ra_im : same for ra \n");
  fprintf(matchstream,"# dec_im : same for dec \n");
  fprintf(matchstream,"# flux_im : flux on image \n");
  fprintf(matchstream,"# ra : ra from catalog (degrees)\n");
  fprintf(matchstream,"# dec : dec from catalog\n");
  fprintf(matchstream,"# MAG_cat : magnitude from the usno catalog\n");
  fprintf(matchstream,"# apflux : aperture flux\n");
  fprintf(matchstream,"# distsec : distance beween catalog and image in arcsec\n");
  fprintf(matchstream,"# end\n");


  for (SEStarCIterator it = imageList.begin(); it != imageList.end();++it  )
    {
      const SEStar &star = **it;
      Point ra_dec_im = Pix2RaDec(star);
      fprintf(matchstream,"%8.2f %8.2f %12.6f %12.6f %12.1f %12.6f %12.6f %6.2f %12.1f %8.3f\n",
	      //             x     y     ra     dec   flux   ra      dec   mag  apflux
	      star.x, star.y, ra_dec_im.x, ra_dec_im.y, star.flux,ra_dec_im.x, ra_dec_im.y, -10., star.flux, 0.);
    } 

  fclose(matchstream);
}

/********************************************************************/
/*************         USNO PROCESS STUFF                *************/
/********************************************************************/

// To match an image with the usno catalogue
bool UsnoProcess(DbImage &dbimage, const bool write)
{
  if (!dbimage.IsValid())
    {
      cerr << " Be careful ! " << dbimage.Name() << " must be an image name ! " << endl;
    }
  
  string catalogName = dbimage.ImageCatalogName(SExtractor);
  if (!FileExists(catalogName.c_str()))
    {
      cerr << "The SExtractor Catalogue associated to image " << dbimage.Name() << " doesn't exist !! " << endl;
      return false;
    }  
  string fitsFileName = dbimage.FitsImageName(Calibrated);
  if (!FileExists(fitsFileName))
    {
      cerr << "The  calibrated fits image " << dbimage.Name() << " doesn't exist !! " << endl;
      return false;
    }

  return UsnoProcess(fitsFileName, catalogName, &dbimage, write);
}


bool UsnoProcess(const string &fitsFileName, const string &catalogName, 
		 DbImage *dbimage, const bool write)
{  

  FitsHeader header(fitsFileName,(write)? RW : RO);

  cout << endl << "Matching with USNO catalog..."<< endl;
  cout << fitsFileName << endl;
  cout << "INSTRUMENT = " << header.KeyVal("TOADINST") << endl;
  cout << "CCD no = " << header.KeyVal("TOADCHIP") << endl;
  cout << "FILTER = " << header.KeyVal("TOADFILT") << endl;
  cout << "OBJECT = " << header.KeyVal("TOADOBJE") << endl;
  cout << "REMEMBER: matching with USNO + Zeropoint made in R-USNO !" << endl;
  cout << endl;
   
  SEStarList sestarlist(catalogName);

  if (sestarlist.size()==0) {cerr << " UsnoProcess: SExtractor catalog empty " << endl; return false;}

  cout << " Read " << sestarlist.size() << " objects from SexCat " << endl;   


  BaseStarList usnoList;
  Point tangentPoint;
  StarMatchList *guessedMatch = 
    FindMatchUsno(fitsFileName, sestarlist, tangentPoint, usnoList);

  if (usnoList.size()==0) 
    {cerr << " UsnoProcess: USNO list empty: give up " << endl;return false;}

  // The list we get here is in units of pixels using the pixel scale from header.
  // the transformation that was used to transform the catalogList is usnoToPix


  // Quand ca ne match pas, j'ai toujours vu un nombre de paires <= 4, so...
  if (guessedMatch->Nused() < 5)
    {
      cerr << " MATCH FAILED... " << endl;
      cerr << "Guessed Transfo = " << endl;
      cerr << *(guessedMatch->Transfo()) << endl;
      cerr << "Chi2 " << guessedMatch->Chi2() << " Nused " << guessedMatch->Nused() << " Residue " 
	   << guessedMatch->Residual() << endl << endl;

      cerr << endl << " ****************************** " << endl;
      return false;
    }

  const Gtransfo *guess = guessedMatch->Transfo();


  // Recollect a better match with a controlled accuracy
  double toldiff = 1.5;
  StarMatchList *accurateMatch = ListMatchCollect(*SE2Base(&sestarlist), 
						  usnoList, 
						  guess, 
						  toldiff/3600.);
  int nMatches = accurateMatch->size();
  delete guessedMatch;

  cout << " collected " << nMatches << " matches with tolerance " 
       << toldiff << " arcsec " << endl;

  if (nMatches < 0.2*usnoList.size())
    {
      cout << " refined failed, giving up " << endl;
      return false;
    }

  int match_order = 3;
  if (getenv("MATCH_ORDER")) match_order = atoi(getenv("MATCH_ORDER"));
  // check that the problem is enough constrained
  while (((match_order+1)*(match_order+2))/2 > nMatches)
    {
      match_order--;
      cout << " reducing match_order to " << match_order << endl;
    }
  accurateMatch->SetTransfoOrder(match_order);
  match_clean(*accurateMatch);
  cout << " nmatches after clean " << accurateMatch->size() << endl;
  // return if not enough stars to produce a linear transfo (oderwize: seg fault soon)
  if (accurateMatch->size() < 5)
    { cerr << " not enough matches. giving up " << endl;return false;}

  // argument is number of sigmas for outlier rejection.
  accurateMatch->RefineTransfo(3);
  
  // some printouts
  cout << "Match number after refine = " << accurateMatch->size() << endl;
  cout << " Nused "  << accurateMatch->Nused() << " Residual " 
       << accurateMatch->Residual()*3600. << " arcsec" << endl;

  //DEBUG
  cout << *accurateMatch->Transfo();
  cout << " 2000 2000 " << accurateMatch->Transfo()->apply(Point(2000,2000))<< endl;

  // dump a match file?
  char *matchFileName;
  if ((matchFileName = getenv("DUMP_MATCH")))
    {
      accurateMatch->write(*accurateMatch->Transfo(),matchFileName);
      BaseStarList imageReducedList;
      for (StarMatchCIterator i=accurateMatch->begin(); i != accurateMatch->end(); ++i)
	{
	  imageReducedList.push_back((i->s1));
	}
      imageReducedList.write((string(matchFileName)+".reduced").c_str());
    }


  /* We split the transfo into 2 parts : a linear one for use by 
     official tools (e.g. saoimage), and a cubic correction 
     to the it to be used by our stuff (for astrometry).
     Something that makes the code messy is that the catalog 
     we use is approximately in pixel scale: (it was converted 
     from alpha,delta to what we get using usnoToPix.
     So to go from image to (alpha,delta) the transfo reads 
     usnoToPix.invert() * FoundTransfo, if this makes sense.
     The choosen implementation changed around may 03, to follow
     the anticipated normalization of optical distortions 
     handling. 

  */

  // get the best linear description of our match
  GtransfoLin linFit;
  linFit.fit(*accurateMatch);

  GtransfoCub cubCorr;
  GtransfoQuad quadCorr;
  GtransfoQuad *correction = NULL;
  if (match_order>=2)
    {
      if (match_order == 3) correction = &cubCorr;
      if (match_order == 2) correction = &quadCorr;
      if (correction) correction->fit(*accurateMatch, &linFit, NULL);
      else
	{
	  cerr << " ERROR in " << __FILE__ << endl;
	  cerr << " match_order is wrong " << match_order 
	       << " no corrections for non-linear distotions " << endl;
	}
    }

  // write the WCS in the fits header :
  if (write)
    {
      TanPix2RaDec pix2RaDec(linFit, tangentPoint);
      if (correction) pix2RaDec.SetCorrections(correction);
      TanWCS2Header(fitsFileName, pix2RaDec);
      if (dbimage) 
	{
	  FillMatchFile(*dbimage, pix2RaDec);
	  FillAllStarsFile(*dbimage, pix2RaDec);
	}
      else // write the match_usno.dat file
	{
	  FillMatchFile(FitsHeader(fitsFileName), sestarlist, pix2RaDec,
			DirName(fitsFileName)+"/match_usno1.dat");
	} 
    }

  // photometric zero point   
  double zeropoint, errzero;
  if (GetUsnoZeroPoint(accurateMatch, RColor, zeropoint, errzero))
    {

      cout << "Residue of the match with usno catalog " << accurateMatch->Residual() << endl
	   << " Zero point " << zeropoint << endl
	   << "Dispersion on Zero Point (in mag)" << errzero << endl;
      if (write)
	{
	  header.AddOrModKey("RMATCHUS", accurateMatch->Residual(), "Residue of the match with usno catalog");
	  header.AddOrModKey("NMATCHUS", int(accurateMatch->size()), " number of objects matched to USNO catalog");
	  header.AddOrModKey("ZEROUSNO", zeropoint, "Zero point : Mag(usno)=zero+Mag(imag)");
	  header.AddOrModKey("DZEROUSN", errzero, "Dispersion on Zero Point (in mag)");
	}

      double skysig = header.KeyVal("SEXSIGMA");
      double seeing = header.KeyVal("SESEEING");
      double flux21 = pow(10., -(21.-zeropoint)/2.5);
      double flux23 = pow(10., -(23.-zeropoint)/2.5);
      double flux25 = pow(10., -(25.-zeropoint)/2.5);
      double aine = M_PI * (2.36 * 2. * seeing) * (2.36 * 2. * seeing);
      double ssn21 = flux21 / sqrt(flux21 + aine * skysig * skysig);
      double ssn23 = flux23 / sqrt(flux23 + aine * skysig * skysig);
      double ssn25 = flux25 / sqrt(flux25 + aine * skysig * skysig);


      cout << "Signal/Noise (at mag 21) = " << ssn21 << " (at mag 25) = " << ssn25 << endl << endl;
      if (write)
	{
	  header.AddOrModKey("USNOSB21", ssn21, "Signal / Noise at mag = 21");
	  header.AddOrModKey("USNOSB23", ssn23, "Signal / Noise at mag = 23");
	  header.AddOrModKey("USNOSB25", ssn25, "Signal / Noise at mag = 25");
	}
    }
  return true;
}




