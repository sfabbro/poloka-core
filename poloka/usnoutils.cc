#define USNOUTILS__CC
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdio>

#include <poloka/fitsimage.h>
#include <poloka/gtransfo.h>
#include <poloka/sestar.h>
#include <poloka/fileutils.h>
#include <poloka/vutils.h>
#include <poloka/frame.h>
#include <poloka/dbimage.h>
#include <poloka/starmatch.h>
#include <poloka/astroutils.h>
#include <poloka/usnoutils.h>
#include <poloka/wcsutils.h>
#include <poloka/listmatch.h>
#include <poloka/fastfinder.h>
#include <poloka/imageutils.h>
#include <poloka/matchexception.h>
#include <poloka/datacards.h>

#ifndef M_PI
#define     M_PI            3.14159265358979323846  /* pi */
#endif


/* TO DO
- Enforce UsnoCollect and UsnoRead to provide some projection transfo 
so that we would apply the transfo and define errors in the tangent plane 
(ie. and orthonormal coordinate system).

- Detect if the match is ok only on a fraction of the CCD (this
  happened on the grid data. To do this, compute the rms of the x and
  y coordinates of the matches and compare them with the CCD size /
  sqrt(12)...
*/

static double sq(const double &x) { return x*x;}


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





/****************************************************************
  Does the actual work of reading a USNO file. The format of the USNO
  files (picked up in the USNO distribution) is documented at the end
  of this source file
*****************************************************************/

//#define DEBUG_READ /* to actually debug the reading of a file */

static void readusno(const string &filebase,double minra,double maxra,
		     double mindec,double maxdec, UsnoColor Color,
		     BaseStarList &ApmList, 
		     const GtransfoLinShift& T= GtransfoLinShift(0,0))
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
	  /* We have to set errors here, because there is no single place 
	     below where it could easily be done. I wish we had a projection 
	     transfo in hand to define errors in the tangent plane.
	     For the value, we chose 0.3" r.m.s
	  */
	  s->vx = sq(0.3/3600/cos(dec*M_PI/180.));
	  s->vy = sq(0.3/3600);
	  s->vxy = 0;
	  if (s->x <minra || s->x > maxra) 
	    cout << " problem with usno star  : " << s ;
          *s = T.apply(*s); // no need to transfor errors : the transfo is a shift
	  ApmList.push_back(s);
        }        
    }
}
  
 fclose(ifp);
}


// read an ascii file containing alpha,delta mag, 1 item per line.
/* \page usno_file_format Astrometric file format
    You can provide your own match catalog to matchusno, either on the 
    command line, or via the USNOFILE environment variable. The code
    expects and ascii file and interprets the 3 first items of
    every line as alpha, delta (equatorial) and mag (where the mag 
    is mainly used to 
    isolate the brightest objects). Sexagesimal coordinates are accepted.
*/
  
static void read_ascii_astrom_file(const string &FileName, 
				  double MinRa,  double MaxRa,
				  double MinDec, double MaxDec, 
				  BaseStarList &ApmList, 
				  const GtransfoLinShift& T= GtransfoLinShift(0,0))
{
  FILE *file = fopen(FileName.c_str(),"r");
  if (!file)
    {
      std::cerr << " cannot open (supposedly ascii USNO file ) " 
		<< FileName << std::endl;
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
  std::cout << " collected " << ApmList.size() 
	    << " objects from " << FileName << std::endl;
  fclose(file);
}

void ConvertMagToFlux(BaseStarList *List, const double Zp)
{
  
  for (BaseStarIterator si = List->begin(); si != List->end(); ++si)
    {
      BaseStar *s = (*si);
      if (s->flux < 40)  s->flux = pow(10., -(s->flux-Zp)*0.4);
      else s->flux = 0;
    }
}

static void actual_usno_read(const string &usnodir,
			     double minra, double maxra, 
			     double mindec, double maxdec, UsnoColor Color,
			     BaseStarList &ApmList)
{
  int cat0,cat1;
  
  cout << " reading usno in window (" << minra << ',' << maxra << ") (" << mindec << ',' << maxdec << ")" << endl;

  /* Read parameters */
  
#ifdef DEBUG_READ
  fprintf(stderr,"Going to read USNO from directory %s\n",usnodir.c_str());
  fprintf(stderr,"minra=%lf, maxra=%lf\n",minra,maxra);
  fprintf(stderr,"mindec=%lf, maxdec=%lf\n",mindec,maxdec);
#endif

  /* Make sure the region is reasonable */

  if (minra == maxra || mindec == maxdec) return;
  if (minra > maxra) swap (minra,maxra);
  if (mindec > maxdec) swap(mindec,maxdec);
  // the following checks are questionnable.... update them (P.A, 01/07)
  //  if (maxra>=360.) return;
  //  if (mindec<-90. || maxdec>90.) return;
  if (mindec < -90.) mindec = -90;
  if (maxdec > 90.) maxdec = 90;

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
  int count = ApmList.size();
  cout << " collected " << count  << " objects" << endl;
}


int UsnoRead(const Frame &W, UsnoColor Color, BaseStarList &ApmList)
{
  return UsnoRead(W.xMin, W.xMax, W.yMin, W.yMax, Color, ApmList);
}



//! ra's and dec's in degrees. Handles limits across alpha=0. (not with USNOFILE)
int UsnoRead(double minra, double maxra, 
	     double mindec, double maxdec, UsnoColor Color, 
	     BaseStarList &ApmList)
{
  const char *ascii_source = NULL;

  if (MatchPrefs.astromCatalogName != "" ) 
    ascii_source = MatchPrefs.astromCatalogName.c_str();
  else
    {
      char *env_var = getenv("USNOFILE");
      if (env_var) ascii_source = env_var;
    }
  if (ascii_source)
    read_ascii_astrom_file(ascii_source, minra,maxra,mindec,maxdec, 
			   ApmList);
  else
    {
      char *usno_dir = getenv("USNODIR");
      if (usno_dir)
	{
	  actual_usno_read(usno_dir, 
			   minra, maxra, mindec, maxdec, Color, 
			   ApmList);
	}
      else
	{
	  std::cerr << " ERROR : You should define USNODIR or USNOFILE env var, " 
		    << std::endl
		    << " or provide ASTROM_CATALOG_NAME via datacards to run this code " << std::endl;
	  return 0;
	}
    }
  return ApmList.size();
}

static bool check_guess(const GtransfoLin &Guess)
{
  //check that the matrix corresponds roughly to a rotation or a flip
  double mod12 = sq(Guess.A11())+ sq(Guess.A12());
  double mod22 = sq(Guess.A21())+ sq(Guess.A22());
  // vectors of same size...
  //  cout << (mod12-mod22)/(mod12+mod22) << endl;
  if (fabs((mod12-mod22)/(mod12+mod22)) > 0.2) return false;
  double cos_angle = (Guess.A11()*Guess.A21()+ Guess.A12()*Guess.A22())
    /sqrt(mod12*mod22);
  // ... and  roughly orthogonal
  if (fabs(cos_angle)> 0.2) return false;
  return true;
}
  

static GtransfoLin *FindShift(const BaseStarList &ImageList, 
			      const BaseStarList &UsnoCat, 
			      const unsigned &MinMatches)
{
  GtransfoLin *shift = NULL;
  // try to guess a simple shift, in case the guessed WCS is good:
  const double dsmax = 10;
  cout << " trying to guess a shift (tolerance " << dsmax << ")\n";
  shift = ListMatchupShift(ImageList, UsnoCat, 
			   GtransfoIdentity(), 200., dsmax);
  StarMatchList *match = ListMatchCollect(ImageList,
					  UsnoCat, shift, dsmax);
  unsigned shiftMatches = match->size();
  delete match;
  if (shiftMatches >= MinMatches)
    {
      // considered as good :
      cout << " found " << shiftMatches << " matches with a simple shift" 
	   << endl;
      return shift;
    }
  else
    {
      cout << " only " << shiftMatches << " matches found "<< endl;
      return NULL;
    }
  return shift;
}


//! Guesses the transfo from sestarlist to catalog.  
/*! TangentPoint denote the tangent point used to project the USNO
  catalog on a tangent plane, with coordinates expressed in degrees.
  Hence, UsnoCat constains coordinates in degrees in the tangent
  plane, and incidentally fluxes (in dummy unit) instead of
  magnitudes.  The tangent point can be extracted from UsnoToPix.  The
  found matches and the guessed transfo derived from them are in the
  returned StarMatchList, to be deleted by the caller.*/


static GtransfoLin *FindCombiMatch(BaseStarList &ImageList, BaseStarList &UsnoCat,const MatchConditions &conditions)
{

  cout << " combinatorial search with lists of " << conditions.NStarsL1 << ' ' <<  conditions.NStarsL2 << " stars " << endl;
  
  StarMatchList *match = MatchSearchRotShiftFlip(ImageList, UsnoCat, conditions);
  if (!match) return NULL;

  GtransfoLin guessCorr = *dynamic_cast<const GtransfoLin *>(match->Transfo());
  cout << " correction to guessed WCS, found using " << match->size() 
       << " pairs" << endl<< guessCorr << endl;
  if (!check_guess(guessCorr))
    {
      cout << " the guess is really too far from a rotation " << endl
	   << " giving up" << endl;
      delete match;
      throw(MatchException(" Found Guess is really too far from a rotation "));
      return NULL;
    }
  return new GtransfoLin(guessCorr);

}


static void match_clean(StarMatchList &List)
{

  cout << " count before satur removal " << List.size() << endl;
  for (StarMatchIterator si=List.begin(); si!= List.end();)
    {
      SEStar *object = (SEStar *) (BaseStar *) si->s1;
      if (object->IsSaturated() || object->FlagBad()) si = List.erase(si);
      else si++;
    }
  cout << " count after satur & bad removal " << List.size() << endl;
}


      

//! fills a file to be used when building finding charts, or to check the match.
static void FillMatchFile(const FitsHeader &Header, const SEStarList &ImageList, 
			  const StarMatchList &MatchList, const Gtransfo &Pix2RaDec, const string &MatchFileName)
{

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

  BaseStarList usedUsno;

  for (StarMatchCIterator it = MatchList.begin(); it != MatchList.end(); ++it)
    {
      const SEStar *s1 = (const SEStar *) (const BaseStar *) it->s1; // since *s1 belongs to imageList ...
      const BaseStar *s2 = it->s2;
      usedUsno.push_back(s2);
      double cosdec = cos(s2->y*M_PI/180.);
      Point ra_dec_im = Pix2RaDec.apply(*s1);
      double distsec = sqrt(sq((ra_dec_im.x-s2->x)*cosdec)+sq(ra_dec_im.y-s2->y))*3600.;

      fprintf(matchstream,"%8.2f %8.2f %12.6f %12.6f %12.1f %12.6f %12.6f %6.2f %12.1f %8.3f\n",
	      //             x     y     ra     dec   flux   ra      dec   mag  apflux
	      s1->x , s1->y , ra_dec_im.x, ra_dec_im.y, s1->flux, s2->x , s2->y, s2->flux, s1->Flux_auto(), distsec);
    }

  //also write distances to unmatched USNO stars

  // collect the whole catalog
  
  Frame pixFrame(Header); // boundaries of the image
  Frame raDecFrame = ApplyTransfo(pixFrame,Pix2RaDec);
  Frame bigraDecFrame = raDecFrame.Rescale(1.1);
  BaseStarList unusedUsno;
  UsnoRead( bigraDecFrame, RColor, unusedUsno);

  FastFinder finder(*SE2Base(&ImageList));
  Gtransfo *raDec2Pix = Pix2RaDec.InverseTransfo(0.01,pixFrame);


  // ignore matched stuff (collected 15 line above)
  for (BaseStarIterator i=unusedUsno.begin(); i!= unusedUsno.end(); ++i)
    {
      BaseStar *s2 = *i;
      BaseStar *b = usedUsno.FindClosest(*s2);
      if (b->Distance(*s2)<0.0001) continue;
      Point imageCoords = raDec2Pix->apply(*s2);
      if (!pixFrame.InFrame(imageCoords)) continue;
      double maxDist = 20; // pixels
      const SEStar *s1 = (const SEStar *) finder.FindClosest(imageCoords,maxDist);
      if (s1)
	{
	  Point ra_dec_im = Pix2RaDec.apply(*s1);
	  double cosdec = cos(s2->y*M_PI/180.);
	  double distsec = sqrt(sq((ra_dec_im.x-s2->x)*cosdec)+sq(ra_dec_im.y-s2->y))*3600.;
	  fprintf(matchstream,"%8.2f %8.2f %12.6f %12.6f %12.1f %12.6f %12.6f %6.2f %12.1f %8.3f\n",
	      //             x     y     ra     dec   flux   ra      dec   mag  apflux
	      s1->x , s1->y , ra_dec_im.x, ra_dec_im.y, s1->flux, s2->x , s2->y, s2->flux, s1->Flux_auto(), distsec);
	}
    }
  fclose(matchstream);
}

//! fills a file to be used when building finding charts, or to check the match.
static void FillMatchFile(const DbImage &Image, const Gtransfo &Pix2RaDec, const StarMatchList &MatchList)
{
  FitsHeader head(Image.FitsImageName(Calibrated));
   FillMatchFile(head,
		 SEStarList(Image.ImageCatalogName()), MatchList,
		 Pix2RaDec, Image.Dir()+"/match_usno.dat");
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
      Point ra_dec_im = Pix2RaDec.apply(star);
      fprintf(matchstream,"%8.2f %8.2f %12.6f %12.6f %12.1f %12.6f %12.6f %6.2f %12.1f %8.3f\n",
	      //             x     y     ra     dec   flux   ra      dec   mag  apflux
	      star.x, star.y, ra_dec_im.x, ra_dec_im.y, star.flux,ra_dec_im.x, ra_dec_im.y, -10., star.flux, 0.);
    } 

  fclose(matchstream);
}

/********************************************************************/
/*************         USNO PROCESS STUFF                *************/
/********************************************************************/




MatchCards::MatchCards()
{
   linMatchCut = 1.5;
   linMatchMinCount = 10;
   distortionDegree = 3;
   secondMatchCut = 1.0;
   writeWCS = true;
   asciiWCS = false;
   wcsFileName = "";
   astromCatalogName = "";
   dumpMatches=true;
   ignoreSatur = false;
   ignoreBad = false;
   cards_read = false;
   minSigToNoise = 0.;
}

bool MatchCards::ReadCards(const string &FileName)
{
  if (!FileExists(FileName))
    {
      std::cerr << " cannot open datacards filename " << FileName << std::endl;
      return false;
    }
  if (cards_read) return false;
  if (getenv("USNOFILE")) astromCatalogName = getenv("USNOFILE");
  DataCards cards(FileName);
#define READ_IF_EXISTS(VAR,TAG,TYPE) \
if (cards.HasKey(TAG)) VAR=cards.TYPE(TAG)
  READ_IF_EXISTS(linMatchCut,"LIN_MATCH_CUT",DParam);
  READ_IF_EXISTS(linMatchMinCount,"LIN_MATCH_MIN_COUNT",IParam);
  READ_IF_EXISTS(distortionDegree,"DISTORTION_DEGREE",IParam);
  READ_IF_EXISTS(secondMatchCut,"SECOND_MATCH_CUT",DParam);
  READ_IF_EXISTS(writeWCS,"WRITE_WCS",IParam);
  READ_IF_EXISTS(asciiWCS,"ASCII_WCS",IParam);
  READ_IF_EXISTS(wcsFileName,"WCS_FILE_NAME",SParam);
  READ_IF_EXISTS(astromCatalogName,"ASTROM_CATALOG_NAME",SParam);
  READ_IF_EXISTS(dumpMatches,"DUMP_MATCHES",IParam);
  READ_IF_EXISTS(ignoreBad,"IGNORE_BAD",IParam);
  READ_IF_EXISTS(ignoreSatur,"IGNORE_SATUR",IParam);
  READ_IF_EXISTS(minSigToNoise, "MIN_SIG_TO_NOISE", DParam);

  astromCatalogName = DbConfigFindCatalog(astromCatalogName);

  if (asciiWCS && wcsFileName == "")
    {
      std::cout << " you have requested an ascii WCS file without providing a WCS_FILE_NAME" << std::endl;
      wcsFileName = "ascii_wcs.head";
      std::cout << " using \"" << wcsFileName << "\"" << std::endl;
    }
  cards_read = true;
  return true;
}


static bool UsnoCollect(const Frame &usnoFrame, const TanPix2RaDec &Wcs, BaseStarList &UsnoCat)
{
  UsnoCat.clear();
  UsnoRead(usnoFrame, RColor, UsnoCat);
  if (UsnoCat.size() == 0)
    {
      std::cerr << " Could not collect anything from a ref catalog : giving up" << std::endl;
      throw(MatchException("No objects grabbed in the reference catalog.")); 
      return false;
    }

  ConvertMagToFlux( &UsnoCat);

  TanRaDec2Pix UsnoToPix = Wcs.invert();

  // convert the ( ra, dec) to pixel scale
  UsnoCat.ApplyTransfo(UsnoToPix);
  return true;
}



/* Yes, this routine is awfully long. It does:

   1) get both object lists, and needed info from image header
   2) try 2 matching algorithms and check the result
   3) from the found match, fit distorted relation, recollect matches, and refit
   4) split the result into a linear term (the CDi_j in WCS's and distortion terms
   5) output all that, lists, Fitsheader
   6) figure out a zero point and write it into the header
   
*/
bool UsnoProcess(const string &fitsFileName, const string &catalogName, 
		 DbImage *dbimage)
{

  // Part 1.
    TanPix2RaDec guessWcs;
    Frame skyRegion, pixFrame;
    double pixSize = 0;

    // Part 1.1 Collect needed stuff from header and provide some printouts

    {  // cfitsio does not accept to reopen RW a file already opened RO: 
      FitsHeader header(fitsFileName);
      const string cat_name((MatchPrefs.astromCatalogName != "" )? MatchPrefs.astromCatalogName : "USNO-R");
      
      cout << endl << "Matching " << fitsFileName << " with " << cat_name  << endl;
      cout << "INSTRUMENT = " << header.KeyVal("TOADINST") << endl;
      cout << "CCD no = " << header.KeyVal("TOADCHIP") << endl;
      cout << "FILTER = " << header.KeyVal("TOADFILT") << endl;
      cout << "OBJECT = " << header.KeyVal("TOADOBJE") << endl;
      cout << endl;
      if (!GuessLinWCS(header, guessWcs))
	{
	  throw(MatchException(" UsnoProcess: no WCS in file="+fitsFileName));
	  return false;
	}
      pixFrame= Frame(header);
      skyRegion = ApplyTransfo(pixFrame, guessWcs, LargeFrame).Rescale(1.2); // add 20% , i.e, 10 % on all sides
      pixSize = header.KeyVal("TOADPIXS");
    }

  // Part 1.2 load image catalog
  SEStarList sestarlist(catalogName);
  int initialSize = sestarlist.size();
  cout << " Read " << initialSize << " objects from SexCat " << std::endl; 

  /* remove objects flagged as "bad" (Ccd defects and cosmics mainly
  but also some bright stars that may be useful for matching). 
  This is  a good idea when using a "dense catalog".   */
  if (MatchPrefs.ignoreBad)
    {
      for (SEStarIterator i = sestarlist.begin(); i !=sestarlist.end(); )
	{
	  SEStar &s = **i;
	  if (s.FlagBad()) i=sestarlist.erase(i);
	  else ++i;
	}

       std::cout << " left with " << sestarlist.size() 
		 << " after removing bad objects" << std::endl;
    }

  
  // 
  // If a min sig-to-noise is requested, remove the objects with a
  // sig-to-noise exceeding this value. 
  // --nrl  [2010-04-19 Mon]
  // 
  if( MatchPrefs.minSigToNoise > 0. )
    {
      unsigned int nb_killed = 0;
      for( SEStarIterator I = sestarlist.begin(); I!=sestarlist.end(); )
	{
	  SEStar & s = **I;
	  if( s.EFlux() < 0. ) 
	    {
	      I = sestarlist.erase(I);
	      nb_killed++;
	      continue;
	    }
	  if( s.flux / s.EFlux() < MatchPrefs.minSigToNoise ) 
	    {
	      I = sestarlist.erase(I);
	      nb_killed++;
	      continue;
	    }
	  I++;
	}
      cout << "[UsnoProcess] removed " << nb_killed << " objects of SNR lower than "
	   << MatchPrefs.minSigToNoise << endl
	   << "              left with " << sestarlist.size() 
	   << " objects."
	   << endl;
    }
  
  
  
  // now count "good objects" for judging latter the match quality (from the fraction of matched objects)
  unsigned ngoodImageObjects = 0;
  for (SEStarCIterator i = sestarlist.begin(); i !=sestarlist.end(); ++i )
    {
      const SEStar &s = **i;
      if (s.FlagBad() || s.Flag()) continue;
      if (s.EFlux() > 0.2*s.flux) continue;
      ngoodImageObjects++;
    }
  cout << "number of image objects with flag=0 && flagbad=0 && S/N>5 : " 
       << ngoodImageObjects << endl;
  

  // if requested, remove saturated stars
  if (MatchPrefs.ignoreSatur)
    {
      for (SEStarIterator i = sestarlist.begin(); i !=sestarlist.end(); )
	{
	  SEStar &s = **i;
	  if (s.IsSaturated()) i = sestarlist.erase(i);
	  else ++i;
	}
      std::cout << " left with " << sestarlist.size() 
		<< " after satur removal" << std::endl;
    }

  // basic check
  if (sestarlist.size()==0) 
    {
      throw(MatchException(" UsnoProcess: empty image catalog file="+catalogName));
      return false;
    }

  // Part 1.3 collect reference catalog sub part.
  TanRaDec2Pix RaDecToPix = guessWcs.invert();
	
  BaseStarList usnoList;
  UsnoCollect(skyRegion, guessWcs, usnoList);
  //  DEBUG
  //  usnoList.write("usno.list");

  // Part 2 : find an initial match.
  StarMatchList *initialMatch = NULL;
  BaseStarList usnoListCopy(usnoList);
  BaseStarList *basestarlist = SE2Base(&sestarlist);

  // 2 methods for finding initial match
  for (int method = 1; method <=2; ++method)
    {
      GtransfoLin *guessCorr = NULL;
      if (method == 2) 
	{ // usnoList may have been corrupted, re-initialize it:
	  usnoListCopy.CopyTo(usnoList);
	}

      // Part 2.1 : guess

      bool can_try_shift = (getenv("NO_SHIFT") == NULL);
      if (method == 1 && can_try_shift)   // try a shift
	{
	  unsigned minShiftMatches = min(min(ngoodImageObjects*2/3, unsigned(25)),
					 unsigned(usnoList.size()/3));
	  guessCorr = FindShift(*basestarlist, usnoList, minShiftMatches);
	}
    

      if (!guessCorr && method == 1) method++;
      if (method == 2)
	{
	  MatchConditions conditions;
	  // areas 
	  double usnoWindowSizeInPix = skyRegion.Area()*fabs(RaDecToPix.Jacobian(RaDecToPix.TangentPoint()));
	  double imageSize = pixFrame.Area();
	  
	  conditions.NStarsL1 = 80;
	  conditions.NStarsL2 = int (conditions.NStarsL1 * (usnoWindowSizeInPix/imageSize));
	  //	  conditions.PrintLevel = 1;
	  conditions.MaxTrialCount = 16;
	  guessCorr = FindCombiMatch(*basestarlist, usnoList, conditions);

	}
      if (guessCorr)
	cout << " guessed correction (in pixel space) to WCS " << endl << *guessCorr << endl;

      if (!check_guess(*guessCorr))
	{
	  cout << " the guess is really too far from a rotation " << endl;
	  delete guessCorr;
	  guessCorr = NULL;
	  if (method == 1) continue;
	}	  

      if (!guessCorr)
	throw(MatchException(" UsnoProcess : could not guess an initial transformation "));

      // Part 2.2 : check that the initialy collected ref catalog roughly covers the whole image
      Frame newSkyRegion = ApplyTransfo(pixFrame,guessWcs*(*guessCorr),LargeFrame).Rescale(1.1);
      Frame overlap = newSkyRegion*skyRegion;
      if (overlap.Area() < newSkyRegion.Area()*0.9) // poor initial guess
	{
	  cout << " FindMatchUsno: bad initial guess of window on USNO, recollecting" << endl;
	  cout << " shifting by " 
	       << newSkyRegion.Center()-skyRegion.Center() << endl;
	  // keep the convention that usnoList is expressed in coordinates defined by guessWcs
	  UsnoCollect(newSkyRegion, guessWcs, usnoList); 
	}

      // Part 2.3 : Recollect a better match with a controlled accuracy ...
      StarMatchList *match = ListMatchCollect(*basestarlist, 
					      usnoList, 
					      guessCorr, 
					      MatchPrefs.linMatchCut/pixSize);
      // add here a clean, because it will be applied later and might kill the guess
      match_clean(*match);
      // ... and check if it is acceptable
      int nMatches = match->size();
      cout << " collected " << nMatches << " matches with tolerance " 
	   << MatchPrefs.linMatchCut << " arcsec " << endl;
      int minMatchCount = min(100, int(min(0.2*usnoList.size(), ngoodImageObjects*0.5)));
      if (nMatches > minMatchCount)
	{
	  initialMatch = match;
	  initialMatch->RefineTransfo(5.);
	  break; // we have a guess which provides an adequate number of matches, we now build a WCS.
	}
      delete match;
      if (method == 1) cout << " poor guess : trying something else " << endl;

    } // end loop on matching methods

  if (!initialMatch)
    {
      cout << " refined failed, giving up for " << fitsFileName << endl;
      throw(MatchException("Matching transfo refined failed" ));
      return false;
    }

  // up to here, the coordinates of both lists were pixels
  // we convert both the current guess wcs correction and the current 
  // usnoList to degrees in the tangent plane. The choosen tangent point
  // is the one of  guessWcs

  // Part 2.4 : transform outputs to degrees in the tangent plane

  GtransfoLin guessCorr = *dynamic_cast<const GtransfoLin *>(initialMatch->Transfo());
  GtransfoLin pix2TP = guessWcs.LinPart() * guessCorr;
  usnoList.ApplyTransfo(guessWcs.LinPart());
  Point tangentPoint = guessWcs.TangentPoint();
  // we have changed coordinate system, rather than hacking the StarMatchList, we just recollect:
  delete initialMatch;
  initialMatch = ListMatchCollect(*basestarlist, 
				  usnoList, 
				  &pix2TP, 
				  MatchPrefs.linMatchCut/3600.); // linMatchCut (& co)  are provided in arcsec in TP 

  // Done, we can proceed towards fitting a distorted transformation, recollecting matches and refitting

  // Part 3 :  account for  distortions

  // Part 3.1 : fit with a higher degree
  // which distortion degree?
  int match_order = MatchPrefs.distortionDegree;
  int nMatches = initialMatch->size();
  // check that the problem is enough constrained
  while (((match_order+1)*(match_order+2))/2 > nMatches)
    {
      match_order--;
      cout << " reducing match_order to " << match_order << endl;
    }
  initialMatch->SetTransfoOrder(match_order);
  match_clean(*initialMatch);
  cout << " Number of matches after clean " << initialMatch->size() << endl;
  initialMatch->RefineTransfo(3);
  // recheck list size
  
  int new_match_order = match_order;
  nMatches = initialMatch->size();
  while (((new_match_order+1)*(new_match_order+2))/2 > nMatches)
    {
      new_match_order--;
      cout << " number of matches too short for the initially selected order " 
	   << " order = " << new_match_order << endl;
      
      initialMatch->SetTransfoOrder(new_match_order);
      initialMatch->RefineTransfo(3);
      match_order = new_match_order;
      if (match_order == 1) break;
    }

  std::cout << ' ' << initialMatch->size() << " matches after refine," 
	    << " (1d res= " 
	    << initialMatch->Residual()*3600. << "\"," 
	    << "order = " << match_order << ')'	<< std::endl;

  // Part 3.2 : recollect matches using this better transformation
  StarMatchList *accurateMatch = ListMatchCollect(*basestarlist, 
					     usnoList, 
					     initialMatch->Transfo(), 
					     MatchPrefs.secondMatchCut/3600.);
  std::cout << " collected " << accurateMatch->size() 
	    << " matches (cut = " << MatchPrefs.secondMatchCut 
	    << "\", transfo order = "  << initialMatch->TransfoOrder() 
	    << ")" <<  std::endl;
  accurateMatch->SetTransfoOrder(match_order);
  accurateMatch->RefineTransfo(3);
  delete initialMatch; initialMatch = NULL;


  new_match_order = match_order;
  nMatches = accurateMatch->size();
  while (((new_match_order+1)*(new_match_order+2))/2 > nMatches)
    {
      new_match_order--;
      cout << " number of matches too short for the initially selected order " 
	   << " order = " << new_match_order << endl;
      
      accurateMatch->SetTransfoOrder(new_match_order);
      accurateMatch->RefineTransfo(3);
      match_order = new_match_order;
      if (match_order == 1) break;
    }


  std::cout << " number of matches after refine "  << accurateMatch->size() 
	    << ", " << "1d residual " 
	    << FitResidual(accurateMatch->Dist2(), *accurateMatch, 
			   *accurateMatch->Transfo())*3600 
	    << std::endl;

  // OK, Bizarre to see that here.... never mind
  char *matchFileName;
  if ((matchFileName = getenv("DUMP_MATCH")))
    {
      accurateMatch->write(matchFileName, &pix2TP);
      BaseStarList imageReducedList;
      for (StarMatchCIterator i=accurateMatch->begin(); i != accurateMatch->end(); ++i)
	{
	  imageReducedList.push_back((i->s1));
	}
      imageReducedList.write((string(matchFileName)+".reduced").c_str());
    }

  // Part 4.
  /* We split the transfo into 2 parts : a linear one for use by 
     official tools (e.g. saoimage), and a cubic correction 
     to the fit to be used by our stuff (for astrometry).
     

  */
  // Part 4.1 : get the best linear description of our match, so that tools that ignore distortion terms
  // in WCS (such as ds9) are not too wrong
  GtransfoLin linFit;
  linFit.fit(*accurateMatch);
  double residual = 3600.*FitResidual(*accurateMatch, linFit);
  cout << " 1d residual before fitting corrections " << residual 
       << " arcsec " << endl;
  
  // Part 4.2 : fit the distortions. In principle the linear part is essentialy unity.

  GtransfoPoly *correction = NULL;

  if (match_order>=2)
    {
      if (match_order <= 3) correction = new GtransfoPoly(match_order);
      if (correction) 
	{
	  StarMatchList linFitApplied;
	  accurateMatch->ApplyTransfo(linFitApplied, &linFit);
	  correction->fit(linFitApplied);
	  residual = FitResidual(linFitApplied, *correction)*3600;
	  cout << " 1d residual after fitting corrections " 
	       << residual << " arcsec " << endl;
	}
      else
	{
	  cerr << " ERROR in " << __FILE__ << endl;
	  cerr << " match_order is wrong " << match_order 
	       << " no corrections for non-linear distortions " << endl;
	}
      cout << " fitted correction :  "  << endl << *correction << endl;
    }

  // Part 5 : outputs
  // wether we write the WCS into the actual fits or in some ASCII fits-like header (A la scamp) 
  string outFitsFileName = fitsFileName;
  if (MatchPrefs.wcsFileName != "")
    {
      if (dbimage)
	{
	  if (strstr(MatchPrefs.wcsFileName.c_str(),"%s"))
	    {
	      char toto[128];
	      sprintf(toto, ("%s/"+MatchPrefs.wcsFileName).c_str(),
		      dbimage->Dir().c_str(), dbimage->Name().c_str());
	      outFitsFileName = toto;
	    }
	  else
	    outFitsFileName = dbimage->Dir()+"/"+MatchPrefs.wcsFileName;
	}
      else outFitsFileName = MatchPrefs.wcsFileName;
    }

  if (outFitsFileName == fitsFileName && MatchPrefs.asciiWCS)
    {
      MatchPrefs.writeWCS = false;
      std::cout << " I refuse to overwrite a fits file by an ascii file " << std::endl;
      std::cout << " no output " << std::endl;
    }
	
  string asciiWcsFileName;
  if (MatchPrefs.asciiWCS && MatchPrefs.writeWCS)
    {
      asciiWcsFileName = outFitsFileName;
      outFitsFileName = tmpnam(NULL);
    }

  // Part 5.1 : write the WCS in the fits header :
  if (MatchPrefs.writeWCS)
    {
      TanPix2RaDec pix2RaDec(linFit, tangentPoint);
      if (correction) pix2RaDec.SetCorrections(correction);
      cout << " writing WCS to " << outFitsFileName << endl; 
      TanWCS2Header(outFitsFileName, pix2RaDec);
      FitsHeader header(outFitsFileName, RW);
      header.AddOrModKey("RMATCHUS", residual, 
			 " 1d geom residual to ref catalog (arcsec)");
      header.AddOrModKey("NMATCHUS", int(accurateMatch->size()), 
			 " number of objects matched to ref catalog");
      string astromRef;
      if (MatchPrefs.astromCatalogName != "") 
	astromRef = MatchPrefs.astromCatalogName;
      else 
	{
	  char *usno_dir = getenv("USNODIR");
	  if (usno_dir) astromRef = BaseName(usno_dir); // strip path
	  else 
	    {
	      std::cerr << " ERROR: Cannot figure out which astrom catalog I used...., contact maintainer " << std::endl;
	    }
	}
      header.AddOrModKey("REFCAT",BaseName(astromRef), " Name of the ref. cat. astrom. and photom.");
      
      // Part 5.2 write to disk some extra info: the actual matches
      if (MatchPrefs.dumpMatches)
	{
      /* before writing lists to disk, we transform the usnoList back to ra
	 and dec, i.e we transform it from the tangent plane to
	 sideral sphere.  This is needed because accurateMatch has
	 references to UsnoList elements. The Transfo reads (using the
	 fact that GtransfoLin() == id)) :
      */

	  GtransfoLin id;
	  TanPix2RaDec tan2RaDec(id, pix2RaDec.TangentPoint());
	  usnoList.ApplyTransfo(tan2RaDec);
	  // write the match_usno.dat file
	  if (dbimage) 
	    {
	      FillMatchFile(*dbimage, pix2RaDec, *accurateMatch);
	      FillAllStarsFile(*dbimage, pix2RaDec);
	    }
	  else 
	    {
	      FitsHeader head(fitsFileName);
	      FillMatchFile(head, sestarlist, 
			    *accurateMatch, pix2RaDec,
			    DirName(fitsFileName)+"/match_usno1.dat");
	    } 
	}
    }

  // Part 6 : photometric zero point   
  double zeropoint, errzero;
  if (GetUsnoZeroPoint(accurateMatch, RColor, zeropoint, errzero))
    {

      cout << "1d Residue of the match with usno catalog " << residual << endl
	   << " Zero point " << zeropoint << endl
	   << "Dispersion on Zero Point (in mag)" << errzero << endl;

      if (MatchPrefs.writeWCS)
	{
	  FitsHeader header(outFitsFileName, RW);
	  header.AddOrModKey("ZEROUSNO", zeropoint, 
			     "Zero point : Mag(usno)=zero+Mag(imag)");
	  header.AddOrModKey("DZEROUSN", errzero, 
			     "Dispersion on Zero Point (in mag)");
	  /* ZP key is considered as a serious zero point:
	     write it if we use a "serious" catalog
	  */
	  if (MatchPrefs.astromCatalogName != "" )
	    header.AddOrModKey("ZP", zeropoint,
			       " zp w.r.t. "+
			       BaseName(MatchPrefs.astromCatalogName));
	}

      double skysig;
      double seeing;
      if (true)
	/* still this cfitsio limitation: cannot reopen RW a file
	   already opened RO */
	{
	  FitsHeader header(outFitsFileName);
	  seeing = header.KeyVal("SESEEING");
	  skysig = header.KeyVal("SEXSIGMA");
	}
      double flux21 = pow(10., -(21.-zeropoint)/2.5);
      double flux23 = pow(10., -(23.-zeropoint)/2.5);
      double flux25 = pow(10., -(25.-zeropoint)/2.5);
      double aine = M_PI * (2.36 * 2. * seeing) * (2.36 * 2. * seeing);
      double ssn21 = flux21 / sqrt(flux21 + aine * skysig * skysig);
      double ssn23 = flux23 / sqrt(flux23 + aine * skysig * skysig);
      double ssn25 = flux25 / sqrt(flux25 + aine * skysig * skysig);


      std::cout << "Signal/Noise (at mag 21) = " 
		<< ssn21 << " (at mag 25) = " << ssn25 << std::endl;

      if (MatchPrefs.writeWCS)
	{
	  FitsHeader header(outFitsFileName, RW);
	  header.AddOrModKey("USNOSB21", ssn21, "Signal / Noise at mag = 21");
	  header.AddOrModKey("USNOSB23", ssn23, "Signal / Noise at mag = 23");
	  header.AddOrModKey("USNOSB25", ssn25, "Signal / Noise at mag = 25");
	}
    }

  // if ascii output requested, convert now
    if (MatchPrefs.asciiWCS && MatchPrefs.writeWCS)
      {
	FitsHeader head(outFitsFileName);
	ofstream s(asciiWcsFileName.c_str());
	s << head << "END    " << std::endl;
	s.close();
      }
    if (MatchPrefs.asciiWCS && MatchPrefs.writeWCS)
      {
	// remove the tmp file (it was open in the previous block!)
	remove(outFitsFileName.c_str());
      }
  return true;
}


#ifdef COMMENT

/* found in read.use in the USNO I catalog directory */

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
