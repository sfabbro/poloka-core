#include <iostream>
#include <cmath>
#ifndef M_PI
#define     M_PI            3.14159265358979323846  /* pi */
#endif

#include "basestar.h"
#include "starlist.h"
#include "histo2d.h"
#include "gtransfo.h"
#include "starmatch.h"

#include "listmatch.h"
#include "starmatch.h"
#include "fastfinder.h"


// cuts.. limits, etc for combinatorial match


#include "datacards.h"

static void read_card(DataCards &Cards, const std::string &Key, double &D)
{
  if (Cards.HasKey(Key)) D= Cards.DParam(Key);
}

static void read_card(DataCards &Cards, const std::string &Key, int &I)
{
  if (Cards.HasKey(Key)) I = Cards.IParam(Key);
}

static double sqr(const double& x) { return x*x; }

#include "fileutils.h"

MatchConditions::MatchConditions(const std::string &DatacardsName)
{  /* default values */
  NStarsL1 = 70; NStarsL2=70;
  MaxTrialCount = 4;
  NSigmas = 3.;
  MaxShiftX = 50; MaxShiftY = 50; 
  SizeRatio = 1; DeltaSizeRatio = 0.1*SizeRatio;
  MinMatchRatio = 1./3.;
  PrintLevel = 0;
  Algorithm = 2;
  if (DatacardsName != "")
    {
      if (!FileExists(DatacardsName))
	{
	  cerr << " cannot find file " << DatacardsName << " to read MatchConditions " << endl;
	  return;
	}
      DataCards cards(DatacardsName);
      read_card(cards,"DATMATCH_NL1", NStarsL1);
      read_card(cards,"DATMATCH_NL2", NStarsL2);
      read_card(cards,"DATMATCH_MAXTRIAL", MaxTrialCount);
      read_card(cards,"DATMATCH_NSIG_CUT", NSigmas);
      read_card(cards,"DATMATCH_MAXSHIFTX", MaxShiftX);
      read_card(cards,"DATMATCH_MAXSHIFTY", MaxShiftY);
      read_card(cards,"DATMATCH_SIZERATIO", SizeRatio);
      read_card(cards,"DATMATCH_DELTA_SIZERATIO", DeltaSizeRatio);
      read_card(cards,"DATMATCH_MINMATCHRATIO", MinMatchRatio);
      read_card(cards,"DATMATCH_PRINTLEVEL", PrintLevel);
      read_card(cards,"DATMATCH_ALGO", Algorithm);
    }
  if (getenv("MATCH_PRINT")) PrintLevel = atoi(getenv("MATCH_PRINT"));
}

/* a Segment is a pair of stars form the same image. it is used for matching starlists */

class Segment{
   public:
  /* data */
  double r, dx,dy;
  const BaseStar *s1,*s2;
  int s1rank;

  /* constructor (could set last argument to identity by default)  */
  Segment(const BaseStar *S1, const BaseStar *S2, const int S1Rank, const Gtransfo &Tin)
      {s1rank = S1Rank; s1=S1; s2=S2;  
      Point P1 = Tin.apply(*S1); Point P2= Tin.apply(*S2); dx = P2.x - P1.x; dy = P2.y - P1.y; r = sqrt(dx*dx+dy*dy);}

  /* arg(Seg2/(*this)) if considered as complex(dx,dy) */
  double relative_angle(Segment *Seg2) 
           { return atan2( Seg2->dx*dy - dx*Seg2->dy,  dx*Seg2->dx + dy*Seg2->dy);}

  friend ostream& operator <<(ostream & stream, const Segment &S)
           { stream << " dx " << S.dx << " dy " << S.dy << " r " << S.r << endl; return stream;}
};


class SegmentList : public list<Segment> {
  public:
  //  SegmentList(const BaseStarList &L, const int NStar);
  SegmentList(const BaseStarList &L, const int NStar, const Gtransfo &Tin = GtransfoIdentity());
};

typedef list<Segment>::iterator SegmentIterator;
typedef list<Segment>::const_iterator SegmentCIterator;



static bool DecreasingLength(const Segment &first, const Segment &second)
{
return (first.r > second.r);
}


SegmentList::SegmentList(const BaseStarList &L, const int NStars, const Gtransfo &Tin) 
{
  BaseStarCIterator si1,si2, siStop;
  const BaseStar *s1, *s2;

/* find the fence */
  siStop = L.begin();
  int limit = min(NStars, int(L.size())) - 1; // -1 because test happens after incrementation
//cout << " DEBUG NStars " << NStars << " L.Size() " << L.size() << endl;
  for (int count = 0; count < limit ; count++) ++siStop;


  // iterate on star pairs
  int rank = 0;
  for (si1 = L.begin(); si1 != siStop; ++si1, rank++) for (si2 = siStop; si2 != si1; --si2)
    {
      s1 = *si1; s2 = *si2;
      if (s1 == s2) cerr << " Catastrophe in SegmentList::SegmentList" << endl;
      push_back(Segment(s1,s2, rank, Tin));
      
    }
this->sort(DecreasingLength); /* allows a break in loops */
// cout << " DEBUG : size " << size() << endl;
}


//#include <pair>

struct SegmentPair : public pair<Segment*, Segment*> {
SegmentPair (Segment *f, Segment *s) : pair<Segment*, Segment*>(f,s) {};
};

typedef list<SegmentPair> SegmentPairList;
typedef SegmentPairList::iterator SegmentPairListIterator;
typedef SegmentPairList::const_iterator SegmentPairListCIterator;


static StarMatchList* MatchListExtract(const SegmentPairList &PairList, int Rank1, int Rank2, const Gtransfo &Tin)
{
  /* first Select in the segment pairs list the ones which make use of star rank1 in segment1
and star s2 in segment2 */

StarMatchList *matchList = new StarMatchList;

for (SegmentPairListCIterator spi = PairList.begin(); spi != PairList.end(); spi++)
  {
  const SegmentPair& a_pair = *spi;
  if (a_pair.first->s1rank != Rank1 || a_pair.second->s1rank != Rank2) continue;
  /* now we store as star matches both ends of segment pairs , 
     but only once the beginning of segments because they all have the same, 
     given the selection 3 lines above  */
  if (matchList->size() == 0) 
     matchList->push_back(StarMatch(Tin.apply(*(a_pair.first->s1)), *(a_pair.second->s1), 
                            a_pair.first->s1, a_pair.second->s1));
  /* always store the match at end */
  matchList->push_back(StarMatch(Tin.apply(*(a_pair.first->s2)), *(a_pair.second->s2), 
                            a_pair.first->s2, a_pair.second->s2));
  }
// cout << " matchList size " << matchList->size() << endl;
return matchList;
}


static bool DecreasingQuality(const StarMatchList* first, const StarMatchList *second)
{
int idiff = first->size() - second->size();
if   (idiff != 0) return ( idiff > 0); else return(first->Dist2() < second->Dist2());
}

/* many matching solutions (StarMatchList) will be compared. Store them in a SolList : */

class SolList : public list<StarMatchList*> {  /* building a class just to have an " active destructor" */
public :
typedef SolList::iterator SolIterator;
  ~SolList() {  for (SolIterator s = begin(); s !=end(); ++s) delete *s;}
};


static void dump_input_list(const BaseStarList &L, const int NStars, 
			      const Gtransfo &Tin, const std::string FileName)
{
  BaseStarList list;
  L.ExtractHead(list,NStars);
  list.ApplyTransfo(Tin);
  list.write(FileName);
}
  



/* This one searches a general transformation by histogramming the relative size and orientation
of star pairs ( Segment's) built from the 2 lists */


static StarMatchList *ListMatchupRotShift_Old(BaseStarList &L1, BaseStarList &L2, 
                                          const Gtransfo &Tin, const MatchConditions &Conditions)
{
SegmentList sList1(L1, Conditions.NStarsL1, Tin);
SegmentList sList2(L2, Conditions.NStarsL2, GtransfoIdentity());

  /* choose the binning of the histogram so that 
     1: ratio = 1 and rotation angle = n * (pi/2) are bin centers. since 
     the angle is computed using atan2, its range is [-pi,pi],
     and the histogram range is [-pi-eps, pi-eps], so
     if (angle>pi- angleOffset) angle -= 2*pi before filling.   */
int nBinsR = 21;
int nBinsAngle = 180; /* can be divided by 4 */
double angleOffset = M_PI/nBinsAngle;
double minRatio = Conditions.MinSizeRatio();
double maxRatio = Conditions.MaxSizeRatio();
Histo2d histo(nBinsR, minRatio, maxRatio,  
                  nBinsAngle, -M_PI- angleOffset, M_PI - angleOffset);

SegmentIterator segi1,segi2;
Segment *seg1,*seg2; 
double ratio, angle;
for (segi1 = sList1.begin(); segi1 != sList1.end(); ++segi1)
  {
  seg1 = &(*segi1);
  if (seg1->r == 0) continue;
  for (segi2 = sList2.begin(); segi2 != sList2.end(); ++segi2)
    {
    seg2 = &(*segi2);
      /* if one considers the 2 segments as complex numbers z1 and z2, ratio=mod(z1/z2) and angle = arg(z1/z2) */
      /* I did not put a member function in Segment to compute both because we apply a cut on ratio before actually
         computing the angle (which involves a call to atan2 (expensive)) */
    ratio = seg2->r/seg1->r;
    if (ratio > maxRatio) continue;
    if (ratio < minRatio) break;  /* use th fact that segment lists are sorted by decresing length */
    angle = seg1->relative_angle(seg2);
    if (angle > M_PI - angleOffset) angle -= 2.*M_PI;
    histo.Fill(ratio,angle);
    }
  }
double binr, bina;
histo.BinWidth( binr, bina);

SolList Solutions;
/* now we want to find in the (r,theta) bins that have the highest counts, the star pair
   (one in l1, one in L2) that contribute to the largest number of segment pairs in this bin :
   so, we histogram a couple of integer that uniquely defines the stars, for the segment pairs
   that contribute to the maximum bin. We choose to histogram the rank of s1 of segment 1
   versus the rank of s1 for segment 2 */

for (int i = 0; i<Conditions.MaxTrialCount; ++i)
  {
  double ratioMax, angleMax;
  double maxContent = histo.MaxBin(ratioMax,angleMax);
  histo.Fill(ratioMax,angleMax, - maxContent);
  if (Conditions.PrintLevel >= 1) 
    {
      cout << " valMax " << maxContent
	   << " ratio " << ratioMax << " angle " << angleMax << endl;
    }
  minRatio = ratioMax - binr/2; maxRatio = ratioMax+binr/2;
  double minAngle = angleMax - bina/2; double maxAngle = angleMax+bina/2;
  SegmentPairList pairList;
  Histo2d historank(Conditions.NStarsL1, 0., Conditions.NStarsL1, 
                    Conditions.NStarsL2, 0., Conditions.NStarsL2); 
  /* reloop on segment pairs to select the ones in this specific bin */
  
  for (segi1 = sList1.begin(); segi1 != sList1.end(); ++segi1)
    {
    seg1 = &(*segi1);
    if (seg1->r == 0) continue;
    for (segi2 = sList2.begin(); segi2 != sList2.end(); ++segi2)
      {
      seg2 = &(*segi2);
      ratio = seg2->r/seg1->r;
      if (ratio > maxRatio) continue;
      if (ratio < minRatio) break;  /* use the fact that segment lists are sorted by decresing length */
      angle = seg1->relative_angle(seg2);
      if (angle > M_PI - angleOffset) angle -= 2.*M_PI;
      if (angle < minAngle || angle > maxAngle) continue;
      pairList.push_back(SegmentPair(seg1, seg2)); /* store the match */
      historank.Fill(seg1->s1rank + 0.5 , seg2->s1rank+0.5);      
      }
    }
  for (int iteration=0; iteration < Conditions.MaxTrialCount; iteration++)
    {
    double dr1,dr2;
    double maxval = historank.MaxBin(dr1,dr2);
    /* set this bin to zero so that next iteration will find next maximum */
    historank.Fill(dr1,dr2,-maxval);
    // cout << maxval << ' ' ; cout << dr1 << ' ' << dr2 << endl;
    StarMatchList *a_list = MatchListExtract(pairList, int(dr1), int(dr2), GtransfoIdentity());
    a_list->RefineTransfo(Conditions.NSigmas); // mandatory for the sorting fields to be filled
    // cout << " list size after Refine " << a_list->size() << endl;
    Solutions.push_back(a_list);
    }
  }/* end of loop on (r,theta) bins */
  Solutions.sort(DecreasingQuality);
  StarMatchList *best = *Solutions.begin();
  /* remove the first one from the list */
  Solutions.pop_front();
  if (Conditions.PrintLevel >=1) 
    {
      cout << " best solution " << best->Residual() << " npairs " << best->size() << endl << *(best->Transfo());
      cout << " Chi2 " << best->Chi2() << ','
	   << " Number of solutions " << Solutions.size() << endl;
    }
  return best;
}



#include "histo4d.h"

/* this matching routine searches brutally a match between lists in
    the 4 parameter space: size ratio, rotation angle, x and y
    shifts. This is done by histogramming where combinations of four
    objets (2 on each list) fall in this 4 parameter space.

    One trick is that rather than using actual offsets, we histogram
    object indices of the combination: 
*/

static StarMatchList *ListMatchupRotShift_New(BaseStarList &L1, BaseStarList &L2, 
					      const Gtransfo &Tin, 
					      const MatchConditions &Conditions) 
{
  if (L1.size() <= 4 || L2.size() <= 4)
    {
      cout << " ListMatchupRotShift_New : (at least) one of the lists is too short " << endl;
      return NULL;
    }

  SegmentList sList1(L1, Conditions.NStarsL1, Tin); 
  SegmentList sList2(L2, Conditions.NStarsL2, GtransfoIdentity());

 if (getenv("DUMP_INPUT_LISTS"))
   {
     dump_input_list(L1, Conditions.NStarsL1, Tin, "input1.list");
     dump_input_list(L2, Conditions.NStarsL2, GtransfoIdentity(), "input2.list");
   }



  /* choose the binning of the histogram so that 
     1: ratio = 1 and rotation angle = n * (pi/2) are bin centers. since 
     the angle is computed using atan2, its range is [-pi,pi],
     and the histogram range is [-pi-eps, pi-eps], so
     if (angle>pi- angleOffset) angle -= 2*pi before filling.   */
int nBinsR = 21;
int nBinsAngle = 180; /* can be divided by 4 */
double angleOffset = M_PI/nBinsAngle;
double minRatio = Conditions.MinSizeRatio();
double maxRatio = Conditions.MaxSizeRatio();
 SparseHisto4d histo(nBinsR, minRatio, maxRatio,
		     nBinsAngle, -M_PI- angleOffset, M_PI - angleOffset,
		     Conditions.NStarsL1, 0., Conditions.NStarsL1, 
		     Conditions.NStarsL2, 0., Conditions.NStarsL2,
		     sList1.size()*sList2.size()); 

SegmentIterator segi1,segi2;
Segment *seg1,*seg2; 
double ratio, angle;

for (segi1 = sList1.begin(); segi1 != sList1.end(); ++segi1)
  {
  seg1 = &(*segi1);
  if (seg1->r == 0) continue;
  for (segi2 = sList2.begin(); segi2 != sList2.end(); ++segi2)
    {
    seg2 = &(*segi2);
      /* if one considers the 2 segments as complex numbers z1 and z2, ratio=mod(z1/z2) and angle = arg(z1/z2) */
      /* I did not put a member function in Segment to compute both because we apply a cut on ratio before actually
         computing the angle (which involves a call to atan2 (expensive)) */
    ratio = seg2->r/seg1->r;
    if (ratio > maxRatio) continue;
    if (ratio < minRatio) break;  /* use th fact that segment lists are sorted by decresing length */
    angle = seg1->relative_angle(seg2);
    if (angle > M_PI - angleOffset) angle -= 2.*M_PI;
    histo.Fill(ratio,angle,seg1->s1rank + 0.5 , seg2->s1rank+0.5);
    }
  }


SolList Solutions;
/* now we find the highest bins of the histogram, and recover the original objects.
   This involves actually re-looping on the combinations, but it is much 
   faster that the original histogram filling loop, since we only compute 
   angle and ratio for Segments that have the right first object
*/

 int oldMaxContent = 0; 

 for (int i = 0; i<4*Conditions.MaxTrialCount; ++i) // leave a limit to make avoid (almost)  infinite loops
   {
     double pars[4];
     int maxContent = histo.MaxBin(pars);
     if (maxContent == 0) break;
     if (Conditions.PrintLevel >=1) 
       {
	 cout << " ValMax " << maxContent 
	      <<" ratio " << pars[0] << " angle " << pars[1] << endl;
       }
     histo.ZeroBin(pars);
     if (i>0)
       { /* the match possibilities come out in a random order when they have the same content.
	    so, we stop investigating guesses when the content goes down AND the requested search depth 
            (MaxTrialCount) is reached */
	 if (maxContent < oldMaxContent && i>=Conditions.MaxTrialCount) break;

       }
     oldMaxContent = maxContent;
     /* reloop on segment pairs to select the ones in this specific bin */
     int rank1L1 = int(pars[2]);
     int rank1L2 = int(pars[3]);
     double minAngle, maxAngle;
     histo.BinLimits(pars,0, minRatio, maxRatio);
     histo.BinLimits(pars,1, minAngle, maxAngle);

     StarMatchList *a_list = new StarMatchList;
     
     for (segi1 = sList1.begin(); segi1 != sList1.end(); ++segi1)
       {
	 seg1 = &(*segi1);
	 if (seg1->s1rank != rank1L1) continue;
	 if (seg1->r == 0) continue;
	 for (segi2 = sList2.begin(); segi2 != sList2.end(); ++segi2)
	   {
	     seg2 = &(*segi2);
	     if (seg2->s1rank != rank1L2) continue;
	     // push in the list the match corresponding to end number 1 of segments
	     if (a_list->size() == 0) 
	       a_list->push_back(StarMatch(*(seg1->s1), *(seg2->s1), seg1->s1, seg2->s1));
	     ratio = seg2->r/seg1->r;
	     if (ratio > maxRatio) continue;
	     if (ratio < minRatio) break;  /* use the fact that segment lists are sorted by decresing length */
	     angle = seg1->relative_angle(seg2);
	     if (angle > M_PI - angleOffset) angle -= 2.*M_PI;
	     if (angle < minAngle || angle > maxAngle) continue;
	     /* here we have 2 segments which have the right
		- length ratio
		- relative angle
		- first objects (objects on the end number 1).
		The objects on the end number 2 are the actual matches : */
	     a_list->push_back(StarMatch(*(seg1->s2), *(seg2->s2), seg1->s2, seg2->s2));
	   }
       }

     // a basic check for sanity of the algorithm :
     
     if (int(a_list->size() ) != maxContent+1 )
       {
	 cerr << " There is an internal inconsistency in ListMatchupRotShift " << endl
	      << " maxContent  = " << maxContent << endl
	      << " matches->size() = " << a_list->size() << endl;
	 cerr << "please store the involved images and contact the developpers" << endl;
       }
     a_list->RefineTransfo(Conditions.NSigmas);
     Solutions.push_back(a_list);
    }

  if (Solutions.size() == 0)
    {
      cout << " error In ListMatchup : not a single pair match " << endl;
      cout << " Probably, the relative scale of lists is not within bounds" << endl;
      cout << " here : " << minRatio << ' ' << maxRatio << endl;
      return NULL;
    }


  Solutions.sort(DecreasingQuality);
  StarMatchList *best = *Solutions.begin();
  /* remove the first one from the list */
  Solutions.pop_front();
  if (Conditions.PrintLevel >=1) 
    {
      cout << " best solution " << best->Residual() << " npairs " << best->size() << endl << *(best->Transfo());
      cout << " Chi2 " << best->Chi2() << ','
	   << " Number of solutions " << Solutions.size() << endl;
    }
  return best;
}


static StarMatchList *ListMatchupRotShift(BaseStarList &L1, BaseStarList &L2, 
                                          const Gtransfo &Tin, const MatchConditions &Conditions)
{
  if (Conditions.Algorithm == 1) return ListMatchupRotShift_Old(L1, L2, Tin, Conditions);
  else return ListMatchupRotShift_New(L1, L2, Tin, Conditions);
}


StarMatchList *MatchSearchRotShift(BaseStarList &L1, BaseStarList &L2, const MatchConditions &Conditions)
{
L1.FluxSort();
L2.FluxSort();

return ListMatchupRotShift(L1, L2, GtransfoIdentity() , Conditions);
}

StarMatchList *MatchSearchRotShiftFlip(BaseStarList &L1, BaseStarList &L2, const MatchConditions &Conditions)
{
L1.FluxSort();
L2.FluxSort();


GtransfoLin flip(0,0,1,0,0,-1);
StarMatchList *flipped   =  ListMatchupRotShift(L1,L2,flip, Conditions);
StarMatchList *unflipped =  ListMatchupRotShift(L1,L2, GtransfoIdentity(), Conditions);
 if (! flipped  || !unflipped) return NULL;
if (Conditions.PrintLevel >=1)
  {
  cout << " unflipped  Residual " << unflipped->Residual() << " nused " << unflipped->size() << endl;
  cout << "   flipped  Residual " << flipped->Residual() << " nused " << flipped->size() << endl;
  }
if (DecreasingQuality(flipped,unflipped))
  {
   if (Conditions.PrintLevel >=1) 
       cout << " keeping flipped solution" << endl;
    delete unflipped;
    // One should NOT apply the flip to the result because the matchlist
    // (even the flipped one) contains the actual coordinates of stars.
    // MatchListExtract is always called with GtransfoIdentity() as last parameter
    return flipped;
  }
else;
  {
  if (Conditions.PrintLevel >=1) 
     cout << " keeping unflipped solution" << endl;
  delete flipped;
  return unflipped;
  }
}


#ifdef STORAGE
// timing : 2.5 s for l1 of 1862 objects  and l2 of 2617 objects
GtransfoLin *ListMatchupShift(const BaseStarList &L1, const BaseStarList &L2, const Gtransfo &Tin, double MaxShift)
{
  int ncomb = L1.size() * L2.size();
  if (!ncomb) return NULL;
  int nx;
  if (ncomb > 10000) nx = 100; else nx = (int) sqrt(ncomb);

  Histo2d histo(nx, -MaxShift, MaxShift, nx, -MaxShift, MaxShift);
  
  BaseStarCIterator s1,s2;
  double x1,y1;
  for (s1 = L1.begin(); s1 != L1.end(); ++s1)
    {
      Tin.apply((*s1)->x, (*s1)->y, x1, y1);
      for (s2 = L2.begin(); s2 != L2.end(); ++s2)
	{
	  histo.Fill((*s2)->x - x1, (*s2)->y - y1);
	}
    }
  double dx=0, dy=0;
  histo.MaxBin(dx,dy);
  return new GtransfoLinShift(dx,dy);
}
#endif /*STORAGE*/





// timing : 140 ms for l1 of 1862 objects  and l2 of 2617 objects (450 MHz, "-O4") MaxShift = 200.
GtransfoLin *ListMatchupShift(const BaseStarList &L1, const BaseStarList &L2, const Gtransfo &Tin, double MaxShift, double BinSize)
{
  int nx;
  if (BinSize == 0)
    {
    int ncomb = L1.size() * L2.size();
    if (ncomb > 10000) nx = 100; else nx = (int) sqrt(double(ncomb));
    if (!ncomb) return NULL;
    }
  else nx = int(2*MaxShift/BinSize+0.5);

  Histo2d histo(nx, -MaxShift, MaxShift, nx, -MaxShift, MaxShift);
  double binSize = 2*MaxShift/nx;
  
  BaseStarCIterator s1;
  FastFinder finder(L2);
  double x1,y1;
  for (s1 = L1.begin(); s1 != L1.end(); ++s1)
    {
      Tin.apply((*s1)->x, (*s1)->y,x1,y1);
      FastFinder::Iterator it = finder.begin_scan(Point(x1,y1), MaxShift);
      while (*it)
	{
	  const BaseStar *s2 = *it;
	  histo.Fill(s2->x - x1, s2->y - y1);
	  ++it;
	}
    }
  //  double dx=0, dy=0;
  //  histo.MaxBin(dx,dy);
  //  return new GtransfoLinShift(dx,dy);
  SolList Solutions;
  for (int i=0; i<4; ++i)
    {      
      double dx = 0, dy = 0;
      double count = histo.MaxBin(dx,dy);
      histo.Fill(dx,dy,-count); // zero the maxbin
      GtransfoLinShift shift(dx,dy);
      Gtransfo *newGuess = GtransfoCompose(&shift, &Tin);
      StarMatchList *raw_matches = ListMatchCollect(L1, L2, newGuess, binSize);
      delete newGuess;
      StarMatchList *matches = new StarMatchList;
      raw_matches->ApplyTransfo(*matches, &Tin);
      delete raw_matches;
      matches->SetTransfoOrder(1);
      matches->RefineTransfo(3.);
      //      cout << *matches->Transfo() << endl;
      Solutions.push_back(matches);
    }
  Solutions.sort(DecreasingQuality);
  GtransfoLin *best = new GtransfoLin(* const_cast<GtransfoLin*>(dynamic_cast<const GtransfoLin*> (Solutions.front()->Transfo())));
  //  cout << " best \"shift\" found " << endl << *best; 
  return best;
}



#ifdef STORAGE

// this is the old fashioned way...

StarMatchList *ListMatchCollect_Slow(const BaseStarList &L1, const BaseStarList &L2,const Gtransfo *Guess, const double MaxDist)
{
  StarMatchList *matches = new StarMatchList;
  /****** Collect ***********/
  for (BaseStarCIterator si = L1.begin(); si != L1.end(); ++si)
    {
      const Point *p1 = (*si);
      const Point p2 = Guess->apply(*p1);
      const BaseStar *neighbour = L2.FindClosest(p2);
      //   BaseStar *neighbour = finder.FindClosest(p2,MaxDist);
      if (!neighbour) continue;
      double distance =p2.Distance(*neighbour); 
      if (distance < MaxDist)
	{
	  matches->push_back(StarMatch(*p1,*neighbour,*si,neighbour));
	  // assign the distance, since we have it in hand: 
	  matches->back().distance = distance;
	}
    }
  return matches;
}
#endif

// here is the real active routine:

StarMatchList *ListMatchCollect(const BaseStarList &L1, 
				const BaseStarList &L2,
				const Gtransfo *Guess, const double MaxDist)
{
  StarMatchList *matches = new StarMatchList;
  /****** Collect ***********/
  FastFinder finder(L2);
  for (BaseStarCIterator si = L1.begin(); si != L1.end(); ++si)
    {
      const BaseStarRef &p1 = (*si);
      Point p2 = Guess->apply(*p1);
      const BaseStar *neighbour = finder.FindClosest(p2,MaxDist);
      if (!neighbour) continue;
      double distance =p2.Distance(*neighbour); 
      if (distance < MaxDist)
	{
	  matches->push_back(StarMatch(*p1,*neighbour,*si,neighbour));
	  // assign the distance, since we have it in hand: 
	  matches->back().distance = distance;
	}

    }
  matches->SetTransfo(Guess);

  return matches;
}

#ifdef STORAGE
// unused
//! iteratively collect and fits, with the same transfo kind, until the residual increases
StarMatchList *CollectAndFit(const BaseStarList &L1, const BaseStarList &L2,
			     const Gtransfo *Guess, const double MaxDist)
{
  const Gtransfo *bestTransfo = Guess;
  StarMatchList *prevMatch = NULL;
  while (true)
    {
      StarMatchList *m = ListMatchCollect(L1,L2,bestTransfo, MaxDist);
      m->SetTransfo(bestTransfo);
      m->RefineTransfo(3.);
      cout << " iterating : resid " << m->Residual() 
	   << ' ' << " size " << m->size() << endl;
      if (!prevMatch || 
	  (prevMatch 
	   && m->Residual() < prevMatch->Residual()*0.999
	   && m->Chi2()>0)
	  )
	{
	  if (prevMatch) delete prevMatch;
	  prevMatch = m;
	  bestTransfo = m->Transfo();
	}
      else 
	{
	  delete m;
	  break;
	}
    }
  return prevMatch;
}
#endif 



StarMatchList *ListMatchCollect(const BaseStarList &L1, const BaseStarList &L2, const double MaxDist)
{
  StarMatchList *matches = new StarMatchList;
  FastFinder finder(L2);
  for (BaseStarCIterator si = L1.begin(); si != L1.end(); ++si)
    {
      const BaseStarRef &p1 = (*si);
      const BaseStar *neighbour = finder.FindClosest(*p1,MaxDist);
      if (!neighbour) continue;
      double distance =p1->Distance(*neighbour); 
      if (distance < MaxDist)
	{
	  matches->push_back(StarMatch(*p1,*neighbour,*si,neighbour));
	  // assign the distance, since we have it in hand: 
	  matches->back().distance = distance;
	}
    }

  matches->SetTransfo(new GtransfoIdentity);

  return matches;
}


#ifdef STORAGE
StarMatchList *ListMatchCollectWU(const BaseStarList &L1, const BaseStarList &L2,const Gtransfo *Guess, const double MaxDist, BaseStarList &unmatches)
{
  int compteur =0;
  StarMatchList *matches = new StarMatchList;
  /****** Collect ***********/
  FastFinder finder(L2);
  for (BaseStarCIterator si = L1.begin(); si != L1.end(); ++si)
    {
      BaseStar *star = (*si);
      Point *p1 = (*si);
      Point p2 = Guess->apply(*p1);
      BaseStar *neighbour = finder.FindClosest(p2,MaxDist);
      if (!neighbour) 
	{
	  unmatches.push_back(star);
	  continue;
	}
      double distance =p2.Distance(*neighbour); 
      if (distance < MaxDist)
	{
	  ++compteur;
	  matches->push_back(StarMatch(*p1,*neighbour,*si,neighbour));
	  // assign the distance, since we have it in hand: 
	  matches->back().distance = distance;
	}
      else
	{
	  unmatches.push_back(star);
	}
    }
  
  // tentatively: put the call to the slow strightforward code in the new fast one, and compare
  // the matches count. just in case....
  StarMatchList *matches_slow = ListMatchCollect_Slow(L1, L2, Guess, MaxDist);
  if (matches_slow->size() != matches->size()) 
    {
      cerr << " CA ne va pas du tout avec le FastFinder " << endl;
      cerr <<" writing match_fast.list and match_slow.list" << endl;
      matches->write("match_fast.list"); matches_slow->write("matches_slow.list");
    }
  delete matches_slow;
  // end of the tentative check.

  return matches;
}

#endif


static bool is_transfo_ok(const StarMatchList* match, const double& pixSizeRatio2, const size_t nmin) {

  if ((fabs(fabs(dynamic_cast<const GtransfoLin*>(match->Transfo())->Determinant())-pixSizeRatio2)/pixSizeRatio2 < 0.2) && (match->size() > nmin))
    return true;
  cout << " is_transfo_ok: no\n";
  match->DumpTransfo();  
  return false;
}

// utility to check current transfo difference
static double transfo_diff(const BaseStarList &List, const Gtransfo *T1, const Gtransfo*T2) {
  double diff2 = 0;
  FatPoint tf1;
  Point tf2;
  int count = 0;
  for (BaseStarCIterator it = List.begin(); it != List.end(); ++it) {
    const BaseStar *s = *it;
    T1->TransformPosAndErrors(*s, tf1);
    T2->apply(*s, tf2);
    double dx = tf1.x - tf2.x;
    double dy = tf1.y - tf2.y;
    diff2 += (tf1.vy*dx*dx + tf1.vx*dy*dy - 2*tf1.vxy*dx*dy) / (tf1.vx*tf1.vy-tf1.vxy*tf1.vxy);
    count++;
  }
  if (count) 
    return diff2 / double(count);
  return 0;
}

static double median_distance(const StarMatchList* match, const Gtransfo* transfo) {
  size_t nstars = match->size();
  vector<double> resid(nstars);
  vector<double>::iterator ir = resid.begin();
  for (StarMatchCIterator it = match->begin(); it != match->end(); ++it, ++ir)
    *ir = sqrt(transfo->apply(it->point1).Dist2(it->point2));
  sort(resid.begin(), resid.end());
  return (nstars & 1) ? resid[nstars/2] : (resid[nstars/2-1] + resid[nstars/2])*0.5;
}


Gtransfo* ListMatchCombinatorial(const BaseStarList &List1, const BaseStarList &List2, const MatchConditions& Conditions) {
  BaseStarList L1, L2;
  List1.CopyTo(L1); L1.FluxSort();
  List2.CopyTo(L2); L2.FluxSort();

  cout << " ListMatchCombinatorial: find match between " << L1.size() << " and " << L2.size() << " stars...";
  StarMatchList *match = MatchSearchRotShiftFlip(L1, L2, Conditions);
  Gtransfo *transfo = 0;
  double pixSizeRatio2 = sqr(Conditions.SizeRatio);
  size_t nmin = min(size_t(10), size_t(min(List1.size(), List2.size())*Conditions.MinMatchRatio));

  if (is_transfo_ok(match, pixSizeRatio2, nmin))
    transfo = match->Transfo()->Clone();
  else {
    delete match;
    cout << "FAILED\n ListMatchCombinatorial: direct transfo failed, trying reverse";
    match = MatchSearchRotShiftFlip(L2, L1, Conditions);
    if (is_transfo_ok(match, pixSizeRatio2, nmin))
      transfo = match->InverseTransfo();
    else {
      cout << "FAILED\n";
      if (transfo) delete transfo;
    }
  }
  delete match;

  if (transfo) {
    cout << "FOUND\n";
    if (Conditions.PrintLevel >= 1)
      cout << " ListMatchCombinatorial: found the following transfo\n"
	   << *transfo << endl;
  } else
    cerr << "FAILED\n. Error: ListMatchCombinatorial: failed to find a transfo\n";  
  return transfo;
}

Gtransfo* ListMatchRefine(const BaseStarList& List1, const BaseStarList& List2, Gtransfo* transfo, const int maxOrder) {

  if (!transfo) { return 0; }

  // some hard-coded constants that could go in a param file
  const double brightDist = 2.;  // distance in pixels in a match
  const double fullDist = 4.;    // distance in pixels in a match between entire lists
  const double nSigmas = 3.;     // k-sigma clipping on residuals
  const size_t nStars  = 500;    // max number of bright stars to fit

  int order = 1;
  size_t nstarmin = 3;

  BaseStarList L1, L2;
  List1.CopyTo(L1); L1.FluxSort(); L1.CutTail(nStars);
  List2.CopyTo(L2); L2.FluxSort(); L2.CutTail(nStars);

  StarMatchList *fullMatch = ListMatchCollect(List1, List2, transfo, fullDist);
  StarMatchList *brightMatch = ListMatchCollect(L1, L2, transfo, brightDist);
  double curChi2 = ComputeChi2(*brightMatch, *transfo) / brightMatch->size();

  cout << " ListMatchRefine: start  "
       << " med.resid "  << median_distance(fullMatch, transfo) 
       << " #match " << fullMatch->size() 
       << endl;

  do { // loop on transfo order on full list of stars
    Gtransfo* curTransfo = brightMatch->Transfo()->Clone();
    unsigned iter = 0;
    double transDiff;
    do { // loop on transfo diff only on bright stars
      brightMatch->SetTransfoOrder(order);
      brightMatch->RefineTransfo(nSigmas);
      transDiff = transfo_diff(L1, brightMatch->Transfo(), curTransfo);
      curTransfo = brightMatch->Transfo()->Clone();
      delete brightMatch;
      brightMatch = ListMatchCollect(L1, L2, curTransfo, brightDist);
    } while (brightMatch->size() > nstarmin && transDiff > 0.05 && ++iter < 5);
    
    double prevChi2 = curChi2;
    curChi2 = ComputeChi2(*brightMatch, *curTransfo) / brightMatch->size();

    delete fullMatch;
    fullMatch = ListMatchCollect(List1, List2, curTransfo, fullDist);
    cout << " ListMatchRefine: order " << order
	 << " med.resid "  << median_distance(fullMatch, curTransfo) 
	 << " #match " << fullMatch->size()
	 << endl;
    if (((prevChi2 - curChi2) > 0.01*curChi2) && curChi2 > 0) {
      cout << " ListMatchRefine: order " << order << " was a better guess\n";
      delete transfo;
      transfo = brightMatch->Transfo()->Clone();
    }
    nstarmin = brightMatch->Transfo()->Npar();
  } while (++order <= maxOrder);

  delete brightMatch;
  delete fullMatch;

  return transfo;
}

/*
RUN_ROOTCINT
LINKDEF_CONTENT : #pragma link C++ class MatchConditions-;
LINKDEF_CONTENT : #pragma link C++ function MatchSearchRotShift(BaseStarList&, BaseStarList&, const MatchConditions&)
LINKDEF_CONTENT : #pragma link C++ function MatchSearchRotShiftFlip(BaseStarList&, BaseStarList&, const MatchConditions&)
LINKDEF_CONTENT : #pragma link C++ function ListMatchCollect(const BaseStarList&, const BaseStarList&, const Gtransfo *Guess, const double);
*/

#ifdef USE_ROOT
#include "root_dict/listmatchdict.cc"
#endif /* USE_ROOT */

