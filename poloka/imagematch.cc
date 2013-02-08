#include <iostream>
#include <iomanip>

// CHANGE MATCH TO KEEP
#include "fitsimage.h"
#include "wcsutils.h"
#include "dbimage.h"
#include "listmatch.h"
#include "starmatch.h"
#include "sestar.h"
#include "fileutils.h"
#include "fitstoad.h"
#include "imagematch.h"
#include "frame.h"

/* TO DO:
   1 - We should really weight the match. Right now, we give equal weight to faint stars,
   bringing sometimes bad transfos when little stars are in common, or bright stars are saturated. 
   Weight should be present, but we should make sure there's a star on every corner of matched 
   common area. 
   2 - If both images have a match to the usno catalog, we could provide a transfo 
   by identifying the same USNO stars in the images, since we keep the image-catalog/
   Usno-catalog when we match to usno.
   3 - Do a likelihood ratio test to see whether going higher order in the fit is significant
*/

// remove bad zones of a list: not in bad zone, keep satur
static void KeepGoodForMatch(SEStarList &List) 
{
  for (SEStarIterator it= List.begin(); it!=List.end(); )
    {
      SEStar *pstar =  *it ;
      if (!(isnan(pstar->vx)) && !(isnan(pstar->vy)) 
	  && (pstar->FlagBad() == 0 && pstar->Flag()<8 )  // not in bad zone, keep satur
	  && (pstar->flux/pstar->EFlux() >10)
	  ) ++it ;
      else it = List.erase(it); // with counted references no need to delete;
    }  
}

// utility to check existence of list and fits image and cut the list for a max of 500 brightest objects
bool ListAndFitsCheckForMatch(const DbImage &Im, SEStarList &Sl)
{
  string listName = Im.AperCatalogName();
  
  if (!FileExists(listName))
    {
      cerr << " ListAndFitsCheck : WARNING no aperse catalog for " << Im.Name() << ", using sextractor" << endl;
      listName = Im.ImageCatalogName();
      if (!FileExists(listName))
	{
	  cerr << " ListAndFitsCheck : FAILURE no sextractor catalog for " << Im.Name() << endl;
	  return false;
	}
    }
  cout << " ListAndFitsCheck : MESSAGE, using catalog " << listName << endl;
  Sl.read(listName);
  KeepGoodForMatch(Sl);
  if (Sl.size() <= 3)
    {
      cerr << " ListAndFitsCheck : FAILURE not enough good stars for " << Im.Name() << endl;
      return false;
    }

  string FitsName = Im.FitsImageName(Calibrated);
  if (!FileExists(FitsName))
    {
      cerr << " ListAndFitsCheck : FAILURE no calibrated image for " <<  Im.Name() << endl;
      return false;
    }
  return true;
}

 
void ShiftGuess(const BaseStarList &List1, const BaseStarList &List2, CountedRef<Gtransfo> &Guess)
{
  cout << " ShiftGuess : doing a quick shift guess " << endl;
  GtransfoLin *shift = ListMatchupShift(List1, List2, *Guess, 400.);  // 400 = maxshift.
  GtransfoLin newtransfo = (*shift) * (dynamic_cast<GtransfoLin&>(*Guess));
  Guess = dynamic_cast<GtransfoLin*>(newtransfo.Clone());
  delete shift;
}

void BrutalGuess(const BaseStarList &List1, const BaseStarList &List2,  CountedRef<Gtransfo> &Guess, 
		 const FitsHeader &Head1, const FitsHeader &Head2, int N_STAR_L)
{
  cout << " BrutalGuess : starting a combinatorial match " << endl;
  double pixscale1 = Head1.KeyVal("TOADPIXS");
  double pixscale2 = Head2.KeyVal("TOADPIXS");
  double pixSizeRatio = pixscale1/pixscale2;

  MatchConditions conditions;
  conditions.SizeRatio = pixSizeRatio;
  conditions.DeltaSizeRatio = 0.1 * conditions.SizeRatio;
  conditions.NStarsL1 =  N_STAR_L ;
  conditions.NStarsL2 =  N_STAR_L ;
  conditions.PrintLevel = 0;
  BaseStarList L1, L2;
  List1.CopyTo(L1); List2.CopyTo(L2);

  StarMatchList *matchingList = 0;
  if (TelInstName(Head1) == TelInstName(Head2)) 
    matchingList = MatchSearchRotShift(L1, L2, conditions);
  else 
    matchingList = MatchSearchRotShiftFlip(L1, L2, conditions);

  Guess = dynamic_cast<GtransfoLin*>( matchingList->Transfo()->Clone());
  delete matchingList; 
}

static bool testNewTransfo(const StarMatchList* matchList, const GtransfoLin* guess, 
			   const double pixSizeRatio2, size_t &nmatch, const size_t nmin)
{

  size_t nmatch_new = matchList->size();  
  
  if ((fabs(fabs(guess->Determinant())-pixSizeRatio2)/pixSizeRatio2 < 0.2)  
      && (nmatch_new > nmin)) return true;
  if (nmatch_new > nmatch) {
    nmatch = nmatch_new;
  }
  //cout << "nmatch_new, nmin = " <<  nmatch_new << " " << nmin << endl;
  //cout << "size ratio " << fabs(fabs(guess->Determinant())-pixSizeRatio2)/pixSizeRatio2 << endl;
  cout << " MatchGuess : guess was bad: \n";
  matchList->DumpTransfo();
  
  return false;
}


static void CutList(const SEStarList &SList, BaseStarList &BList, int N_CUT_LIST)
{
  BList.ClearList();
  for (SEStarCIterator it=SList.begin(); it != SList.end(); ++it)
    BList.push_back(new BaseStar(**it));
  BList.FluxSort();
  BList.CutTail(N_CUT_LIST);
}

bool MatchGuess(const BaseStarList &List1, const BaseStarList &List2,
		const FitsHeader &Head1, const FitsHeader &Head2,
		CountedRef<Gtransfo> &One2Two, CountedRef<Gtransfo> &Two2One, float min_match_ratio)
{

  cout << " MatchGuess : starting guess with "  
    << List1.size() << " " << List2.size() << " bright objects\n";
  
  // tolerance used to match stars in pixels through all routines below
  const double init_toldist = 4; 
  

  One2Two = new GtransfoLin();// the identity by default

  double pixscale1 = Head1.KeyVal("TOADPIXS");
  double pixscale2 = Head2.KeyVal("TOADPIXS");

  int Nx = Head1.KeyVal("NAXIS1");
  int Ny = Head1.KeyVal("NAXIS2");
  int N_STAR_L = 70 ;
  // patch for big frames
  if (Nx > 10000 && Ny > 10000)
    N_STAR_L = 2000 ;


  double pixSizeRatio2 = pixscale1*pixscale1/pixscale2/pixscale2;

  StarMatchList *matchList = ListMatchCollect(List1, List2, One2Two, init_toldist);
  size_t nmatch = matchList->size();
  size_t nmin = min(List1.size(), List2.size())/3;
  if(min_match_ratio>0)
    nmin = size_t(min(List1.size(), List2.size())*min_match_ratio);
  if(nmin<10) {
    cout << "MatchGuess, warning : request nmin=" << nmin << " (from min_match_ratio=" << min_match_ratio << ") but use nmin=10" << endl;
    nmin = 10;
  }
  cout << "MatchGuess, minimum allowed number of matches = " << nmin << endl;
  

  bool debug_match = true;
  
  
  if(debug_match) cout << " 1 - Try WCS " << endl;
  cout << " MatchGuess : doing a quick guess with WCS" << endl;
  if (!WCSTransfoBetweenHeader(Head1, Head2, dynamic_cast<GtransfoLin&>(*One2Two))) 
    cout << " MatchGuess : one of images does not have WCS \n";
  delete matchList;
  matchList = ListMatchCollect(List1, List2, One2Two , init_toldist);
  if (!testNewTransfo(matchList,One2Two , pixSizeRatio2, nmatch, nmin))
    {
      if(debug_match) cout << " 2 - WCS failed: try a quick shift with this WCS" << endl; 
      ShiftGuess(List1, List2, One2Two);
      delete matchList;
      matchList = ListMatchCollect(List1, List2, One2Two, init_toldist);
      if (!testNewTransfo(matchList, One2Two, pixSizeRatio2, nmatch, nmin))
	{
	  if(debug_match) cout << " 3 - Both above failed: try brutal combinatorial approach" << endl; 
	  BrutalGuess(List1, List2, One2Two, Head1, Head2, N_STAR_L);
	  delete matchList;
	  matchList = ListMatchCollect(List1, List2,One2Two , init_toldist);
	  if (!testNewTransfo(matchList, One2Two, pixSizeRatio2, nmatch, nmin))
	    {
	      cout << " MatchGuess : bad initial brutal match \n";
	      delete matchList;
	      return false;
	    }
	  else cout << " MatchGuess : found a guess with BrutalGuess \n";
	}
      else cout << " MatchGuess : found a guess with ShiftGuess \n";
    }
  else cout << " MatchGuess : found a guess with WCS \n";
  
  cout << " MatchGuess : guessed transfo : \n";
  cout << *(matchList->Transfo()) << endl << " list size : " << matchList->size() << endl;

  // no need to delete with countedref
  // this is included in countedref operator =
  //if (One2Two) delete One2Two;
  //if (Two2One) delete Two2One;
 
  Two2One = (dynamic_cast<GtransfoLin&>(*One2Two)).invert().Clone(); // clone need to get a new gtransfo
  delete matchList; 
  //delete guess; // we don't use clone so no delete here

  return true;
}


// utility to check current transfo difference
static double transfo_diff(const BaseStarList &List, const Gtransfo *T1, const Gtransfo*T2)
{
  double diff2 = 0;
  FatPoint tf1;
  Point tf2;
  int count = 0;
  for (BaseStarCIterator si = List.begin(); si != List.end(); ++si)
    {
      const BaseStar *s = *si;
      T1->TransformPosAndErrors(*s, tf1);
      T2->apply(*s, tf2);
      double dx = tf1.x-tf2.x;
      double dy = tf1.y-tf2.y;
      diff2 += (tf1.vy*dx*dx + tf1.vx*dy*dy - 2*tf1.vxy*dx*dy)/(tf1.vx*tf1.vy-tf1.vxy*tf1.vxy);
      count ++;
    }
  if (count) return diff2/double(count);
  return 0;
}

// iterative collect / fit until transfo stabilizes.
// increase order if transfo looks ok
int RefineGuess(const BaseStarList &List1, const BaseStarList &List2, 
		CountedRef<Gtransfo> &One2Two, CountedRef<Gtransfo> &Two2One,
		const string& image1_name, const string& image2_name  /* names are just for dumping results (facultative)*/, 
		int max_order 
		)
{
  cout << " RefineGuess : starting guess with "  
    << List1.size() << " " << List2.size() << " bright objects\n";

  double refine_toldist = 4;
  int bestorder = 1;
  

  cout << " RefineGuess : Will try until order " << max_order << endl;
  int order = 1;
  double prevChi2Red;
  size_t minToMatch = 3;

  Gtransfo *usedToCollect = One2Two->Clone();
  StarMatchList *refinedMatch = ListMatchCollect(List1, List2, usedToCollect, refine_toldist);  
  double curChi2Red = ComputeChi2(*refinedMatch, *usedToCollect) / refinedMatch->Dof();
  
  int oldprec = cout.precision();
  cout << setprecision(3);
  cout << " RefineGuess : Order " << order 
       << " Resid " << FitResidual(*refinedMatch, *usedToCollect) 
       << " Nused " << refinedMatch->size() 
       << " Chi2/DOF " << curChi2Red << endl;
  cout << " RefineGuess : Start looping \n";

  int bestNused=0;
  double bestChi2PDF=0.;
  double bestResid=0.;
  do // loop on transfo order
    {
      prevChi2Red = curChi2Red;
      int loops = 0;
      double change;

      do // loop on transfo change
	{ 
	  if (refinedMatch) delete refinedMatch;
	  refinedMatch = ListMatchCollect(List1, List2, usedToCollect, refine_toldist);
	  refinedMatch->SetTransfoOrder(order);
	  refinedMatch->RefineTransfo(3.);
	  change = transfo_diff(List1, refinedMatch->Transfo(), usedToCollect);
	  ++loops;
	  int dof = refinedMatch->Dof();
	  if (dof == 0)
	    {
	      curChi2Red = -1;
	      break;
	    }
	  curChi2Red = refinedMatch->Chi2() / dof;
	  cout << " RefineGuess : Order " << order 
	       << " Resid " << refinedMatch->Residual() 
	       << " Npair " << refinedMatch->size() 
	       << " Chi2/DOF " << curChi2Red << endl;
	  usedToCollect = refinedMatch->Transfo()->Clone(); 
	} 
      while (change > 0.05 && loops < 5);

      // the short version of likelihood ratio test
      if (((prevChi2Red-curChi2Red) > 0.01*curChi2Red) && curChi2Red > 0 )
	{
	  cout << " RefineGuess : Order " << order << " was a better guess \n";
	  One2Two = refinedMatch->Transfo()->Clone();
	  Two2One = refinedMatch->InverseTransfo();
	  bestorder = order;
	  bestChi2PDF = curChi2Red;
	  bestNused = refinedMatch->size();
	  bestResid = refinedMatch->Residual();
	}
      //else break;
      minToMatch = usedToCollect->Npar();
      order++;
      // refine_toldist -= 0.5;
    } 
  while (refinedMatch->size() > minToMatch && order <= max_order);
  cout << setprecision(oldprec);
  // this dump is usefull to check if everything is ok in the log
  cout << " RefineGuess_SUMMARY_image1_image2_order_resid_nused_chi2pdf " 
       << image1_name << " "
       << image2_name << " "
       << bestorder << " " 
       << bestResid << " " 
       << bestNused << " "
       << bestChi2PDF << endl;
  
  delete refinedMatch;
  delete usedToCollect;
  return bestorder;
}


static void FrameList(const SEStarList &SList, const FitsHeader &Head, 
		      const Gtransfo *SListToHead, BaseStarList &BList, int N_FRAME_LIST )
{//HACK
  if (N_FRAME_LIST<=0 ) N_FRAME_LIST = 300 ;
  Frame frame(Head);
  BList.ClearList();
  for (SEStarCIterator it = SList.begin(); it != SList.end(); ++it)
    {
      if (frame.InFrame(SListToHead->apply(**it)))
	BList.push_back(new BaseStar(**it));
    }
  BList.FluxSort();
  BList.CutTail(N_FRAME_LIST);
}



// N_CUT_LIST, = 500, 50000 pour big image
// N_FRAME_LIST = 300, 10000 pour big image
// max_order = 3, 9 pour big image


bool ImageListMatch(const DbImage &DbImage1, const DbImage &DbImage2, 
		    CountedRef<Gtransfo> &One2Two, CountedRef<Gtransfo>  &Two2One, float min_match_ratio,  int max_order, int N_CUT_LIST, int N_FRAME_LIST)
{
  cout << "\n ImageListMatch : Matching " << DbImage1.Name() 
       << " and " << DbImage2.Name() << endl;

  One2Two = new GtransfoIdentity;
  Two2One = new GtransfoIdentity;

  if (DbImage1 == DbImage2)
    {
      cout << " ImageListMatch : same image : transfo = identity ! for " 
	   << DbImage1.Name() << endl;
      return true;
    }

  cout << " ImageListMatch : preparing lists and images \n" ;
  SEStarList sl1, sl2;

  if (!ListAndFitsCheckForMatch(DbImage1, sl1) || 
      !ListAndFitsCheckForMatch(DbImage2, sl2))
    return false;

  return ImageListMatch(DbImage1, sl1, DbImage2, sl2, One2Two, Two2One, min_match_ratio, max_order, N_CUT_LIST, N_FRAME_LIST);
}


// N_CUT_LIST, = 500, 50000 pour big image
// N_FRAME_LIST = 300, 10000 pour big image
// max_order = 3, 9 pour big image

bool ImageListMatch(const DbImage &DbImage1, const SEStarList& SL1,
		    const DbImage &DbImage2, const SEStarList& SL2,
		    CountedRef<Gtransfo> &One2Two, CountedRef<Gtransfo>  &Two2One, float min_match_ratio, int max_order, int N_CUT_LIST, int N_FRAME_LIST )
{
  string fitsName1 = DbImage1.FitsImageName(Calibrated);
  string fitsName2 = DbImage2.FitsImageName(Calibrated);
  FitsHeader head1(fitsName1), head2(fitsName2);
  BaseStarList bl1,bl2;

  cout << " ImageListMatch : doing a first guess \n" ;
  CutList(SL1, bl1, N_CUT_LIST);
  CutList(SL2, bl2, N_CUT_LIST);

  // 1 - Try a bunch of method to get an initial match
  if (!MatchGuess(bl1, bl2, head1, head2, One2Two, Two2One, min_match_ratio)) 
    {
      // in some cases the inverse transfo works better
      if (!MatchGuess(bl2, bl1, head2, head1, Two2One, One2Two, min_match_ratio)) 
	{
	  cerr << " ImageListMatch : FAILURE : Unable to get an initial match\n";
	  return false;
	}
    }
  
  // keep only objects in common transformed frame
  FrameList(SL1, head2, One2Two, bl1, N_FRAME_LIST);
  FrameList(SL2, head1, Two2One, bl2, N_FRAME_LIST);
  
  // 2 - Refine the initial match with initial full list in framed match
  cout << " ImageListMatch : refining initial guess \n" ;
  
  RefineGuess(bl1, bl2, One2Two, Two2One, DbImage1.Name(),DbImage2.Name(), max_order);
  cout << " ImageListMatch : final refined transfo : \n";
  cout << *One2Two << endl;

  // We do not necessarily have in refinedMatch the exact sample used for the fit
  if (getenv("DUMP_MATCH"))
    {
      StarMatchList *allMatch = ListMatchCollect(bl1, bl2, One2Two, 2.5);
      string name = DbImage1.Name()+"."+DbImage2.Name()+".match.list";
      cout << " writing match list " << name << endl;
      allMatch->write(name, One2Two);
    }

  return true;
}


#ifdef USE_ROOT

/*
RUN_ROOTCINT
LINKDEF_CONTENT : #pragma link C++ function MatchGuess(const BaseStarList &, const BaseStarList &,const FitsHeader &, const FitsHeader &,Gtransfo* &, Gtransfo* &);
LINKDEF_CONTENT : #pragma link C++ function  RefineGuess(const BaseStarList &, const BaseStarList &, 	Gtransfo* &, Gtransfo* &); 
LINKDEF_CONTENT : #pragma link C++ function ImageListMatch(const DbImage &, const DbImage &, Gtransfo* &, Gtransfo* &T);
*/

#include "root_dict/imagematchdict.cc"
#endif /* USE_ROOT */
