#ifndef LISTMATCH__H
#define LISTMATCH__H

#include <string>

#include <poloka/basestar.h>
class Gtransfo; class GtransfoLin;
#include <poloka/starmatch.h>


struct MatchConditions 
{
  int NStarsL1, NStarsL2;
  int MaxStarsL1, MaxStarsL2;
  int MaxTrialCount;
  double NSigmas;
  double MaxShiftX, MaxShiftY;
  double SizeRatio, DeltaSizeRatio, MinMatchRatio;
  double Dist, MaxDist;
  int PrintLevel;
  int Algorithm;
  int MaxOrder;
  
  MatchConditions(const std::string &DatacardsName = "");
  void read(const std::string &DatacardsName = "");
  void init();

  double MinSizeRatio() const { return SizeRatio - DeltaSizeRatio;}  
  double MaxSizeRatio() const { return SizeRatio + DeltaSizeRatio;}  

} ;


/*! \file
    \brief Combinatorial searches for linear transformations to go from
           list1 to list2.

    The following routines search a geometrical transformation that make
two lists of stars to match geometrically as well as possible. They are used
either to match two images of the same sky area, or an image with a catalogue.
They assume that fluxes assigned to stars are actual fluxes, i.e. the brighter
the star, the higher the flux. They however only rely on flux ordering, 
not values.
 */


//! searches a geometrical transformation that goes from List1 to List2.
/*!  The found transformation is a field of the returned object, as well as the star pairs
(the matches) that were constructed.  (see StarMatchList class definition for more details).
The various cuts are contained in Conditions (see listmatch.h) for its contents.
This routine searches a transformation that involves a shift and a rotation. */

StarMatchList* MatchSearchRotShift(BaseStarList &L1, BaseStarList &L2, const MatchConditions &Conditions);

//!same as above but searches also a flipped solution.

StarMatchList* MatchSearchRotShiftFlip(BaseStarList &L1, BaseStarList &L2, const MatchConditions &Conditions);

//! assembles star matches.
/*! It picks stars in L1, transforms them through Guess, and collects
closest star in L2, and builds a match if closer than MaxDist). */

StarMatchList *ListMatchCollect(const BaseStarList &L1, const BaseStarList &L2,const Gtransfo *Guess, const double MaxDist);

//! same as before except that the transfo is the identity

StarMatchList *ListMatchCollect(const BaseStarList &L1, const BaseStarList &L2, const double MaxDist);

//! searches for a 2 dimensional shift using a very crude histogram method.



GtransfoLin *ListMatchupShift(const BaseStarList &L1, 
			      const BaseStarList &L2, 
			      const Gtransfo &Tin, 
			      double MaxShift, double BinSize = 0);

//! computes number of matches in a StarMatchList within 3 sigmas of chi2
size_t StarMatchCheck(const StarMatchList* Match);

//! wrapper around all possible combinatorial matching routines above
GtransfoRef ListMatchCombinatorial(const BaseStarList &List1,
				   const BaseStarList &List2,
				   const MatchConditions& Conditions=MatchConditions());

//! refine (fit) a Gtransfo up to maxOrder polynomial. maxOrder=0 not doing anything
GtransfoRef ListMatchRefine(const BaseStarList& List1,
			    const BaseStarList& List2,
			    GtransfoRef transfo, 
			    const MatchConditions& Conditions=MatchConditions());

//! wrapper of the two above routines
GtransfoRef ListMatch(const BaseStarList& List1,
		      const BaseStarList& List2,
		      MatchConditions Conditions=MatchConditions());

#endif /* LISTMATCH__H */
