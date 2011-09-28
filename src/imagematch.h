// This may look like C code, but it is really -*- C++ -*-
#ifndef IMAGEMATCH__H
#define IMAGEMATCH__H


#include "sestar.h"
class FitsHeader;
class DbImage;
class Gtransfo;
class GtransfoLin;


/*! \file 
    \brief This file contains wrappers to call match guessing routines.
*/

//! checks that necessary files exists, loads the catalog and filters it.
bool ListAndFitsCheckForMatch(const DbImage &Im, SEStarList &Sl);

//! tries to guess the Linear geometric transformation to go from List1 to List2. 
/*! It uses
the image headers to check if there is any WCS info that may avoid the combinatorial search 
(from listmatch.cc) for the guess. If images come from different telescope/instruments, and the combinatorial search 
is to be run, a flip is allowed for. */

bool MatchGuess(const BaseStarList &List1, const BaseStarList &List2,
		const FitsHeader &Head1, const FitsHeader &Head2,
		CountedRef<Gtransfo> &One2Two, CountedRef<Gtransfo> &Two2One, float min_match_ratio=0);


//! From an initial first order guess, fit and refine it
int RefineGuess(const BaseStarList &List1, const BaseStarList &List2, 
		CountedRef<Gtransfo> &One2Two, CountedRef<Gtransfo> &Two2One,
		const string& image1_name="", const string& image2_name="", int max_order=3); //names are just for dumping results

//! calls MatchGuess, refines the guess with RefineGuess
bool ImageListMatch(const DbImage &DbImage1, const DbImage &DbImage2, 
                    CountedRef<Gtransfo> &One2Two, CountedRef<Gtransfo> &Two2One, float min_match_ratio=0, int max_order=3,int N_CUT_LIST=500, int N_FRAME_LIST=300 );

//! same as above but with already loaded catalogs. Use ListAndFitsCheckForMatch to lad and filter the catalog.
/*! this routine was introduced to avoid reloading one of the catalogs when matching many images to the same reference */ 
bool ImageListMatch(const DbImage &DbImage1, const SEStarList& SL1,
		    const DbImage &DbImage2, const SEStarList& SL2,
		    CountedRef<Gtransfo> &One2Two, 
		    CountedRef<Gtransfo>  &Two2One, float min_match_ratio=0, int max_order=3,int N_CUT_LIST=500, int N_FRAME_LIST=300 );
#endif /* IMAGEMATCH__H */
