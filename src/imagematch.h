// This may look like C code, but it is really -*- C++ -*-
#ifndef IMAGEMATCH__H
#define IMAGEMATCH__H


#include "basestar.h"
class FitsHeader;
class DbImage;
class Gtransfo;
class GtransfoLin;


/*! \file 
    \brief This file contains wrappers to call match guessing routines.
*/

//! tries to guess the Linear geometric transformation to go from List1 to List2. 
/*! It uses
the image headers to check if there is any WCS info that may avoid the combinatorial search 
(from listmatch.cc) for the guess. If images come from different telescope/instruments, and the combinatorial search 
is to be run, a flip is allowed for. */

bool MatchGuess(const BaseStarList &List1, const BaseStarList &List2,
		const FitsHeader &Head1, const FitsHeader &Head2,
		CountedRef<Gtransfo> &One2Two, CountedRef<Gtransfo> &Two2One);

//! From an initial first order guess, fit and refine it
int RefineGuess(const BaseStarList &List1, const BaseStarList &List2, 
		CountedRef<Gtransfo> &One2Two, CountedRef<Gtransfo> &Two2One);

//! calls MatchGuess, refines the guess with RefineGuess
bool ImageListMatch(const DbImage &DbImage1, const DbImage &DbImage2, 
                    CountedRef<Gtransfo> &One2Two, CountedRef<Gtransfo> &Two2One);

#endif /* IMAGEMATCH__H */
