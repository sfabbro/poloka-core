// This may look like C code, but it is really -*- C++ -*-
//  daophotutils.h
//
// Last change: Time-stamp: <05 Mar 2003 at 15:01:54 by Sebastien Fabbro>
//
//

#ifndef DAOPHOTUTILS__H
#define DAOPHOTUTILS__H

#include "reducedimage.h"

//! Build a DAOPHOT PSF
void MakeDaoPsf(ReducedImage &Rim, const bool Redo=false);
void MakePrecisePsf(const ReducedImage &Rim);

//! Fit a DAOPHOT PSF on SExtractor catalog
void MakeAllStar(const ReducedImage &Rim, const bool Redo=false);

//! produces a PSF and fit stars with it
void MakeDaoPsfCat(ReducedImage &Rim, const bool DoPsf=true, const bool DoCat=false, 
		   const bool Manually=false, const bool Combine=false, const bool Redo=false);

//! Build a fake stars using DAOPHOT routines and PSF
void MakeFakeStarImage(const ReducedImage &Rim, const bool Redo=false);


#endif // DAOPHOTUTILS__H
