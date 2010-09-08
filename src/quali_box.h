// This may look like C code, but it is really -*- C++ -*-
#ifndef QUALI_BOX__H
#define QUALI_BOX__H

#include <iostream>

#include "gtransfo.h"
#include "fitsimage.h"
#include "image.h"
#include "starmatch.h"
#include "sestar.h"
#include "align_box.h"
#include "listmatch.h"



void Test_Qual_N_same(int N, const FitsHeader & header,  
		 SEStarList & stl1, SEStarList & stl2, 
		 Gtransfo *tf, const Image **pimage, double *Fond,double * SigFond, 
		 int dtaille, ostream & pr);
//   j'ai vire le dernier argument parce que passe sans reference, il sert a rien....
//		 int dtaille, ostream & pr, StarMatchList *lmatchphot);

void
Test_Qual_N_same(int N, StarMatchList & liste, const Image **pimage,
		 double *Fond, double * SigFond, int dtaille, ostream & pr);


#endif /* QUALI_BOX__H */
