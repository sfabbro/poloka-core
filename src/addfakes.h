#ifndef ADDFAKES__H
#define ADDFAKES__H

#include "reducedimage.h"
class SENearStarList;


void 
MakeListSn( const ReducedImage *AnImage, ReducedImage *New, SENearStarList &NearStarRef, bool AssociateGal=false);

#endif
