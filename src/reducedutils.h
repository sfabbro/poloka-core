// This may look like C code, but it is really -*- C++ -*-
#ifndef REDUCEDUTILS__H
#define REDUCEDUTILS__H

#include "reducedimage.h"
#include "sestar.h"

class Gtransfo;

//! arrange and copy a list of ReducedImage into different sets of images
void ArrangeByInstBandNight(const ReducedImageList &ImList, vector<ReducedImageList> &ImageSets);

//! build a name such CFHT_MEGACAM_20030201_R from a ReducedImage
string ImageSetName(const ReducedImage& AnImage);

//! simple routine to find a geometric reference (best seeing among smallest pixel sizes)
const ReducedImage* BestResolutionReference(const ReducedImageList &ImList);

//! remove elements of a ReducedImageList which don't overlap with a given reference
void FilterWithOverlap(const ReducedImageList &Images, const ReducedImage &aReference,
		       ReducedImageList &overlapImages, const double MinOverlap=1);


//! align photometrically a night to a reference with weigthed mean on the star flux ratios
//! (for a more accurate photometric ratio, you should use KernelFit or PsfMatch)
//! The routine assumes by default that images are already aligned( cur2ref == 0 : no need to geom 
//! transform the lists. Provide the actual transfo (use ImageListMatch to find it) if needed.

double QuickPhotomRatio(const ReducedImage &CurImage, const ReducedImage &RefImage, 
			double &error, const Gtransfo* Cur2Ref = 0);

//! assumes that the Gtransfo is known.
double QuickPhotomRatio(const SEStarList &CurList, const SEStarList &RefList, 
			double &error, const Gtransfo *transfo);

//! quick and robust photometric ratio between two BaseStarList's
double MedianPhotomRatio(const BaseStarList &CurList, const BaseStarList &RefList,  const Gtransfo *Transfo);

//! quick and robust photometric ratio between two BaseStarList's
double MedianPhotomRatio(const BaseStarList &CurList, const BaseStarList &RefList);

#endif // REDUCEDUTILS__H
