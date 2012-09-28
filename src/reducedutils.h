// This may look like C code, but it is really -*- C++ -*-
#ifndef REDUCEDUTILS__H
#define REDUCEDUTILS__H

#include "reducedimage.h"
#include "sestar.h"
#include "gtransfo.h"
#include "listmatch.h"

string GtransfoName(const DbImage& Ref, const DbImage& Src);
string ShiftedName(const string &ToShift, const string &Ref);

//! arrange and copy a list of ReducedImage into different sets of images
void ArrangeByInstBandNight(const ReducedImageList &ImList, vector<ReducedImageList> &ImageSets);

//! build a name such CFHT_MEGACAM_20030201_R from a ReducedImage
string ImageSetName(const ReducedImage& AnImage);

//! simple routine to find a geometric reference (best seeing among smallest pixel sizes)
const ReducedImage* BestResolutionReference(const ReducedImageList &ImList);

//! remove elements of a ReducedImageList which don't overlap with a given reference
void FilterWithOverlap(const ReducedImageList &Images, const ReducedImage &aReference,
		       ReducedImageList &overlapImages, const double MinOverlap=1);

//! returns a frame containing all images
Frame UnionFrame(const ReducedImageList& ImList, const ReducedImage* Reference=0);

//! produce an image with margins large enough to include all images to align
void MakeUnionRef(const ReducedImageList& ToAlign, const ReducedImage& Reference, const string& unionName);

//! loads a decent basestarlist to use with matching routines
//! preferences: 1. PSFStarList, 2. AperSEStarList, 3. SEStarList
void LoadForMatch(const ReducedImage& Im, BaseStarList& BList, const double& MinSN=15);

//! convenient wrapper to find a Gtransfo between 2 images only based on WCS information
GtransfoRef FindTransfoFromWCS(const ReducedImage& Src, const ReducedImage& Dest);

//! wrapper for lists to call matching routines
//! tries: 1. WCS composition 2. Combinatorics then refines
GtransfoRef FindTransfo(const BaseStarList& SrcList, const BaseStarList& DestList,
			const ReducedImage& Src, const ReducedImage& Dest);

//! convenient wrapper to find a Gtransfo between 2 images using the above routine
GtransfoRef FindTransfo(const ReducedImage& Src, const ReducedImage& Dest);

string ImageResample(const ReducedImage& Im, const ReducedImage& Ref, const GtransfoRef ImToRef=GtransfoRef(), const GtransfoRef RefToIm=GtransfoRef());
string ImageResample(const ReducedImage& Im, const GtransfoRef ImToRef, const GtransfoRef RefToIm);
string ImageIntegerShift(const ReducedImage& Im, const ReducedImage& Ref, const GtransfoRef ImToRef=GtransfoRef());

#endif // REDUCEDUTILS__H
