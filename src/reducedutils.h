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


//! align photometrically a night to a reference with weigthed mean on the star flux ratios
//! (for a more accurate photometric ratio, you should use KernelFit or PsfMatch)
//! The routine assumes by default that images are already aligned( cur2ref == 0 : no need to geom 
//! transform the lists. Provide the actual transfo (use ImageListMatch to find it) if needed.

double QuickPhotomRatio(const ReducedImage &CurImage, const ReducedImage &RefImage, 
			double &error, const Gtransfo* Cur2Ref = 0);

//! assumes that the Gtransfo is known. Computes "cur"/"ref"
double QuickPhotomRatio(const SEStarList &CurList, const SEStarList &RefList, 
			double &error, const GtransfoRef transfo);

//! quick and robust photometric ratio given a StarMatchList (s1->flux/s2->flux).
double MedianPhotomRatio(const StarMatchList &MatchList);

//! quick and robust photometric ratio between two BaseStarList's. (Cur/Ref)
double MedianPhotomRatio(const BaseStarList &CurList, const BaseStarList &RefList,  const Gtransfo *Transfo=NULL);



struct FluxPair {
  double f1, sig1, f2, sig2;
  FluxPair(const double &F1, const double &S1, 
	   const double &F2, const double &S2) : 
    f1(F1), sig1(S1), f2(F2), sig2(S2) {};
};

typedef std::list<FluxPair> FluxPairList;

//! compute R such that f1=R*f2 on average.
bool SlowPhotomRatio(const FluxPairList &L, const double NSigChi2Cut, double &R, double &Var);

//! wrapper to SlowPhotomRatio
bool PhotomRatio(const DbImage &Rim, const DbImage &Ref, double& Ratio, double &Error, GtransfoRef &Im2Ref);

inline bool PhotomRatio(const DbImage &Rim, const DbImage &Ref, double& Ratio, double &Error) {
  GtransfoRef g;
  return PhotomRatio(Rim, Ref, Ratio, Error, g);
}

inline double PhotomRatio(const DbImage &Rim, const DbImage &Ref) {
  double e, r;
  if (PhotomRatio(Rim,Ref,r,e))
    return r;
  return -1; // bad idea, but you should not use this wrapper routine anyway.
}

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
string ImageIntegerShift(const ReducedImage& Im, const ReducedImage& Ref, const GtransfoRef ImToRef=GtransfoRef());

#endif // REDUCEDUTILS__H
