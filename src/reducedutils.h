// This may look like C code, but it is really -*- C++ -*-
#ifndef REDUCEDUTILS__H
#define REDUCEDUTILS__H

#include "reducedimage.h"
#include "sestar.h"

class Gtransfo;

//! a class representing a set of overlapping images of the same night, 
//! same instrument and same filter
class NightSet : public ReducedImageList {
public:
  NightSet(const ReducedImage &anImage);
  NightSet(const NightSet &Other);
  ~NightSet(){};
  //! true if belongs to the same set
  bool IsFromSameSet(const ReducedImage &Another) const;
  //! returns the generic string name of the set
  string GenericName() const;
  //! allows cout << NighSet;
  friend ostream& operator << (ostream & Stream, const NightSet& MySet) 
       {MySet.dump(Stream); return Stream;}
  void dump(ostream &Stream = cout) const;
  double Seeing() const;
  double Exposure() const;
  double SignalToNoise23() const;
};

void SelectBestImages(ReducedImageList &Images);
const ReducedImage *BuildReference(const string &RefName, 
				   const ReducedImageList &Images, 
				   const bool ToAlign=true);

//! arrange and copy a list of ReducedImage into different NightSet
void ArrangeByNight(const ReducedImageList &ImList, vector<NightSet> &AllNights);

//! simple routine to find a geometric reference (best seeing among smallest pixel sizes)
const ReducedImage* BestResolutionReference(const ReducedImageList &ImList);

//! remove elements of a ReducedImageList which don't overlap with a given reference
void FilterWithOverlap(const ReducedImageList &Images, const ReducedImage &aReference,
		       ReducedImageList &overlapImages, const double &MinOverlap=1);


//! align photometrically a night to a reference with weigthed mean on the star flux ratios
//! (for a more accurate photometric ratio, you should use KernelFit or PsfMatch)
//! The routine assumes by default that images are already aligned( cur2ref == NULL : no need to geom 
//! transform the lists. Provide the actual transfo (use ImageListMatch to find it) if needed.


double QuickPhotomRatio(const ReducedImage &CurImage, const ReducedImage &RefImage, 
			double &error, const Gtransfo* Cur2Ref = NULL);

//! assumes that the Gtransfo is known.
double QuickPhotomRatio(const SEStarList &CurList, const SEStarList &RefList, 
			double &error, const Gtransfo *transfo);


#endif // REDUCEDUTILS__H
