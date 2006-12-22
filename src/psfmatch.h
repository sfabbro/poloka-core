// This may look like C code, but it is really -*- C++ -*-
#ifndef PSFMATCH__H
#define PSFMATCH__H

#include "reducedimage.h"
#include "frame.h"
#include "basestar.h"
#include "kernelfit.h"


class Kernel;
class Stamp;


//! A class that wraps calls to KernelFit. Used both to carry out subtractions (ImageSubtraction) and just kernel fitting (for the light curve).

class PsfMatch {
private:
  bool ref_is_best;
  ReducedImageRef best,worst;
  string refName, newName;
  Frame intersection;
  double seeing,sigmaBack,photomRatio;
  KernelFitRef fit;

public:
  string NotFilteredStarListName();

  //!
  PsfMatch(const ReducedImageRef Ref, const ReducedImageRef New, const PsfMatch *APreviousMatch = NULL, bool noswap=false);
  ~PsfMatch();
  PsfMatch();

  BaseStarList objectsUsedToFit;
  //! Carry out kernel fit. argument to enable to keep the images
  bool FitKernel(const bool KeepImages = false);
  int FilterStarList(const double MaxDist=1);
  double PhotomRatio() const {return photomRatio;}
  double Chi2() const;
  int Nstars() const;
  int Nparams() const;
  double SigmaBack() const {return sigmaBack;}
  Frame Intersection() const {return intersection;}
  void KernelToWorst(Kernel &Result, const double &x, const double &y) const;
  void BackKernel(Kernel &Diffback,const double &xc, const double &yc);
  bool VarianceConvolve(Image &Variance);
  bool FitExists() const {return ((fit == NULL)? false : true); }
  void ConvolveBest(ReducedImage &ConvImage);
  bool Subtraction(ReducedImage &RImage, bool KeepConvolvedBest = false);
  bool RefIsBest() const { return ref_is_best;}
  ReducedImage* Ref() { return ((ref_is_best)? best : worst); }
  ReducedImage* New() { return ((ref_is_best)? worst : best); }
  ReducedImage* Best() { return best;}
  ReducedImage* Worst() { return worst;}
  KernelFitRef const GetKernelFit() const {return fit;}
  void SetKernelFit(KernelFit *kernel);

  void dump(std::ostream &stream = std::cout) const{;
  stream << "refName, newName;" << refName << " " <<  newName << endl ;
  stream << "Frame " ; intersection.dump();
  stream << endl << "seeing,sigmaBack,photomRatio : " 
	 << seeing << " " << sigmaBack<< " " 
	 << photomRatio << endl ;
  stream << "KernelFit";
  fit->dump();  
  }


};

#endif //PSFMATCH__H
