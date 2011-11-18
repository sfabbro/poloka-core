#ifndef KERNELFITTER__H
#define KERNELFITTER__H

#include "kernelfit.h"
#include "reducedimage.h"
#include "frame.h"



class KernelFitter : public KernelFit
{
 private :
  ReducedImageRef best, worst;
  bool refIsBest;
  Frame commonFrame;
  double largestSeeing;

 public :
  KernelFitter(const ReducedImageRef &Ref, const ReducedImageRef &New, const bool NoSwap = false);
  KernelFitter(const string& FileName);

  int DoTheFit();
  bool RefIsBest() const { return refIsBest;}
  ReducedImage *Best() {return best;}
  ReducedImage *Worst() {return worst;}

  Frame CommonFrame() const { return commonFrame;}
  double LargestSeeing() const { return largestSeeing;}
  void WriteKernel(bool overwrite=false) const;
  bool ReadKernel();
};




#endif /* KERNELFITTER__H */
