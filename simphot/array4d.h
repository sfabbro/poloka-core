#ifndef ARRAY4D__H
#define ARRAY4D__H

#include "array2d.h"


typedef float CoeffType;


typedef Array2D<CoeffType> CoeffBlock;



class Array4D : public Array2D<CoeffBlock>
{
  public :
  Array4D(const int Xmin=0, const int Ymin=0, const int Xmax=0, const int Ymax=0):
  Array2D<CoeffBlock>(Xmin, Ymin, Xmax, Ymax) {};

  Array4D(const IntFrame &Frame) : Array2D<CoeffBlock>(Frame) {};

  IntFrame InternalFrame() const;

    void Transpose(Array4D &Result) const ;

    void Convolve(const Array2D<double> &Kernel, Array4D &Result) const;

    size_t MemSize() const;

    // ! to call under gdb
    const CoeffBlock& GetElement(const int i, const int j) const;
    
};



#endif /* ARRAY4D__H */
