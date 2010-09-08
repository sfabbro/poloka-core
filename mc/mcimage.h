// This may look like C code, but it is really -*- C++ -*-
#ifndef MCIMAGE__H
#define MCIMAGE__H 

#include "gtransfo.h"
#include "basestar.h"
#include "reducedimage.h"

typedef enum AddMethod{WModel =0, WGaussian, WDaoPsf};


class SimSNStarList;

class MCImage : public ReducedImage 
{
  ReducedImageRef Source ;
  GtransfoRef Transfo_fromref ;

  SimSNStarList* SNList;

  AddMethod Methode ;
  
SimSNStarList*  PrepareList(double zp=-99.0);
public :
  MCImage(string name, const ReducedImageRef source, const GtransfoRef Transfo, SimSNStarList* List, AddMethod a_method  ); 
  MCImage(string name,const ReducedImageRef source, SimSNStarList* List, AddMethod a_method  );
  MCImage(string name, const ReducedImageRef Ref, const ReducedImageRef source, SimSNStarList* List, AddMethod a_method  );

bool MakeFits();
bool MakeCatalog();
bool MakeWeight() ;
bool MakeBad() ;
bool MakeSatur(); 
};


#endif
