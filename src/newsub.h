// This may look like C code, but it is really -*- C++ -*-
#ifndef NEWSUB__H
#define NEWSUB__H

#include "reducedimage.h"
#include "candidatestar.h"
#include "candstar.h"
#include "yquemstar.h"
#include <list>

typedef enum TypeSub{TypeRef =0, TypeNew1, TypeNew2, TypeNew};

typedef list<ReducedImageList> ReducedNewList;
typedef list<ReducedImageList>::const_iterator ReducedNewCIterator;
typedef list<ReducedImageList>::iterator ReducedNewIterator;


//class DataCards;

struct DatSim {
  int numberOfFakes;
  double minMag, maxMag;
  
  DatSim() { numberOfFakes = 100; minMag = 22; maxMag = 26;}
  void LitDataCards(DataCards &);
  DatSim(const string &FileName);
  void Print();
};


class SENearStarList;
class ImageSum;
class ImageSubtraction;


typedef list<ImageSubtraction*> ImageSubtractionList;
typedef list<ImageSubtraction*>::const_iterator ImageSubtractionCIterator;
typedef list<ImageSubtraction*>::iterator ImageSubtractionIterator;


/* \file */

//! Handling of an actual subtraction (shift-coadd-subtract-detect). See \ref subfile for the way to drive it.
struct NewSub {
  
  ImageSubtraction *Sub;
  
  ImageSubtractionList ListOfSub;
  
  int NumberOfSub;
  ReducedImageList Ref;
  ReducedNewList ListNewList;
  //ReducedImageList New2;
  ReducedImageList AllNew;
  ReducedImageList AllImages;
  bool overwrite;
  bool onlyDet;
  bool onlyOneSub;
  SENearStarList *FakeList;

  //string DatacardsName;
  bool AddFakes;
  bool AssociateGal;
  bool FixRef;
  
  
  const ReducedImage *GeometricReference;
  //! the constructor. see \ref subfile for the syntax of the "subfile"
  NewSub(const string &SubFileName, const bool Overwrite = false, const bool OnlyDet = false);
  ~NewSub();
  int  DoIt();
  void FindGeometricReference();
  int DoOneSub(const ReducedImage *RefStack, const ImageSum *NewStack, const string &SubName, ImageSubtraction *&Sub);


  
  void MatchDetectionsWithFakes(SEStarList *Detections, const string &MatchListName);
  void ApplyCutsAndWrite(ImageSubtraction &ASubtraction);
  void Construction(ImageSubtraction &ASub, YquemStarList & stlcand);
  void ConstructMatchAndCut();
  void Cut_Write(YquemStarList & stlcand, string & cutscanname);
  void Cut_Write(CandStarList & stlcand);
  
  
  
};






#endif /* NEWSUB__H */
