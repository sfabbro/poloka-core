#ifndef ROOTSTUFF__H
#define ROOTSTUFF__H

#ifdef __CINT__
#define USE_ROOT
#endif /*__CINT__ */

class TObject;
//! a routine that loads the first object of a file
int read_single_object_file(const char *FileName, TObject *obj);

#ifdef USE_ROOT
#include <TObject.h>


#else
// undefine Root macros
#define ClassDef(a,b)
#define ClassImp(a)
#define ClassDefT(a,b)
#define ClassDefT2(a,b)
#endif



#endif /* ROOTSTUFF__H */
