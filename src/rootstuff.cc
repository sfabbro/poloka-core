#include <iostream>
#include "rootstuff.h"


#ifdef USE_ROOT
#include <TFile.h>
#include <TKey.h>
#endif /* USE_ROOT */

using namespace std;

int read_single_object_file(const char *FileName, TObject *obj)
{
#ifdef USE_ROOT
  TFile tfile(FileName);
  TIter nextkey(tfile.GetListOfKeys());
  TKey *key = (TKey*)nextkey();
  if (key) 
    {
      key->Read(obj);
    }
  else 
    {
      cerr << " could not read : " << FileName << endl;
      return 0;
    }
  tfile.Close(); 
  return 1;
#else
  cerr << " Error : Cannot load Root file \""<< FileName 
       << "\" without USE_ROOT defined at compilation! " << endl;
  return 0;
#endif
}

