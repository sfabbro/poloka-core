#include <iostream>

#include "fileutils.h"
#include "fitsimage.h"
#include "dbimage.h"
#include "superflat.h"



int main(int argc, char ** argv)
{
if (argc == 1)
  {
    cerr << " flatfield [-o (overwrite existing flatfielded images)] <dbimage(s)> " << endl;
    exit(1);
  }
bool overwrite = false;
for (int i=1;  i<argc; ++i) /* argv[0] is the executable name */
  {
  
  if (strcmp(argv[i],"-o") == 0) {overwrite = true; continue;}  
  char *name = argv[i];
  DbImage dbimage(name);

  if (!dbimage.IsValid())
    {
      cerr << " flatfield : cannot find "<< name  << endl;
      continue;
    }

  string outName = dbimage.FitsImageName(Calibrated);
  if (!overwrite && FileExists(outName))
    {
      cerr << " flatfielded image already exists for image " << name << endl;
      continue;
    }

  string inName = dbimage.FitsImageName(Raw);
  
  if (FileExists(dbimage.FitsFlatName()))
    {
      cout << " starting flatfielding of " << name << endl;
      FlatFieldImage(inName, dbimage.FitsFlatName(),
                    dbimage.FitsBiasName(), dbimage.FitsFringeName(),
                    outName);
    }
  else
    {
      Image flat = Image(); /* empty image -> no flatfielding */
      cerr << " no flat for image " << name << endl;
      cerr << " Ou alors, c'est une image qui vient de Cambridge : on modife juste le header ! " << endl;
      ImageAlreadyFlatFielded(inName, flat, outName);      
    }
  }
}
