#include <iostream>
#include <string>
#include <cstdio>
#include "fitsimagearray.h"
#include "fileutils.h"
#include "fitstoad.h"

static void usage()
{
  cout << " usage : split_fits2 file [file ...] (input files cannot be gzipped)" << endl;
}

int main(int nargs, char **args)
{
string outDirName = "./";
int nIn;
 if (nargs <=1) {usage(); exit(-1);}

for (int i=1; i<nargs; ++i)
  {
    string fileName = args[i];
    FitsImageArray inFile(fileName);
    if (!inFile.IsValid())
      {
	cerr << " file " << fileName << " is not a valid fits file " << endl;
        usage();
        continue;
      }
    cout << " splitting " << fileName << " from " << TelInstName(inFile) << endl;
    inFile.SplitAndWrite();
  }
return 1;
}
