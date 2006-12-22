#include <cstdio>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "fileutils.h"
#include "fitsimage.h"
#include "superflat.h"
#include "fitsset.h"



static void usage(const char *pg_name, const char * message = NULL)
{
  if (message) cerr << message << endl;
  cerr << pg_name ;
  cerr << " \t-l <file list> : Build a master flat frame from <file list>" << endl;
  cerr << "\t\t -b <file name> : bias name to use to subtract bias. no such arg implies no bias subtraction" << endl;
  cerr << "\t\t -o <output fits> " << endl;
  cerr << "\t\t -Bias : do a MasterBias from <file list>" << endl;
  cerr << "\t\t - s <skyflat name> : to make a Master in 2 steps and then a better tagging of dead pixels " << endl;
  exit(-1);
}


int main(int argc, char**argv)
{
char *fileList = NULL;
char *biasName = NULL;
char *outName = NULL;
char *skyName = NULL;
int masterBias = 0;

for (int i=1; i< argc; i++)
  {
  char *arg = argv[i];
  if (arg[0] != '-') argc = 1;
  else 
  switch (arg[1])
    {
    case 'l' : i++ ; fileList = argv[i]; break;
    case 'b' : i++ ; biasName = argv[i]; break;
    case 'o' : i++ ; outName = argv[i] ; break;
    case 's' : i++ ; skyName = argv[i] ; break;
    case 'B' :
      if (arg[2] == 'i' && arg[3] == 'a' && arg[4] == 's')
      {
	masterBias = 1;
	break;
      }
      else argc = 1;
    default : argc = 1;
    }
  }

if (argc <= 1) usage (argv[0]);
if (!fileList)  usage(argv[0], " no file list given : abort " );
if (!outName) usage(argv[0], " provide an output file name !");

if (masterBias == 1)
  {
    FitsSet fitsFileSet(fileList, false);
    Image *MasterBias = MakeRawMedian(fitsFileSet);
    FitsHeader first(fitsFileSet[0]);
    FitsImage masterBiasFits(outName, first, *MasterBias);
    delete MasterBias;
    masterBiasFits.AddOrModKey("OBJECT","MasterBias","see history for more informations");
    masterBiasFits.AddCommentLine("masterbias done using " + fitsFileSet.AllNames());
    
    return EXIT_SUCCESS;
  }


FitsSet fitsFileSet(fileList);

Image *meanBias = NULL;
if (biasName != NULL)
  {
  meanBias = new FitsImage(biasName);
  }

Image *skyFlat = NULL;
if (skyName != NULL)
  {
  skyFlat = new FitsImage(skyName);
  }

Image *Flat = MakeSuperFlat(fitsFileSet, meanBias, skyFlat);
if (meanBias) delete meanBias;
if (skyFlat) delete skyFlat;

if (!Flat) return 0;

Image *DeadPixels = DeadPixelImageAndFlatSmoothing(*Flat,0.6,1.4);

// write them (and release memory)
string deadName = AddPrefix("dead_",outName);


FitsHeader first(fitsFileSet[0]);
FitsImage flatFits(outName, first, *Flat);
delete Flat;
FitsImage deadFits(deadName, first, *DeadPixels);
deadFits.ModKey("BITPIX",8);
delete DeadPixels;

flatFits.AddOrModKey("OBJECT","MasterFlat","see history for more informations");
deadFits.AddOrModKey("OBJECT","Dead Pix Mask","see history for more informations");
flatFits.AddCommentLine("masterflat done using " + fitsFileSet.AllNames());
deadFits.AddCommentLine("dead pixels mask deduced from " + string(outName));

return EXIT_SUCCESS;
}
