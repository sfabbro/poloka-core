#include <iostream>
#include <vector>

#include "fitsimage.h"
#include "imagesum.h"

static void usage(const char *progName)
{
  cerr << "Usage: " << progName << " [OPTIONS] <DbImage(s)>\n"
       << "  register, resample and stack images with a given method into dbimage named sum. \n"
       << "   OPTIONS are \n"
       << "     -stack  <n> specify the stacking method among: (default is 1)\n"
       << "              1  Weighted average\n"
       << "              2  Clipped weighted average\n"
       << "              3  Median\n"
       << "     -weight <n> specify the weight method among: (default is 2)\n"
       << "              1  Point source optimal (inverse sky variance and seeing)\n"
       << "              2  Extended source optimal (inverse sky variance)\n"
       << "              3  No global weighting  (just local weighting)\n"
       << "              4  No weights at all\n"
       << "     -sum <name> give a name to the stacked dbimage (default is sum)\n"
       << "     -geo <name> specify the geometric reference name (default is first image)\n"
       << "     -pho <name> specify the photometric reference name (default is first image)\n"
       << "     -j : just stack, do not register and resample\n";
  exit(-1);
}

int main(int nargs, char **args)
{
  // if nothing is given
  if (nargs < 2){usage(args[0]);}
  if (nargs == 2){cerr << " Sum at least 2 images !!\n\n"; usage(args[0]);}

  // take care of arguments
  bool justSum = false;
  StackingMethod stackMethod = WeightedAverage;
  WeightingMethod weightMethod = ExtendedSourceOptimal;
  string sumName("sum"), geoName("NOGEO"), phoName("NOPHO");

  ReducedImageList toSum;
  for (int i=1; i<nargs; ++i)
    {
      char *arg = args[i];
      // images
      if (arg[0] != '-') 
	{
	  ReducedImage *im = new ReducedImage(arg);
	  toSum.push_back(im);
	  continue;
	}

      // options
      arg++;
      sscanf(arg, "%s", arg);
      if (strcmp(arg,"stack")==0)   { ++i; stackMethod =  (StackingMethod) atoi(args[i]); continue;}
      if (strcmp(arg,"weight")==0) { ++i; weightMethod = (WeightingMethod) atoi(args[i]); continue;}
      if (strcmp(arg,"geo")==0) { ++i; geoName = args[i]; continue;}
      if (strcmp(arg,"sum")==0) { ++i; sumName = args[i]; continue;}
      if (strcmp(arg,"pho")==0) { ++i; phoName = args[i]; continue;}
      if (strcmp(arg,"j")==0) { justSum = true; continue;}

      // unrecognized option
      usage(args[0]);      
    }

  if (geoName == "NOGEO") geoName = toSum.front()->Name();
  ReducedImage geoRef(geoName);

  if (phoName == "NOPHO") phoName = toSum.front()->Name();
  ReducedImage phoRef(phoName);


  cout << " Will stack " << toSum.size() << " images with:\n"
       << "     - stack name      : " << sumName << endl
       << "     - geometric ref   : " << geoName << endl
       << "     - photometric ref : " << phoName << endl
       << "     - stack method    : " << name_of_stackingMethod(stackMethod) << endl
       << "     - weight method   : " << name_of_weightingMethod(weightMethod) << endl;


  if (!justSum) ImagesAlignAndSum(toSum, geoRef, sumName, 
				  DoFits | DoCatalog | DoWeight | DoSatur, 
				  &phoRef, weightMethod, stackMethod);
  else
    {
      ImageSum sum(sumName, toSum, &phoRef, weightMethod, stackMethod);
      sum.Execute(DoFits | DoCatalog | DoWeight | DoSatur);
    }
  
  return EXIT_SUCCESS;
}

