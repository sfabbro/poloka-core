#include <iostream>
#include <vector>

#include "imagesubtraction.h"
#include "fitsimage.h"

static void usage(const char *progName)
{
  cerr << "Usage: " << progName << " [OPTIONS] <DbImage(s)> " << endl;
  cerr << "  match PSF of <DbImage(s)> with a reference image \n" 
       << "   OPTIONS:\n"
       << "    -ref <name>: indicate the reference image with <name>. Default is first one. \n"
       << "    -sub  : subtract and perform detection. \n"
       << "    -conv : convolve the best seeing image of the two. \n"
       << "   Note: All the images have to be on the same pixel grid ! \n";

  exit(-1);
}

int main(int nargs, char **args)
{
  // if nothing is given
  if (nargs < 2){usage(args[0]);}
  if (nargs == 2){cerr << " PSF match at least 2 images !!\n\n"; usage(args[0]);}
  
  // take care of arguments
  string refName("NOREF");
  bool dosub = false, doconv=false;
  
  ReducedImageList toMatch;

  for (int i=1; i<nargs; ++i)
    {
      char *arg = args[i];
      // images
      if (arg[0] != '-') 
	{
	  ReducedImage *im = new ReducedImage(arg);
	  toMatch.push_back(im);
	  continue;
	}

      // options
      arg++;
      if (strcmp(arg,"ref")==0) { ++i; refName = args[i]; continue;}
      if (strcmp(arg,"sub")==0) { dosub=true; continue;}
      if (strcmp(arg,"conv")==0) { doconv=true; continue;}

      // unrecognized option
      usage(args[0]);
    }

  if (refName == "NOREF") refName = toMatch.front()->Name();
  ReducedImage ref(refName);
  for (ReducedImageCIterator it = toMatch.begin(); it != toMatch.end(); ++it)
    {
      const ReducedImage *current = *it;
      PsfMatch match(ref, *current);
      match.FitKernel(true);      
      if (doconv) 
	{
	  string convname = "C_"+match.Best()->Name()+match.Worst()->Name();
	  ReducedImage conv(convname);
	  match.ConvolveBest(conv);
	}
      if (dosub) 
	{
	  string subname = SubtractedName(refName, current->Name());
	  ImageSubtraction sub(subname, ref, *current, &match);
	  sub.Execute(DoFits | DoWeight | DoCatalog);
	}
    }
 
  return EXIT_SUCCESS;
}

