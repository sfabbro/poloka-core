#include <cstdio>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include "fileutils.h"
#include "fitsimage.h"
#include "imageback.h"
#include "image.h"

int main(int argc, char **argv)
{

  int backStep = 64;
  bool outflag = false;
  bool justback = false;
  bool writeflag = false;
  
  string in_name;
  string out_name;

  // take care of options  
  for (int i=1; i< argc; i++)
    {
      char *arg = argv[i];
      if (arg[0] != '-') {in_name=argv[i];}
      else switch (arg[1])
	{
	case 'o' : i++; out_name=argv[i];outflag=true;break;
	case 's' : i++; backStep=atoi(argv[i]); break;
	case 'b' : justback=true;break;
	case 'w' : writeflag=true;break;
	default : argc = 2;
	}
    }
  if (argc <= 1)
    {
      cout << "\n" <<argv[0]<< " builds and removes a background of a FITS image" << endl;
      cout << " usage : removeback <image in> <options>"<<endl;
      cout << "  <options>:"<<endl;
      cout << "   -o <file name> : Redirect background subtracted image in <file name>" << endl;
      cout << "   -s <number> : Box size <number> X <number> to step backimage. Default is 64." << endl;
      cout << "   -b : Remove large scale sky variations only, not sky level."<<endl;
      cout << "   -w : Write the back image as a FITS file. "<<endl;
      cout << "\n";
      return 0;
    }


  FitsImage inFits(in_name);
      

  // remove temporary sky level if  -b option is set
  float skylev=0.0;
  if (justback)
    {
      if (!inFits.HasKey("SKYLEV"))
	{
	  float sigexp=0.0;
	  inFits.MeanSigmaValue(&skylev,&sigexp);
	  float sigth=0.0;
	  if (skylev > 0.0) sigth = sqrt(skylev);
	  inFits.AddKey("SKYLEV",skylev,"Sky Level in photoelectrons");
	  inFits.AddOrModKey("SKYSIGTH",sigth,"Sigma of Sky Level obtained from sqrt(SKYLEV)");
	  inFits.AddOrModKey("SKYSIGEX",sigexp,"Sigma of Sky Level obtained from Image");
	}
      else skylev = inFits.KeyVal("SKYLEV");
      inFits -= skylev;
    }
  
  // compute and subtract back image 
  Image outImage(inFits.Nx(),inFits.Ny());
  inFits.Surface(backStep,outImage);
  float sigsky=0.0,dumb=0.0;
  outImage.MeanSigmaValue(&dumb,&sigsky);
  if (justback) outImage += skylev;
  
  // write back image is -w option is set
  if (writeflag) 
    {
      Image backImage = inFits - outImage;
      string back_name = "back_"+in_name;
      FitsImage backFits(back_name,backImage);
      backFits.AddKey("MESHSTEP",backStep,"Mesh Size used to create the surface background");
    }

  // write subtracted FITS image
  if (!outflag) out_name = "backsub_"+in_name; 
  FitsImage outFits(out_name,inFits,outImage);
  outFits.AddKey("SKYSIGBK",sigsky,"Sigma of Sky Level obtained from Image after background removal");
  char comment[64];
  sprintf(comment,"Background subtracted using %d X %d step mesh",backStep,backStep);
  outFits.AddCommentLine(comment);
  
}
