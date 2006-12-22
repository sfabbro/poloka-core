
#include <iostream>


#include "fitsimage.h"
#include "wcsutils.h"
#include "reducedimage.h"
#include "gtransfo.h"
#include "ccd.h"


static void usage(const char * prog)
{
  std::cerr << "usage : " << std::endl
	    << prog << " -r <refName> -s <Nx>x<Ny> <new DbImages> " << std::endl;
  exit(EXIT_FAILURE);
}

int main(int nargs, const char**args)
{
  std::string refName;
  StringList newNames;
  int nx=0;
  int ny=0;
  for (int i=1; i< nargs; ++i)
    {
      const char *arg = args[i];
      if (arg[0] == '-')
	switch (arg[1])
	  {
	  case 'r' : ++i; refName = args[i]; break;
	  case 's' : 
	    ++i; 
	    if (sscanf(args[i], "%dx%d", &nx,&ny) != 2)
	      {
		std::cerr << " cannot decode NXxNY in " << args[i] 
			  << std::endl;
		usage(args[0]);
	      }
	    break;
	  default: usage(args[0]);
	  }
      else
	{
	  newNames.push_back(string(arg));
	}
    }
  if (refName == "")
    {
      std::cerr << " ERROR: no ref provided " << std::endl;
      usage(args[0]);
    }
  ReducedImage ref(refName);
  if (!ref.IsValid())
    {
      std::cerr << " ERROR : cannot find " << refName << std::endl;
      exit(EXIT_FAILURE);
    }
  if (nx ==0 || ny == 0)
    {
      std::cerr << " you did not provide the split you want " << std::endl;
      usage(args[0]);
    }
  Toads_ccd::Ccds ccds(refName);
  for (StringCIterator i= newNames.begin(); i != newNames.end(); ++i)
    if (!ccds.AddImage(*i))
      {
	std::cerr << " drop image " << *i << std::endl;
      }
  
  int overlap = 20;
  const Frame &refFrame = ccds.WholeRefFrame();
  int xChunk = (int(refFrame.Nx())+(nx-1)*overlap)/nx+1;
  int yChunk = (int(refFrame.Ny())+(ny-1)*overlap)/ny+1;
  int xStep = xChunk-overlap;
  int yStep = yChunk-overlap;
    
  int resample_count = 0;

  int count = 0;
  for (int j=0; j<ny; ++j)
    for (int i=0; i< nx; ++i)
      {
	double ymax = refFrame.Ny() - j*yStep;
	double ymin = ymax - yChunk;;
	double xmin = i*xStep;
	double xmax = xmin+xChunk;
	Frame subFrame(Point(xmin,ymin),Point(xmax,ymax));
	subFrame *= refFrame; // on sides, it goes beyond large ref limits
	StringList insideList;
	ccds.SelectInsiders(subFrame, insideList);
	// format : makesubfile <count> <refName> <subrRefframe> '<listOfNew>'
	char pseudoccd[8];
	sprintf(pseudoccd,"%02d", count);
	std::cout << " makesubfile "
		  << pseudoccd << ' ' 
		  << ref.Name() << ' ' 
		  << '\'' 
		  << '[' << subFrame.xMin << ':' << subFrame.xMax << ','
		  << subFrame.yMin << ':' << subFrame.yMax << "] "
		  << "\' "
		  << '\'' << insideList << '\'' << std::endl;
	count++;
	resample_count += insideList.size();
      }
  std::cerr << " nresamplings/nimages " << resample_count <<'/' << newNames.size() << '=' << double(resample_count)/double(newNames.size()) << std::endl;
  return EXIT_SUCCESS;
}

	  

