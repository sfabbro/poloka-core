#include <iostream>


#include "transformedimage.h"
#include "usnoutils.h"
#include "fitsimage.h"
#include "imagesum.h"

#include <vector>

void usage(char *progName)
{
  cerr << progName << " -n <sumName> -r <refName> <images> " << endl;
}

int main(int nargs, char **args)
{
  string refName,sumName;
  vector<string> toSum;
  for (int i=1; i<nargs; ++i)
    {
      char *arg = args[i];
      if (arg[0] != '-') 
	{
	  toSum.push_back(arg);
	  continue;
	}
      switch(arg[1])
	{
	case 'r' : ++i; refName = args[i]; break;
	case 'n' : ++i; sumName = args[i]; break;
	default : usage(args[0]); exit(1);
	}
    }
  //ReducedImage sum(sumName);
  // cout << sum.Saturation() << endl;
  ImagesAlignAndSum(toSum, refName, sumName, DoFits | DoCatalog | DoSatur);
  //DbImage dbimage(sumName);
  string listName =  (DbImage(sumName)).ImageCatalogName(SExtractor);
  string Name =  (DbImage(sumName)).FitsImageName(Calibrated);
  bool write = true;
  UsnoProcess(Name,listName, NULL);
  return EXIT_SUCCESS;
}





#ifdef STORAGE
DbImage ref(args[1]);
DbImage to_align(args[2]);

ImageGtransfo imtransfo(ref,to_align);

TransformedImage timage("T"+to_align.Name(), 
			to_align,
			&imtransfo);
cout << timage.FitsName() << timage.CatalogName() << endl;
FitsImage toto(timage.FitsName(),RW);
timage.SetSeeing(12.);
#endif





