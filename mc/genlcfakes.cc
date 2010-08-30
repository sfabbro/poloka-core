#include "simulation.h"
#include "toadscards.h" // for SetDatacardsFileName()

#include "reducedimage.h"
#include "wcsutils.h"
#include "gtransfo.h"
#include "fitsimage.h"
#include "addlcfakes.h"

#include <algorithm> // max min sort
using namespace std;

/* ingredients :
   list of stars
   list of galaxies

   in radec or image coordinates

   standard ingredients of the fake generation
   as defined in DatSim


in datacards, expect:
# for other method names see in datacards/sub.datacard
@GENERATION_METHOD FAKES_ADAPTED 
@NUMBER_OF_FAKES
# mandatory : for the band you are considering:
@MIN_FAKE_MAG_I
@MAX_FAKE_MAG_I
@MIN_FAKE_MAG_R
@MAX_FAKE_MAG_R
@MIN_FAKE_MAG_G
@MAX_FAKE_MAG_G


# difference between host and SN
@MIN_DELTA_MAG
@MAX_DELTA_MAG
@GEOM_REF_NAME <dbimage>
# optional :
@MODEL_STARS filename <RADEC (if coordinates in ra,dec)>
@GALAXIES filename <RADEC ....>

*/


bool  read_list_from_file(StringList &imageList, 
			  const string &ImageListFileName)
{
  FILE *f = fopen(ImageListFileName.c_str(),"r");
  if (!f)
    {
      cout << " cannot open " << ImageListFileName << endl;
      return false;
    }
  char line[512];
  while (fgets(line,512,f))
    {
      if (line[0] == '#') continue;
      char imagename[512];
      sscanf(line,"%s",imagename);
      imageList.push_back(string(imagename));
    }
  fclose(f);
  return true;
}


#include "dicstar.h"
#include "listmatch.h"

bool read_star_list(DataCards &cards, const string &Key,
		    const ReducedImage &geomRef, SEStarList &SelectedStars)
{
  SEStarList catalog(geomRef.CatalogName());


  if (cards.HasKey("MODEL_STARS"))
    {
      Gtransfo *transfo = NULL;
      string modelFileName ;
      if (cards.NbParam("MODEL_STARS")==2)
	{
	  modelFileName = cards.SParam("MODEL_STARS",1);
	  string list_coords = cards.SParam("MODEL_STARS",2);	
	  if (strstr(list_coords.c_str(),"RADEC"))
	    {
	      Gtransfo *wcs = NULL;
	      FitsHeader  head(geomRef.FitsName());
	      if (!WCSFromHeader(head, wcs))
		exit(EXIT_FAILURE);
	      transfo = wcs->InverseTransfo(0.01, Frame(head));
	      delete wcs;
	      cout << "lists are assumed to be provided in sideral coordinates" << endl;
	    }
	}
      else // default : coordinates in geom_ref coordinates
	{
	  modelFileName = cards.SParam("MODEL_STARS");
	}
      DicStarList models(modelFileName);
      if (transfo) models.ApplyTransfo(*transfo);
      BaseStarList &baseModels = (BaseStarList &) models;
      StarMatchList* matches = ListMatchCollect(baseModels, 
						*SE2Base(&catalog),
						2./*pixels*/);
      // should remove ambiguities....
      // now extract the matches
      for (StarMatchIterator i = matches->begin(); i != matches->end(); ++i)
	{
	  SelectedStars.push_back(dynamic_cast<SEStar*>((BaseStar*)(i->s2)));
	}
      delete matches;
    }
  // remove objects too close form the edges and not within the image
  SelectedStars.CutEdges(geomRef.UsablePart(), 10. /*pixels */);
  return true;
}


int main(int nargs, const char **args)
{
  if (nargs<=1)
    {
      cerr << " usage : "<< endl
	   << args[0] << " <datacards file name> " << endl;
      exit(EXIT_FAILURE);
    }
  string cardsFileName(args[1]);
  SetDatacardsFileName(cardsFileName);// for DatSim to read the right file
  DataCards cards(cardsFileName);
  if (!cards.HasKey("IMAGE_LIST"))
    {
      cout << " The fake simulator needs an image list in the IMAGE_LIST key"  << endl;
      exit(EXIT_FAILURE);
    }
  StringList imageList;
  read_list_from_file(imageList, cards.SParam("IMAGE_LIST"));
  string geomRefName;
  if (cards.HasKey("GEOM_REF_NAME"))
    geomRefName = cards.SParam("GEOM_REF_NAME");
  else
    geomRefName = imageList.front();

  string photRef;
  if (!cards.HasKey("PHOT_REF"))
    {
      cout << " The fake simulator needs an image list in the PHOT_REF key"  
	   << endl;
      exit(EXIT_FAILURE);
    }      
  photRef = cards.SParam("PHOT_REF");


  double mmjd_start;
  double mmjd_stop;
  if (cards.HasKey("MMJD_START") && cards.HasKey("MMJD_STOP"))
    {
      mmjd_start = cards.DParam("MMJD_START");
      mmjd_stop = cards.DParam("MMJD_STOP");
    }
  else
    {
      double *mmjds = new double[imageList.size()];
      int count = 0;
      for (StringCIterator i = imageList.begin(); i != imageList.end(); ++i)
	{
	  ReducedImage im(*i);
	  double mmjd = im.ModifiedModifiedJulianDate();
	  mmjds[count++] = mmjd;
	}
      sort(mmjds, mmjds+count);
      mmjd_start = mmjds[count/4]-0.001;
      mmjd_stop = mmjds[(3*count)/4]+0.001;
      delete [] mmjds;
    }

  cout << " mmjd_{start,stop} " <<  mmjd_start << ' ' <<  mmjd_stop << endl;
 

  ReducedImage geomRef(geomRefName);
  if (!geomRef.IsValid())
    {
      cout << " could not locate geomref " << geomRefName << endl;
      exit(EXIT_FAILURE);
    }

  SEStarList modelStars;
  if (cards.HasKey("MODEL_STARS"))
    {
      if (!read_star_list(cards, "MODEL_STARS", geomRef, modelStars))
	{
	  cout << " could not collect model stars" << endl;
	  exit (EXIT_FAILURE);
	}
    }
  else
    {
      SEStarList catalog(geomRef.CatalogName());
      FitsImage image(geomRef.FitsName());
      PreSelectModelStars(catalog, modelStars, image);
    }

  SEStarList galaxies;
  if (cards.HasKey("GALAXIES") && 
      !read_star_list(cards, "GALAXIES   ", geomRef, galaxies  ))
    {
      cout << " could not collect galaxies   " << endl;
      exit (EXIT_FAILURE);
    }
  // no else : FroSimWModel will collect galaxies if the 
  // provided list is empty



  ForSimWModel simulator(geomRef, modelStars, galaxies, NULL);
  SimSNWModelStarList *sns = simulator.MakeListSNWModel();

  SNAdder adder;
  
  ofstream stream("generated.list");
  cout << "########## generated ################" << endl;
  sns->dump();
  cout << "#####################################" << endl;
  stream << "# xModel : " << endl
       << "# yModel : " << endl
       << "# xgen : " << endl
       << "# ygen : " << endl
       << "# photFactor : " << endl
       << "# mmjd_start : " << endl
       << "# mmjd_end : " << endl
       << "#end " << endl;
  for (SimSNWModelStarCIterator i = sns->begin(); i != sns->end(); ++i)
    {
      const SimSNWModelStar &s = **i;
      double photFactor = s.flux/s.model_on_ref.flux;
      double xModel = s.model_on_ref.x;
      double yModel = s.model_on_ref.y;
      int xShift = int(s.x-xModel+0.5);
      int yShift = int(s.y-yModel+0.5);
      stream << xModel << ' ' 
	   << yModel << ' ' 
	   << s.x << ' ' 
	   << s.y << ' ' 
	   << photFactor << ' ' 
	   << mmjd_start << ' '
	   << mmjd_stop << ' ' 
	   << endl;
      BaseStar model(xModel,yModel,0);
      SNObservations observations(model, xShift, yShift);
            for (StringCIterator i = imageList.begin(); i != imageList.end(); ++i)
	{
	  const string &imageName = *i;
	  ReducedImage ri(imageName);
	  double mmjd = ri.ModifiedModifiedJulianDate();
	  double this_pf = photFactor;
	  if (mmjd<mmjd_start || mmjd > mmjd_stop) this_pf = 0;
	  observations.AddObservation(imageName,this_pf);
	}
      adder.AddObject(observations);
    }
  stream.close();
  adder.GenerateImages();
  adder.GenerateLightCurveFile("config_fakesns.txt",photRef);
  return EXIT_SUCCESS;
}
	  

  
  
