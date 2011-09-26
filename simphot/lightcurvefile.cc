#include "lightcurvefile.h"
#include "lightcurvesyntaxerror.h"
#include "simphotfit.h"


// routine to read a DbImage in a lightcurve file. 
// so far kept local to this program
static ReducedImage* read_image(const string& line)
{
  vector<string> words;
  DecomposeString(words, line);
  RemovePattern(words.front()," ");

  ReducedImage *red = new ReducedImage(words.front());
  /*
  if (!red->ActuallyReduced())
    {
      throw(LightCurveSyntaxError(line, " image is not reduced " ));
      return NULL;
    }
  //cout << "mmjd " << words.front() << " " << red->ModifiedModifiedJulianDate() << endl;
  */
  return red;
}




#include <fstream> // for ifstream

LightCurveFile::LightCurveFile(const string &LCFileName)
{
  geomRef = (ReducedImage *)NULL;
  writeVignettes = false;
  writeMatrices = false ; 
  subDirPerObject = false;
  useProperMotions = false;

  ifstream stream(LCFileName.c_str());
  if (!stream)
    throw(PolokaException(" lightfile \""+LCFileName+"\" could not be opened"));

  bool objin=false, rimin=false, phoin=false, pmlistin=false;
  string line;
 
  while (getline(stream, line))
    {
      if (line[0] == '#' || line.length() == 0) continue; 
      
      if (line.find("OBJECTS") != line.npos)
	{ rimin=false; objin=true; phoin=false;  pmlistin = false; continue; }
      
      if (line.find("IMAGES") != line.npos)   
	{ rimin=true;  objin=false; phoin=false;  pmlistin = false; continue; }

      if (line.find("PHOREF") != line.npos)   
	{ rimin=false;  objin=false; phoin=true; pmlistin = false; continue; }

      if (line.find("PMLIST") != line.npos)
	{ rimin=false;  objin=false; phoin=false; pmlistin = true; continue;  }

      if (objin) 
	{
	  ObjectToFit *obj = new ObjectToFit(line);
	  DecodeFitType(obj->FitType()); // may throw
	  objects.push_back(obj);
	  continue;
	}

      if (rimin) 
	{ 
	  ReducedImage *im = read_image(line); 
	  if (im && images.Find(im->Name()) == images.end()) 
	      images.push_back(im); 
	  else
	    if (im)
	      {
		delete im;
		throw(LightCurveSyntaxError(line, "image appears twice in "+LCFileName ));
	      }
	  continue; 
	}
      if (phoin) // PHOREF 
	{ 
	  ReducedImage *im = read_image(line); 
	  if (im) 
	    {
	      geomRef = im; 
	      if (images.Find(im->Name()) == images.end()) 
		images.push_back(im);
	    }
	  continue; 
	}
      if (pmlistin)
	{
	  pmStarListFileName = line;
	  RemovePattern(pmStarListFileName," ");
	  useProperMotions = true;
	  continue;
	}

      // if we get here, we have a problem
      throw(LightCurveSyntaxError(line, "could not interpret line"));
    }// end of parsing
  stream.close();

  if (!geomRef)
    {
      cout << " No Geom REf provided (PHOREF card...) : cannot fit " << endl;
      throw (LightCurveSyntaxError(""," lightfile \""+LCFileName+"\" does not contain a PHOREF line"));
    }
}



#include "calibratedstar.h"
#include <fstream>

bool  LightCurveFile::SimPhotFitAllCalib(const string &CalibCatalogName, const string &OutputCatalog, int itype, int Nmax, double vignette_size_n_seeing) const
{
  bool status = true;
  string band = geomRef->Band();
  CalibratedStarList cls(CalibCatalogName, geomRef);
  ofstream tuple;
  int star_count = 0;
  for (CalibratedStarCIterator i = cls.begin(); i != cls.end(); ++i, star_count++)
    {

      if (Nmax > 0 && star_count > Nmax) break ;

      const CalibratedStar &cs = **i;
      char line[128];
      sprintf(line," %f %f DATE_MIN=-1e30 DATE_MAX=1e30 NAME=calib_%d TYPE=%d BAND=%s", cs.x,cs.y,cs.id, itype, band.c_str());
      
      ObjectToFit object(line);
      try
	{
	  SimPhotFit simPhotFit(object, *this);

	  if(vignette_size_n_seeing>0) 
	    simPhotFit.vignette_size_n_seeing = vignette_size_n_seeing ;

	  bool thisStatus = simPhotFit.DoTheFit();
	  cout << "Fit Done for star " << star_count << endl ;
	  status &= thisStatus;
	  if (thisStatus)
	    {
	      if (i == cls.begin())// first of the list
		{
		  tuple.open(OutputCatalog.c_str());
		  tuple <<"@CALIBCATALOG " << CalibCatalogName << endl;
		  simPhotFit.WriteTupleHeader(tuple, cls.size());
		}
	      simPhotFit.WriteTupleEntries(tuple, cs);

	    }
	  else
	    {
	      cout << " ERROR : The fit went wrong for object " 
		   << object << endl;
	      cout << " no output for this objects "<< endl;
	    }
	}
      catch (SimPhotFitException e)
	{
	  cout << " Problem in the fit of object \"" 
	       << object.Name() << "\":" << endl;
	  cout << e.message() << endl;
	  status = false;	    
	}
    }
  return status;
}


bool  LightCurveFile::SimPhotFitAll(double vignette_size_n_seeing) const
{
  bool status = true;
  for (ObjectToFitCIterator i = objects.begin(); i != objects.end(); ++i)
    {
      const ObjectToFit &object = **i;
      try
	{
	  SimPhotFit simPhotFit(object, *this);

	  if(vignette_size_n_seeing>0) 
	  simPhotFit.vignette_size_n_seeing = vignette_size_n_seeing ;

	  bool thisStatus = simPhotFit.DoTheFit();
	  status &= thisStatus;
	  if (thisStatus) // output results
	    {
	      string dir = "./";
	      if (subDirPerObject)
		{
		  dir = object.Name();
		  if (object.Band() != "") dir +="_"+object.Band();
		  dir +="/";
		  MKDir(dir.c_str());

		}
	      simPhotFit.Write(dir,writeVignettes,writeMatrices);
	    }
	  else // the fit went wrong:
	    {
	      cout << " ERROR : The fit went wrong for object " 
		   << object << endl;
	      cout << " no output for this objects "<< endl;
	    }
	}
      catch (SimPhotFitException e)
	{
	  cout << " Problem in the fit of object \"" 
	       << object.Name() << "\":" << endl;
	  cout << e.message() << endl;
	  status = false;
	}
    }// end loop on objects
  return status;
}



const PmStar* LightCurveFile::FindNearestPmStar(const Point &Where) const
{
  if (!useProperMotions) return NULL;
  if (pmStarList.size() == 0) 
    {
      cout << " INFO : reading " << pmStarListFileName << " for proper motions " << endl;
      ((PmStarList &) pmStarList).read(pmStarListFileName);  // will throw if file is missing
    }
  return pmStarList.FindClosest(Where.x, Where.y);
}
