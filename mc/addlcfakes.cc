#include "addlcfakes.h"

#include "simulation.h"
#include "sestar.h"


 

#include "fitsimage.h"
#include "reducedimage.h"
#include "fileutils.h"

static bool copy_file(const string &From, 
		      const string &ToDir)
{
  string ToName = AddSlash(ToDir)+BaseName(From);
  if (FileExists(ToName))
    {
      cout << " file " << ToName << " already exists, nothing done " << endl;
      return true;
    }
  string com = "cp -f "+From+" "+ToDir;
  bool ok = (0==system(com.c_str()));
  if (!ok)
    cout << " could not execute :" << com << endl;
  return ok;
}

static bool link_file(const string &From, const string &ToDir)
{
  bool ok = true;
  StringList e;
  ExpandPath( From.c_str(), e);
  for (StringIterator i= e.begin(); i != e.end(); ++i)
    {
      string fileName = BaseName(*i);
      ok &=MakeRelativeLink(*i, AddSlash(ToDir)+fileName);
    }
  return ok;
}

static bool PseudoCopyDbImage(const string &SourceName, const string &DestName)
{
  ReducedImage source(SourceName);
  ReducedImage dest(DestName);  
  if (!dest.IsValid()) dest.Create("here");
  // create a link inside dest that points to source
  string original_pointer = AddSlash(dest.Dir())+"original";
  MakeRelativeLink(FullFileName(source.Dir()), original_pointer);
  // link (or copy) the files themselves
  return (
  link_file(source.CatalogName(), DestName) &&
  link_file(source.FitsWeightName(), DestName) &&
  copy_file(source.FitsName(), DestName) &&
  link_file(source.Dir()+"*.xml", DestName)
  );
}


static string image_with_fakes_name(const string &without_fakes_name)
{
  return without_fakes_name+"_addfakes";
}



Observation::Observation(const string &S, const double &pf) :
    inputName(S), photFactor(pf) 
{
  mmjd = ReducedImage(S).ModifiedModifiedJulianDate();
};


void SNObservations::AddObservation(const string &ImageName, 
				    const double PhotFactor)
{
  push_back(Observation(ImageName, PhotFactor));
}
	    


void SNObservations::GetMinMaxDates(double &mindate, double &maxdate) const
{
  mindate = 1e30;
  maxdate = -1e30;
  for (ObsCIterator o=begin(); o != end(); ++o)
    {
      if (o->photFactor == 0) continue;      
      mindate = min(o->mmjd,mindate);
      maxdate = max(o->mmjd,maxdate);
    }
  mindate -= 0.0001;
  maxdate += 0.0001;
}





/**************************************************************/

/* auxillary classes just used for SNAdder::GenerateImages for
   grouping addition of fakes on a given image. Could easily be made
   without those classes, which remain from a first design.
*/

#include "simsnstar.h" // for ModelStar

   
struct Epoch {
  string inputName; // dbimage name.
  double mmjd; // jd -jd(01-01-2003).
  double seeing;
  list<ModelStar> additions;
  
  Epoch( const string &ImageName);
} ;


Epoch::Epoch(const string &ImageName) : inputName(ImageName)
{
  ReducedImage ri(inputName);
  mmjd = ri. ModifiedModifiedJulianDate();
  seeing = ri.Seeing();
}


struct Epochs : public vector<Epoch> 
{
  Epoch* Find(const string &Name)
  { for (iterator i=begin(); i!=end(); ++i) 
    if (i->inputName == Name) return &(*i);
  return NULL;
  }
  
  Epoch* FindOrAdd(const string &Name)
  {
    Epoch *e = Find(Name);
    if (!e) push_back(Epoch(Name));
    return Find(Name);
  }
  
  void AddObject(const SNObservations &LightCurve);

  Epochs(const SNList &L);

};
  

Epochs::Epochs(const SNList &L)
{
  for (SNCIterator i = L.begin(); i != L.end(); ++i) AddObject(*i);
}

void Epochs::AddObject(const SNObservations &LightCurve)
{
  for (ObsCIterator i = LightCurve.begin();
       i != LightCurve.end(); ++i)
    {
      const Observation &obs = *i;
      Epoch *epoch = FindOrAdd(obs.inputName);
      epoch->additions.push_back(ModelStar(LightCurve.modelStar,
					   int(5*epoch->seeing),
					   LightCurve.xShift,
					   LightCurve.yShift,
					   obs.photFactor));
    }
}

/**************************************************************/





void SNAdder::AddObject(const SNObservations &LightCurve)
{
  snList.push_back(LightCurve);
}



void SNAdder::GenerateImages() const
{
  Epochs epochs(snList);
  GtransfoIdentity identity;
  for (unsigned i=0; i <epochs.size(); ++i)
    {
      const Epoch &epoch = epochs[i];
      string inputName = epoch.inputName;
      string outputName  = image_with_fakes_name(inputName);
      PseudoCopyDbImage(inputName, outputName);
      ReducedImage out(outputName);
      FitsImage outFits(out.FitsName(),RW);
      outFits.SetWriteAsFloat();
      list<ModelStar>::const_iterator s;
      for (s = epoch.additions.begin(); s!= epoch.additions.end(); ++s)
	{
	  s->AddToImage(outFits, outFits, &identity,0,-1,out.Gain());
	  Point where = *s;
	  cout << " adding at " << where << " in " << outFits.FileName() 
	       << " with phot factor " << s->PhotFactor() << " (gain=" << out.Gain() << ")" << endl;
	}
    }
} 


#include <stdio.h> /* for C I/O's */

bool SNAdder::GenerateLightCurveFile(const string &ConfigFileName,
				     const string &PhoRefName) const
{
  Epochs epochs(snList);
  FILE *f = fopen(ConfigFileName.c_str(), "w");
  if (!f) 
    {
      cout << " SNAdder::GenerateLightCurveFile() : could not generate " 
	   << ConfigFileName  << endl;
      return false;
    }
  fprintf(f,"OBJECTS\n");
  int count = 1;
  string band = ReducedImage(epochs[0].inputName).Band();
  for (SNCIterator sn=snList.begin(); sn != snList.end(); ++sn,count++)
    {
      double mindate,maxdate;
      sn->GetMinMaxDates(mindate,maxdate);
      fprintf(f,"%f %f DATE_MIN=%f DATE_MAX=%f NAME=fake_%d TYPE=0 BAND=%s\n",
	      sn->modelStar.x + sn->xShift, 
	      sn->modelStar.y + sn->yShift ,
	      mindate,
	      maxdate,
	      count,
	      band.c_str());
      fprintf(f,"%f %f DATE_MIN=%f DATE_MAX=%f NAME=fixed_fake_%d TYPE=-1 BAND=%s\n",
	      sn->modelStar.x + sn->xShift, 
	      sn->modelStar.y + sn->yShift ,
	      mindate,
	      maxdate,
	      count,
	      band.c_str());
      fprintf(f,"%f %f DATE_MIN=-100 DATE_MAX=50000 NAME=model_%d TYPE=1 BAND=%s\n",
	      sn->modelStar.x , 
	      sn->modelStar.y,
	      count,
	      band.c_str());

    }
  fprintf(f,"IMAGES\n");
  for (unsigned i=0; i < epochs.size(); ++i)
    fprintf(f,"%s\n",image_with_fakes_name(epochs[i].inputName).c_str());
  fprintf(f,"PHOREF\n");
  fprintf(f,"%s\n",PhoRefName.c_str());
  
  fclose(f);
  cout << " written the input file for make_lightcurve " << ConfigFileName 
       << endl;
  return true;
}
