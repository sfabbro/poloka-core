#include <iostream>
#include <fstream>
#include <fileutils.h>
#include <dictfile.h>
#include <basestar.h>
#include <frame.h>
#include <wcsutils.h>
#include <reducedimage.h>
#include <fitsimage.h>
#include <gtransfo.h>
#include <photstar.h>
#include <lightcurve.h>
#include <simfitphot.h>
#include <vutils.h>

#include <map>
#include <iomanip>

using namespace std;

static void usage(const char *pgname)
{
  cerr << pgname << " <dbimage1> <dbimage2> <dbimage3> ... -r <referenceimage> -c <catalog>" << endl ;
  cerr << "options:"<< endl;
  cerr << "     -o <catalog> : output catalog name (default is calibration.list)" << endl;
  exit(1);
}

static string getkey(const string& prefix, const DictFile& catalog) {
  if(!catalog.HasKey(prefix)) {
    cerr << "catalog does not not have info about " << prefix << endl;
    exit(2);
    return "absent";
  }
  return prefix;
}

class CalibratedStar : public BaseStar {
 public:
  CalibratedStar() {};
  CalibratedStar(BaseStar& toto) : BaseStar(toto) {};
  double ra,dec;
  double u,g,r,i,z,x,y;
  double ue,ge,re,ie,ze;
  int id;
};

int main(int argc, char **argv)
{
  string referencedbimage = "";
  string catalogname = "";
  string matchedcatalogname = "calibration.list";
  vector<string> dbimages;
  
  if (argc < 7)  {usage(argv[0]);}
  for (int i=1; i<argc; ++i)
    {
      char *arg = argv[i];
      if (arg[0] != '-')
	{
	  dbimages.push_back(arg);continue;
	}
      switch (arg[1])
	{
	case 'r' : referencedbimage = argv[++i]; break;
	case 'c' : catalogname = argv[++i]; break;
	case 'o' : matchedcatalogname = argv[++i]; break;
	default : 
	  cerr << "unknown option " << arg << endl;
	  usage(argv[0]);
	}
    }
  

  if(!FileExists(catalogname)) {
    cerr << "cant find catalog " << catalogname << endl;
    usage(argv[0]);
  }
  cout << "catalog          = " << catalogname << endl;
  cout << "referencedbimage = " << referencedbimage << endl;
  cout << "n. dbimages      = " << dbimages.size() << endl;
  

  // put all of this info in a LightCurveList which is the food of the photometric fitter
  LightCurveList lclist;
  lclist.RefImage = new ReducedImage(referencedbimage);
  for (unsigned int im=0;im<dbimages.size();++im) {
    lclist.Images.push_back(new ReducedImage(dbimages[im]));
  }
  
  // we know want to put new objects in the list
  lclist.Objects.clear();
  lclist.clear();

  
  FitsHeader header(lclist.RefImage->FitsName());
  // prepare transfo and frames for stars' selection
  Frame W = lclist.RefImage->UsablePart();
  // W = W.Rescale(1.); // remove boundaries
  Gtransfo* Pix2RaDec=0;
  WCSFromHeader(header, Pix2RaDec);
  Gtransfo *RaDec2Pix = Pix2RaDec->InverseTransfo(0.01,W);
  Frame radecW = (W.ApplyTransfo(*Pix2RaDec)).Rescale(1.1);

  DictFile catalog(catalogname);
  
  int requiredlevel=2;
  // get keys for mag
  string band = header.KeyVal("TOADBAND");
  string mag_key=getkey("m"+band,catalog);
  
  BaseStar star;
  int count_total=0;
  int count_total_stars=0;
  int count_ok=0;
  
  double mag;
  char name[100];
  double mag_med;
  double mag_rms;
  map<RefStar*,CalibratedStar> assocs;
  
  for(DictFileCIterator entry=catalog.begin();entry!=catalog.end();++entry) {

    count_total++;
    
    if(int(entry->Value("level"))<requiredlevel)  continue; // not a star with correct level 
    mag=entry->Value(mag_key);
    
    count_total_stars++;
    
    star.x=entry->Value("x"); // ra (deg)
    star.y=entry->Value("y"); // dec (deg)
    
    // now check if star is in image
    if(!radecW.InFrame(star)) continue; // bye bye
    // apply transfo to this star (x,y)=pixels
    RaDec2Pix->apply(star,star);
    // check again
    if (!W.InFrame(star)) continue; // bye bye
    
    count_ok++;
    
    
    // set flux
    star.flux=pow(10.,-0.4*mag);
    
    // ok now we copy this star in a refstar and put it in the lclist.Objects
    CountedRef<RefStar> rstar = new RefStar(lclist.RefImage);
    sprintf(name,"calibstar%d",count_ok);
    rstar->name = name;
    rstar->type = 1; // a star
    rstar->band = band[0];
    rstar->x = star.x;
    rstar->y = star.y;
    rstar->ra = entry->Value("x"); // ra (deg)
    rstar->dec = entry->Value("y"); // dec (deg)
    rstar->jdmin = -1.e30; // always bright
    rstar->jdmax = 1.e30;
    lclist.Objects.push_back(rstar);
    
    // and also creat a lc (something stupid in the design)
    LightCurve lc(rstar);
    for (ReducedImageCIterator im=lclist.Images.begin(); im != lclist.Images.end(); ++im) {
      lc.push_back(*im, new PhotStar(star)); // add one image and one PhotStar
    }
    lclist.push_back(lc);
    
    // we also want to keep calibration info
    CalibratedStar cstar(star);
    cstar.ra=entry->Value("x");
    cstar.dec=entry->Value("y");
    cstar.u=entry->Value("mu");
    cstar.g=entry->Value("mg");
    cstar.r=entry->Value("mr");
    cstar.i=entry->Value("mi");
    cstar.z=entry->Value("mz");
    cstar.ue=entry->Value("emu");
    cstar.ge=entry->Value("emg");
    cstar.re=entry->Value("emr");
    cstar.ie=entry->Value("emi");
    cstar.ze=entry->Value("emz");
    cstar.flux=star.flux;
    cstar.id=count_total;
    cstar.x=star.x;
    cstar.y=star.y;
    
    //keep a link between the two of them
    assocs[(RefStar*)rstar] = cstar;
    //if(count_ok>0)
    //break;
  }
  
  cout << count_total << " objects in the catalog" << endl;
  cout << count_total_stars << " correct stars (with mag in band " << band << ") in the catalog" << endl;
  cout << count_ok << " stars in this image" << endl;
  
  // ok now let's do the fit
  SimFitPhot doFit(lclist,false);
  doFit.dowrite=false; // don't write anything before all is done
  
  // does everything
  for_each(lclist.begin(), lclist.end(), doFit);
  
  // now we want to write many many things, let's make a list
  ofstream stream(matchedcatalogname.c_str());
  stream << "@NSTARS " << count_ok << endl;
  stream << "@NIMAGES " << dbimages.size() << endl;
  stream << "#x :" << endl;
  stream << "#y :" << endl;
  stream << "#flux :" << endl;
  stream << "#error :" << endl;
  stream << "#sky :" << endl;
  stream << "#skyerror :" << endl;
  stream << "#mag :" << endl;
  stream << "#mage :" << endl;
  stream << "#ra : initial " << endl;
  stream << "#dec : initial " << endl;
  stream << "#ix : initial x" << endl;
  stream << "#iy : initial y" << endl;
  stream << "#u : from catalog" << endl;
  stream << "#g : from catalog" << endl;
  stream << "#r : from catalog" << endl;
  stream << "#i : from catalog" << endl;
  stream << "#z : from catalog" << endl;
  stream << "#ue : from catalog" << endl;
  stream << "#ge : from catalog" << endl;
  stream << "#re : from catalog" << endl;
  stream << "#ie : from catalog" << endl;
  stream << "#ze : from catalog" << endl;
  stream << "#img : image number" << endl;
  stream << "#star : start number in the catalog" << endl;
  stream << "#chi2pdf : chi2 of PSF photometry" << endl;
  stream << "#end" <<endl;
  stream << setprecision(12);
  
  
  for(LightCurveList::iterator ilc = lclist.begin(); ilc!= lclist.end() ; ++ilc) { // loop on lc
    CalibratedStar cstar=assocs.find(ilc->Ref)->second;
    //cout << "=== " << cstar.r << " " << cstar.flux << " ===" << endl;
    int count_img=0;
    double chi2pdf=ilc->chi2ndf();
    for (LightCurve::const_iterator it = ilc->begin(); it != ilc->end(); ++it) { // loop on points
      count_img++;
      const Fiducial<PhotStar> *fs = *it;
      if(fabs(fs->flux)<0.001) // do not print unfitted fluxes
	continue;
      stream << fs->x << " ";
      stream << fs->y << " ";
      stream << fs->flux << " ";
      if(fs->varflux>0)
	stream << sqrt(fs->varflux) << " ";
      else
	stream << 0 << " ";
      stream << fs->sky << " ";
      if(fs->varsky>0 && fs->varsky<10000000.)
	stream << sqrt(fs->varsky) << " ";
      else
	stream << 0 << " ";
      
      // mag
      if (band=="u") stream << cstar.u << " " << cstar.ue << " ";
      if (band=="g") stream << cstar.g << " " << cstar.ge << " ";
      if (band=="r") stream << cstar.r << " " << cstar.re << " ";
      if (band=="i") stream << cstar.i << " " << cstar.ie << " ";
      if (band=="z") stream << cstar.z << " " << cstar.ze << " ";
      
      stream << cstar.ra << " ";
      stream << cstar.dec << " ";
      stream << cstar.x << " ";
      stream << cstar.y << " ";
      stream << cstar.u << " ";
      stream << cstar.g << " ";
      stream << cstar.r << " ";
      stream << cstar.i << " ";
      stream << cstar.z << " ";
      stream << cstar.ue << " ";
      stream << cstar.ge << " ";
      stream << cstar.re << " ";
      stream << cstar.ie << " ";
      stream << cstar.ze << " ";
      stream << count_img << " "; 
      stream << cstar.id << " "; 
      stream << chi2pdf << " ";      
      stream << endl;
    }
  }
  stream.close();
  return EXIT_SUCCESS;
}

