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
#include <map>

using namespace std;

static void usage(const char *pgname)
{
  cerr << pgname << " <config_file> <catalog>" << endl ;
  exit(1);
}


class CalibratedStar : public BaseStar {
 public:
  CalibratedStar() {};
  CalibratedStar(BaseStar& toto) : BaseStar(toto) {};
  double ra,dec;
  double u,g,r,i,z;
  int id;
};

int main(int argc, char **argv)
{
  if (argc < 3)  {usage(argv[0]);}

  string config_filename = argv[1];
  string catalog_name = argv[2];
  if(!FileExists(config_filename)) {
    cerr << "cant find " << config_filename  << endl;
    usage(argv[0]);
  }
  
  if(!FileExists(catalog_name)) {
    cerr << "cant find catalog " << catalog_name << endl;
    usage(argv[0]);
  }
  
  ifstream str(config_filename.c_str());
  LightCurveList lclist(str); // read a config file
  
  // we know want to put new objects in the list
  lclist.Objects.clear();
  lclist.clear();

  
  FitsHeader header(lclist.RefImage->FitsName());
  // prepare transfo and frames for stars' selection
  Frame W = lclist.RefImage->UsablePart();
  W = W.Rescale(0.8); // remove boundaries
  Gtransfo* Pix2RaDec=0;
  WCSFromHeader(header, Pix2RaDec);
  Gtransfo *RaDec2Pix = Pix2RaDec->InverseTransfo(0.01,W);
  Frame radecW = (W.ApplyTransfo(*Pix2RaDec)).Rescale(1.2);

  DictFile catalog(catalog_name);
  
  //get the band of this image
  string band = header.KeyVal("TOADBAND");
  if(!catalog.HasKey(band)) {
    cerr << "catalog does not not have info about band " << band << endl;
    exit(2);
  }

  
  BaseStar star;
  int count_total=0;
  int count_total_stars=0;
  int count_ok=0;
  
  double mag;
  char name[100];
  
  map<RefStar*,CalibratedStar> assocs;
  
  for(DictFileCIterator entry=catalog.begin();entry!=catalog.end();++entry) {

    count_total++;
    
    // apply a cut on cgal
    if(double(entry->Value("cgalcat"))>0.2) continue; // not a nice star
    
    // now check mag
    mag=entry->Value(band);
    if(mag<=0 || mag>50) continue; // crazy mag
    
    count_total_stars++;
    
    star.x=entry->Value("ra");
    star.y=entry->Value("dec");
    // now check if star is in image
    if(!radecW.InFrame(star)) continue; // bye bye
    // apply transfo to this star
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
    rstar->ra = entry->Value("ra");
    rstar->dec = entry->Value("dec");
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
    cstar.ra=entry->Value("ra");
    cstar.dec=entry->Value("dec");
    cstar.u=entry->Value("u");
    cstar.g=entry->Value("g");
    cstar.r=entry->Value("r");
    cstar.i=entry->Value("i");
    cstar.z=entry->Value("z");
    cstar.flux=star.flux;
    cstar.id=count_total;

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
  ofstream stream("calibration.list");
  stream << "#star :" << endl;
  stream << "#ra :" << endl;
  stream << "#dec :" << endl;
  stream << "#u :" << endl;
  stream << "#g :" << endl;
  stream << "#r :" << endl;
  stream << "#i :" << endl;
  stream << "#z :" << endl;
  stream << "#img :" << endl;
  stream << "#flux :" << endl;
  stream << "#error :" << endl;
  stream << "#end" <<endl;
  // first let's try to compute the ZP
  double weight;
  int count=0;
  double sumzp=0;
  double sumzp2=0;
  double sumweight=0;
  double zp;
  for(LightCurveList::iterator ilc = lclist.begin(); ilc!= lclist.end() ; ++ilc) { // loop on lc
    CalibratedStar cstar=assocs.find(ilc->Ref)->second;
    //cout << "=== " << cstar.r << " " << cstar.flux << " ===" << endl;
    int count_img=0;
    for (LightCurve::const_iterator it = ilc->begin(); it != ilc->end(); ++it) { // loop on points
      count_img++;
      const Fiducial<PhotStar> *fs = *it;
      if(fabs(fs->flux)<0.001) // do not print unfitted fluxes
	continue;
      
      stream << cstar.id << " ";
      stream << cstar.ra << " ";
      stream << cstar.dec << " ";
      stream << cstar.u << " ";
      stream << cstar.g << " ";
      stream << cstar.r << " ";
      stream << cstar.i << " ";
      stream << cstar.z << " ";
      stream << count_img << " ";
      stream << fs->flux << " ";
      if(fs->varflux>0)
	stream << sqrt(fs->varflux) << " ";
      else
	stream << 0 << " ";
      stream << endl;
      
      weight = 1./fs->varflux;
      count ++;
      sumweight += weight;
      zp = 2.5*log10(fs->flux/cstar.flux);
      sumzp += zp*weight;
      sumzp2 += zp*zp*weight;
    }
  }
  zp = sumzp/sumweight;
  double rms = sqrt(sumzp2/sumweight -zp*zp);
  printf("RESULT_zp_rms_error= %6.6f %6.6f %6.6f\n",zp,rms,rms/sqrt(float(count)));
  stream.close();
  return EXIT_SUCCESS;
}

