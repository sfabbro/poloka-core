// -*- C++ -*-
#include <getopt.h>
#include <math.h>

#include <map>
#include <vector>
#include <iostream>

#include <reducedimage.h>
#include <fileutils.h>
#include <dicstar.h>
#include <gtransfo.h>

#include <listmatch.h>
#include <starmatch.h>


using namespace std;


bool cmp_aperse_zp(string const& dbimname, map<string,string> const&,
		   double zp, double ezp);
bool cmp_psf_zp(string const& dbimname, map<string,string> const&,
		double zp, double ezp);


void usage()
{
  cerr << "usage: calibratedbim [OPTIONS] <dbimnames>" << endl;
  cerr << " OPTIONS:" << endl;
  cerr << "   -h    print this message" << endl;
  exit(-1);
}



int main(int argc, char** argv)
{
  int i;
  string dbimname = "";
  
  char c;
  while( (c=getopt(argc, argv, "aph" )) != -1 )
    switch(c) {
    case 'h':
      usage();
      break;
    default:
      usage();
    }
  if(argc == optind) usage();
  
  vector<string>  dbims;
  for(i=optind;i<argc;i++)
    dbims.push_back(argv[i]);
  
  // read the reference catalogs
  //  map<string,DicStarList*> catalogs;
  //  catalogs["D1"] = new DicStarList("/sps/snls/CALIBRATION/secondarycatalogs_3/D1.list");
  //  catalogs["D2"] = new DicStarList("/sps/snls/CALIBRATION/secondarycatalogs_3/D2.list");
  //  catalogs["D3"] = new DicStarList("/sps/snls/CALIBRATION/secondarycatalogs_3/D3.list");
  //  catalogs["D4"] = new DicStarList("/sps/snls/CALIBRATION/secondarycatalogs_3/D4.list");
  map<string,string> catalog_names;
  catalog_names["D1"] = "/sps/snls/CALIBRATION/secondarycatalogs_3/D1.list";
  catalog_names["D2"] = "/sps/snls/CALIBRATION/secondarycatalogs_3/D2.list";
  catalog_names["D3"] = "/sps/snls/CALIBRATION/secondarycatalogs_3/D3.list";
  catalog_names["D4"] = "/sps/snls/CALIBRATION/secondarycatalogs_3/D4.list";
  
  double aper_zp, eaper_zp;
  double psf_zp, epsf_zp;
  
  vector<string>::iterator I;
  for(I=dbims.begin();I!=dbims.end();I++) {
    cmp_aperse_zp(*I, catalog_names, aper_zp, eaper_zp);
    cmp_psf_zp(*I,    catalog_names, psf_zp,  epsf_zp);
  }
}




struct Measurement {
  double val;
  double eval;
};


Measurement wmean(vector<Measurement>& vx, double& xi2, int& ndf, bool renorm=false)
{
  double sum=0, w2sum=0;
  vector<Measurement>::iterator I;
  for(I=vx.begin();I!=vx.end();I++) {
    double w2 = I->eval * I->eval;
    sum   += I->val / w2 ;
    w2sum += 1. / w2;
  }
  
  Measurement ret;
  ret.val  = sum / w2sum;
  ret.eval = sqrt(1./w2sum);
  
  xi2=0;
  for(I=vx.begin();I!=vx.end();I++) {
    double r = (I->val - ret.val) / I->eval;
    xi2 += r*r;
  }
  ndf = vx.size() - 1;
  if(ndf<=0) xi2 = 0;
  else xi2 /= ndf;
  
  if(renorm)
    ret.eval *= sqrt(xi2);
  
  return ret;
}




bool cmp_aperse_zp(string const& dbimname, map<string,string> const& catalogs,
		   double zp, double ezp) 
{
  ReducedImage redim(dbimname);
  if(!redim.IsValid()) {
    cout << "invalid reduced image: " << dbimname << endl;
    return false;
  }
  
  string apercat = dbimname + "/aperse.list";
  
  if(!FileExists(apercat)) {
    cout << "invalid apercat: " << apercat << endl;
    return false;
  }
  
  DicStarList stl(apercat);
  
  string target = redim.Target();
  map<string,string>::const_iterator Icat;
  Icat = catalogs.find(target);
  if(Icat==catalogs.end()) {
    cout << "unable to find catalog for target: " << target << endl;
    return false;
  }
  DicStarList catstl(Icat->second);
  
  string band = redim.Band();
  if(band!="g" && band!="r" && 
     band!="i" && band!="z") {
    cout << "invalid band: " << band << endl;
    return false;
  }
  
  // OK. We should be all set now...
  
  // so, we apply the WCS transfo to the catalog
  Gtransfo* pix2radec = redim.PixelsToRaDec();
  if(!pix2radec) return false;
  stl.ApplyTransfo(*pix2radec);
  
  // then, we attempt to find a transformation
  StarMatchList* mlist = ListMatchCollect((BaseStarList&)stl,
					  (BaseStarList&)catstl,
					  0.00055); // arcsec
  
  // and compute the ZP, in the right band
  vector<Measurement> mvec;
  StarMatchIterator Istar;
  for(Istar=mlist->begin();Istar!=mlist->end();Istar++) {
    DicStar* s1 = (DicStar*)(Istar->s1);
    DicStar* s2 = (DicStar*)(Istar->s2);
    double flx  = s1->getval("apfl6");
    double eflx = s1->getval("eapfl6");
    double mag  = s2->getval("m" + band);
    double emag = s2->getval("em" + band);
    if(mag<=0 || emag <=0) continue;
    
    Measurement M;
    M.val  = mag + 2.5*log10(flx);
    M.eval = 1.0857362047581294 * eflx / flx;
    mvec.push_back(M);
  }
  
  int ndf;
  double xi2;
  Measurement ret = wmean(mvec, xi2, ndf, true);
  cout << " ZP_APER " << ret.val << " " << ret.eval << " " << xi2 << " " << ndf << " " << endl;
  
  return true;
}




bool cmp_psf_zp(string const& dbimname, map<string,string> const& catalogs,
		double zp, double ezp)
{
  ReducedImage redim(dbimname);
  if(!redim.IsValid()) {
    cout << "invalid reduced image: " << dbimname << endl;
    return false;
  }
  
  string psfcat = dbimname + "/psfstars.list";
  if(!FileExists(psfcat)) {
    cout << "invalid psfcat: " << psfcat << endl;
    return false;
  }
  
  DicStarList stl(psfcat);
  
  string target = redim.Target();
  map<string,string>::const_iterator Icat;
  Icat = catalogs.find(target);
  if(Icat==catalogs.end()) {
    cout << "unable to find catalog for target: " << target << endl;
    return false;
  }
  DicStarList catstl(Icat->second);
  
  string band = redim.Band();
  if(band!="g" && band!="r" && 
     band!="i" && band!="z") {
    cout << "invalid band: " << band << endl;
    return false;
  }
  
  StarMatchList* mlist = ListMatchCollect((BaseStarList&)stl,
					  (BaseStarList&)catstl,
					  0.00055); // arcsec
  
  vector<Measurement> mvec;
  StarMatchIterator Istar;
  for(Istar=mlist->begin();Istar!=mlist->end();Istar++) {
    DicStar* s1 = (DicStar*)(Istar->s1);
    DicStar* s2 = (DicStar*)(Istar->s2);
    double flx  = s1->getval("flux");
    double eflx = s1->getval("eflux");
    double mag  = s2->getval("m" + band);
    double emag = s2->getval("em" + band);
    if(mag<=0 || emag <=0) continue;
    
    Measurement M;
    M.val  = mag + 2.5*log10(flx);
    M.eval = 1.0857362047581294 * eflx / flx;
    mvec.push_back(M);
  }
  
  int ndf;
  double xi2;
  Measurement ret = wmean(mvec, xi2, ndf, true);
  cout << " ZP_PSF " << ret.val << " " << ret.eval << " " << xi2 << " " << ndf << " " << endl;
  
  return true;
}
