#include <iostream>
#include <fstream>
#include <fileutils.h>
#include <dicstar.h>
#include <basestar.h>
#include <gaussianfit.h>

#include <map>
#include <iomanip>

using namespace std;

static void usage(const char *pgname)
{
  cerr << pgname << " -c <dicstarlist>" << endl;
  cerr << "options:  -m <keyofreferencemagnitude> (default is \"mag\")"<< endl;
  cerr << "          -f <keyofflux>               (default is \"flux\")" << endl ;
  cerr << "          -o <output ASCII file>       (default is zpfitter.dat)" << endl ;
  cerr << "          -d mjd_min mjd_max       (date range)" << endl ;
  
  
  exit(1);
}

static double sqr(const double& x) {return x*x;};

//////////////////////////////////////////////////////
// ZPSTAR
//////////////////////////////////////////////////////
class zpstar {
private :
  double sum_1;
  double sum_w;
  double sum_fw;
  double sum_s;
  double sum_f2w;
  double sum_w01;
  double sum_fw01;
  double sum_f2w01;
  
  
public :
  double init_ra,init_dec;
  double init_x,init_y;  
  double fitted_ra,fitted_dec;
  double fitted_x,fitted_y;
  double mag,error;
  double psfchi2,zpchi2,zpchi2_01;
  double nmag;
  double ndist; 
  double neic;
  double mag_min,mag_max;
  double g,ge,r,re,i,ie,z,ze;
  double catalog_mag;
  double catalog_mage;
  double zp;
  double sky;
  double flux;
  double fluxrms;
  int satur;
  int contam;
  int var;
  int var_3;
  int var_rms;
  int nostar;
  int id;
  
  // not saved in file:
  double flux_min,flux_max;
  
  int nmeas() const {
    return int(sum_1);
  }

  void reset() {
    id = -1;
    sum_1 = 0;
    sum_w = 0;
    sum_fw = 0;
    sum_f2w = 0;
    sum_w01 = 0;
    sum_fw01 = 0;
    sum_f2w01 = 0;
    sum_s = 0;
    mag = -1;
    error = -1;
    init_ra = 0;
    init_dec = 0;
    fitted_ra = 0;
    fitted_dec = 0;
    init_x = 0;
    init_y = 0;
    fitted_x = 0;
    fitted_y = 0;
    satur = 0;
    contam = 0;
    var = 0;
    var_3 = 0;
    var_rms = 0;
    nostar = 0;
    flux_min = 1.e30;
    flux_max = -1.e30;
    
  }

  void load(const DicStar* istar) {

    id = int(istar->getval("star"));

    init_ra = istar->getval("ra");
    init_dec = istar->getval("dec");
    init_x = istar->getval("ix");
    init_y = istar->getval("iy");
    fitted_x = istar->getval("x");
    fitted_y = istar->getval("y");
    
    psfchi2 = istar->getval("chi2pdf");
    if (istar->HasKey("neic")) {neic = istar->getval("neic");}
    else {neic=0;}
    if (istar->HasKey("neid")) {ndist = istar->getval("neid");}
    else {ndist=0;}
    g = istar->getval("g");
    ge = istar->getval("ge");
    r = istar->getval("r");
    re = istar->getval("re");
    i = istar->getval("i");
    ie = istar->getval("ie");
    z = istar->getval("z");
    ze = istar->getval("ze");  
    catalog_mag = istar->getval("mag");
    catalog_mage= istar->getval("mage");
    
    double w = 1./sqr(istar->getval("error"));
    flux = istar->getval("flux");
		     
    sum_1 += 1;
    sum_w += w;
    sum_fw += flux*w;
    sum_f2w += flux*flux*w;
    sum_s += istar->getval("sky");
    
    if(flux>0) { // why requiring (flux>0) ?? P.A.
      double w01 = 1./(sqr(istar->getval("error"))+sqr(0.01*istar->getval("flux")));
      sum_w01 += w01;
      sum_fw01 += flux*w01;
      sum_f2w01 += flux*flux*w01;

   
    }
    
    if(istar->getval("satur")>0.5) satur = 1;
    
    if(flux<flux_min)
      flux_min = flux;
    if(flux>flux_max)
      flux_max = flux;
    
  }
  
  void process() {
    flux = sum_fw/sum_w;
    error = 1./sqrt(sum_w); // flux
    fluxrms = sum_f2w/sum_w-flux*flux;
    fluxrms = (fluxrms>0) ? sqrt(fluxrms) : -1;
    error = error/flux*2.5/log(10.); // mag
    mag   = -2.5*log10(flux)+zp;
    mag_min = -2.5*log10(flux_max)+zp;
    mag_max = -2.5*log10(flux_min)+zp;
    sky = sum_s/sum_1;
    
    zpchi2 = (sum_f2w + sum_w*sqr(flux) - 2.*sum_fw*flux)/(sum_1-1.);
    if(zpchi2>1)
      error *= sqrt(zpchi2);
    
    zpchi2_01 = (sum_f2w01 + sum_w01*sqr(flux) - 2.*sum_fw01*flux)/(sum_1-1.);
    
    // scores (for WNR cat neic=0, no star will be contam flagged)
    nmag = neic/flux*2.5/log(10.);
    if(nmag>0.005)
      contam = 1;
    
    if(fluxrms>0 && flux>0 && fluxrms/flux>0.05)
      var_rms = 1;

   if(zpchi2_01>5)
      var = 1;
   if(zpchi2_01>3)
      var_3 = 1;
    
    if(psfchi2>100)
      nostar = 1;
  }
  
  zpstar() {
    zp=0;
    reset();
  }
  
  void write(ostream &s) const {
    s << init_ra << " ";
    s << init_dec << " ";
    
    s << fitted_x << " ";
    s << fitted_y << " ";
    
    
    
    s << init_x << " ";
    s << init_y << " ";
    
    s << id << " ";
    s << mag << " ";
    s << error << " ";
    s << flux << " ";
    s << fluxrms << " ";
    s << sum_1 << " ";
    s << zpchi2 << " ";
    s << zpchi2_01 << " ";
    s << nmag << " ";
    s << ndist << " ";
    s << mag_min << " ";
    s << mag_max << " ";
    s << sky << " ";
    s << psfchi2 << " ";
    s << satur << " ";
    s << contam << " ";
    s << var << " ";
    s << var_3 << " ";
    s << var_rms << " ";
    s << nostar << " ";
    s << g << " ";
    s << ge << " ";
    s << r << " ";
    s << re << " ";
    s << i << " ";
    s << ie << " ";
    s << z << " ";
    s << ze << " ";
  }
  
  friend ostream& operator << (ostream &stream, const zpstar &s)
  { s.write(stream); return stream;}
};


//////////////////////////////////////////////////////
// MAIN
//////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  
  string catalogname = "";
  string magkey = "mag";
  string fluxkey = "flux";
  string outputfilename = "zpfitter.dat";
  double mjd_min = 0;
  double mjd_max = 1.e12;
  
  if (argc < 3)  {usage(argv[0]);}
  for (int i=1; i<argc; ++i)
    {
      char *arg = argv[i];
      if (arg[0] != '-')
	{
	  cerr << "unexpected parameter " << arg << endl;
	  usage(argv[0]);
	}
      switch (arg[1])
	{
	case 'm' : magkey = argv[++i]; break;
	case 'c' : catalogname = argv[++i]; break;
	case 'f' : fluxkey = argv[++i]; break;
	case 'o' : outputfilename = argv[++i]; break;
	case 'd' : mjd_min = atof(argv[++i]); mjd_max = atof(argv[++i]);  break;
	default : 
	  cerr << "unknown option " << arg << endl;
	  usage(argv[0]);
	}
    }
  
  
  if(!FileExists(catalogname)) {
    cerr << "cant find catalog " << catalogname << endl;
    usage(argv[0]);
  }
  DicStarList catalog(catalogname);
  if ( catalog.empty()) {
    cerr << "catalog is empty" << endl;
    usage(argv[0]);
  }
  int nstars=-1;
  if(catalog.GlobVal().HasKey("NSTARS")) nstars=int(catalog.GlobVal().getDoubleValue("NSTARS"));
  int nimages=-1;
  if(catalog.GlobVal().HasKey("NIMAGES")) nimages=int(catalog.GlobVal().getDoubleValue("NIMAGES"));
  
  // fill a vector of ZPs
  unsigned int nentries = catalog.size();
  double* values = new double[nentries];
  double flux,mag;
  
  


  float contamination_cut = 0.005;
  double mjd;
  int count = 0;
  
  DicStarIterator entry=catalog.begin();
  
  while(entry!=catalog.end()) {
    
    bool bad = false;

    if(mjd_min>0) {
      mjd = (*entry)->getval("mjd");
      bad |= (mjd<mjd_min);
      bad |= (mjd>mjd_max);
    }

    bad |= ( (*entry)->getval("satur") >0.5 ); // cause has saturated pixels
    
    if(fluxkey=="flux")
      flux = (*entry)->getval("flux"); // warning WNR cats are not well ordered for disctar cause not x y flux (no time a this point for reprocessing). So ask explicitely for "flux" key 
    else
      flux = (*entry)->getval(fluxkey);

    mag = (*entry)->getval(magkey);
    bad |= (mag<10);
    bad |= ((*entry)->getval("error")<1.e-6);
    bad |= (flux<=1);
    double contamination_of_neighbour ; 
    if ((*entry)->HasKey("neic"))
	contamination_of_neighbour = (*entry)->getval("neic");
    else
	contamination_of_neighbour = 0;

    bad |= (contamination_of_neighbour/flux>contamination_cut); // cause contamination of neighbour
    if(bad) {
      entry = catalog.erase(entry);
    }else{
      values[count++]=2.5*log10(flux)+mag;
      ++entry;
    }
  }
  if (count==0) {
    cerr << "no valid entries in catalog!" << endl;
    delete [] values;
    return EXIT_FAILURE;
  }
  double zp,rms;
  FILE *file = fopen(outputfilename.c_str(),"w");
  
  // gaussian
  zp = gaussianfit(values,count,zp,rms,3.,true);
  cout << zp << " " << rms << endl;
  if(rms>0) { 
    zp = gaussianfit(values,count,zp,rms,2.,false);
    cout << zp << " " << rms << endl;
    zp = gaussianfit(values,count,zp,rms,1.5,false);
    cout << zp << " " << rms << endl;
    zp = gaussianfit(values,count,zp,rms,1.,false);
    cout << zp << " " << rms << endl;
  }else{
    rms=-rms;
  }
  fprintf(file,"@PSFZP %6.6f\n",zp);
  fprintf(file,"@PSFZPERROR %6.6f\n",rms);
  fprintf(file,"@NMEASUREMENTS %d\n",count);
  
  printf("@PSFZP %6.6f %6.6f %d\n",zp,rms,count);
  if(nstars>0) {
    fprintf(file,"@NSTARS %d\n",nstars);
  }
  if(nimages>0) {
    fprintf(file,"@NIMAGES %d\n",nimages);
  }
  
  delete [] values;
  
  
  // build a processed star list
  
  
  ofstream stream("zpstars.list");
  stream << setprecision(10);
  stream << "# using catalog " << catalogname << endl;
  stream << "# using mag key '" << magkey << "'" << endl;
  stream << "# ra : as in calib. catalog" << endl;
  stream << "# dec : as in calib. catalog" << endl;
  stream << "# x : fitted x" << endl;
  stream << "# y : fitted y" << endl;
  stream << "# ix : input x" << endl;
  stream << "# iy : input y" << endl;
  stream << "# id : number in calib. catalog" << endl;
  stream << "# mag : mag derived from match to catalog." << endl;
  stream << "# error : error on mag." << endl;
  stream << "# flux : flux in ref. image units" << endl;
  stream << "# fluxrms : scatter around the above" << endl;
  stream << "# nmeas : number of measurements involved in the average" << endl;
  stream << "# zpchi2 : chi2pdf ass. stable star without sys. effect" << endl;
  stream << "# zpchi2_01 : chi2pdf ass. stable star with 0.01 sys. effect" << endl;
  stream << "# nmag : level of contamination from neighbours (mag.)" << endl;
  stream << "# ndist : distance of nearest neighbour (=neid, pixels)" << endl;
  stream << "# mag_min : min. meas. magnitude" << endl;
  stream << "# mag_max : max. meas. magnitude" << endl;
  stream << "# sky : av. fitted residual sky level" << endl;
  stream << "# psfchi2 : chi2 pdf of PSF photometry" << endl;
  stream << "# satur : 0=ok" << endl;
  stream << "# contam : 0=ok, else contamination" << endl;
  stream << "# var : 0=ok else variable star" << endl;
  stream << "# var_3 : 0=ok else variable star" << endl;
  stream << "# var_rms : 0=ok else variable star" << endl;
  stream << "# nostar : 0=ok else probably a galaxy" << endl;   
  stream << "# g : " << endl;
  stream << "# ge : " << endl;
  stream << "# r : " << endl;
  stream << "# re : " << endl;
  stream << "# i : " << endl;
  stream << "# ie : " << endl;
  stream << "# z : " << endl;
  stream << "# ze : " << endl;
  stream << "# end " << endl;
  
  zpstar current_star;
  current_star.zp = zp;
  
  double sum_chi2_of_images = 0; // compute total chi2 with zp and star magnitudes as free parameters
  double sum_ndf_of_images  = 0; // 
  double sum_chi2_of_images_01 = 0; // same, assuming 0.01 sys. effect on, say, flatfielding or psf
  
  double sum_chi2_of_stars = 0; // compute total chi2 once dispersion between measurements of a star propagated into error
  double sum_ndf_of_stars  = 0; 
  
  int nstar_var = 0, nstar_contam=0, nstar_nostar=0;
  int nstar_var_3 = 0 ;
  int nstar_var_rms = 0 ;
  int nstar_not_var_et_rms = 0 ;
  int nstar_not_var_3_et_rms = 0 ;
  int nstar_tot=0;
  for(DicStarCIterator entry=catalog.begin();entry!=catalog.end();++entry) {
    
    // if new star, dump results here
    if(current_star.id != (*entry)->getval("star")) {
      current_star.process();
      nstar_tot++;
      if (current_star.contam> 0 ) nstar_contam++;
      if (current_star.var> 0 ) nstar_var++;
      if (current_star.var_3> 0 ) nstar_var_3++;
      if (current_star.var_rms> 0 ) nstar_var_rms++;
      if (current_star.var_rms> 0 && current_star.var <=0 ) nstar_not_var_et_rms++;
      if (current_star.var_rms> 0 && current_star.var_3 <=0 ) nstar_not_var_3_et_rms++;
      
      if (current_star.nostar> 0 ) nstar_nostar++;

      if(current_star.id>0) {
	stream << current_star << endl;
	sum_chi2_of_images += current_star.zpchi2*current_star.nmeas();
	sum_chi2_of_images_01 += current_star.zpchi2_01*current_star.nmeas();
	sum_ndf_of_images += current_star.nmeas() - 1; // -1 for the magnitude free parameter	
	sum_chi2_of_stars += sqr(current_star.mag-current_star.catalog_mag)/(sqr(current_star.catalog_mage)+sqr(current_star.error));
	sum_ndf_of_stars +=1;
      }
      current_star.reset();
    }
    current_star.load(*entry);            
  }
  stream.close();
  
  sum_ndf_of_images -= 1 ; // for the free zp parameter
  sum_ndf_of_stars -= 1 ; // for the free zp parameter


  cout << "@NSTAR_INIT " << nstar_tot << endl ;
  cout << "@NSTAR_VAR " << nstar_var << endl ; 
  cout << "@NSTAR_VAR_PCT " << 100.*nstar_var/nstar_tot  << endl ;
  cout << "@NSTAR_VAR_3 " << nstar_var_3 << endl ; 
  cout << "@NSTAR_VAR_3_PCT " << 100.*nstar_var_3/nstar_tot  << endl ;
  cout << "@NSTAR_VAR_RMS " << nstar_var_rms << endl ; 
  cout << "@NSTAR_VAR_RMS_PCT " << 100.*nstar_var_rms/nstar_tot  << endl ;
  cout << "@NSTAR_NOT_VAR_AND_RMS " << nstar_not_var_et_rms << endl ; 
  cout << "@NSTAR_NOT_VAR_AND_RMS_PCT " << 100.*nstar_not_var_et_rms/nstar_tot  << endl ;
  cout << "@NSTAR_NOT_VAR_3_AND_RMS " << nstar_not_var_3_et_rms << endl ; 
  cout << "@NSTAR_NOT_VAR_3_AND_RMS_PCT " << 100.*nstar_not_var_3_et_rms/nstar_tot  << endl ;

  cout << "@NSTAR_CONTAM " << nstar_contam << endl ;
  cout << "@NSTAR_CONTAM_PCT " << 100.*nstar_contam/nstar_tot  << endl ;
  cout << "@NSTAR_NOSTAR " << nstar_nostar << endl ;


  fprintf(file,"@NSTAR_INIT %d\n",nstar_tot);
  fprintf(file,"@NSTAR_VAR %d \n",nstar_var); 
  fprintf(file,"@NSTAR_VAR_PCT %2.1f \n",100.*nstar_var/nstar_tot );
  fprintf(file,"@NSTAR_VAR_3 %d \n",nstar_var_3); 
  fprintf(file,"@NSTAR_VAR_3_PCT %2.1f \n",100.*nstar_var_3/nstar_tot );
  fprintf(file,"@NSTAR_VAR_RMS %d \n",nstar_var_rms); 
  fprintf(file,"@NSTAR_VAR_RMS_PCT %2.1f \n",100.*nstar_var_rms/nstar_tot );
  fprintf(file,"@NSTAR_CONTAM %d \n",nstar_contam);
  fprintf(file,"@NSTAR_CONTAM_PCT %2.1f \n",100.*nstar_contam/nstar_tot );
  fprintf(file,"@NSTAR_NOSTAR %d \n",nstar_nostar);


  fprintf(file,"@NSTAR_NOSTAR_PCT %2.1f \n",100.*nstar_nostar/nstar_tot );
  fprintf(file,"@CHI2PDF_OF_IMAGES %6.6f\n",sum_chi2_of_images/sum_ndf_of_images);
  fprintf(file,"@CHI2PDF_OF_IMAGES_01 %6.6f\n",sum_chi2_of_images_01/sum_ndf_of_images);
  fprintf(file,"@NDF_OF_IMAGES %6.6f\n",sum_ndf_of_images);   
  fprintf(file,"@CHI2PDF_OF_STARS %6.6f\n",sum_chi2_of_stars/sum_ndf_of_stars);
  fprintf(file,"@NDF_OF_STARS %6.6f\n",sum_ndf_of_stars);  
  fclose(file);
  return EXIT_SUCCESS;
}


