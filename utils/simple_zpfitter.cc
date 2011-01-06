#include <iostream>
#include <fstream>
#include <fileutils.h>
#include <dicstar.h>
#include <basestar.h>
#include <vutils.h>


using namespace std;

static void usage(const char *pgname)
{
  cerr << pgname << " <zpstarlist> <keyofreferencemagnitude> ( <output ASCII file> )" << endl;
  cerr << "           (default output file name is simple_zpfitter.dat)" << endl ;
  exit(1);
}

double weightedclipmean(double *values, double *weights, int &nval, double &sigma, const double &k=3.5, const int niter=4) {
  
  double mean,median,sigmat,sumw;
  Dmean_median_sigma(values,nval,mean,median,sigmat);
  double clip = k * sigmat;
  double low = median - clip;
  double high = median + clip;
  int j,n=0;
  int nold = nval;
  if (nval==1)
    {
      sigma = 0.0;
      return values[0];
    }
  for (int i=0; i<niter;i++)
    {
      mean = 0;
      sigma = 0;
      sumw = 0;
      n = 0;
      for (j=0;j<nval;++j)
	{
	  if (values[j] > low && values[j] < high ) 
	    {
	      n++;
	      mean  += values[j]*weights[j];
	      sigma += values[j]*values[j]*weights[j];
	      sumw  += weights[j];
	    }
	}
      mean /= sumw;
      sigma = sigma/sumw - mean*mean;
      if (sigma >0 ) sigma = sqrt(sigma);
      else {sigma = 0; break;}
      clip = k * sigma;
      if (nold == n) break;
      low = mean - clip;
      high = mean + clip;
      nold = n;
    }
  nval = n;
  return mean;
}



//////////////////////////////////////////////////////
// MAIN
//////////////////////////////////////////////////////
int main(int argc, char **argv)
{

  

  if (argc < 3)  {usage(argv[0]);}

  string zpstarlist_filename = argv[1];
  string key = argv[2];

  // leger hack : si band = y -> key = 1 dans zpstarslist
  if (key == "y") {key = "i" ;} 
  
  
  double mag_min = 10; // just to remove dummy stuff
  int  nmeas_min = 10;  // min. number of measurements
  double mag_max = 0;
  double rms_max = 0.02;
  double nsig_clip = 2.5;

  

  switch(key.c_str()[0]) {
  case 'g' :
    mag_max = 20.;
    break;
  case 'r' :
    mag_max = 20.;
    break;
  case 'i' :
    mag_max = 20.;
    break;
  case 'y' :
    mag_max = 20.;
    break;
  case 'z' :
    mag_max = 19.;
    break;
  default :
    cerr << "cannot deal with mag key '" << key << "'" << endl;
    exit(EXIT_FAILURE);
    
  }
  
  cout << "Cuts applied : " << endl;
  cout << "mag_min   " << mag_min   << endl;
  cout << "mag_max   " << mag_max   <<  endl;
  cout << "nmeas_min " << nmeas_min <<  endl;
  cout << "rms_max   " << rms_max   <<  endl;
  


  string outputfilename = "simple_zpfitter.dat";
  if(argc>=4) outputfilename = argv[3];
  
  if(!FileExists(zpstarlist_filename)) {
    cerr << "cant find list " << zpstarlist_filename << endl;
    usage(argv[0]);
  }
  
  DicStarList zpstarlist(zpstarlist_filename);
  int n_stars = zpstarlist.size();
  
  if ( n_stars == 0 ) {
    cerr << "zpstarlist is empty" << endl;
    usage(argv[0]);
  }
  
  double* zps = new double[n_stars];
  double* weights = new double[n_stars];
  
  int n=0;
  for(DicStarCIterator entry=zpstarlist.begin();entry!=zpstarlist.end();++entry) {
    double mag  = (*entry)->getval(key);
    double flux = (*entry)->getval("flux");
    double fluxrms = (*entry)->getval("fluxrms");
    double error = (*entry)->getval("error");
    if(flux<=0) continue;
    if(mag<mag_min) continue;
    if(mag>mag_max) continue;
    if(error<=0) continue;
    if(fluxrms/flux>rms_max) continue; // keep only high S/N SNe
    if((*entry)->getval("nmeas")<nmeas_min) continue;
    zps[n] = mag+2.5*log10(flux);
    weights[n] = 1./(error*error);
    n++;
  }
  int n_selected = n;
  double sigma=0.01;
  double zp;
  zp = weightedclipmean(zps,weights,n,sigma,nsig_clip,10);
  
  /*
  zp = gaussianfit(zps,n,zp,sigma,3.,true);
  cout << zp << " " << sigma << endl;
  if(sigma>0) { 
    zp = gaussianfit(zps,n,zp,sigma,2.,false);
    cout << zp << " " << sigma << endl;
    zp = gaussianfit(zps,n,zp,sigma,1.5,false);
    cout << zp << " " << sigma << endl;
    zp = gaussianfit(zps,n,zp,sigma,1.,false);
    cout << zp << " " << sigma << endl;
  }else{
    sigma=-sigma;
  }
  */
  FILE *file = fopen(outputfilename.c_str(),"w");
  
  fprintf(file,"@PSFZP %6.6f\n",zp);
  fprintf(file,"@PSFZPERROR %6.6f\n",sigma);
  fprintf(file,"@NSTARS_WITH_PHOTOMETRY %d\n",n_stars);
  fprintf(file,"@NSTARS_SELECTED %d\n",n_selected);
  fprintf(file,"@NSTARS_FOR_ZP %d\n",n);
  fprintf(file,"@INPUTSTARLIST %s\n",zpstarlist_filename.c_str());
  fprintf(file,"@MAGMAX %f\n",mag_max);
  fprintf(file,"@FLUXRMSMAX %f\n",rms_max);
  fprintf(file,"@NMEASMIN %d\n",nmeas_min);
  fprintf(file,"@NSIGMACLIP %f\n",nsig_clip);
  fclose(file);
  delete [] zps;
  delete [] weights;
  printf("@PSFZP %6.6f %6.6f %d\n",zp,sigma,n_stars);
  return EXIT_SUCCESS;
}


