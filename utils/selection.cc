/* 
 * $Source: /cvs/snovae/toads/poloka/utils/Attic/selection.cc,v $
 * $Revision: 1.1 $
 * $Author: nrl $
 * $Date: 2004/02/20 10:48:45 $
 * $Name:  $
 */

#include <vector>
#include <string>

#include "fitsimage.h"
#include "fringeutils.h"


// selection : select a list of fits files according to some criteria UNDER DEVELOPPEMENT !!

// brutal check whether pixel is in a dead column or not
bool PixIsInDeadColumn(const Image &image, int x, int y, int nx, int ny, Pixel mean, Pixel sigma) {
  //Pixel nodiff = sigma;
  Pixel maxdiff1 = 3*sigma;
  Pixel maxdiff2 = 4*sigma;
  
  int jump=2;
  if(jump==2) {// go to pair values
    x-=(x%2);
    y-=(y%2);
  }
  Pixel val = fabs(image(x,y)-mean);
  if(val<maxdiff1)
    return false;
  bool is_overnsigma2=0;
  if(val>maxdiff2) {
    if(y==0 || y==(ny-jump) || x==0 || x==(nx-jump))
      return true;
    is_overnsigma2=true;
  }
  
  //int n_the_same = 0;
  //int n_very_different = 0;
  int n_overnsigma1 = 0;
  int n_overnsigma2 = 0;
  Pixel val2=0; 
  for(int j=-1;j<2;j+=1) {
    for(int i=-1;i<2;i+=1) {
      if(i==0 && j==0)
	continue;
      if(i!=0 && j!=0)
	continue;
      
      val2=fabs(image(x+i*jump,y+j*jump)-mean);
      if(val2>maxdiff1)
	n_overnsigma1++;
      if(val2>maxdiff2)
	n_overnsigma2++;
    }
  }
  return ( (n_overnsigma2>1) || (is_overnsigma2 && n_overnsigma1>1));
  //return ( (n_overnsigma2>0) || (n_overnsigma1>1));
  //return false;
}

int MakeDeadImage(string &imagename,FitsImage &deadreference, bool save_diff) {
  FitsImage image(imagename,RO);
  if(! image.IsValid()) {
    cout << "ERROR in MakeDeadImage : " << imagename.c_str() << " is not valid" << endl;
    return -1;
  }
  
  Pixel mean,sigma;
  image.MeanSigmaValue(&mean,&sigma);
  cout << "mean = " << mean << " sigma = " << sigma << endl;
  Pixel min,max;
  image.MinMaxValue(&min,&max);
  cout << "min = " << min << " max = " << max << endl;

  int nx = image.Nx();
  int ny = image.Ny();
  float bscale = (float)image.KeyVal("BSCALE");
  cout << "bscale = " << bscale << endl;
  
  
  Image dead(image.Nx(),image.Ny());
  
  int ndead = 0;
  int ndead_added= 0;
  int ndead_missing= 0;
  int ndead_ref= 0;
  Pixel min_deadref,max_deadref;
  deadreference.MinMaxValue(&min_deadref,&max_deadref);
  cout << "dead ref min max = " << min_deadref << " " <<  max_deadref << endl;  
  float zeroval = (max_deadref+min_deadref)/2.;
  bool isdead = false;
  bool refisdead = false;
  for (int j =2; j<ny-2;j++) { // don't take ito account borders
    for (int i =2; i<nx-2;i++) {
      isdead = PixIsInDeadColumn(image,i,j,nx,ny,mean,sigma);
      refisdead = deadreference(i,j)>zeroval;
      
      if(!isdead && !refisdead) dead(i,j)=0;
      if(isdead && refisdead) dead(i,j)=1;
      if(isdead && !refisdead) dead(i,j)=2;
      if(!isdead && refisdead) dead(i,j)=-1;
      if(isdead) ndead++;
      if(refisdead) ndead_ref++;
      if(isdead && !refisdead) ndead_added++;
      if(!isdead && refisdead) ndead_missing++;
      
    }
  }
  
  //cout << "ndead columns = " << ndead << " (" << (100.*ndead)/(nx*ny) << "%)" << endl;
  cout << "ndead_found= " << ndead << endl;
  cout << "ndead_added= " << ndead_added << endl;
  cout << "ndead_missing= " << ndead_missing << endl;
  cout << "ndead_ref= " << ndead_ref << endl;
  if(save_diff) {
    string deadname = imagename;
    deadname += ".dead";
    FitsImage(deadname,dead);
  }
  return 0;
}



void DumpHelp(const char *programname) {
  cout << programname << "  <FITS Image(s)> " << endl 
       << " option:  -size                   : Checks image size" << endl
       << "          -filter                 : Checks filter" << endl
       << "          -dead <dead reference>  : Dead pixels/columns" << endl
       << "          -savediff               : save difference" << endl
       << endl;
  exit(-1);
}


int main(int argc,char **argv)
{
  if(argc <=1)
    DumpHelp(argv[0]);
  
  bool check_size = false;
  bool check_filter = false;
  bool check_dead = false;
  bool save_diff = false;
  string deadreferencename = "";
  vector<string> filenames;

  cout << "try to do only dead for the moment" << endl;
  for (int i=1; i< argc; i++) {
    if (!strcmp(argv[i],"--help") || !strcmp(argv[i],"-h")) {
      DumpHelp(argv[0]);
    }
    else if (!strcmp(argv[i],"-size")) {
      check_size = true;
      continue;
    }
    else if (!strcmp(argv[i],"-filter")) {
      check_filter = true;
      continue;
    }
    else if (!strcmp(argv[i],"-dead")) {
      check_dead = true;
      deadreferencename = argv[++i];
      continue;
    }
    else if (!strcmp(argv[i],"-savediff")) {
      save_diff = true;
      continue;
    }
    else {
      string filename = argv[i];
      filenames.push_back(filename);
    } 
    
  }
  
  int nimages = filenames.size();
  if( nimages <= 0 ) {
    cout << " No input images !!" << endl;
    DumpHelp(argv[0]);
  }
  if(check_dead) {
    cout << "Checking deads compared with " << deadreferencename.c_str() << endl;
  }
  FitsImage deadreference(deadreferencename,RO);
  if(check_dead && ! deadreference.IsValid())
    return -1;
  for(int im = 0; im< nimages ;im++) {
    if(check_dead) {
      MakeDeadImage(filenames[im],deadreference,save_diff);
    }
  }
  
  return 0;
}
