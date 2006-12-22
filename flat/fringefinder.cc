/* 
 * $Source: /cvs/snovae/toads/poloka/flat/fringefinder.cc,v $
 * $Revision: 1.2 $
 * $Author: guy $
 * $Date: 2006/12/22 13:35:40 $
 * $Name:  $
 */
#include <iostream>

// TOADS
#include "fitsimage.h"
#include "matvect.h"
#include "fringeutils.h"
#include "fitsimagearray.h"


// get cvs version of the code
#define CVSVERSION "$Revision: 1.2 $"


void DumpHelp(const char* programname) {
  cout << " syntax  : " << programname << " -p <image1> <image2> -d <dead1> <dead2>...." << endl;
  cout << " options : " << endl;
  cout << "            -f <fringefilename> : name of fringe pattern (default is fringe.fits)" << endl;
  cout << "            -nvec #             : max number of eigenvectors to write (default is the number of input images)" << endl;
  cout << "            -v                  : verbose" << endl;
  cout << "            --help (-h)         : this help" << endl; 
  cout << endl;
  cout << " fringefinder will make a PCA on n input images (with flag -p) using optionnally n dead patterns" << endl;
  cout << " if there are dead patterns, a file fringefilename_dead will be written in addition to the PCA vectors fringefilename" << endl;
  cout << " this dead image is a OR of all input dead patterns " << endl;
  //cout << " " << endl;
  
}

int main(int argc, char **argv) {

  if (argc < 2 ) {
    DumpHelp(argv[0]);
    exit(1);
  }
  
  vector<string> filelist;
  vector<string> deadlist;
  vector<string> fileoklist;
  string fringesname = "fringe.fits";
  
  bool verbose = false;

  int nvec = 1000;

  for (int i=1; i < argc; i++) {
    if (!strcmp(argv[i],"--help") || !strcmp(argv[i],"-h")) {
      DumpHelp(argv[0]);
      exit(1);
    }
    else if (!strcmp(argv[i],"-f")) {
      fringesname = argv[++i];
      continue;
    }
    else if (!strcmp(argv[i],"-nvec")) {
      nvec = atoi(argv[++i]);
      continue;
    }
    else if (!strcmp(argv[i],"-v")) {
      verbose = true;
      continue;
    }
    else if (!strcmp(argv[i],"-p")) {
      int j = i+1;
      while(true) {
	if(j>=argc) {
	  i=j; // this will end the whole loop
	  cout << "end of loop" << endl;
	  break;
	}
	if((argv[j][0]=='-')) {
	  i=j-1;
	  cout << "an option " << argv[j] << endl;
	  break; // we continue the main loop
	}
	string filename = argv[j];
	cout << "adding file " << filename << endl;
	filelist.push_back(filename);
	j++;
      }
      continue;
    }
    else if (!strcmp(argv[i],"-d")) {
      int j = i+1;
      while(true) {
	if(j>=argc) {
	  i=j; // this will end the whole loop
	  break;
	}
	if((argv[j][0]=='-')) {
	  i=j-1;
	  cout << "an option " << argv[j] << endl;
	  break; // we continue the main loop
	}
	string filename = argv[j];
	cout << "adding dead " << filename << endl;
	deadlist.push_back(filename);
	j++;
      }
      continue;
    }
    else{
      cerr << "unknown parameter sequence" << endl;
      DumpHelp(argv[0]);
      exit(1);
      
    }
  } 
  
  int nimages = (int)filelist.size();
  if(nimages<2) {
    cout << "You need at least 2 files to build the fringe patterns (" << nimages << ")" << endl;
    DumpHelp(argv[0]);
    exit(1);
  }
  

  nvec = min(nvec,nimages);
  
  if(verbose) {
    cout << "Creates a FitsImageArray using the first header for reference" << endl;
  }

  // just dump images
//   cout << "Input fringe patterns : " << endl;
//   for(int i=0;i<filelist.size();i++)
//     cout << filelist[i] << endl;
//   cout << "dead patterns : " << endl;
//   for(int i=0;i<deadlist.size();i++)
//     cout << deadlist[i] << endl;
//   cout << "Output Fringe file name = " <<  fringesname << endl;
  

  

  // Creates a FitsImageArray using the first header for reference
  // Uses the constructor with FitsHeader in order to keep all information about instrument, filter ...etc
  // Writes the first image in the first HDU
  FitsHeader first(filelist[0],RO);
  if(!first.IsValid()) {
    cout << "Something was wrong with the first file" << endl;
    cout << "(Try " << argv[0] << " --help)" << endl;
    exit(-1);
  }
  FitsImageArray fringes(fringesname,first,true);
  if(!fringes.IsValid()) {
    cout << "Something was wrong with the first file" << endl;
    cout << "(Try " << argv[0] << " --help)" << endl;
    exit(-1);
  }

  // make the dead image
  string deadname = fringesname+"_dead.gz";
  Image* deadimage = 0;
  if(deadlist.size()>0) {
    cout << "making dead image " << deadname << " = OR of :" << endl;
    for(unsigned int i=0;i<deadlist.size();i++)
      cout << deadlist[i] << endl;
    
    
    deadimage = new Image(fringes.Nx(),fringes.Ny());
    
    Pixel *pix = deadimage->begin();
    Pixel *endpix = deadimage->end();
    
    for(unsigned int i=0;i<deadlist.size();i++) {
      FitsImage image(deadlist[i]);
      if(!image.IsValid()) {
	cout << "dead image " << deadlist[i] << " is not valid!" << endl;
	continue;
      }
      if(deadimage->Nx()!=image.Nx() || deadimage->Ny()!=image.Ny()) {
	cout << "dead image size of " << deadlist[i] << " does not match that of first image" << endl;
	continue;
      }
      pix = deadimage->begin();
      Pixel *inputpix = image.begin();
      for(;pix!=endpix;++pix,++inputpix) {
	if(*inputpix>0.5)
	  *pix=1;
      }
    }
    // writing it as binary
    FitsImage deadfits(deadname,fringes,*deadimage);
    deadfits.ModKey("BITPIX",8);
  }

  
  // Checks all of these Keys to add a file
  fringes.SetCriterion(FitsImageArray::FILTER|FitsImageArray::SIZE|FitsImageArray::CHIP);
  
  if(verbose) {
    cout << "Scanning the list of files to see if everything is ok ..." << endl;
  }
  
  fileoklist.push_back(filelist[0]);
  for(int i=1;i<nimages;i++) {
    if(verbose)
      cout << "Checking " << filelist[i].c_str() << " ..."<< endl;
    FitsHeader tmp_header(filelist[i],RO);
    if(!tmp_header.IsValid())
      continue;
    if(!fringes.CheckHeader(tmp_header)) {
      cerr << filelist[i] << " is not a valid image" << endl;
      continue;
    }
    fileoklist.push_back(filelist[i]);
  }
  
  nimages = (int)fileoklist.size();
  if(nimages<2) {
    cout << "You need at least 2 VALID files to build the fringe patterns (" << nimages << ")" << endl;
    DumpHelp(argv[0]);
    exit(1);
  }
  cout << nvec << " " << nimages << endl;
  nvec = min(nvec,nimages);
  

  // Now add stuff in the header specific to this program
  char comment[256];
  
  fringes.AddCommentLine("----------------------------------------");
  
  sprintf(comment,"Images produced by fringefinder version %s",CVSVERSION);
  cout << comment << endl; fringes.AddCommentLine(comment);
  
  sprintf(comment,"Fringes output file name is %s",fringesname.c_str());
  cout << comment << endl; fringes.AddCommentLine(comment);
  
  sprintf(comment,"%d fringe vectors are produced using the following %d initial patterns:",nvec,nimages);
  cout << comment << endl; fringes.AddCommentLine(comment);
  
  fringes.AddOrModKey("NEXTEND",nvec,"Number of fringe vectors");

  for(int i=0;i<nimages;i++) {
    cout << fileoklist[i].c_str() << endl;
    fringes.AddCommentLine(fileoklist[i].c_str());
  }
  
  if(verbose) {
    cout << "Compute Scalar products Matrix ... " << endl;
  }
  
  Mat scalar_product_matrix = FringeUtils::ScalarProductMatrix(fileoklist,deadimage);
  if(scalar_product_matrix.SizeX()==0 || (scalar_product_matrix.SizeX()!=scalar_product_matrix.SizeY()))
    return 0;
  
  if(verbose) {
    cout << "Scalar products Matrix :" << endl;
    cout << scalar_product_matrix << endl;
    cout << endl;
    cout << "Diagonalizing ... " << endl;
  }
  
  fringes.AddCommentLine("Scalar products Matrix :");
  for(int i=0;i<nimages;i++) {
    sprintf(comment," ");
    for(int j=0;j<nimages;j++) {
      sprintf(comment+strlen(comment),"%f ",scalar_product_matrix(i,j));
    }
    fringes.AddCommentLine(comment);
  }
  
  
  
  Vect eigenvalues;
  Mat eigenvectors;
  int ierr = DiagonalizeRealSymmetricMatrix(scalar_product_matrix,eigenvectors,eigenvalues);
  if(ierr!=0) {
    cout << "!!! ERROR in DiagonalizeRealSymmetricMatrix (routine eisrs1 of CERNLIB) !!! error code = " << ierr << endl;
    return 0;
  }
  
  int signe = 1;
  for(int i=0;i<nimages;i++)
    if(eigenvectors(i,0)<0)
      signe = -1;
  

  if(verbose) {  
    cout << "Eigen Values:" << endl;
    cout << eigenvalues << endl;
    cout << "Eigen Vectors:" << endl;
    cout << eigenvectors << endl;
  }
  
  fringes.AddCommentLine("Eigen Values:");
  sprintf(comment," ");
  for(int i=0;i<nimages;i++) {
    sprintf(comment+strlen(comment),"%f ",eigenvalues(nimages-1-i));
  }
  fringes.AddCommentLine(comment);
  fringes.AddCommentLine("Eigen Vectors:");
  for(int i=0;i<nimages;i++) {
    sprintf(comment," ");
    for(int j=0;j<nimages;j++) {
      sprintf(comment+strlen(comment),"%f ",signe*eigenvectors(i,(nimages-1-j)));
    }
    fringes.AddCommentLine(comment);
  }
  fringes.AddCommentLine("----------------------------------------");
  for(int ivec = 0; ivec< nvec ; ivec++) {
    
    if(verbose) {
      cout << "Building vector " << ivec << " ..." << endl;
    }
    
    double scale = 1./sqrt(eigenvalues(nimages-1-ivec));
    Image vector = FitsImage(fileoklist[0])*(signe*eigenvectors(0,nimages-1-ivec))*scale;
    for(int im=1;im<nimages;im++) {
      vector += FitsImage(fileoklist[im])*(signe*eigenvectors(im,nimages-1-ivec))*scale;
    }
    if(verbose) {
      cout << "normalize this vector" << endl;
    }
    Pixel mean,sigma;
    vector.SkyLevel(&mean,&sigma);
    vector -= mean;
    vector /= sigma;
    if(deadimage) {
      vector *= (*deadimage-1); // set to zero dead pixels
    }
    if(verbose) {
      cout << "write it" << endl;
    }
    fringes.Append(vector);
    fringes.AddKey("EIGENVAL",eigenvalues(nimages-1-ivec),"Eigen value of this Eigen vector");  
    sprintf(comment,"fringe%03d",ivec);  
    fringes.AddOrModKey("EXTNAME",comment,"Id of this Eigen vector");
  }
  
  if(verbose) {
    cout << "done" << endl;
  }    
  
  return 0;
}
