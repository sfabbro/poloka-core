/* 
 * $Source: /cvs/snovae/toads/poloka/utils/Attic/fringefinder.cc,v $
 * $Revision: 1.1 $
 * $Author: nrl $
 * $Date: 2004/02/20 10:48:46 $
 * $Name:  $
 */
#include <iostream>

// TOADS
#include "fitsimage.h"
#include "vutils.h"
#include "fringeutils.h"
#include "fitsimagearray.h"


// get cvs version of the code
#define CVSVERSION "$Revision: 1.1 $"


void DumpHelp(const char* programname) {
  cout << " syntax  : " << programname << " <image1> <image2> ...." << endl;
  cout << " options : " << endl;
  cout << "            -f <fringefilename> : name of fringe pattern (default is fringe.fits)" << endl;
  cout << "            -nvec #             : max number of eigenvectors to write (default is the number of input images)" << endl;
  cout << "            -v                  : verbose" << endl;
  cout << "            --help (-h)         : this help" << endl; 
  cout << endl;
}

int main(int argc, char **argv) {

  if (argc < 2 ) {
    DumpHelp(argv[0]);
    exit(1);
  }
  
  vector<string> filelist;
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
    else {
      string filename = argv[i];
      filelist.push_back(filename);
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
    if(!fringes.CheckHeader(tmp_header))
      continue;
    fileoklist.push_back(filelist[i]);
  }
  
  nimages = (int)fileoklist.size();
  if(nimages<2) {
    cout << "You need at least 2 VALID files to build the fringe patterns (" << nimages << ")" << endl;
    DumpHelp(argv[0]);
    exit(1);
  }
  
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
  
  double* spm = FringeUtils::ScalarProductMatrix(fileoklist);
  
  if(verbose) {
    cout << "Scalar products Matrix :" << endl;
    for(int i=0;i<nimages;i++) {
      for(int j=0;j<nimages;j++) {
	cout << spm[i+nimages*j] << " ";
      }
      cout << endl;
    }
    cout << "Diagonalizing ... " << endl;
  }
  
  fringes.AddCommentLine("Scalar products Matrix :");
  for(int i=0;i<nimages;i++) {
    sprintf(comment," ");
    for(int j=0;j<nimages;j++) {
      sprintf(comment+strlen(comment),"%f ",spm[i+nimages*j]);
    }
    fringes.AddCommentLine(comment);
  }
  
  
  // Allocating memory
  double *eigenvalues = new double[nimages];
  double *eigenvectors = new double[nimages*nimages];
  int ierr = DiagonalizeRealSymmetricMatrix(nimages,spm,eigenvectors,eigenvalues);
  if(ierr!=0) {
    cout << "!!! ERROR in DiagonalizeRealSymmetricMatrix (routine eisrs1 of CERNLIB) !!! error code = " << ierr << endl;
    return 0;
  }
  
  int signe = 1;
  for(int i=0;i<nimages;i++)
  if(eigenvectors[i]<0)
    signe = -1;
  

  if(verbose) {  
    cout << "Eigen Values:" << endl;
    for(int im=0;im<nimages;im++) {
      cout << im << " " << eigenvalues[nimages-1-im] << endl;
    }
    cout << "Eigen Vectors:" << endl;
    for(int i=0;i<nimages;i++) {
      for(int j=0;j<nimages;j++) {
	cout << signe*eigenvectors[i+nimages*(nimages-1-j)] << " ";
      }
      cout << endl;
    }
  }
  
  fringes.AddCommentLine("Eigen Values:");
  sprintf(comment," ");
  for(int i=0;i<nimages;i++) {
    sprintf(comment+strlen(comment),"%f ",eigenvalues[nimages-1-i]);
  }
  fringes.AddCommentLine(comment);
  fringes.AddCommentLine("Eigen Vectors:");
  for(int i=0;i<nimages;i++) {
    sprintf(comment," ");
    for(int j=0;j<nimages;j++) {
      sprintf(comment+strlen(comment),"%f ",signe*eigenvectors[i+nimages*(nimages-1-j)]);
    }
    fringes.AddCommentLine(comment);
  }
  fringes.AddCommentLine("----------------------------------------");
  for(int ivec = 0; ivec< nvec ; ivec++) {
    
    if(verbose) {
      cout << "Building vector " << ivec << " ..." << endl;
    }
    
    double scale = 1./sqrt(eigenvalues[nimages-1-ivec]);
    Image vector = FitsImage(fileoklist[0])*(signe*eigenvectors[0+nimages*(nimages-1-ivec)])*scale;
    for(int im=1;im<nimages;im++) {
      vector += FitsImage(fileoklist[im])*(signe*eigenvectors[im+nimages*(nimages-1-ivec)])*scale;
    }
    fringes.Append(vector);
    fringes.AddKey("EIGENVAL",eigenvalues[nimages-1-ivec],"Eigen value of this Eigen vector");  
    sprintf(comment,"fringe%03d",ivec);  
    fringes.AddOrModKey("EXTNAME",comment,"Id of this Eigen vector");
  }
  
  if(verbose) {
    cout << "done" << endl;
  }    

  delete [] spm;
  delete [] eigenvectors;
  delete [] eigenvalues;
  return 0;
}
