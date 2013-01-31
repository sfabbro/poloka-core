#include <iostream>
#include <fstream>

#include <fitsimage.h>

using namespace std;

// this is a very simple program

int usage(char* pg) {
  cout << pg << ": performs simple image manipulation" << endl; 
  cout << endl;
  cout << "usage : " << pg << " <im1> {+,*,/}  <scalar>  <im2> <result> " << endl;
  cout << "ex : " <<  pg  << " toto.fits + -12 tata.fits tutu.fits" << endl;
  cout << " result = im1 {+,*,/} scalar*im2" << endl;
  cout << endl;
  cout << " or : " << pg << " <scalar> <im1> <result> " << endl;
  cout << "ex : " <<  pg  << " 12 toto.fits tutu.fits" << endl;
  cout << " result = scalar*im1" << endl;
  cout << " -f  : to get a Float coded output image " << endl;
  return 0;
} 

int main(int argc0, char** argv0) {

  if (argc0 > 50)
    {
      cout << " wrong number of arguments... maybe you forgot to quote a * ??" << endl;
      return EXIT_FAILURE;
    }

  char* argv[100];
  argv[0] = argv0[0];
  int argc = 1; 

  bool floatOutput = false;

  bool printlist = false ;

  for (int i=1; i<argc0; ++i)
    {
      if (string(argv0[i]) == "-f")
        {
	  floatOutput = true; continue;
	}
      if (string(argv0[i]) == "-P")
        {
	  printlist = true; continue;
	}
      argv[argc] = argv0[i];
      argc++;
    }

  if(argc != 6 && argc !=4) {
    cout << " wrong argc " << argc << endl;
    usage(argv[0]);
    return EXIT_FAILURE;
  }
  if(argc == 4) {
    float scale = atof(argv[1]);
    FitsImage image1(argv[2]);
    if(!image1.IsValid())
      return usage(argv[0]);
    
    Image image2=image1;
    image2*=scale;
    FitsImage im(argv[3],image1,image2);
  }
  if(argc == 6) {
    FitsImage image1(argv[1]);
    if(!image1.IsValid())
      return usage(argv[0]);
    char operation = argv[2][0];
    float scale = atof(argv[3]);
    FitsImage image2(argv[4]);
    if(!image2.IsValid())
      return usage(argv[0]);
    Image image3;
    
    image2*=scale;
    switch(operation){
    case '+':
      image3 = image1+image2;
      break;
    case '*':
      image3 = image1*image2;
      break;
    case '/':
      image3 = image1/image2;
      break;
    default:
      cout << "unknown operator " << operation << endl;
      break;
    }
    FitsImage im(argv[5],image1,image3);
    if (floatOutput) im.SetWriteAsFloat();

    if(printlist)
      {
	ofstream pr("xyf.list");
	pr << "#x : " << endl << "#y : " << endl << "#f : " << endl << "#end " << endl ;
	for(int i = 0 ; i < im.Nx() ; i++)
	  for(int j = 0 ; j < im.Ny() ; j++)
	    pr << i << " " << j << " " << im(i,j) << endl ;
	pr.close() ;
      }
  }
  
  
  return 0;
}
