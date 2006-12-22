#include <iostream>

#include <fitsimage.h>

using namespace std;

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
  return 0;
} 

int main(int argc, char** argv) {
  if(argc != 6 && argc !=4) {
    return usage(argv[0]);
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
  }
  
  
  return 0;
}
