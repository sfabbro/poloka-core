#include <iostream>

#include <fitsimage.h>

using namespace std;

int usage(char* pg) {
  cout << pg << ": performs simple image manipulation" << endl; 
  cout << "usage : " << pg << " <image1> {+,-,*,/} <image2> = <result> " << endl;
  return 0;
}

int main(int argc, char** argv) {
  if(argc < 6) {
    return usage(argv[0]);
  }
  
  FitsImage image1(argv[1]);
  if(!image1.IsValid())
    return usage(argv[0]);
  FitsImage image2(argv[3]);
  if(!image2.IsValid())
    return usage(argv[0]);
  Image image3;
  char operation = argv[2][0];

  switch(operation){
  case '+':
    image3 = image1+image2;
    break;
  case '-':
    image3 = image1-image2;
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

  {FitsImage im(argv[5],image1,image3);}
  
  return 0;
}
