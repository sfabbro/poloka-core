#include <iostream>
#include <fstream>
#include <matvect.h>

using namespace std;


int main(int argc, char **argv)
{
  
  Mat m1(4,4);
  m1.Identity();

  Mat m2=m1;
  m2.Zero();
  
  
  for(int i=0;i<4;i++)
    m2(0,i)=1;
  
  m2 -= (2.*m1);
  m1 *= 2;
  Mat m3 = m1*m2;
  cout << m1 << endl;
  cout << m2 << endl;
  cout << m3 << endl;
  
  Vect v1(3);
  v1(0)=1;
  v1(1)=2;
  v1(2)=3;

  cout << v1 << endl;
  Vect v2 = v1*2.;
  cout << v2 << endl;
  cout << v1*v1 << endl;

  Mat m4 = v1;
  cout << m4 << endl;
  
  Mat m5(3,3);
  m5.Identity();
  Vect toto = m5*v1;
  cout << m5*v1 << endl;
  cout << toto << endl;
  cout << "==========" << endl;
  cout << double(v1.transposed()*m5*v1) << endl;
  cout << double(v1.transposed()*v1) << endl;
  cout << Mat(v1)*v1.transposed() << endl;

  Mat m6(10,10);
  for(unsigned int i=0;i<m6.SizeX();++i)
    m6(i,0) = i;
  cout << m6 << endl;
  m6.writeFits("mat.fits");
  Mat m7;
  m7.readFits("mat.fits");
  cout << m7 << endl;

  Mat m8 = m6.SubBlock(2,8,0,3);
  cout << m8 << endl;
  
  return 0;
}

