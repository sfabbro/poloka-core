#include <iostream>
#include <matvect.h>

using namespace std;

Mat::Mat(const unsigned int NX, const unsigned int NY) { 
 data=NULL;
 nx=ny=0;
 if(NX<=0 || NY<=0) {
   cout << "Mat::Mat ERROR NX,NY =" << NX << "," << NY <<  endl;
 }else{
   allocate(NX,NY);
 }
}

Mat::Mat(const Mat& other) {
  data=NULL;
  nx=ny=0;
  allocate(other.SizeX(),other.SizeY());
  memcpy(data,other.Data(),sizeof(double)*nx*ny); 
}

void Mat::allocate(const unsigned int NX, const unsigned int NY) {
  if(NX<=0 || NY<=0) {
    cout << "Mat::allocate ERROR NX,NY =" << NX << "," << NY <<  endl;
  }
  if(nx!=NX || ny!=NY) {
    nx=NX; 
    ny=NY; 
    if (data) 
      delete [] data;
    data = new double[nx*ny];
  } 
  Zero();
}

double Mat::operator () (const unsigned int i, const unsigned int j) const {
#ifdef MATVECT_CHECK_BOUNDS
  if (i>=nx || j>=ny || i<0 || j<0) { 
    cout << "Mat::operator () overflow i,j nx,ny " 
	 << i << "," << j << " "
	 << nx << "," << ny << " "
	 << endl;
    abort();
  }
#endif
  return data[i*ny+j];
}

double& Mat::operator () (const unsigned int i, const unsigned int j) {
#ifdef MATVECT_CHECK_BOUNDS
  if (i>=nx || j>=ny || i<0 || j<0) { 
    cout << "Mat::operator () overflow i,j nx,ny " 
	 << i << "," << j << " "
	 << nx << "," << ny << " "
	 << endl;
    abort();
  }
#endif
  return data[i*ny+j];
}

void Mat::dump(ostream& Stream) const {
    for(unsigned int j=0;j<ny;j++) {
      Stream << "0.." << nx-1 << "," << j;
      for(unsigned int i=0;i<nx;i++) {
	Stream << " " << float((*this)(i,j));
      }
      Stream << endl;
    }
  }

//=================================================================

Vect::Vect(const unsigned int N) {
  data = NULL; 
  n=0;
  if(N<=0) {
    cout << "Vect::Vect ERROR N = " << N <<  endl;
  }else{
    allocate(N);
  }
}

Vect::Vect(const Vect& other) {
  data = NULL; 
  n=0;
  allocate(other.Size());
  memcpy(data,other.Data(),sizeof(double)*n); 
}

void Vect::allocate(const unsigned int N) {
  if(N<=0) {
    cout << "Vect::allocate ERROR N = " << N <<  endl;
  }
  if(n!=N) {
    n=N;
    if (data) 
      delete [] data;
    data = new double[n];
  }
  Zero();
};


double Vect::operator () (const unsigned int i) const {
#ifdef MATVECT_CHECK_BOUNDS
  if (i>=n || i<0) {
    cout << "Vec::operator () overflow i,n " 
	 << i << "," << n << endl;
    abort();
  }
#endif
  return data[i];
}

double& Vect::operator () (const unsigned int i) {
#ifdef MATVECT_CHECK_BOUNDS
  if (i>=n || i<0) {
    cout << "Vec::operator () overflow i,n " 
	 << i << "," << n << endl;
    abort();
  }
#endif
  return data[i];
}

void Vect::dump(ostream& Stream) const {
  for(unsigned int i=0;i<n;i++) {
    Stream << i << " " << float((*this)(i)) << endl;
  }
}
