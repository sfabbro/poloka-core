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
      Stream << "0.." << nx-1 << "," << j << ": ";
      for(unsigned int i=0;i<nx;i++) {
	Stream << " " << float((*this)(i,j));
      }
      Stream << endl;
    }
  }


void Mat::Identity() {
  if(nx!=ny) {
    cout << "Mat::Identity ERROR nx!=ny" <<endl;
    abort();
  }
  Zero();
  for(unsigned int i=0;i<nx;++i)
    (*this)(i,i)=1.;
}

static bool same_size(const Mat& m1, const Mat& m2)
{
  if ((m1.SizeX() == m2.SizeX()) && (m1.SizeX() == m2.SizeY())) return true;
  cout << " matrices have different sizes" << endl;
  abort(); // we should in fact throw an exception.
  return false;
}

static bool same_size(const Vect& v1, const Vect& v2)
{
  if (v1.Size() == v2.Size()) return true;
  cout << " vectors have different sizes" << endl;
  abort(); // we should in fact throw an exception.
  return false;
}

Mat Mat::operator +(const Mat& Right) const
{
  same_size(*this,Right);
  Mat res = (*this);
  res += Right;
  return res;
}

Mat Mat::operator -(const Mat& Right) const
{
  same_size(*this,Right);
  Mat res = (*this);
  res -= Right;
  return res;
}

void Mat::operator +=(const Mat& Right)
{
  same_size(*this,Right);
  double *a = data;
  const double *b = Right.Data();
  unsigned int size = nx*ny;
  for(unsigned int i=0;i<size;++i, ++a, ++b)
    *a += *b;
}

void Mat::operator -=(const Mat& Right)
{
  same_size(*this,Right);
  double *a = data;
  const double *b = Right.Data();
  unsigned int size = nx*ny;
  for(unsigned int i=0;i<size;++i, ++a, ++b)
    *a -= *b;
}

Mat Mat::operator *(const double Right) const 
{
  Mat res = (*this);
  res *= Right;
  return res;
}

Mat operator *(const double Left, const Mat &Right)
{
  Mat res = Right;
  res *= Left;
  return res;
}
 
void Mat::operator *=(const double Right)
{
  unsigned int size = nx*ny;
  double *a = data;
  for(unsigned int i=0;i<size;++i, ++a)
    *a *= Right;
}

Mat Mat::operator *(const Mat& Right) const
{
  if(nx != Right.SizeY()) {
    cout << "Mat::operator *= ERROR nx != Right.SizeY()" << endl;
    abort();
  }
  Mat res(Right.SizeX(),ny);
  
  for(unsigned int i=0;i<res.SizeX();++i) {
    for(unsigned int j=0;j<res.SizeY();++j) {
      for(unsigned int k=0;k<nx;++k) {
	res(i,j) += (*this)(k,j)*Right(i,k);
      }
    }
  }
  return res;
}

Mat Mat::operator *(const Vect& Right) const
{
  return (*this)*Mat(Right);
}

void Mat::operator *=(const Mat& Right)
{
  Mat res = (*this)*Right;
  (*this) = res;
}

Mat::operator double() const
{
  if(nx!=1 || ny !=1) {
    cout << "Mat::operator double() error, nx=ny=1 needed, you have nx=" 
	 << nx <<" ny=" << ny << endl;
    abort();
  }
  return (*this)(0,0);
}

Mat::operator Vect() const
{
  if(nx!=1) {
    cout << "Mat::operator Vect() error, nx=1 needed, you have nx=" 
	 << nx << endl;
    abort();
  }
  Vect res(ny);
  for(unsigned int i=0;i<ny;i++) {
    res(i) = (*this)(0,i);
  }
  return res;
}

Mat Mat::transposed() const {
  Mat res(ny,nx);
  
  for(unsigned int i=0;i<nx;i++) {
    for(unsigned int j=0;j<ny;j++) {
      res(j,i)=(*this)(i,j);
    }
  }
  return res;
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
    Stream << i << ":  " << float((*this)(i)) << endl;
  }
}

Vect Vect::operator *(const double Right) const 
{
  Vect res = (*this);
  res *= Right;
  return res;
}

Vect operator *(const double Left, const Vect &Right)
{
  Vect res = Right;
  res *= Left;
  return res;
}
 
void Vect::operator *=(const double Right)
{
  double *a = data;
  for(unsigned int i=0;i<n;++i, ++a)
    *a *= Right;
}


Vect Vect::operator +(const Vect& Right) const
{
  same_size(*this,Right);
  Vect res = (*this);
  res += Right;
  return res;
}

Vect Vect::operator -(const Vect& Right) const
{
  same_size(*this,Right);
  Vect res = (*this);
  res -= Right;
  return res;
}

void Vect::operator +=(const Vect& Right)
{
  same_size(*this,Right);
  double *a = data;
  const double *b = Right.Data();
  for(unsigned int i=0;i<n;++i, ++a, ++b)
    *a += *b;
}

void Vect::operator -=(const Vect& Right)
{
  same_size(*this,Right);
  double *a = data;
  const double *b = Right.Data();
  for(unsigned int i=0;i<n;++i, ++a, ++b)
    *a -= *b;
}

double Vect::operator *(const Vect& Right) const
{
  same_size(*this,Right);
  double res = 0;
  const double *a = data;
  const double *b = Right.Data();
  for(unsigned int i=0;i<n;++i, ++a, ++b)
    res += (*a)*(*b);
  return res;
}


Mat Vect::transposed() const {
  Mat trans(n,1);
  for(unsigned int i=0;i<n;++i) {
    trans(i,0)=(*this)(i);
  }
  return trans;
}

Vect::operator double() const
{
  if(n!=1) {
    cout << "Vect::operator double() error, n=1 needed, you have n=" 
	 << n << endl;
    abort();
  }
  return (*this)(0);
}

Vect::operator Mat() const {
//Mat Vect::asMat() const {
  Mat mat(1,n);
  for(unsigned int i=0;i<n;++i) {
    mat(0,i)=(*this)(i);
  }
  return mat;
}

