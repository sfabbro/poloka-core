#include <iostream>
#include <matvect.h>
#include <fitsio.h>

using namespace std;

#define dsinv dsinv_
#define dfact dfact_
#define dfinv dfinv_
#define dfeqn dfeqn_
#define eisrs1 eisrs1_

// using cernstuff (from cernlib)
extern "C" 
{
  void dsinv(int *N, double *A, int *IDIM, int *IFAIL);
  void dfact(int *N, double *A, int *idim, double *r, int *ifail, double *Det, int *jfail);
  void dfinv(int *n, double *A, int *idim, double *r);
  void dfeqn(int *n, double *a, int *idim, double *r, int *k, double *b);
  void eisrs1(int *NM,int *N,double *AR,double *WR,double *ZR,int *IERR,double *WORK);
}

// using lapack 
extern "C" {
  void dposv_(char *, int *, int *, double *, int *, double *, int *, int *);
  void dpotri_(char *, int *, double *, int *, int *);
};

int cholesky_solve(Mat &A, Vect &B, char* UorL)
{
#ifdef FNAME
  cout << " >  cholesky_solve" << endl;
#endif  

  if(A.SizeX() != A.SizeY() || A.SizeX()<=0) {
    cout << "error in matvect.cc, cholesky_solve Matrix A must be symmetric and you have nx,ny = "
	 << A.SizeX() << "," << A.SizeY() << endl;
    abort();
  }
  if(B.Size() != A.SizeX()) {
    cout << "error in matvect.cc, cholesky_solve Vector B must have a dimention B.Size()=A.SizeY() and you have B.Size(),A.SizeY() = "
	 << B.Size() << "," << A.SizeY() << endl;
    abort();
  }

  double *a = A.NonConstData();
  double *b = B.NonConstData();
  
  int n = B.Size();
  int nhrs = 1, info = 0;

  dposv_(UorL, &n, &nhrs, a, &n, b, &n, &info);

  if (info != 0) 
    cerr << " cholesky_solve(" << a << "," << b << "," << n
	 << ") : Warning: Cholesky factorization failure . info =" 
	 << info <<  " (>0 is not pos.def)" << endl;

  return info;
}

int cholesky_invert(Mat &A, char* UorL)
{  
  if(A.SizeX() != A.SizeY() || A.SizeX()<=0) {
    cout << "error in matvect.cc, cholesky_invert Matrix A must be symmetric and you have nx,ny = "
	 << A.SizeX() << "," << A.SizeY() << endl;
    abort();
  }
  
  int info = 0;
  double *a = A.NonConstData();
  int n = A.SizeX();

  //  Now invert using the factorization done in dposv_

  dpotri_(UorL, &n, a, &n, &info);

  if (info != 0) 
    cerr << " cholesky_invert(" << a << "," << n
	 << ") : Warning: inversion failure . info =" 
	 << info <<  " (>0 is not pos.def)" << endl;

  return info;
}



//==================================================================




Mat::Mat(const unsigned int NX, const unsigned int NY) { 
 data=NULL;
 nx=ny=0;
 if(NX<0 || NY<0) {
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
  if(NX<0 || NY<0) {
    cout << "Mat::allocate ERROR NX,NY =" << NX << "," << NY <<  endl;
  }
  if(nx!=NX || ny!=NY) {
    nx=NX; 
    ny=NY; 
    if (data) 
      delete [] data;
    if(nx*ny>0)
      data = new double[nx*ny];
    else
      data = 0;
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
  return data[i+j*nx];
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
  return data[i+j*nx];
}

int Mat::writeASCII(ostream& Stream) const {
  Stream << nx << " " << ny << endl;
  for(unsigned int j=0;j<ny;j++) {
    //Stream << "0.." << nx-1 << "," << j << ": ";
    for(unsigned int i=0;i<nx;i++) {
      Stream << " " << float((*this)(i,j));
    }
    Stream << endl;
  }
  return 0;
}

int Mat::readASCII(istream& Stream) {
  unsigned int  fnx,fny;
  
  Stream >> fnx >> fny;
  allocate(fnx,fny);
  
  double val;
  for(unsigned int j=0;j<ny;j++) {
    for(unsigned int i=0;i<nx;i++) {
      Stream >> val;
      (*this)(i,j)=val;
    }
  }
  return 0;
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
  
  for(unsigned int j=0;j<res.SizeY();++j) {
    for(unsigned int i=0;i<res.SizeX();++i) {  
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

Mat & Mat::operator =(const Mat& Right){
  allocate(Right.SizeX(),Right.SizeY());
  memcpy(data,Right.Data(),nx*ny*sizeof(double));
  return (*this);
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

Mat Mat::SubBlock
(unsigned int x_min,unsigned int x_max,unsigned int y_min,unsigned int y_max) const {
  if( x_min<0 || x_max >= nx || y_min <0 || y_max>= ny ) {
    cout << "Mat::SubBlockFromIndexes ERROR, trying to get a sub-matrix with indices" << endl;
    cout << "x_min,x_max,y_min,y_max = "
	 << x_min << ","
	 << x_max << ","
    	 << y_min << ","
	 << y_max << endl;
    cout << "nx,ny = "<< nx << "," << ny << endl;
    abort();
  }
  unsigned int nx_new = (x_max-x_min+1);
  unsigned int ny_new = (y_max-y_min+1);
  Mat res(nx_new,ny_new);
  for(unsigned int j=0;j<ny_new;++j)
    for(unsigned int i=0;i<nx_new;++i)
      res(i,j) = (*this)(i+x_min,j+y_min);
  return res;
}

Mat Mat::WithoutRows(unsigned int y_min,unsigned int y_max) const {
  if( y_min <0 || y_max < y_min || y_max >= ny ) {
    cout << "Mat::WithoutRows ERROR " << endl;
    abort();
  } 
  unsigned int nrows = y_max-y_min+1;
  Mat res(nx,ny-nrows);
  for(unsigned int j = 0 ;j<y_min ;j++)
    for(unsigned int i=0;i<nx;i++)
      res(i,j)=(*this)(i,j);
  for(unsigned int j = y_max+1 ;j<ny ;j++)
    for(unsigned int i=0;i<nx;i++)
      res(i,j-nrows)=(*this)(i,j);
  return res;
}

Mat Mat::WithoutColumns(unsigned int x_min,unsigned int x_max) const {
  if( x_min <0 || x_max < x_min || x_max >= nx ) {
    cout << "Mat::WithoutColumns ERROR " << endl;
    abort();
  } 
  unsigned int ncols = x_max-x_min+1;
  Mat res(nx-ncols,ny);
  for(unsigned int i=0;i<x_min;i++)
    for(unsigned int j = 0 ;j<ny ;j++)
      res(i,j)=(*this)(i,j);
  for(unsigned int i=x_max+1;i<nx;i++)
    for(unsigned int j = 0 ;j<ny ;j++)
      res(i-ncols,j)=(*this)(i,j);
  return res;
}





int Mat::readFits(const string &FitsName) {
  int status = 0;
  fitsfile *fptr = 0;
  fits_open_file(&fptr,FitsName.c_str(),0,&status);
  if (status)
   {
     cerr << " when opening file : " << FitsName ;
     cerr << " we got these successive errors:" << endl;
     fits_report_error(stderr, status);
     return status;
   }
  // first get the size of the image/matrix
  status=0;
  char value[256];
  fits_read_key(fptr, TSTRING, "NAXIS1", &value, NULL, &status);
  if (status)
    {
      cerr << " when reading content of : " << FitsName ;
      cerr << " we got these successive errors:" << endl;
      fits_report_error(stderr, status);
      return status;
    }
  int n1 = atoi(value);
  fits_read_key(fptr, TSTRING, "NAXIS2", &value, NULL, &status);
  if (status)
    {
      cerr << " when reading content of : " << FitsName ;
      cerr << " we got these successive errors:" << endl;
      fits_report_error(stderr, status);
      return status;
    }
  int n2 = atoi(value);
  if(n1<=0 || n2<=0) {
    cout << "Mat::readFits error NAXIS1,NAXIS2 = " << n1 << "," << n2 << endl;
    return -1;
  }
  allocate(n1,n2);

  status = 0;
  float nullval = 0;
  int anynull;
  fits_read_img(fptr, TDOUBLE, 1, nx*ny, &nullval,  data, &anynull, &status);
  if (status)
    {
      cerr << " when reading content of : " << FitsName ;
      cerr << " we got these successive errors:" << endl;
     fits_report_error(stderr, status);
     return status;
    }

  status = 0;
  fits_close_file(fptr, &status);
  if (status)
    {
     cerr << " when closing file : " << FitsName ;
     cerr << " we got these successive errors:" << endl;
     fits_report_error(stderr, status);
   }
  return status;
}


int Mat::writeFits(const string &FitsName) const {
  // we cannot use fitsimage cause we want to write it in double precision
  
  int status = 0;
  fitsfile *fptr = 0;

  remove(FitsName.c_str());

  // enum FitsFileMode {RO = 0, RW = 1};
  fits_create_file(&fptr,FitsName.c_str(),  &status);
  if (status)
   {
     cerr << " when creating file : " << FitsName ;
     cerr << " we got these successive errors:" << endl;
     fits_report_error(stderr, status);
     return status;
   }
  // set a minimal header
  status = 0;
  long naxes[2];
  naxes[0]=nx;
  naxes[1]=ny;
  fits_write_imghdr(fptr,-64,2,naxes,&status);
  if (status)
   {
     cerr << " when writing minial header  ";
     cerr << " we got these successive errors:" << endl;
     fits_report_error(stderr, status);
     return status;
   }

  // say to cfitsio to take into account this new BITPIX
  status = 0;
  fits_flush_file(fptr,&status);
  if (status)
   {
     cerr << " when flushing  ";
     cerr << " we got these successive errors:" << endl;
     fits_report_error(stderr, status);
     return status;
   }
  status = 0;
  fits_write_img(fptr, TDOUBLE, 1, nx*ny, data, &status);
  if (status)
   {
     cerr << " when writing data ";
     cerr << " we got these successive errors:" << endl;
     fits_report_error(stderr, status);
     return status;
   }
  status = 0;
  fits_close_file(fptr, &status);
  if (status)
    {
     cerr << " when closing file : " << FitsName ;
     cerr << " we got these successive errors:" << endl;
     fits_report_error(stderr, status);
   }
  return status;
}

void Mat::Symmetrize(const char* UorL) {
  
  if(nx!=ny) {
    cout << "Mat::Symmetrize ERROR nx!=ny nx,ny = " << nx << "," << ny << endl;
    abort();
  }
  
  
  for(unsigned int j=0;j<ny;j++)
    for(unsigned int i=j+1;i<nx;i++)
      if(UorL[0]=='L') { // x >= y
	(*this)(j,i)=(*this)(i,j);
      }else{
	(*this)(i,j)=(*this)(j,i);
      }
} 

int Mat::SymMatInvert() {
  if(nx!=ny) {
    cout << "Mat::Symmetrize ERROR nx!=ny nx,ny = " << nx << "," << ny << endl;
    abort();
  }
  double *A = data;
  int n = nx;
   int ierr;
  dsinv(&n,A,&n,&ierr);
  return (!ierr);
}

//=================================================================

Vect::Vect(const unsigned int N) {
  data = NULL; 
  n=0;
  if(N<0) {
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
  if(N<0) {
    cout << "Vect::allocate ERROR N = " << N <<  endl;
  }
  if(n!=N) {
    n=N;
    if (data) 
      delete [] data;
    if(n>0)
      data = new double[n];
    else
      data = 0;
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

Vect & Vect::operator =(const Vect& Right){
  allocate(Right.Size());
  memcpy(data,Right.Data(),n*sizeof(double));
  return (*this);
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

