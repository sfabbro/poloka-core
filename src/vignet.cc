#include "vignet.h"
#include "fitsimage.h" //for FitsHeader

Vignet::Vignet()
  : hx(0),hy(0),ic(0),jc(0),dxc(0),dyc(0),istart(0),jstart(0),iend(0),jend(0)
{}

Vignet::Vignet(const double &X, const double &Y, const Image &Source, const int HMaxX, const int HMaxY)
  : Kernel(HMaxX,HMaxY), hx(HMaxX), hy(HMaxY), 
  ic(int(X)), jc(int(Y)),dxc(X-ic), dyc(Y-jc),
  istart(max(0,ic - hSizeX)), jstart(max(0,jc - hSizeY)),
  iend(min(ic+hSizeX+1, Source.Nx())), jend(min(jc+hSizeY+1, Source.Ny()))

{
  if (Nx()*Ny() != (istart-iend)*(jstart-jend)) 
    cerr << " we miss pixels for (" << ic << "," << jc << ")" << endl;
  for (int j=jstart; j<jend; ++j)
    for (int i= istart; i<iend; ++i)
      (*this)(i-ic,j-jc) = Source(i,j);
}

Vignet::Vignet(const double &X, const double &Y, const int HMaxX, const int HMaxY)
  : Kernel(HMaxX,HMaxY), hx(HMaxX), hy(HMaxY), 
  ic(int(X)), jc(int(Y)),dxc(X-floor(X)), dyc(Y-floor(Y)),
  istart(max(0,ic - hSizeX)), jstart(max(0,jc - hSizeY)),
  iend(istart+Nx()), jend(jstart+Ny())
{
}

Vignet::Vignet(const string &FileName)
{
  readFits(FileName);
}

void Vignet::Resize(const int HNewX, const int HNewY)
{
  if ((HNewX > hSizeX) ||(HNewY > hSizeY) || (HNewX < 0) ||(HNewY < 0))
    {
      cerr << " Vignet::Resize can't resize vignet(" 
	   << hSizeX << "," << hSizeY<< ") with HnewX = " 
	   << HNewX << " HNewY = " << HNewY<< endl;
      return;
    }
  hx = HNewX;
  hy = HNewY;
}

void Vignet::Resize(const double &ScaleFactor)
{
  if (ScaleFactor > 1 || ScaleFactor <0)
    {
      cerr << " Vignet::Resize can't resize with scale = " 
	   << ScaleFactor << endl;
      return;
    }
  hx = int(hSizeX*ScaleFactor);
  hy = int(hSizeY*ScaleFactor);
}

void Vignet::Detect(const double &PosThresh, const double &NegThresh)
{
  // sweep a pos. and neg. 3x3 pixels aperture on sub
  Kernel plus(hx,hy);
  Kernel minus(hx,hy);

  // 3 adjacent pixels above threshold
  int npix=3;
  int sump = 0;
  int summ = 0;  

  for (int j=-hy+1; j < hy; ++j)
    for (int i=-hx+1; i < hx; ++i)
      {
	double value = (*this)(i,j);
	if (value > PosThresh) 
	  {
	    int p=0;
	    sump++;
	    for (int k=i-1;k<=i+1;++k)
	      for (int l=j-1;l<=j+1;++l) {if ((*this)(k,l)>PosThresh) ++p;}
	    if (p>=npix){plus(i,j) = value; continue;}
	    }
	if (value < NegThresh) 
	  {
	    summ++;
	    int p=0;
	    for (int k=i-1;k<=i+1;++k)
	      for (int l=j-1;l<=j+1;++l) {if ((*this)(k,l)<NegThresh) ++p;}
	    if (p>=npix){minus(i,j) = value;}
	  }
      }
  if (sump>summ) {plus.MaxPixel(dxc,dyc);}
  else if (summ>sump) {minus.MinPixel(dxc,dyc);}  
}

double Vignet::Aperture(double &VarAper, const double &Radius, 
			const double &VarPix, const double &Sky) const
{
  // bounds check
  if (floor(dxc-Radius)< -hx || floor(dyc-Radius) < -hy || 
      ceil(dyc+Radius) >= hy || ceil(dxc+Radius)  >= hy)
    {
      double minrad =  min(hx, hy);
      cerr << " Vignet::Aperture radius too large :" 
	   << ceil(Radius+max(fabs(dxc),fabs(dyc))) <<" > " << minrad << endl;
      return 0;
    }

  // sum up
  double aper = 0;
  double rad0 = (Radius-0.5)*(Radius-0.5);
  double rad1 = (Radius+0.5)*(Radius+0.5);
  //double rad2 = Radius * Radius;
  for(int j = -hy; j<=hy; ++j)
    {
      double dy2 = (j-dyc)*(j-dyc);
      for(int i = -hx; i<=hx; ++i)
	{
	  double dist2 = (i-dxc)*(i-dxc) + dy2;
	  double mask = 0;
	  if (dist2 < rad1) 
	    {
	      mask = 1;
	      if (dist2 > rad0) mask = 0.5 - sqrt(dist2) + Radius;	      
	    }
	  double pixVal = ((*this)(i,j)-Sky) * mask;
	  VarAper += mask*mask*VarPix + pixVal*mask;
	  aper += pixVal;
	}
    }

  return aper;
}

double Vignet::WeightedAperture(double &VarAper, const Kernel &Model, 
				const double &VarPix, const double &Sky) const
{
  if (Model.HSizeX() < hSizeX || Model.HSizeY() < hSizeY ) 
    {
      cerr << " Vignet::WeigthedAperture : model too small " << endl;
      return -1;
    }

  double norm = Model.sum();
  if (norm == 0) return 0;
  double aper = 0;
  double sumweight = 0;
  for(int j = -hy; j<=hy; ++j)
    for(int i = -hx; i<=hx; ++i)
      {
	double normModel = Model(i,j) / norm;
	double weight = normModel / (VarPix + normModel*norm);
	aper += weight * ((*this)(i,j) - Sky);
	sumweight += weight * normModel;
      }
  VarAper = 1. / sumweight;
  return aper*VarAper;
}


void Vignet::WeightedRecentroid(double &Flux, const double &Sky,
				const double &SigX, const double &SigY)
{
  int iter = 0;
  int maxiter = 10;
  double eps = 0.005;
  double xprev,yprev;
  double alphax = -0.5/(SigX*SigX);
  double alphay = -0.5/(SigY*SigY);

  do
    {
      xprev = dxc;
      yprev = dyc;
      double sumw = 0;
      for (int j=-hx; j<=hy; ++j)
	for (int i=-hx; i<=hx; ++i)
	  {
	    double w = exp((i-dxc)*(i-dxc)*alphax + (j-dyc)*(j-dyc)*alphay);
	    double pixelVal = (*this)(i,j) - Sky;
	    Flux += pixelVal*w;
	    dxc += double(i)*w*pixelVal;
	    dyc += double(j)*w*pixelVal;
	    sumw += w;
	  }

      if (Flux>0)
	{
	  dxc /= Flux;
	  dyc /= Flux;
	}
      Flux /= sumw;
      iter++;
    } 
  while (((fabs(xprev-dxc)>eps) || (fabs(yprev-dyc)>eps)) && (iter < maxiter));
}

void Vignet::readFits(const string &FileName)
{
  Kernel::readFits(FileName);
  FitsHeader head(FileName);
  dxc = head.KeyVal("DXC");
  dyc = head.KeyVal("DYC");
  ic = head.KeyVal("IC");
  jc = head.KeyVal("JC");
  istart = head.KeyVal("ISTART");
  jstart = head.KeyVal("JSTART");
  iend = head.KeyVal("IEND");
  jend = head.KeyVal("JEND");
}

void Vignet::writeFits(const string &FileName) const
{
  Kernel::writeFits(FileName);
  FitsHeader head(FileName,RW);
  head.AddOrModKey("DXC",dxc,"x coordinate from (0,0) on middle");
  head.AddOrModKey("DYC",dyc,"y coordinate from (0,0) on middle");
  head.AddOrModKey("IC",ic,"Center pixel where the vignet was grabbed");
  head.AddOrModKey("JC",jc,"Center pixel where the vignet was grabbed");
  head.AddOrModKey("ISTART",istart,"Down-Left pixel where the vignet was grabbed");
  head.AddOrModKey("JSTART",jstart,"Down-Left pixel where the vignet was grabbed");
  head.AddOrModKey("IEND",iend,"Up-Right pixel where the vignet was grabbed");
  head.AddOrModKey("JEND",jend,"Up-Right pixel where the vignet was grabbed");
}

