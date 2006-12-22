#include "imagepsf.h"
#include "model.h"

class Accumulator
{
private :
  double sum;
  double sum2;
  int count;

public :
  Accumulator() : sum(0), sum2(0), count(0) {};

  void Clear() {sum = 0; sum2 = 0; count = 0;}
  void AddVal(const double &Val) { sum+= Val; sum2+=Val*Val; count++;}
  int Count() const {return count;}
  double Average() const { return sum/count;}
  double Rms() const 
  {
    if (count == 0 || count == 1) return 0;
    return sqrt((sum2 - sum*sum/double(count))/double(count-1));
  }

  void print(const string &Message, ostream & s=cout) const
  {
    s << Message << " = " << Average() << " (rms = " << Rms() << ")" << endl;
  }
};



#include "reducedimage.h"

int main(int nargs, char **args)
{
  string name(args[1]);
  ReducedImage* ri = new ReducedImage(name);
  ImagePSF imagepsf(*ri,false);
  Frame frame(ri->UsablePart());

  // first check the intrapixel variability
  IntPoint posInImage(frame.Center().x, frame.Center().y);
  Accumulator xc,yc,sum;


  for (int j = -9; j < 10; ++j)
  for (int i=-9; i <10; ++i)
    {
      PixelBlock psf;
      Point delta(i/10.,j/10.);
      ComputePSFPixels(imagepsf, delta+posInImage, posInImage, psf);
      // if (i==0 && j== 0) psf.WriteFits(name+".psf.fits");
      double mx,my,mx2,my2,mxy;
      psf.Moments(mx,my,mx2,my2,mxy);
      xc.AddVal(mx-delta.x);
      yc.AddVal(my-delta.y);
      sum.AddVal(psf.Sum());
    }
  cout << " variations intra-pixel " << endl;
  xc.print("xc");
  yc.print("yc");
  sum.print("integrale");

  xc.Clear();
  yc.Clear();
  sum.Clear();

  int npoints = 20;
  double xstep = frame.Nx()/double(npoints+1);
  double ystep = frame.Ny()/double(npoints+1);
  for (double x = xstep/2 ; x< frame.Nx(); x += xstep)
  for (double y = ystep/2 ; y< frame.Ny(); y += ystep)
    {
      PixelBlock psf;
      Point pos(x,y);
      int ix = int(floor(x+0.5));
      int iy = int(floor(y+0.5));
      IntPoint intPos(ix,iy);
      Point delta(pos-intPos);

      ComputePSFPixels(imagepsf, pos, intPos, psf);
      double mx,my,mx2,my2,mxy;
      double integral = psf.Sum();
      psf.Moments(mx,my,mx2,my2,mxy);
      xc.AddVal(mx*integral-delta.x);
      yc.AddVal(my*integral-delta.y);
      sum.AddVal(psf.Sum());
    }
  cout << " variations dans l'image " << endl;
  xc.print("xc");
  yc.print("yc");
  sum.print("integrale");
      


  return EXIT_SUCCESS;
}

  
