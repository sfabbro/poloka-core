#include "cdaophot.h"
#include "daophot.h"
#include "daophotpsf.h"
#include "daopsfphot.h"

DaoPsfPhot::DaoPsfPhot(const ReducedImage& Rim) : dao(Rim), psf(Rim)
{
  data = new float[MAXBOX*MAXBOX];

  perr   = dao.opt[PercError].Value() * 0.01;
  pkerr  = dao.opt[ProfError].Value() * 0.01;
  radius = min(dao.opt[PsfRadius].Value(), float(((psf.npsf-1.)/2. - 1.)/2.) );
  
  dao.Attach(Rim.FitsName());
}

void DaoPsfPhot::PeakFit(PhotStar *Star)
{
  x      = Star->x - DAOPHOT_TOADS_SHIFT;
  y      = Star->y - DAOPHOT_TOADS_SHIFT;
  deltax = (x-1.) / psf.Xpsf()-1.;
  deltay = (y-1.) / psf.Ypsf()-1.;
  scale  = Star->flux;
  sky    = Star->sky;
  lx     = max(1, int(x-radius)+1);
  ly     = max(1, int(y-radius)+1);
  nx     = min(SIZE.ncol, int(x-radius)+1) - lx +1;
  ny     = min(SIZE.nrow, int(y-radius)+1) - ly +1;
 
  x += -lx+1;
  y += -ly+1;
 
  RDARAY("DATA", &lx, &ly, &nx, &ny, &MAXBOX, data, &status);

  PKFIT(f, &nx, &ny, &MAXBOX, &x, &y, &scale, &sky, &radius,
	&lowbad, &dao.opt[AduHighDatum].Value(), &dao.opt[Gain].Value(), &dao.opt[ReadNoise].Value(), &perr, &pkerr,
	&psf.bright, &psf.type, psf.param,
	&MAXPAR, &psf.npar, psf.table, 
	&MAXPSF, &MAXEXP, &psf.npsf, &psf.nexp,
	&psf.nfrac, &deltax, &deltay, &errmag, &chi,
	&sharp, &niter, &global_sky);

  Star->x       = x + DAOPHOT_TOADS_SHIFT;
  Star->y       = y + DAOPHOT_TOADS_SHIFT;
  Star->flux    = scale;
  Star->varflux = errmag*errmag;

}

void DaoPsfFit::operator () (PhotStar *Star)
{
  PeakFit(Star);
}
