
DaoPsfFit::DaoPsfFit()
{
  f = new float[MAXBOX*MAXBOX];
  perr   = opt[PercError].Value()*0.01;
  pkerr  = opt[ProfError].Value()*0.01;
  radius = min(opt[PsfRadius].Value(), float(((psf->npsf-1.)/2. - 1.)/2.) );
}


~DaoPsfFit::DaoPsfFit() 
{
  if (f) delete [] f;
}

//! fit only one SEStar with a DaoPsf
void DaoPsfFit::PeakFit(PhotStar *Star)
{
  x      = Star.x - DAOPHOT_TOADS_SHIFT;
  y      = Star.y - DAOPHOT_TOADS_SHIFT;
  deltax = (x-1) / psf.Xpsf()-1.;
  deltay = (y-1) / psf.Ypsf()-1.;
  scale  = Star->flux;
  sky    = Star->sky;
  lx     = max(1, int(x-radius)+1);
  ly     = max(1, int(y-radius)+1);
  nx     = min(SIZE.ncol, int(x-radius)+1) - lx +1;
  ny     = min(SIZE.nrow, int(y-radius)+1) - ly +1;
 
  x += -lx+1;
  y += -ly+1;
 
  RDARAY("DATA", &lx, &ly, &nx, &ny, &maxbox, f, &status);

  PKFIT(f, &nx, &ny, &MAXBOX, &x, &y, &scale, &sky, &radius,
	&lowbad, &opt[AduHighDatum].Value(), &opt[Gain].Value(), &opt[ReadNoise].Value(), &perr, &pkerr,
	&psf.bright, &psf.type, psf.param,
	&MAXPAR, &psf.npar, psf.table, 
	&MAXPSF, &MAXEXP, &psf.npsf, &psf.nexp,
	&psf.nfrac, &deltax, &deltay, &errmag, &chi,
	&sharp, &niter, &global_sky);

  Star.x       = x + DAOPHOT_TOADS_SHIFT;
  Star.y       = y + DAOPHOT_TOADS_SHIFT;
  Star.flux    = scale;
  Star.varflux = errmag*errmag;
  //Star.chi     = chi;

}

