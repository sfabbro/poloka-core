#include <fstream>
#include <iomanip>
#include "fileutils.h"
#include "datacards.h"
#include "datdetec.h"

#ifndef M_PI
#define     M_PI            3.14159265358979323846  /* pi */
#endif


DatDetec::DatDetec(const string &DatacardsFileName)
{
  if (FileExists(DatacardsFileName))
    {
      DataCards cards(DatacardsFileName);
      LitDataCard(cards);
    }
  else
    {
      cerr << " Cannot read " << DatacardsFileName << endl;
      Default();
    }
}


DatDetec::DatDetec()
{
  prlevel = 0 ;
  type_det = 0 ;
  // detection avec SExtractor
  scale_sigma = -1 ;
  // detection avec simple detection
  n_rad_bad  = 0. ;
  rad_bad  = 0. ;
  n_rad_locmax = 0. ;
  rad_locmax = 0 ;
  n_rayon = 0 ;
  rayon = 0 ;
  n_seuil_amp = 0. ;
  seuil_amp = 0. ;
  seuil_flux= 0. ;
  n_seuil_flux= 0. ;


  dist_min = 0 ;

  // ct signal-to-noise
  sigtonoise = 0. ;
  sigtonoise_1 = 0. ;
  sigtonoise_2 = 0. ;

  //cut comparaison avec ref
  dass_ref_min = 0. ;
  frac_fref = 0. ;
  frac_faper_ref = 0. ;
  frac_aref = 0. ;

  // cut flux anneau
  N_anneau = 0. ;
  N_anneau_frac = 0. ;


  // donnees images et sous-images
  fd_ref = 0.  ;
  fd_new1 = 0.  ;
  fd_new2 = 0.  ;
  fd_new = 0.  ;
  sig_fd_ref = 0. ;
  sig_fd_new1 = 0. ;
  sig_fd_new2 = 0. ;
  sig_fd_new = 0. ;
}

void
DatDetec::Print(ostream & s) const 
{
  s  << "*****DatDetec******" << endl ;
  s << "prlevel: " <<  prlevel  << endl ;
  s << "Detecteur (0=SE, 1=SD 2=CvDetection): " <<   type_det  << endl ;
    s  << "---- Detection SExtractor" << endl ;
s << "facteur par lequel on x le sigma pr baisser seuil detection: " <<  scale_sigma  << endl ;
    s  << "---- Simple Detection " << endl ;
    s  << "Distance min. de pixel bad (n seeing): " << n_rad_bad << endl ;
    s  << "Distance min. de pixel bad : " << rad_bad << endl ;    
    s  << "Rayon ds lequel max local (n seeing): " << n_rad_locmax << endl ;
    s  << "Rayon ds lequel max local : " << rad_locmax << endl ;   
    s  << "Rayon flux & bary-centre  (n seeing): " << n_rayon << endl ;
    s  << "Rayon flux & bary-centre : " << rayon  << endl ;
    s  << "Seuil sur amplitude (n sigma): " << n_seuil_amp  << endl ;
    s  << "Seuil sur amplitude : " << seuil_amp  << endl ;
    s  << "Seuil sur flux (n sigma): " << n_seuil_flux  << endl ;
    s  << "Seuil sur flux : " << seuil_flux  << endl ;

    s  << "---- Cuts" << endl ;    
    s << "dist min. entre les 3 detections  : " <<  dist_min  << endl ;
    s << "rayon calcul noise / seeing   : " << n_rayon  << endl ;
    s << "rayon calcul noise   : " << rayon << endl ;
    s << "S/N sous-image 1 : " << sigtonoise_1 << endl ;
    s << "S/N sous-image 2 : " << sigtonoise_2 << endl ;
    s << "S/N image : " << sigtonoise << endl ;
    s << "dist a objet ref min   : " << dass_ref_min << endl ;
    s << "flux / flux_ref   : " << frac_fref << endl<< endl ;
    s << "pix max obj. / pix max obj. ref   : " << frac_aref << endl<< endl ;
    s << "flux aper / flux aper ref   : " << frac_faper_ref << endl<< endl ;
    s << "cut sur flux ds un anneau / bruit    : " << N_anneau << endl<< endl ;
    s << "cut sur flux ds un anneau / flux   : " << N_anneau_frac << endl<< endl ;
    s << " Images " << endl ;
    s << " fond ref    : " << fd_ref << endl  ;
    s << " sigma ref   : " << sig_fd_ref<< endl  ;
    s << " fond new 1  : " << fd_new1 << endl  ;
    s << " sigma new 1 : " <<  sig_fd_new1 << endl ;
    s << " fond new 2  : " << fd_new2 << endl  ;
    s << " sigma new 2 : " <<  sig_fd_new2 << endl ;
    s << " fond new    : " << fd_new << endl  ;
    s << " sigma new   : " <<  sig_fd_new << endl ;
}



void
DatDetec::Default()
{

  prlevel = 0 ;
  type_det = 1 ;
  scale_sigma = 0.4 ;
  dist_min = 2. ;
  n_rad_bad = 1. ;
  rad_bad = 2. ;
  n_rad_locmax= 1. ;
  rad_locmax = 1 ;
  n_rayon = 2.5  ;
  rayon = 5  ;
  n_seuil_amp = 1.5 ;
  seuil_amp = 0. ;
  n_seuil_flux = 0. ;
  seuil_flux = 0. ;
  
  sigtonoise = 3.  ;
  sigtonoise_1 = 3.  ;
  sigtonoise_2 = 3.  ;
  dass_ref_min = 3. ;

  N_anneau = 200 ;
  N_anneau_frac = 200 ;


  frac_fref = 0. ;
  frac_aref = 0.15 ;
  frac_faper_ref = 0. ;
}


void
DatDetec::LitDataCard(DataCards & data)
{
  prlevel = data.IParam("DETEC_PRLEVEL");
  type_det = data.IParam("DETEC_TYPE");
  scale_sigma = data.DParam("DETEC_SCALESIG");
  dist_min = data.DParam("DETEC_DMIN");
  n_rad_bad = data.DParam("DETEC_RBAD_NSEEING");
  rad_bad = data.DParam("DETEC_RBAD");
  n_rad_locmax = data.DParam("DETEC_RMAXLOC_NSEEING");
  rad_locmax = data.IParam("DETEC_RMAXLOC");
  n_rayon = data.DParam("DETEC_RAD_NSEEING");
  rayon = data.DParam("DETEC_RAD");
  n_seuil_amp  = data.DParam("DETEC_SAMP_NSIG");
  seuil_amp  = data.DParam("DETEC_SAMP");
  seuil_flux  = data.DParam("DETEC_SFLUX");
  n_seuil_flux  = data.DParam("DETEC_SFLUX_NSIG");

  sigtonoise_1  = data.DParam("DETEC_SN_1");
  sigtonoise_2  = data.DParam("DETEC_SN_2");
  sigtonoise  = data.DParam("DETEC_SN");
  dass_ref_min  = data.DParam("DETEC_DREF");

  N_anneau   = data.DParam("DETEC_N_RING");
  N_anneau_frac   = data.DParam("DETEC_N_RING_PRCT");

  frac_fref  = data.DParam("DETEC_FRAC_FREF");
  frac_aref  = data.DParam("DETEC_FRAC_AREF");
  frac_faper_ref  = data.DParam("DETEC_FRAC_FAREF");
}



void
DatDetec::ComputeRadius(double seeing)
{
  if ( n_rayon > 0 )
    rayon = n_rayon * seeing  ;
  if ( n_rad_bad > 0 )
    rad_bad = n_rad_bad * seeing  ;
  //  if ( n_rad_locmax > 0 )   rad_locmax = n_rad_locmax * seeing  ;
}

#include <math.h>

void
DatDetec::Compute_Seuil(double seeing, double sigma_fd)
{
  
  if ( n_seuil_amp > 0 )
     seuil_amp = n_seuil_amp *  sigma_fd ;
  seuil_flux = 2. * M_PI * seeing * seeing * seuil_amp ;
  double NN = 0 ;
  if ( n_rayon > 0 )
    NN = n_rayon ;
  else
    if ( seeing > 1.e-10)
      NN = rayon / seeing ;
  double correc = 1. - exp ( -0.5 * NN * NN ) ;
  seuil_flux *=  correc ;
} 

void
DatDetec::ComputeRadius_Seuil(double seeing, double sigma_fd)
{
  ComputeRadius(seeing);
  Compute_Seuil(seeing, sigma_fd);
}




