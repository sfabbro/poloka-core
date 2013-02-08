#ifndef  DATDETEC__H 
#define  DATDETEC__H

#include <string>

using namespace std;
enum {T_SExtractor = 0, T_SimpleDetection,T_CvDetection,T_CvSExtractor };

class DataCards;
class FitsHeader;


struct DatDetec{
  int prlevel ;
  int type_det  ;
  double scale_sigma ;
  double n_rad_locmax ;
  int rad_locmax ;
  double n_rad_bad ;
  double rad_bad ;
  double n_rayon;
  double rayon;
  double n_seuil_amp;
  double seuil_amp;
  double n_seuil_flux;
  double seuil_flux;

  double dist_min;
  double sigtonoise_1;
  double sigtonoise_2;
  double sigtonoise;
  double dass_ref_min ;

  double N_anneau ;
  double N_anneau_frac ;

  double frac_fref ;
  double frac_faper_ref ;
  double frac_aref ;

  double fd_ref   ;
  double fd_new1   ;
  double fd_new2   ;
  double fd_new   ;
  double sig_fd_ref  ;
  double sig_fd_new1  ;
  double sig_fd_new2  ;
  double sig_fd_new  ;

  DatDetec();
  DatDetec(const string &DatacardsFileName); 
  void LitDataCard(DataCards & data);
  void Default();
  void Print(ostream & s=cout) const ;

  void ComputeRadius(double seeing);
  void Compute_Seuil(double seeing, double sigma_fd) ;
  void ComputeRadius_Seuil(double seeing, double sigma_fd) ;
 
};
#endif








