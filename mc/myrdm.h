#ifndef MYRDM__H
#define MYRDM__H



class Image;

// ================ Random Stuff ===================
// random between 0 and 1 using random()
double my_rndm();
// give a seed computed on date to random()
void RandomSeed();
float NormalGaussRand();
// tirage de x et y selon une loi gaussienne centree
// en (0,0) de parametres sig_a, sig_b, theta
// re-pompee de EROS
void NormGauss(double & x, double & y, double sig_a, double sig_b, double theta);
//  returns pixel the value with noise added, pixel in ADU
//(or in photons if Gain = 1)
float Gaussian_Noisy(float pixel, float Gain=1., float ROnoise=0.) ;
void ImgAddNoise(Image & img,float Gain=1., float ROnoise=0.);

#endif /* MYRDM__H */
