#ifndef _HypCalculator_h_
#define _HypCalculator_h_

#include "functions.h"

class HypCalculator
{
 public:
  HypCalculator();
  HypCalculator(const HypCalculator& myHypCalculator);
  ~HypCalculator();
  
  // For the below functions, more details can be found in DLMF: $15.8 at https://dlmf.nist.gov/15.8
  
  // Initializetion of eta, useful if 2F1 is calculated for the same eta for many z values, but at the same eta (from which a, b, c are derived)
  void initialize_eta(const double _eta);

  // Calculation of the final quantity, 2F1(ieta, 1+ieta, 1, x-i0), where x is real
  complex<double> F1_F2_F3(const double x);
  
 private:
  // Calculation of the Gaussian 2F1 hypergeometric function 2F1(a,b;c;z), using the power series around z=0
  complex<double> Gauss_2F1_series0(const complex<double> a, const complex<double> b, const complex<double> c, const double z);

  // Calculation of the Gaussian 2F1 hypergeometric function 2F1(a,b;c;1-z), using Gauss_2F1_series0
  complex<double> Gauss_2F1_1mZ_noint(const complex<double> a, const complex<double> b, const complex<double> c, const double z); 
  
  // Calculation of the Gaussian 2F1 hypergeometric function 2F1(a,b;c;z, using the power series around z=0, at the initialized a, b, c values
  complex<double> Gauss_2F1_series0(const double z);
  
  // Calculation of the Gaussian 2F1 hypergeometric function 2F1(a,b;c;1-z), using Gauss_2F1_series0, at the initialized a, b, c values
  complex<double> Gauss_2F1_1mZ_noint(const double z);
  
  // Calculation of the Gaussian 2F1 hypergeometric function 2F1(a,b;c;1/z), using Gauss_2F1_series0, at the initialized a, b, c values
  complex<double> Gauss_2F1_1perZ_noint(const double z);
  
  // Calculation of the Gaussian 2F1 hypergeometric function 2F1(a,b;c;1/z), using l'Hospital version of the power series, at the initialized a, b, c values
  complex<double> Gauss_2F1_c_1(const double z);
  
  // Internal constants and variables
  double EPSILON;
  double EPSILON_ABC;
  double euler_gamma;
  bool is_eta_initialized;
  complex<double> sinabc;
  complex<double> sinab;
  complex<double> Gamma_c;
  complex<double> Gamma_ab1c;
  complex<double> Gamma_ca;
  complex<double> Gamma_cb;
  complex<double> Gamma_a;
  complex<double> Gamma_b;
  complex<double> Gamma_c1ab;
  complex<double> A;
  complex<double> B;
  complex<double> C;
  double eta;
  double pieta_cth_pieta;
  complex<double> sinhpieta_pieta;
  complex<double> tilde_s0eta;
  complex<double> tilde_nu0;
};

#endif // _HypCalculator_h_
