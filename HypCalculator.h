#ifndef _HypCalculator_h_
#define _HypCalculator_h_

#include "functions.h"

class HypCalculator
{
 public:
  HypCalculator();
  HypCalculator(const HypCalculator& myHypCalculator);
  ~HypCalculator();

  complex<double> Gauss_2F1_series0(const complex<double> a, const complex<double> b, const complex<double> c, const double z); // using power series around z=0
  complex<double> Gauss_2F1_1mZ_noint(const complex<double> a, const complex<double> b, const complex<double> c, const double z); // using power series around z=0
  void initialize_abc(const complex<double> a, const complex<double> b, const complex<double> c);
  void initialize_eta(const double _eta);
  complex<double> Gauss_2F1_series0(const double z); // using power series around z=0
  complex<double> Gauss_2F1_1mZ_noint(const double z); // using power series around z=0
  complex<double> Gauss_2F1_1perZ_noint(const double z);
  complex<double> Gauss_2F1_c_1(const double z);

  complex<double> F1_F2_F3(const double x); // Fplus(x-i0) = 2F1(ieta, 1+ieta, 1, x-i0)
 private:
  double EPSILON;
  double EPSILON_ABC;
  double euler_gamma;

  bool is_abc_initialized;
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
