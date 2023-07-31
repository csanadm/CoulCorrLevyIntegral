#include "HypCalculator.h"

HypCalculator::HypCalculator()
{
  EPSILON = 1e-12;
  EPSILON_ABC = 1e-6;
  euler_gamma = 0.5772156649;

  complex<double> complex_zero(0., 0.);

  is_eta_initialized = false;

  sinabc      = complex_zero;
  sinab       = complex_zero;
  Gamma_c     = complex_zero;
  Gamma_ab1c  = complex_zero;
  Gamma_ca    = complex_zero;
  Gamma_cb    = complex_zero;
  Gamma_a     = complex_zero;
  Gamma_b     = complex_zero;
  Gamma_c1ab  = complex_zero;

  A = complex_zero;
  B = complex_zero;
  C = complex_zero;

  pieta_cth_pieta = 0.;
  tilde_s0eta = complex_zero;
  tilde_nu0 = complex_zero;
  sinhpieta_pieta = complex_zero;

  eta = 0.;
};

HypCalculator::HypCalculator(const HypCalculator& myHypCalculator)
{
  is_eta_initialized = myHypCalculator.is_eta_initialized;
}

HypCalculator::~HypCalculator()
{
}

complex<double> HypCalculator::Gauss_2F1_series0(const complex<double> a, const complex<double> b, const complex<double> c, const double z)
{
  complex<double> result(-9999., 0);
  if(abs(z) < 0.9999)
  {
    complex<double> n(0., 0.);
    complex<double> term(1., 0);
    result = term;
    while(abs(term) > EPSILON)
    {
      term = term * z / (n + 1.) * (a + n) * (b + n) / (c + n);
      result += term;
      n += 1. ;
    }
  }
  return result;
}

complex<double> HypCalculator::Gauss_2F1_1mZ_noint(const complex<double> a, const complex<double> b, const complex<double> c, const double z)
{
  complex<double> result(-9999., 0);
  if(abs(1. - z) < 0.9999)
  {
    complex<double> _sinabc = sin(M_PI * (c-a-b));
    if(abs(_sinabc) > EPSILON_ABC)
      result = M_PI / _sinabc * Gamma(c) * ( Gauss_2F1_series0(a, b, a+b+1.-c, 1.-z) / Gamma(a+b+1.-c) / Gamma(c-a) / Gamma(c-b)
                                           - Gauss_2F1_series0(c-a, c-b, c+1.-a-b, 1.-z) / Gamma(a) / Gamma(b) / Gamma(c+1.-a-b) / pow(1.-z, a+b-c));
  }
  return result;
}

void HypCalculator::initialize_eta(const double _eta)
{
  eta = _eta;
  
  complex<double> a(0., eta);
  complex<double> b(1., eta);
  complex<double> c(1., 0.);
  
  sinabc     = sin(M_PI * (c-a-b));  
  sinab      = sin(M_PI * (a-b));  
  Gamma_c    = Gamma(c);
  Gamma_ab1c = Gamma(a+b+1.-c); 
  Gamma_ca   = Gamma(c-a); 
  Gamma_cb   = Gamma(c-b); 
  Gamma_a    = Gamma(a); 
  Gamma_b    = Gamma(b); 
  Gamma_c1ab = Gamma(c+1.-a-b); 
  A = a;
  B = b;
  C = c;
  
  if(eta == 0.)
  {
    pieta_cth_pieta = 1.;
    sinhpieta_pieta = 1.;
  }
  else
  {
    pieta_cth_pieta = M_PI * eta / tanh(M_PI * eta);
    sinhpieta_pieta = sinh(M_PI * eta) / (M_PI * eta);
  }
  
  complex<double> s0(2. * eta * euler_gamma - eta , -1. * pieta_cth_pieta);
  tilde_s0eta = s0;
  complex<double> ieta1(1., eta);
  tilde_s0eta += 2. * eta * digamma(ieta1);
  tilde_nu0 = -1. * eta;
  
  is_eta_initialized = true;
}


complex<double> HypCalculator::Gauss_2F1_series0(const double z)
{
  complex<double> result(-9999., 0);
  if(is_eta_initialized)
  {
    if(abs(z) < 0.9999)
    {
      complex<double> n(0., 0.);
      complex<double> term(1., 0);
      result = term;
      while(abs(term) > EPSILON)
      {
        term = term * z / (n + 1.) * (A + n) * (B + n) / (C + n);
        result += term;
        n += 1. ;
      }
    }
  }
  return result;
}

complex<double> HypCalculator::Gauss_2F1_1mZ_noint(const double z)
{
  complex<double> result(-9999., 0);
  if(is_eta_initialized)
  {
    if(abs(1. - z) < 0.9999)
    {
      if(abs(sinabc) > EPSILON_ABC)
        result = M_PI / sinabc * Gamma_c * ( Gauss_2F1_series0(A, B, A+B+1.-C, 1.-z) / Gamma_ab1c / Gamma_ca / Gamma_cb
                                            -Gauss_2F1_series0(C-A, C-B, C+1.-A-B, 1.-z) / Gamma_a / Gamma_b / Gamma_c1ab / pow(1.-z, A+B-C));
    }
  }
  return result;
}

complex<double> HypCalculator::Gauss_2F1_1perZ_noint(const double z)
{
  complex<double> result(-9999., 0);
  if(is_eta_initialized)
  {
    if(abs(1./z) < 0.9999)
    {
      if(abs(sinab) > EPSILON_ABC)
      {
        result = M_PI / sinab * ( Gauss_2F1_series0(B + 1. -C, B, B + 1. - A, 1./z) / pow(-z, B) / Gamma_a / Gamma_cb - 
                                  Gauss_2F1_series0(A + 1. -C, A, A + 1. - B, 1./z) / pow(-z, A) / Gamma_b / Gamma_ca);
      }
      else
        ;
    }
  }
  return result;
}

complex<double> HypCalculator::Gauss_2F1_c_1(const double z)
{
  complex<double> result(-9999., 0.);
  complex<double> I(0.,1.);
  complex<double> aux(-z, 0.);
  if(is_eta_initialized)
  {
    if(abs(1./z) < 0.9999)
    {
      complex<double> s = - eta * log(aux) + tilde_s0eta; // ide beírtam egy mínuszt
      complex<double> nu = tilde_nu0 / z;
      double n = 0.;
      complex<double> sum = nu * s;
      while(abs(nu * s) > EPSILON)
//      while(abs(sum) > EPSILON)
      {
        nu = nu * (1. + I * eta / (n + 1.)) * (1. + (I * eta - 1.)/(n + 2.)) / z;
        s +=  2. * eta /(n + 1. + I * eta) - eta * (2. * n + 3.)/(n * n + 3. * n + 2.);
        sum += (nu * s);
        n += 1.;
      }   
      result = sinhpieta_pieta * pow(aux, -I * eta) * (1. + sum);
    }
  }
  return result;
}

complex<double> HypCalculator::F1_F2_F3(const double x)
{
  complex<double> result(-9999., 0.);
  if(x >= 0. && is_eta_initialized)
  {
    if(x <= 0.6)
      result = Gauss_2F1_series0(x);
    else if(x <= 1.6)
      result = Gauss_2F1_1mZ_noint(x);
    else
      result = Gauss_2F1_c_1(x);
  }
  return result;

}

