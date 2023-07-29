#include "CoulCorrCalc.h"

CoulCorrCalc::CoulCorrCalc()
{
  HypCalculatorInstance = new HypCalculator();
}

CoulCorrCalc::~CoulCorrCalc()
{
  delete HypCalculatorInstance;
}

double CoulCorrCalc::f_s(const double q, const double Rcc, const double alpha)
{
  return exp(-0.5 * pow(fabs(q*Rcc / (HBARC*1000.)), alpha));
}

double CoulCorrCalc::A_1_s_wo_int(const double x, const double k, const double Rcc, const double alpha, const double eta)
{
  double fs1_1 = (f_s(2.*k*x, Rcc, alpha) - f_s(0., Rcc, alpha)) / x;
  double fs1_2 = (f_s(2.*k/x, Rcc, alpha) - f_s(0., Rcc, alpha)) / x;
  complex<double> func_1 = pow(1. + 1./x, 2.*eta*I) * HypCalculatorInstance->F1_F2_F3(1. / (x*x));
  complex<double> func_2 = pow(1. + x, 2.*eta*I) *  HypCalculatorInstance->F1_F2_F3(x*x);
  return -2./eta * imag(fs1_1 * func_1 + fs1_2 * func_2);
}

double CoulCorrCalc::A_2_s_wo_int(const double x, const double k, const double Rcc, const double alpha, const double eta)
{
 double func_1 = sin(eta * log((1. + x)/(1. - x))) / (x * (x + 1.));
 double func_2 = x * x * exp(M_PI * eta) * (f_s(2. * k * x, Rcc, alpha) - f_s(2. * k, Rcc, alpha)) / (1. - x);
 double func_3 = (f_s(2. * k / x, Rcc, alpha) - f_s(2. * k, Rcc, alpha)) / (1. - x);
 return (2. / eta ) * (func_1 * func_2 - func_1 * func_3);
}

double CoulCorrCalc::A1_int(const double k, const double Rcc, const double alpha, const double eta)
{
  complex<double> a(0., eta);
  complex<double> b(1., eta);
  complex<double> c(1., 0.);
  HypCalculatorInstance->initialize_abc(a, b, c);
  HypCalculatorInstance->initialize_eta(eta);

  double error = 0;
  double lowerlim = 0;
  double upperlim = 1;
  int maxiter = 3;
  double tolerance = 1e-2;
  auto func = bind(&CoulCorrCalc::A_1_s_wo_int, this, placeholders::_1, k, Rcc, alpha, eta);

  double result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(func, lowerlim, upperlim, maxiter, tolerance, &error);
  return result;
}

double CoulCorrCalc::A2_int(const double k, const double Rcc, const double alpha, const double eta)
{
  complex<double> a(0., eta);
  complex<double> b(1., eta);
  complex<double> c(1., 0.);
  HypCalculatorInstance->initialize_abc(a, b, c);
  HypCalculatorInstance->initialize_eta(eta);

  double error = 0;
  double lowerlim = 0;
  double upperlim = 1;
  int maxiter = 3;
  double tolerance = 1e-2;
  auto func = bind(&CoulCorrCalc::A_2_s_wo_int, this, placeholders::_1, k, Rcc, alpha, eta);

  double result =  boost::math::quadrature::gauss_kronrod<double, 15>::integrate(func, lowerlim, upperlim, maxiter, tolerance, &error);
  return result;
}

const double CoulCorrCalc::calc_eta(const double k, const double Mc2)
{
  return FINESTRUCTURE_CONSTANT * Mc2 * 0.5 / k;
}

double CoulCorrCalc::CorrFuncValue(const double alpha, const double R, const double Q)
{
  double k = Q*500; //k = Q/2, but in MeV here
  double Rcc = R*pow(2.,1./alpha);
  double eta = calc_eta(k, Mass_Pi*1000.);
  double Gamow = (2. * M_PI * eta) / (exp(2. * M_PI * eta) - 1.);
  double A1 = A1_int(k, Rcc, alpha, eta);
  double A2 = A2_int(k, Rcc, alpha, eta);
  double value = Gamow * (1. + f_s(2.*k, Rcc, alpha) + (eta / M_PI) * (A1 + A2));
  return value;
}

double CoulCorrCalc::CoulCorrValue(const double alpha, const double R, const double Q)
{
  double k = Q*500; //k = Q/2, but in MeV here
  double Rcc = R*pow(2.,1./alpha);
  double eta = calc_eta(k, Mass_Pi*1000.);
  double Gamow = (2. * M_PI * eta) / (exp(2. * M_PI * eta) - 1.);
  double A1 = A1_int(k, Rcc, alpha, eta);
  double A2 = A2_int(k, Rcc, alpha, eta);
  double C0value = 1. + f_s(2.*k, Rcc, alpha);
  double CCvalue = Gamow * (1. + f_s(2.*k, Rcc, alpha) + (eta / M_PI) * (A1 + A2));
  return CCvalue/C0value;
}

