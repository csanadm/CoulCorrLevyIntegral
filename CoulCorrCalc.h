#ifndef _CoulCorrCalc_h_
#define _CoulCorrCalc_h_

//#include <algorithm>
//#include <functional>
//#include <iostream>
//#include <cmath>
//#include <limits>
#include "basics.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include "HypCalculator.h"
#include "CoulCorrCalc.h"

class CoulCorrCalc
{
 public:
  CoulCorrCalc();
  ~CoulCorrCalc();
  double CorrFuncValue(const double alpha, const double R, const double Q);
  double CoulCorrValue(const double alpha, const double R, const double Q);

private:
  double f_s(const double q, const double R, const double alpha);
  double A_1_s_wo_int(const double x, const double k, const double R, const double alpha, const double eta);
  double A_2_s_wo_int(const double x, const double k, const double R, const double alpha, const double eta);
  double A1_int(const double k, const double R, const double alpha, const double eta);
  double A2_int(const double k, const double R, const double alpha, const double eta);
  const double calc_eta(const double k, const double Mc2); // k: MeV/c, Mc2: MeV, Mc2 is the mass of the particle, not the reduced mass!!!

  HypCalculator* HypCalculatorInstance;
  
  complex<double> I = complex<double>(0., 1.);

  double xlim1 = 0.05;
  double xlim2 = 0.95;

  int NQ = 399;
  double dQ = 0.001;
  double Qmin = 0.001;
  double Qmax = 0.399;

  static const int NA = 15;
  static const int NR = 10;
  static constexpr double alpha_values[NA] = {0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};
  static constexpr double R_values[NR]     = {3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.,11.,12.};
};

#endif // _CoulCorrCalc_h_
