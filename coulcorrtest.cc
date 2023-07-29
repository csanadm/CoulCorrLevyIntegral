#include "CoulCorrCalc.h"

int main()
{
  CoulCorrCalc *ccinstance = new CoulCorrCalc();
  double alpha = 1.2;
  double R = 5.3;
  double Q = 0.05;
  cerr << ccinstance->CoulCorrValue(alpha, R, Q) << endl;
  delete ccinstance;
  return 0;
}

