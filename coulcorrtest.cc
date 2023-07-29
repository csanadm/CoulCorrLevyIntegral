#include "CoulCorrCalc.h"

int main()
{
  CoulCorrCalc *ccinstance = new CoulCorrCalc();
  double alpha = 1.2;
  double R = 5.3;
  double lambda = 0.8;
  for(double Q=0.001; Q<0.2; Q+=0.001)
  {
	double FullCorrFunc = ccinstance->CorrFuncValue(alpha, R, Q);
	double CoulombCorr = ccinstance->CoulCorrValue(alpha, R, Q);
    double PureCorrFunc = 1 + lambda*(FullCorrFunc/CoulombCorr-1);
	double FullCorrFuncLambda = 1 - lambda + lambda*FullCorrFunc;
    cout << Q << "\t" << PureCorrFunc << "\t" << FullCorrFuncLambda << "\t" << CoulombCorr << endl;
  }
  delete ccinstance;
  return 0;
}

