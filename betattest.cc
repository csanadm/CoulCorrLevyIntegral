#include <iostream>
#include "CoulCorrCalc.h"

using namespace std;

int main()
{
  // Create an instance of the CoulCorrCalc class
  CoulCorrCalc *cccinstance = new CoulCorrCalc();
  // Declare parameter variables
  double alpha = 1.2; //1.387; //1.5; //1.2;
  double RLCMS = 5.3; //7.11; //5.0
  double lambda = 0.8; //0.788; //1.0; //0.8;
  double betaT = 0.97;
  // Start loop for calculating the correlation function
  for(double QLCMS=0.001; QLCMS<0.2; QLCMS+=0.002)
  {
    double QPCMS = QLCMS*sqrt(1.-betaT*betaT/3.);
    // Case when the average RPCMS radius is substituted
    cccinstance->SetBetaT(0);
    double RPCMS = RLCMS*sqrt((1.-2.*betaT*betaT/3.)/(1-betaT*betaT));
    double CoulombCorr1 = cccinstance->CoulCorrValue(alpha, RPCMS, QPCMS);
    double PureCorrFunc1 = cccinstance->PureCorrFuncValueLambda(alpha, RLCMS, 1.0, QLCMS);
    double FullCorrFuncLambda1 = 1. - lambda + lambda*CoulombCorr1*PureCorrFunc1;
    // Case when the correlation function is spherically averaged
    cccinstance->SetBetaT(betaT);
    double CoulombCorr2 = cccinstance->CoulCorrValue(alpha, RLCMS, QPCMS);
    double PureCorrFunc2 = PureCorrFunc1;
    double FullCorrFuncLambda2 = 1. - lambda + lambda*CoulombCorr2*PureCorrFunc2;
    // Printout
    cout << QLCMS << "\t" << 1.-lambda+lambda*PureCorrFunc1 << "\t" << FullCorrFuncLambda1 << "\t" << CoulombCorr1;
    cout <<          "\t" << 1.-lambda+lambda*PureCorrFunc2 << "\t" << FullCorrFuncLambda2 << "\t" << CoulombCorr2 << endl;
    //cout << Q << "\t" << CoulombCorr << endl;
  }
  // Delete CoulCorrCalc instance and return
  delete cccinstance;
  return 0;
}

