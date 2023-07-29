#include "CoulCorrCalc.h"

int main()
{
  // Create an instance of the CoulCorrCalc class
  CoulCorrCalc *cccinstance = new CoulCorrCalc();
  // Declare parameter variables
  double alpha = 1.2;
  double R = 5.3;
  double lambda = 0.8;
  // Start loop for calculating the correlation function
  for(double Q=0.001; Q<0.2; Q+=0.001)
  {
    // Coulomb correction
    double CoulombCorr = cccinstance->CoulCorrValue(alpha, R, Q);
    // Pure (no Coulomb) correlation function, with the specified lambda value
    double PureCorrFunc = cccinstance->PureCorrFuncValueLambda(alpha, R, lambda, Q);
    // Full correlation function, including the Coulomb effect, with the specified lambda value
    double FullCorrFuncLambda = cccinstance->FullCorrFuncValueLambda(alpha, R, lambda, Q);
    // Printout
    cout << Q << "\t" << PureCorrFunc << "\t" << FullCorrFuncLambda << "\t" << CoulombCorr << endl;
  }
  // Delete CoulCorrCalc instance and return
  delete cccinstance;
  return 0;
}

