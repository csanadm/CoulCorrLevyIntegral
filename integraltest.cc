#include <iostream>
#include "CoulCorrCalc.h"

using namespace std;

int main()
{
  // Create an instance of the CoulCorrCalc class
  CoulCorrCalc *cccinstance = new CoulCorrCalc();
  // Declare parameter variables
  double alpha = 1.4;
  double R = 6.0;
  double lambda = 1.0;
  double Q = 0.20;
  // Calculate correlation function with various integral settings
  for(int NMaxIter=1; NMaxIter<=12; NMaxIter++)
    for(double epsTolerance=1e-3; epsTolerance>1e-15; epsTolerance/=10)
    {
      cccinstance->SetIntegrationProperties(NMaxIter,epsTolerance);
      // Printout
      cout << NMaxIter << "\t" << -log10(epsTolerance) << "\t" << cccinstance->FullCorrFuncValueLambda(alpha, R, lambda, Q)-1 << "\t" << cccinstance->GetNFuncCalls() << endl;
    }
  // Delete CoulCorrCalc instance and return
  delete cccinstance;
  return 0;
}

