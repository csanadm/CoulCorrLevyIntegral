#include <iostream>
#include "CoulCorrCalc.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

using namespace std;

CoulCorrCalc *cccinstance;
TGraphErrors* gr;

double Qmin;
const int NPARS = 4;
int NDF;

const char *statuses[6] = {
                         "converged",
                         "cov. made pos.def.",
                         "Hesse invalid",
                         "Edm above max",
                         "call lim. reached",
                         "other failure"
                       };
const char *covstatuses[5] = {
                            "not available",
                            "not pos.def.",
                            "approximate",
                            "forced pos.def.",
                            "accurate"
                          };

double FitFunction(const double *x, const double *par)
{
  double N      = par[0];
  double lambda = par[1];
  double R      = par[2];
  double alpha  = par[3];
  double Q      = x[0];
  double corrfunc = cccinstance->FullCorrFuncValueLambda(alpha, R, lambda, Q);
  return N*corrfunc;
}

double MyChi2(const double *par)
{
  double chi2 = 0;
  NDF = 0;
  for(int ix=1;ix<gr->GetN();ix++)
  {
    double Q = gr->GetX()[ix];
    if(Q<Qmin) continue;
    double exp = gr->GetY()[ix];
    double theor = FitFunction(&Q,par);
    double err = gr->GetEY()[ix];
    if(err==0) continue;
    double chi = (exp-theor)/err;
    chi2 += chi*chi;
    NDF++;
  }
  NDF -= NPARS;
  return chi2;
}

int main()
{
  cccinstance = new CoulCorrCalc();
  
  gr = new TGraphErrors("Cqdata.txt","%lg %lg %lg");
  Qmin = 0.017;

  // Choose method upon creation between:
  // kMigrad, kSimplex, kCombined, kScan, kFumili
  ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kCombined );

  min.SetMaxFunctionCalls(1000000);
  min.SetMaxIterations(100000);
  min.SetTolerance(0.001);

  ROOT::Math::Functor f(&MyChi2,NPARS);

  min.SetFunction(f);

  // Set the free variables to be minimized!
  min.SetVariable(0,"N",     1  ,0.01);
  min.SetVariable(1,"lambda",0.6,0.01);
  min.SetVariable(2,"R",     5  ,0.01);
  min.SetVariable(3,"alpha", 1  ,0.01);

  // min.SetFixedVariable(0,"x0",91);

  min.Minimize();
  min.PrintResults();
  min.ProvidesError();
  min.Hesse();
  const double *par = min.X();
  const double *err = min.Errors();
  double chi2 = MyChi2(par);
  cout << "Probability: " << chi2 << "/" << NDF << "->" << TMath::Prob(chi2,NDF) << endl;
  cout << "Parameters:" << endl;
  for(int ipar=0;ipar<NPARS;ipar++) cout << "par" << ipar << "=" << par[ipar] << "+-" << err[ipar] << endl;

  cout << "Fit status:" << endl;
  int fitstatus = min.Status();
  int fitcovstatus = min.CovMatrixStatus();
  if(fitstatus<0 || fitstatus>5) fitstatus=5;
  if(fitcovstatus<-1 || fitcovstatus>3) fitcovstatus=-1;
  if(fitstatus == 0 && fitcovstatus == 3) cout << "Fit converged, full accurate cov. matrix";
  else 
  {
    cout << "fit status: " << statuses[fitstatus] << endl;
    cout << "cov. matrix " << covstatuses[fitcovstatus+1] << endl;
  }
  cout << endl;
  cout << "(fitstatus=" << fitstatus << ",covstatus=" << fitcovstatus << ")" << endl;

  double errup[NPARS] = {0};
  double errdn[NPARS] = {0};

  cout << "Minos errors: " << endl;
  for(unsigned int ipar=0; ipar<NPARS; ipar++)
  {
    min.GetMinosError(ipar,errdn[ipar],errup[ipar]);
    cout << "err" << ipar << ": +" << errup[ipar] << " " << errdn[ipar] << endl;
  }

  return 0;
}

