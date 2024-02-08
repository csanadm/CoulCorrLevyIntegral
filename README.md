# Correlation function calculation with Lévy source and Coulomb FSI

 Bose-Einstein correlation function integral calculation based on a Lévy source, incorporating the Coulomb final-state interaction. Most of the integrals can be performed analytically, only one, very well behaving integral remains, that is performed numerically.

## Description
This package contains a calculation for quantum-statistical correlation functions, including Coulomb-correction, based on analytic results and a final numerical integration. For the calculation to work, the `boost` library is needed, although the numerical integral can be programmed by the user as well. One plotting code and the fit test requires `ROOT` (latter also `Minuit2`), and several `python` libraries are required for the `python` plotters as well.

## File content

### Basics
- [**README.md**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/README.md): This README file
- [**Makefile**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/Makefile): Using `make all`, it will create an executable

### Libraries
- [**CoulCorrCalc.cpp**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/CoulCorrCalc.cpp): The main calculator class, containing the formulas and the final integral (via the Gauss-Kronrod method of `boost`)
- [**CoulCorrCalc.h**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/CoulCorrCalc.h): Heacer file for the `CoulCorrCalc` class
- [**HypCalculator.cpp**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/HypCalculator.cpp): The functions in this class calculate the hypergeometric function 2F1.
- [**HypCalculator.h**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/HypCalculator.h): Heacer file for the `HypCalculator` class
- [**functions.cpp**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/functions.cpp): Auxiliary functions, such as the Gamma function
- [**functions.h**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/functions.h): : Heacer file for `functions.cpp`
- [**basics.h**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/basics.h): Basic constants, needed for all calculations

### Testing of the libraries
- [**coulcorrtest.cc**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/coulcorrtest.cc): An example code for testing the library
- [**coulcorrtestplot.py**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/coulcorrtestplot.py): A `python` plotter for the test result
- [**fitexample.cc**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/fitexample.cc): An example `Minuit2` fit code, fitting the "fake" data provided in [`Cqdata.txt`](https://github.com/csanadm/CoulCorrLevyIntegral/blob/main/Cqdata.txt) (which was generated with [`coulcorrtest.cc`](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/coulcorrtest.cc), with the parameters alpha=1.2, R=5.3 fm, lambda=0.8)

### Testing the precision of the Gauss-Kronrod integral
- [**integraltest.cc**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/integraltest.cc): Calculation of the integral with various precision settings
- [**integraltest.py**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/integraltest.py): A `python` plotter for the result of `integraltest.cc`
- [**integraltest.C**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/integraltest.C): A `ROOT` plotter for the result of `integraltest.cc`

## Running the codes
- A `deps` directory should be created, this is where `make` stores the dependencies.
- A `make all` command creates all executables (they will have the `.exe` suffix), or a specific `make <codename>.exe` command produces just the given executable (here `<codename>` can be for example `coulcorrtest`, `fitexample` or `integraltest`).
- One shall run `coulcorrtest.exe > coulcorrtest.out` and then `python coulcorrtestplot.py` to produce a test plot (can edit `coulcorrtest.cc` to make the plot for different parameters).
- The code `fitexample.cc` is run as `coulcorrtest.exe`, and it produces a fit (and writes out the corresponding messages) to `Cqdata.txt`, as mentioned above.
- The code `integraltest.cc` should also be run as  `integraltest.exe > integraltest.out`, which then can be plotted as `python integraltest.py`, or via ROOT as `root.exe -b -q integraltest.C` (where the variable `const char* inputfile` has the name of the file from which data are to be plotted).

## Example results
The below example has been created using the output from [coulcorrtest.cc](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/coulcorrtest.cc)

![coulcorrtest](https://github.com/csanadm/CoulCorrLevyIntegral/assets/38218165/8ff72bda-34ae-4d04-bdf0-5486dcbdf6f7)

The below examples has been created using the output from [integraltest.cc](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/integraltest.cc)

![integraltest Q040](https://github.com/csanadm/CoulCorrLevyIntegral/assets/38218165/d3362a56-d303-4187-b83c-61817ec2df94)
![integraltest Q200](https://github.com/csanadm/CoulCorrLevyIntegral/assets/38218165/b6cb5c80-b3d2-4925-8336-57806f08608c)

The fit example code [`fitexample.cc`](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/fitexample.cc), using the "fake" data [`Cqdata.txt`](https://github.com/csanadm/CoulCorrLevyIntegral/blob/main/Cqdata.txt) generates the following output:

`Minuit2Minimizer : Valid minimum - status = 0`<br>
`FVAL  = 163.628528208775037`<br>
`Edm   = 2.31624412029060214e-07`<br>
`Nfcn  = 239`<br>
`N         = 0.999426     +/-  0.0016909`<br>
`lambda    = 0.801399     +/-  0.027961`<br>
`R         = 5.3208       +/-  0.159243`<br>
`alpha     = 1.18897      +/-  0.0435182`<br>
`Probability: 163.629/191->0.924927`<br>
`Fit status:`<br>
`Fit converged, full accurate cov. matrix`<br>
`(fitstatus=0,covstatus=3)`<br>
`Minos errors:`<br>
`err0: +0.00166718 -0.00171808`<br>
`err1: +0.0292056 -0.0268804`<br>
`err2: +0.167903 -0.151789`<br>
`err3: +0.0438573 -0.0432895`<br>

The fit output shall be compared to the input parameters of N=1, alpha=1.2, R=5.3 fm, lambda=0.8.

## Publications
- Márton Nagy, Aletta Purzsa, Máté Csanád, Dániel Kincses, <i>A novel method for calculating Bose-Einstein correlation functions with Coulomb final-state interaction</i>, Eur. Phys. J. C (2023) [arXiv:2308.10745](https://arxiv.org/abs/2308.10745)
