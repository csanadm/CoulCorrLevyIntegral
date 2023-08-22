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
- [**fitexample.cc**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/fitexample.cc): An example `Minuit2` fit code, fitting the "fake" data provided in [`Cqdata.txt`](https://github.com/csanadm/CoulCorrLevyIntegral/blob/main/Cqdata.txt)

### Testing the precision of the Gauss-Kronrod integral
- [**integraltest.cc**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/integraltest.cc): Calculation of the integral with various precision settings
- [**integraltest.py**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/integraltest.py): A `python` plotter for the result of `integraltest.cc`
- [**integraltest.C**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/integraltest.C): A `ROOT` plotter for the result of `integraltest.cc`

## Example results
The below example has been created using the output from [coulcorrtest.cc](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/coulcorrtest.cc)

![coulcorrtest](https://github.com/csanadm/CoulCorrLevyIntegral/assets/38218165/8ff72bda-34ae-4d04-bdf0-5486dcbdf6f7)

The below examples has been created using the output from [integraltest.cc](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/integraltest.cc)

![integraltest Q040](https://github.com/csanadm/CoulCorrLevyIntegral/assets/38218165/d3362a56-d303-4187-b83c-61817ec2df94)
![integraltest Q200](https://github.com/csanadm/CoulCorrLevyIntegral/assets/38218165/b6cb5c80-b3d2-4925-8336-57806f08608c)

## Publications
- Márton Nagy, Aletta Purzsa, Máté Csanád, Dániel Kincses, A novel method for calculating Bose-Einstein correlation functions with Coulomb final-state interaction (2023) [arXiv:2308.10745](https://arxiv.org/abs/2308.10745)
