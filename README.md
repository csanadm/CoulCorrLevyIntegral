# Correlation function calculation with Lévy source and Coulomb FSI

 Bose-Einstein correlation function integral calculation based on a Lévy source, incorporating the Coulomb final-state interaction. Most of the integrals can be performed analytically, only one, very well behaving integral remains, that is performed numerically.

## Description
This package contains a calculation for quantum-statistical correlation functions, including Coulomb-correction, based on analytic results and a final numerical integration. For the calculation to work, the `boost` library is needed, although the numerical integral can be programmed by the user as well.

## File content
- [**README.md**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/README.md): This README file
- [**Makefile**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/Makefile): Using `make all`, it will create an executable
- [**coulcorrtest.cc**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/coulcorrtest.cc): An example code for testing the library
- [**coulcorrtestplot.py**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/coulcorrtestplot.py): A python plotter for plotting the test result
- [**coulcorrtest.cc**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/integraltest.cc): An example code for investigating the integral accuracy
- [**coulcorrtestplot.py**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/integraltest.py): A python plotter for the integral testing
- [**CoulCorrCalc.cpp**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/CoulCorrCalc.cpp): The main calculator class, containing the formulas and the final integral (via the Gauss-Kronrod method of `boost`)
- [**CoulCorrCalc.h**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/CoulCorrCalc.h): Heacer file for the `CoulCorrCalc` class
- [**HypCalculator.cpp**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/HypCalculator.cpp): The functions in this class calculate the hypergeometric function 2F1.
- [**HypCalculator.h**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/HypCalculator.h): Heacer file for the `HypCalculator` class
- [**functions.cpp**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/functions.cpp): Auxiliary functions, such as the Gamma function
- [**functions.h**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/functions.h): : Heacer file for `functions.cpp`
- [**basics.h**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/basics.h): Basic constants, needed for all calculations

## Example results
The below example has been created using the output from [coulcorrtest.cc](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/coulcorrtest.cc)

![coulcorrtest](https://github.com/csanadm/CoulCorrLevyIntegral/assets/38218165/8ff72bda-34ae-4d04-bdf0-5486dcbdf6f7)

The below example has been created using the output from [integraltest.cc](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/integraltest.cc)

![integraltest](https://github.com/csanadm/CoulCorrLevyIntegral/assets/38218165/5a0b8ecb-deb0-4c8a-bb91-2cf29b633eb8)

## Publications
- To be listed
