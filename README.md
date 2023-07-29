# CoulCorrLevyIntegral

 Coulomb correction integral calculation based on a Levy source

## Description
This package contains a calculation for quantum-statistical correlation functions, including Coulomb-correction, based on analytic results and a final numerical integration. For the calculation to work, the `boost` library is needed, although the numerical integral can be programmed by the user as well.

## File content
- [**README.md**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/README.md): This README file
- [**Makefile**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/Makefile): Using `make all`, it will create an executable
- [**coulcorrtest.cc**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/coulcorrtest.cc): An example code for testing the library
- [**coulcorrtestplot.py**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/coulcorrtestplot.py): A python plotter for plotting the test result
- [**CoulCorrCalc.cpp**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/CoulCorrCalc.cpp): The main calculator class, containing the formulas and the final integral (via the Gauss-Kronrod method of `boost`)
- [**CoulCorrCalc.h**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/CoulCorrCalc.h): Heacer file for the `CoulCorrCalc` class
- [**HypCalculator.cpp**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/HypCalculator.cpp): The functions in this class calculate the hypergeometric function 2F1.
- [**HypCalculator.h**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/HypCalculator.h): Heacer file for the `HypCalculator` class
- [**functions.cpp**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/functions.cpp): Auxiliary functions, such as the Gamma function
- [**functions.h**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/functions.h): : Heacer file for `functions.cpp`
- [**basics.h**](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/basics.h): 

## Example results
The below example has been created using the output from [coulcorrtest.cc](https://github.com/csanadm/CoulCorrLevyIntegral/blob/master/coulcorrtest.cc)

![coulcorrtest](https://github.com/csanadm/CoulCorrLevyIntegral/assets/38218165/8ff72bda-34ae-4d04-bdf0-5486dcbdf6f7)


## Publications
- To be listed
