
# Mathematica AdditiveGroupPackage Suite For Modular Arithmetic 

A Mathematica package suite for exploring the simplest finite groups Zn, the set of 
integers from 0 to n-1 with group operation addition modulu n. 

The package structure has many purposes like code maintanability and reading, 
modularity for later additions and not least for learnig purposes, which is
prime to this project of course.

The suite consists of four sub packages

  - AdditiveGroupMinimal    The most necessary to start build it yourself as that is
                            best way of learning.
  - AdditiveGroupBasics     Some more than minimal.
  - AdditiveGroup           What is common functionality without using quotient groups.
  - AdditiveGroupQuotients  Functions for exploring quotient groups.
  - AdditiveGroupTheorems   The first famous destinations of algebraic structures namely
                            the correspondence theorem and the three isomorphism theorems.

REMARK: There is a lot of 'in between' functionality that might be very disturbing and 
        these functions may well be without dependencies and can then be removed. 
        Everything will probably be much neater in later releases. They serves as inspiration 
        and can be removed.

To get started after putting the packages in your Mathematica $PackagesPath, eye through
the usage instructions with Doc["Minimal"] for example.

The time predictions are calibrated on a Dell Predision 3250 Intel i7 (x4) 2.9 GHz with 16 GB 
and are probably not right for you.

Enjoy!
