This directory contains algorithms for solving multiobjective optimization problems, described in

M. L. N. Gonçalves, F. S. Lima, and L. F. Prudente, Globally convergent Newton-type methods for multiobjective optimization, technical report, 2022.

- MOPsolverNewton.f90: routine containing the Newton algorithm with safeguards
- MOPsolverNG.f90: routine containing the Newton-Gradient algorithm with safeguards
- MOPsolverNewtonWoSafeg.f90: routine containing the Newton algorithm without safeguards
- MOPsolverSD.f90: routine containing the steepest descent algorithm
- MOPsolverNewtonScalar.f90: routine containing the Scalarized Newton algorithm with safeguards

This folder also contains the third-party free codes: 
1) software Algencan 3.1.1;
    -  E. G. Birgin and J. M. Martı́nez, Practical augmented Lagrangian methods for constrained optimization, SIAM, 2014.
    - https://www.ime.usp.br/~egbirgin/tango/

2) LAPACK, version 3.8.0 (see lapack.f).
   - Copyright (c) 1992-2013 The University of Tennessee and The University of Tennessee Research Foundation.  All rights reserved.
     Copyright (c) 2000-2013 The University of California Berkeley. All rights reserved.
     Copyright (c) 2006-2013 The University of Colorado Denver.  All rights reserved.
   - http://www.netlib.org/lapack/

File main.f90 contains the main program. Modify myproblem.f90 routine to solve your own problem. Alternatively, set a test problem in main.f90 routine - see myproblem.f90.

Instructions:
-------------

File main.f90 contains the main program where you can choose the algorithm to be used.
Modify myproblem.f90 routine to solve your own problem. Alternatively, set a test problem in main.f90 routine; see myproblem.f90.

The codes are written in Fortran 90. Users need to install gfortran.

1) In the terminal, go to the folder and type:

    make

2) Run typing:

    ./MOPsolver

and see the output in the screen.

- out: outer iteration number
- |theta|: optimality measure 
- LS: flag of the line search routine to compute the step size (0 means success)
- IS: flag of the inner solver routine to compute the search direction (0 means success)
- #evalf: number of function evaluations
- #evalg: number of gradient evaluations

