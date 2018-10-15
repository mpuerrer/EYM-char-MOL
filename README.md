# EYM-char-MOL
A C++ code to solve the coupled Einstein-Yang-Mills system as a characteristic IVP in spherical symmetry to study tails due to backscattering of outgoing radiation.

I cast the characteristic initial value problem (in Bondi coordinates) in a form amenable to standard methods for the solution of Cauchy problems. In particular, I solve the nonlinear wave equation for the Yang-Mills field $w$ with a method-of-lines discretization using an upwind-biased 6th order accurate finite difference stencil in space and the classical 4th order explicit Runge-Kutta in time. I use a spatially compactified radius variable $x$ for the evolution. Along with the wave equation which now appears as an advection equation I solve constraint equations which are ordinary differential equations in the slices of constant Bondi time $u$.

This code was used to preoduce the results published in M. Pürrer and P. C. Aichelburg, “Tails for the Einstein-Yang-Mills system,” Class. Quant. Grav. 26 (2009) 035004, arXiv:0810.2648 [gr-qc].


