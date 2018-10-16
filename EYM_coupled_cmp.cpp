/***************************************************************************
 *   Copyright (C) 2008 by Michael Pürrer                                  *
 *   Michael.Puerrer@univie.ac.at                                          *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/*
 * $Id$
 */

/** @file
 *
 *  This is the core numerical evolution code for long-time evolutions of the
 *  Einstein-Yang-Mills system. Its primary purpose is the determination of tails.
 *
 */
 
 
#include <cmath>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>


typedef enum {YM_FLAT, YM_SS, EYM_COUPLED} SCHEME;

std::ofstream h_file, hbar_tails_file;

const int tail_array_sample_every = 100;

inline int floorint(const double& value) {
    return static_cast<int>( std::floor(value) );
}

void cumsimps(double *x, double *y, double *z, int n) {
  double h = 1.0 / (n-1);
  z[0] = 0;

  for (int i=1; i < n-1; i+=2) {
    z[i]   = h/12.0 * (5*y[i-1] + 8*y[i] - y[i+1]);   // odd index
    z[i+1] = h/12.0 * ( -y[i-1] + 8*y[i] + 5*y[i+1]); // even index
  }

  if (n % 2 == 0) // if n is even
    // We can't apply the stencil for odd indices, as it would require y(n).
    // So we resort to the stencil for even indices
    z[n-1] = h/12.0 * (-y[n-3] + 8*y[n-2] + 5*y[n-1]);

  for (int i=1; i < n; i++) // cumulative sum
    z[i] += z[i-1];
}

void cummilne(double *x, double *y, double *z, int n) {
  // cumulative Milne/Boole (Newton-Cotes order 6) quadrature
  assert(n > 4);

  double h = 1.0 / (n-1);
  z[0] = 0;

  for (int i=1; i < n-1; i+=4) {
    z[i]   = h/720.0 * (251*y[i-1] + 646*y[i] - 264*y[i+1] + 106*y[i+2] -  19*y[i+3]);   // index % 4 == 1
    z[i+1] = h/720.0 * (-19*y[i-1] + 346*y[i] + 456*y[i+1] -  74*y[i+2] +  11*y[i+3]);   // index % 4 == 2
    z[i+2] = h/720.0 * ( 11*y[i-1] -  74*y[i] + 456*y[i+1] + 346*y[i+2] -  19*y[i+3]);   // index % 4 == 3
    z[i+3] = h/720.0 * (-19*y[i-1] + 106*y[i] - 264*y[i+1] + 646*y[i+2] + 251*y[i+3]);   // index % 4 == 0
  }

  if (n % 4 == 0) {
    z[n-3] = h/720.0 * (-19*y[n-5] + 346*y[n-4] + 456*y[n-3] -  74*y[n-2] +  11*y[n-1]);
    z[n-2] = h/720.0 * ( 11*y[n-5] -  74*y[n-4] + 456*y[n-3] + 346*y[n-2] -  19*y[n-1]);
    z[n-1] = h/720.0 * (-19*y[n-5] + 106*y[n-4] - 264*y[n-3] + 646*y[n-2] + 251*y[n-1]);
  }
  else if (n % 4 == 2) {
    z[n-1] = h/720.0 * (-19*y[n-5] + 106*y[n-4] - 264*y[n-3] + 646*y[n-2] + 251*y[n-1]);
  }  
  else if (n % 4 == 3) {
    z[n-2] = h/720.0 * ( 11*y[n-5] -  74*y[n-4] + 456*y[n-3] + 346*y[n-2] -  19*y[n-1]);
    z[n-1] = h/720.0 * (-19*y[n-5] + 106*y[n-4] - 264*y[n-3] + 646*y[n-2] + 251*y[n-1]);
  }  

  for (int i=1; i < n; i++) // cumulative sum
    z[i] += z[i-1];
}

inline double mRHS(double x_i, double m_i, double h_i, double hbar_i, int i) {
  if (i==0)
    return 0;
  else {
    double h_i_sq = h_i*h_i;
    double hbar_i_sq = hbar_i*hbar_i;
    return (1 - 2*m_i * (1-x_i) / x_i) * (1-x_i)*(1-x_i) * h_i_sq 
           + hbar_i_sq * (4 + 4*hbar_i + hbar_i_sq) / (2 * x_i*x_i);
  }
}  

// see interpolation_polynomials_equispaced.nb for the computation of the interpolation stencils below:
inline double poly_interpol_left(double* f, int iS, int n) {
  assert(iS >= 0 && iS+3 < n);
  return (5*f[iS] + 15*f[iS+1] - 5*f[iS+2] + f[iS+3]) / 16.0; // 4th order 
}

inline double poly_interpol_right(double* f, int iE, int n) {
  assert(iE-3 >= 0 && iE < n);
  return (f[iE-3] - 5*f[iE-2] + 15*f[iE-1] + 5*f[iE]) / 16.0; // 4th order
}

inline double poly_interpol_center(double* f, int iS, int n) {
  assert(iS >= 0 && iS+3 < n);
  return (-f[iS] + 9*f[iS+1] + 9*f[iS+2] - f[iS+3]) / 16.0; // 4th order
}
  
typedef double (*RHSfunc) (double, double, double, double, int);

void rk4_scalar_i(RHSfunc f, double *x, double f0, int n, double* h, double *hbar, double *psi) {
  double dx, q1, q2, q3, q4;
  psi[0] = f0;
  double h_interp, hbar_interp;
  
  for (int i=0; i < n-1; i++) {
    dx = x[i+1] - x[i]; // let's keep this general
    
    // interpolate h and hbar at i + 1/2
    if (i > 0 && i < n-2) {
      h_interp    = poly_interpol_center(h,    i-1, n);
      hbar_interp = poly_interpol_center(hbar, i-1, n);
    } else if (i == 0) {
      h_interp    = poly_interpol_left(h,    0, n);
      hbar_interp = poly_interpol_left(hbar, 0, n);
    } else { // i == n-1
      h_interp    = poly_interpol_right(h,    n-1, n);
      hbar_interp = poly_interpol_right(hbar, n-1, n);
    }
    
    q1 = dx * f(x[i],          psi[i],          h[i],     hbar[i],     i);
    q2 = dx * f(x[i] + 0.5*dx, psi[i] + 0.5*q1, h_interp, hbar_interp, i);
    q3 = dx * f(x[i] + 0.5*dx, psi[i] + 0.5*q2, h_interp, hbar_interp, i);
    q4 = dx * f(x[i] + dx,     psi[i] + q3,     h[i+1],   hbar[i+1],   i);
    
    psi[i+1] = psi[i] + (q1 + 2*q2 + 2*q3 + q4) / 6.0;
  }
}
  
void YM_flat_cmp_RHS_Q(double* x, double* h, double* Q, int n, double dummy, double dx, bool print_diagnostics = false) {
  double nu = -0.5 / dx;
  
  double* hbar = new double[n];
  
  // integrate constraint
  //cumsimps(x, h, hbar, n);
  cummilne(x, h, hbar, n);
  
  if (print_diagnostics) std::cout << std::endl;  // nothing interesting to print here
  
  // centered 4th order
  double s = 60; // artificial dissipation parameter 
  // stable until s~2000, and convergence order seems to increase!
  Q[0] = -nu/12 * (1-x[0])*(1-x[0]) * (-25*h[0] + 48*h[1] -36*h[2] +16*h[3] -3*h[4]);
  Q[1] = -nu/12 * (1-x[1])*(1-x[1]) * (-25*h[1] + 48*h[2] -36*h[3] +16*h[4] -3*h[5]);
  for (int i=2; i < n-2; i++)
    Q[i] = nu/12.0 * (1-x[i])*(1-x[i]) * (-h[i-2] + 8*h[i-1] - 8*h[i+1] + h[i+2])
         - s*dx * (1-x[i])*(1-x[i]) * (h[i-2] - 4*h[i-1] + 6*h[i] - 4*h[i+1] + h[i+2]);
  
  Q[n-2] = nu/6.0 * (1-x[n-2])*(1-x[n-2]) * (-h[n-4] + 6*h[n-3] -3*h[n-2] - 2*h[n-1]); // 3rd order
  Q[n-1] = 0; // advection outer BC

  double hb, hbs, hbc, xi;
  for (int i=1; i<n; i++) {
    hb = hbar[i];
    hbs = hb*hb;
    hbc = hb*hbs;
    xi = x[i];
    Q[i] += -(1 - xi)*h[i] // sink term
            -0.5 * (2*hb + 3*hbs + hbc) / (xi*xi); // YM self-interaction nonlinearity
  }
  Q[0] = 0; // regularity BC for YM
  // no origin BC for scalar field

  delete [] hbar;
}

void YM_Schwarzschild_cmp_RHS_Q(double* x, double* h, double* Q, int n, double m, double dx, bool print_diagnostics = false) {
  double nu = -0.5 / dx;
  
  double* hbar = new double[n];
  double* H = new double[n];
  
  // integrate constraint
  cummilne(x, h, hbar, n);
  
  // auxiliary field for conservative FD
  for (int i=0; i < n; i++)
    H[i] = (1 - 2*m*(1-x[i])/x[i]) * (1 - x[i])*(1 - x[i]) * h[i];
  
  if (print_diagnostics) std::cout << std::endl;  // nothing interesting to print here
	    
  Q[0] = -nu/6.0 * (-11*H[0] + 18*H[1] - 9*H[2] + 2*H[3]); // 3rd order
  
  for (int i=1; i < n-2; i++)
    Q[i] = nu/6.0 * (2*H[i-1] + 3*H[i] - 6*H[i+1] + H[i+2]);

  Q[n-2] = nu/6.0* (-H[n-4] + 6*H[n-3] -3*H[n-2] - 2*H[n-1]); // 3rd order
  
  // 6th order upwind biased    
  Q[0] = -nu/60.0 * (-147*H[0] + 360*H[1] - 450*H[2] + 400*H[3] - 225*H[4] + 72*H[5] - 10*H[6]);
  Q[1] = -nu/60.0 * (-10*H[0] - 77*H[1] + 150*H[2] - 100*H[3] + 50*H[4] - 15*H[5] + 2*H[6]);
  Q[2] = -nu/60.0 * (2*H[0] - 24*H[1] - 35*H[2] + 80*H[3] - 30*H[4] + 8*H[5] - H[6]);
  
  for (int i=3; i < n-6; i++)
    Q[i] = -nu/60.0 * (2*H[i-2] - 24*H[i-1] - 35*H[i] + 80*H[i+1] - 30*H[i+2] + 8*H[i+3] - H[i+4]); // stable, upwind biased
           
  Q[n-6] = -nu/60.0 * (-10*H[n-7] - 77*H[n-6] + 150*H[n-5] - 100*H[n-4] + 50*H[n-3] - 15*H[n-2] + 2*H[n-1]);
  Q[n-5] = -nu/60.0 * (2*H[n-7] - 24*H[n-6] - 35*H[n-5] + 80*H[n-4] - 30*H[n-3] + 8*H[n-2] - H[n-1]);    
  Q[n-4] = -nu/60.0 * (H[n-8] - 8*H[n-7] + 30*H[n-6] - 80*H[n-5] + 35*H[n-4] + 24*H[n-3] - 2*H[n-2]);
  Q[n-3] = -nu/60.0 * (H[n-7] - 8*H[n-6] + 30*H[n-5] - 80*H[n-4] + 35*H[n-3] + 24*H[n-2] - 2*H[n-1]);
  Q[n-2] = -nu/60.0 * (-2*H[n-7] + 15*H[n-6] - 50*H[n-5] + 100*H[n-4] - 150*H[n-3] + 77*H[n-2] + 10*H[n-1]);

      
  Q[n-1] = 0; // no advection at x=1; this is a characteristic!

  double hb, hbs, hbc, xi;
  for (int i=0; i<n; i++) {
    hb = hbar[i];
    hbs = hb*hb;
    hbc = hb*hbs;
    xi = x[i];
    Q[i] += -0.5 * (2*hb + 3*hbs + hbc) / (xi*xi); // YM self-interaction nonlinearity
  }

  delete [] hbar;
  delete [] H;
}

double max(double *v, int n) {
  double max = v[0];
  for (int i=1; i<n; i++)
    if (v[i] > max)
      max = v[i];
      
  return max;
}

void EYM_coupled_cmp_RHS_Q(double* x, double* h, double* Q, int n, double dummy, double dx, bool print_diagnostics = false) {
  double nu = -0.5 / dx;
  
  double* hbar = new double[n];
  double* H = new double[n];
  double* beta = new double[n];
  double* e2beta = new double[n];
  double* betaRHS = new double[n];
  double* m = new double[n];
  double* two_m_r = new double[n];
  
  // integrate constraints ////////////////////////////////////
  //cumsimps(x, h, hbar, n); // hbar
  cummilne(x, h, hbar, n);
  
  betaRHS[0] = 0; // betaRHS = O(x)
  for (int i=1; i<n; i++) {
    double a = 1-x[i];
    double b = h[i];
    betaRHS[i] = a*a*a / x[i] * b*b;	    
  }
  //cumsimps(x, betaRHS, beta, n); // beta
  cummilne(x, betaRHS, beta, n);
  
  rk4_scalar_i(mRHS, x, 0, n, h, hbar, m); // m

  // initialize auxiliary arrays ///////////////////////////////
  // 2m/r diagnostic
  two_m_r[0] = 0;
  for (int i=1; i<n; i++)
    two_m_r[i] = 2*m[i] * (1-x[i]) / x[i]; 
  
  if (print_diagnostics)
    std::cout << "\tmax 2m/r = " << max(two_m_r, n) << std::endl;
  
  // auxiliary field for conservative FD (thus, we don't need a sink term below)
  H[0] = 0; // H = O(x)
  for (int i=1; i < n; i++) {
    double a = 1-x[i];
    e2beta[i] = exp(2*beta[i]);
    H[i] = e2beta[i] * (1 - 2*m[i]* a/x[i]) * a*a * h[i];
  }
  
  // advection term  ///////////////////////////////////////////
  
  // 6th order upwind biased    
  Q[0] = -nu/60.0 * (-147*H[0] + 360*H[1] - 450*H[2] + 400*H[3] - 225*H[4] + 72*H[5] - 10*H[6]);
  Q[1] = -nu/60.0 * (-10*H[0] - 77*H[1] + 150*H[2] - 100*H[3] + 50*H[4] - 15*H[5] + 2*H[6]);
  Q[2] = -nu/60.0 * (2*H[0] - 24*H[1] - 35*H[2] + 80*H[3] - 30*H[4] + 8*H[5] - H[6]);
  
  for (int i=3; i < n-6; i++)
    Q[i] = -nu/60.0 * (2*H[i-2] - 24*H[i-1] - 35*H[i] + 80*H[i+1] - 30*H[i+2] + 8*H[i+3] - H[i+4]); // stable, upwind biased
           
  Q[n-6] = -nu/60.0 * (-10*H[n-7] - 77*H[n-6] + 150*H[n-5] - 100*H[n-4] + 50*H[n-3] - 15*H[n-2] + 2*H[n-1]);
  Q[n-5] = -nu/60.0 * (2*H[n-7] - 24*H[n-6] - 35*H[n-5] + 80*H[n-4] - 30*H[n-3] + 8*H[n-2] - H[n-1]);    
  Q[n-4] = -nu/60.0 * (H[n-8] - 8*H[n-7] + 30*H[n-6] - 80*H[n-5] + 35*H[n-4] + 24*H[n-3] - 2*H[n-2]);
  Q[n-3] = -nu/60.0 * (H[n-7] - 8*H[n-6] + 30*H[n-5] - 80*H[n-4] + 35*H[n-3] + 24*H[n-2] - 2*H[n-1]);
  Q[n-2] = -nu/60.0 * (-2*H[n-7] + 15*H[n-6] - 50*H[n-5] + 100*H[n-4] - 150*H[n-3] + 77*H[n-2] + 10*H[n-1]);
  
  
  Q[n-1] = 0; // no advection at x=1; this is a characteristic!

  // YM self-interaction nonlinearity  //////////////////////////	  
  double hb, hbs, hbc, xi;
  for (int i=0; i<n; i++) {
    hb = hbar[i];
    hbs = hb*hb;
    hbc = hb*hbs;
    xi = x[i];
    Q[i] += -0.5 * e2beta[i] * (2*hb + 3*hbs + hbc) / (xi*xi);
  }

  Q[0] = 0; // origin regularity BC (it's really the time-derivative of it)

  delete [] hbar;
  delete [] H;
  delete [] beta;
  delete [] e2beta;
  delete [] betaRHS;
  delete [] m;
  delete [] two_m_r;
}

void rk4_Q(double k, int m, const double* psi0, int n, SCHEME scheme) {
  /*
   * Evolve Einstein-Yang-Mills coupled with compactified radial coordinate
   *
   * Input:   k               size of timestep
   *          m               length of t, i.e. the number of timesteps
   *          psi0            initial data for the field hbar
   *          n               length of psi0, i.e. the number of spatial gridpoints
   *          scheme          one of YM_FLAT, YM_SS, EYM_COUPLED
   * Output:  h, hbar_tails   written top files
   */ 
   
  std::cout << "Evolution using " << n << " gridpoints" << std::endl;
  
  double* psi   = new double[n];
  double* x     = new double[n];
  double* x_SS  = new double[n];
  double* h     = new double[n];
  double* temp  = new double[n];
  double* Q     = new double[n];
  double* q1    = new double[n];
  double* q2    = new double[n];
  double* q3    = new double[n];
  double* q4    = new double[n];
  double* hbar  = new double[n];
  
  void (*scheme_RHS_Q) (double*, double*, double*, int, double, double, bool);

  // copy initial data into the psi array
  for (int i=0; i<n; i++) psi[i] = psi0[i]; 

  // regular evolution domain x \in [0, 1]
  double dx = 1.0 / (n-1);
  for (int i=0; i<n; i++) 
    x[i] = i*dx;

  // Schwarzschild evolution domain x \in [2m/(1+2m), 1]
  const double m_SS = 0.5; // Schwarzschild mass parameter      
  double x_SS_origin = 2*m_SS / (1.0 + 2*m_SS);
  double dx_SS = (1.0 - x_SS_origin) / (n-1);
  for (int i=0; i<n; i++) 
    x_SS[i] = x_SS_origin + i*dx_SS;
    
  switch (scheme) {
    case YM_FLAT : 
      scheme_RHS_Q = YM_flat_cmp_RHS_Q;
      std::cout << "Yang-Mills on flat space using compactified radial coordinate x." << std::endl;
      break;            
    case YM_SS : 
      scheme_RHS_Q = YM_Schwarzschild_cmp_RHS_Q;
      std::cout << "Yang-Mills on Schwarzschild using compactified radial coordinate x." << std::endl;
      dx = dx_SS;
      for (int i=0; i<n; i++) x[i] = x_SS[i];
      break;            
    case EYM_COUPLED : 
      std::cout << "Einstein-Yang-Mills coupled using compactified radial coordinate x." << std::endl;
      scheme_RHS_Q = EYM_coupled_cmp_RHS_Q;
      break;    
    default :
      std::cerr << "Invalid scheme specified" << std::endl;
      exit(-1);
  }
    
    
  bool print_diagnostics = false;
  
  for (int i=0; i<m-1; i++) {  // timeloop
    
    for (int j=0; j<n; j++)
      h[j] = psi[j]; // old field vector
     
    if (i % (m/50) == 0) {
      // dump field to file
      h_file << "# u = " << i*k << std::endl;      
      for (int j=0; j<n; j++)
        h_file << x[j] << '\t' << h[j] << '\n';
      h_file << '\n' << std::endl;

      std::cout << "u = " << i*k;
      print_diagnostics = true;
    }

    scheme_RHS_Q(x, h, Q, n, m_SS, dx, print_diagnostics);
    for (int j=0; j<n; j++) {
      q1[j] = k*Q[j];
      temp[j] = h[j] + 0.5*q1[j];
    }

    scheme_RHS_Q(x, temp, Q, n, m_SS, dx, false);  
    for (int j=0; j<n; j++) {
      q2[j] = k*Q[j];
      temp[j] = h[j] + 0.5*q2[j];
    }

    scheme_RHS_Q(x, temp, Q, n, m_SS, dx, false);
    for (int j=0; j<n; j++) {
      q3[j] = k*Q[j];
      temp[j] = h[j] + q3[j];
    }

    scheme_RHS_Q(x, temp, Q, n, m_SS, dx, false);
    for (int j=0; j<n; j++)
      q4[j] = k*Q[j];

    for (int j=0; j<n; j++)
      psi[j] = h[j] + (q1[j] + 2*q2[j] + 2*q3[j] + q4[j]) / 6.0;
    
	  // tails: integrate constraint using the old field -- won't make much of a difference
	  if (scheme == YM_SS)
	    cumsimps(x_SS, h, hbar, n);
    else
      cumsimps(x, h, hbar, n);
    
	  if (i % tail_array_sample_every == 0) {
  	  // tails: locations of observers
      hbar_tails_file << i*k << '\t' 
                      << hbar[floorint(n*0.5)]    << '\t' 
                      << hbar[floorint(n*0.7)]    << '\t' 
                      << hbar[floorint(n*0.8)]    << '\t'
                      << hbar[floorint(n*0.9)]    << '\t' 
                      << hbar[floorint(n*0.95)]   << '\t'
                      << hbar[floorint(n*0.99)]   << '\t' 
                      << hbar[floorint(n*0.995)]  << '\t'
                      << hbar[floorint(n*0.999)]  << '\t' 
                      << hbar[floorint(n*0.9995)] << '\t'
                      << hbar[n-1]      << std::endl;
    }
    
    print_diagnostics = false;
  }

  std::cout << "Finished " << m << " timesteps using " << n << " gridpoints." << std::endl;

  delete [] x;
  delete [] x_SS;
  delete [] h;
  delete [] temp;
  delete [] q1;
  delete [] q2;
  delete [] q3;
  delete [] q4;
  delete [] Q;
  delete [] hbar; 
  delete [] psi;
}

void MOL_Q_driver(int n, double uf, SCHEME scheme) {
  double dx = 1.0 / (n-1);
  double C = 0.9; // du/dr
  double du = C*dx;
  int m = floorint(uf / du);
  
  double *x     = new double[n];
  double *x_SS  = new double[n];
  double *hbar0 = new double[n];
  double *h0    = new double[n];
  
  // regular evolution domain x \in [0, 1]
  for (int i=0; i<n; i++)
    x[i] = i*dx;
  
  // Schwarzschild evolution domain x \in [2m/(1+2m), 1]
  double m_SS = 0.5; // Schwarzschild mass parameter      
  double x_SS_origin = 2*m_SS / (1.0 + 2*m_SS);
  double dx_SS = (1.0 - x_SS_origin) / (n-1);
  for (int i=0; i<n; i++)
    x_SS[i] = x_SS_origin + i*dx_SS;
  
  // initial data
  const double a = 0.28; // EYM strong field max 2m/r ~ 0.6
  // const double a = 0.01;
  const double sigma_inv = 200.0;
  const double x_mid = 0.5;
  for (int i=0; i<n; i++) {
    double xi = x[i] - x_mid;
    hbar0[i] = a * exp(-sigma_inv * xi*xi); // Gaussian for the physical field
    h0[i] = -2*sigma_inv * hbar0[i] * xi; // analytical derivative h = (\bar h)'
  }
  
  // evolve
  rk4_Q(du, m, h0, n, scheme);
  
  delete [] x;
  delete [] x_SS;
  delete [] hbar0;
  delete [] h0;
}


int main(int argc, char* argv[]) {
  int n, s;
  double uf;
  SCHEME scheme;
  
  if (argc != 4 || std::strcmp(argv[1], "-h") == 0) {
    std::cout << 
    "syntax: EYM_coupled_cmp n uf scheme\n\
    n   ... number of spatial gridpoints\n\
    uf  ... evolve for retarded time u in [0,uf]\n\
    scheme: 0 ... YM_FLAT\n\
            1 ... YM_SS\n\
            2 ... EYM_COUPLED\n" << std::endl;
    exit(0);
  }
  
  if ((sscanf(argv[1], "%d", &n) != 1) ||
      (sscanf(argv[2], "%lf", &uf) != 1) ||
      (sscanf(argv[3], "%d", &s) != 1)) {
    std::cerr << "Invalid parameter" << std::endl;
    exit(-1);
  }
    
  switch (s) {
    case 0: 
      scheme = YM_FLAT;
      break;
    case 1: 
      scheme = YM_SS;
      break;
    case 2: 
      scheme = EYM_COUPLED;
      break;
    default:
      std::cerr << "Invalid scheme" << std::endl;
      exit(-1);
  }  
  
  h_file.open("h.dat");
  h_file << "# n = " << n << "\t uf = " << uf << std::endl;
  
  hbar_tails_file.open("hbar_tails.dat");
  hbar_tails_file << "# n = " << n << "\t uf = " << uf << std::endl;
  hbar_tails_file << "# u" << '\t' 
                  << "hbar[n*0.5] \t" 
                  << "hbar[n*0.7]\t" 
                  << "hbar[n*0.8]\t"
                  << "hbar[n*0.9] \t" 
                  << "hbar[n*0.95]\t"
                  << "hbar[n*0.99]\t" 
                  << "hbar[n*0.995]\t"
                  << "hbar[n*0.999]\t" 
                  << "hbar[n*0.9995]\t"
                  << "hbar[n-1]"  << std::endl;                

  MOL_Q_driver(n, uf, scheme);
  
  h_file.close();
  hbar_tails_file.close();
}
