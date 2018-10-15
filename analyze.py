#!/usr/bin/env python

##########################################################################
#  Copyright (C) 2008 by Michael Pürrer                                  #
#  Michael.Puerrer@univie.ac.at                                          #
#                                                                        #
#  This program is free software; you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation; either version 2 of the License, or     #
#  (at your option) any later version.                                   #
#                                                                        #
#  This program is distributed in the hope that it will be useful,       #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#  GNU General Public License for more details.                          #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with this program; if not, write to the                         #
#  Free Software Foundation, Inc.,                                       #
#  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             #
##########################################################################

#
# $Id$
#

# @file
#
#  This script processes the data generated by the core numerical code (EYM_coupled_cmp.cpp).
#  Its primary purpose is to do tail analysis, i.e. subdivide the time interval that has been evolved
#  into a number of windows and use linear least-squares fitting to determine the tail decay in each 
#  of those. Finally, display the resulting tail exponents via time.
#
#

from numpy import *
from scipy import weave
from scipy import stats
from scipy.integrate import cumtrapz
from numpy.linalg import norm
#from scitools.easyviz import *
from pylab import plot, title, show, legend, xlim, ylim, loglog, figure, xlabel, ylabel, draw, ion, hold


def read_h():
  infile = open('h.dat', 'r')
  line = infile.readline() # read header
  n = int(line.split()[3]) # number of spatial gridpoints
  
  h_slice = zeros(0)
  x_slice = zeros(0)
  h = zeros((0,n))
  
  while True:
    line = infile.readline()
    if not line: break
    if len(line.split()) > 0 and line.split()[0]=='#': continue
    values = [float(x) for x in line.split()]
    
    if len(values) == 0:
      if (len(h_slice) > 0): 
        h = append(h, [h_slice], axis=0)
        h_slice = zeros(0)
        x = x_slice
        x_slice = zeros(0)
      continue
      
    x_slice = append(x_slice, values[0])
    h_slice = append(h_slice, values[1])
  
  return x, h 
  
def read_hbar_tails():
  infile = open('hbar_tails.dat', 'r')
  
  s = infile.readline().split() # read header
  n = int(s[3])
  uf = int(s[6])
  line = infile.readline() 
  #hbar_tails_array = zeros((0,6))
  hbar_tails_array = zeros((0,11))
  
  counter = 0
  while True:
    line = infile.readline()
    if not line: break
    values = [float(x) for x in line.split()]
    if len(values) == 0: continue
    if counter % 10 == 0:
      hbar_tails_array = append(hbar_tails_array, [values], axis=0)
      print counter
    counter += 1
    
  return transpose(hbar_tails_array), uf, n
    
def simple_evo():
  x, h = read_h()
  n = len(x)
  
  # turn interactive mode on for dynamic updates.  If you aren't in
  # interactive mode, you'll need to use a GUI event handler/timer.
  ion()
  #line, = plot(x,h[0,:])
  hbar = zeros(n)
  hbar[1:n] = cumtrapz(h[0,:], x)
  line, = plot(x,hbar)
  for i in xrange(0, h.shape[0]): # plot solution
    hbar[1:n] = cumtrapz(h[i,:], x)
    #line.set_ydata(hbar)
    hold(False)
    plot(x,hbar)
    #line.set_ydata(h[i,:])
    draw()

def fit_for_window(u, uStart, uEnd, n, hbar_tails_array, plot_data=0):
  print "Fitting window [%(uStart)g, %(uEnd)g]:" % vars(),
  iStart, iEnd = where((u > uStart) & (u < uEnd))[0][[0,-1]]
  iEnd -= 1 # last entry may be zero, so skip it; otherwise we get nan's in the fits
  
  [dummy, hbar_x9, hbar_x8, hbar_x7, hbar_x6, hbar_x5, hbar_x4, hbar_x3, hbar_x2, hbar_x1, hbar_Scri] = hbar_tails_array
  
  uClip = u[iStart:iEnd+1]
  hbar_ScriClip = hbar_Scri[iStart:iEnd+1]
  hbar_x1Clip = hbar_x1[iStart:iEnd+1]
  hbar_x2Clip = hbar_x2[iStart:iEnd+1]
  hbar_x3Clip = hbar_x3[iStart:iEnd+1]
  hbar_x4Clip = hbar_x4[iStart:iEnd+1]
  
  hbar_x5Clip = hbar_x5[iStart:iEnd+1]
  hbar_x6Clip = hbar_x6[iStart:iEnd+1]
  hbar_x7Clip = hbar_x7[iStart:iEnd+1]
  hbar_x8Clip = hbar_x8[iStart:iEnd+1]
  hbar_x9Clip = hbar_x9[iStart:iEnd+1]
  
  hbar_tails_Clip = (hbar_ScriClip, hbar_x1Clip, hbar_x2Clip, hbar_x3Clip, hbar_x4Clip, 
                     hbar_x5Clip, hbar_x6Clip, hbar_x7Clip, hbar_x8Clip, hbar_x9Clip)
                     
  n10 = n/10
  slope_list = []
  for hbar_x in hbar_tails_Clip:
    (slope, intercept, r, tt, stderr) = stats.linregress(log10(uClip), log10(abs(hbar_x)))
    # Calculates a regression line on two arrays, x and y, corresponding to x,y pairs.
    # Returns: slope, intercept, r, two-tailed prob, stderr-of-the-estimate
    slope_list.append(slope)
    
  for slope in slope_list:
    print "%.2f" % (slope),
  print
  
  if plot_data == 1:
    loglog(uClip, abs(hbar_ScriClip), uClip, abs(hbar_x1Clip), uClip, abs(hbar_x2Clip), 
          uClip, abs(hbar_x3Clip), uClip, abs(hbar_x4Clip), uClip, abs(hbar_x5Clip), 
          uClip, abs(hbar_x6Clip), uClip, abs(hbar_x7Clip), uClip, abs(hbar_x8Clip), uClip, abs(hbar_x9Clip))       
          
  return slope_list  

def tail_analysis():
  print "tail analysis"
  hbar_tails_array, uf, n = read_hbar_tails()
  
  # looking for tails
  # [u_sampled, hbar_x4, hbar_x3, hbar_x2, hbar_x1, hbar_Scri] = hbar_tails_array
  [u_sampled, hbar_x9, hbar_x8, hbar_x7, hbar_x6, hbar_x5, hbar_x4, hbar_x3, hbar_x2, hbar_x1, hbar_Scri] = hbar_tails_array
  
  # find indices for overall fitting window
  winStart, winEnd = 5, uf
  print "Overall fitting window is [%(winStart)g, %(winEnd)g]" % vars()
  
  winSize = 20
  # subdivide into intervals of length winSize
  num_of_windows = (winEnd - winStart) / winSize
  print "Processing a total of ", num_of_windows, " fitting windows."
  # fit for each of those, no plot
  slopes_windows_list = []
  u_window_list = []
  for uStart in range(winStart, winEnd, winSize):
    slopes_windows_list.append(fit_for_window(u_sampled, uStart, uStart + winSize, n, hbar_tails_array))
    u_window_list.append(uStart + winSize/2)
    
    
  figure(1)
  title('Tail Exponents from Fits')
  xlabel('u')
  ylabel('Tail exponents')
  ylim(-6, -2)
  slopes_windows_array = array(slopes_windows_list)
  for i in range(slopes_windows_array.shape[1]):
    plot(u_window_list, slopes_windows_array[:,i])
    
  # write tail-exponents to file
  data = zeros((len(u_window_list), slopes_windows_array.shape[1] + 1))
  data.T[0] = u_window_list
  data.T[1:] = slopes_windows_array.T
  savetxt('tails.dat', data, fmt="%12.6G")
  
  figure(2)
  title('Whole Evolution')
  xlabel('u')
  ylabel('hbar_const_x')    
  print len(u_sampled), len(hbar_Scri), len(hbar_x1), len(hbar_x2), len(hbar_x3), len(hbar_x4), 
  len(hbar_x5), len(hbar_x6), len(hbar_x7), len(hbar_x8), len(hbar_x9)
  
  k = min(len(u_sampled), len(hbar_Scri))
  loglog(u_sampled[0:k], abs(hbar_Scri[0:k]), 
        u_sampled[0:k], abs(hbar_x1[0:k]), 
        u_sampled[0:k], abs(hbar_x2[0:k]), 
        u_sampled[0:k], abs(hbar_x3[0:k]), 
        u_sampled[0:k], abs(hbar_x4[0:k]),
        u_sampled[0:k], abs(hbar_x5[0:k]), 
        u_sampled[0:k], abs(hbar_x6[0:k]), 
        u_sampled[0:k], abs(hbar_x7[0:k]), 
        u_sampled[0:k], abs(hbar_x8[0:k]),
        u_sampled[0:k], abs(hbar_x9[0:k]))

# main
tail_analysis()
 
show()
