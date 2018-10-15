#include <iostream>
#include <cmath>
#include <qd/dd_real.h>

using namespace std;

int main() {
  unsigned int old_cw;
  fpu_fix_start(&old_cw);
  
  cout.precision(30);
  cout.setf(ios_base::scientific);

  int n = 100000;
  dd_real* f = new dd_real[n];
  dd_real* fp = new dd_real[n];
  dd_real dx = 1.0 / (n-1);
  dd_real sigma_inv = "100.0";
  dd_real x_mid = "0.5";
  
  for (int i=0; i<n; i++)
    f[i] = exp( -sigma_inv * (i*dx - x_mid)*(i*dx - x_mid));

  for (int j=0; j<100; j++)
    for (int i=1; i<n-1; i++)
      fp[i] = (f[i+1] - f[i-1]) / (2 * dx);

  cout << fp[1] << endl;
  
  delete [] f;
  delete [] fp;
 
  fpu_fix_end(&old_cw);
  
  return 0;
}