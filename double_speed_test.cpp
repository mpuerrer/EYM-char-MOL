#include <iostream>
#include <cmath>

using namespace std;

int main() {
  cout.precision(16);
  cout.setf(ios_base::scientific);

  int n = 100000;
  double* f = new double[n];
  double* fp = new double[n];
  double dx = 1.0 / (n-1);
  double sigma_inv = 100.0;
  double x_mid = 0.5;
  
  for (int i=0; i<n; i++)
    f[i] = exp( -sigma_inv * (i*dx - x_mid)*(i*dx - x_mid));

  for (int j=0; j<100; j++)
    for (int i=1; i<n-1; i++)
      fp[i] = (f[i+1] - f[i-1]) / (2 * dx);

  cout << fp[1] << endl;
  
  delete [] f;
  delete [] fp;
  
  return 0;
}