# include <cstdlib>
# include <iostream>
# include <cmath>
# include <ctime>

//
//  CHANGE 1) include the OpenMP include file
# include <omp.h>                      //added to include library at 09/06/2016
//
using namespace std;

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for QUAD.
//
//  Discussion:
//
//    QUAD is a C++ program to estimate the integral, from 0 to 10, of
//     f = 50.0 / ( pi * ( 2500.0 * x[i] * x[i] + 1.0 ) )
//
//    If you modify the program in 4 places, it should run correctly
//    under OpenMP.
//
//    If you can do that, then compare the execution time when you
//    use 1, 2, 4 or 8 threads.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  double a = 0.0;
  double b = 10.0;
  double error;
  double exact = 0.49936338107645674464;
  int i;
  int n = 1000000;
  double pi = 3.141592653589793;
  double total;
  double wtime = 0.0;
  double *x;

  cout << "\n";
  cout << "QUAD:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Estimate the integral of f(x) from A to B.\n";
  cout << "  f(x) = 50 / ( pi * ( 2500 * x * x + 1 ) ).\n";
  cout << "  A = " << a << "\n";
  cout << "  B = " << b << "\n";
  cout << "  Exact integral from 0 to 10 is 0.49936338107645674464...\n";
//
//  Load the array X with evenly spaced values between A and B.
//
  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] =  ( ( double ) ( n - i - 1 ) * a   
            + ( double ) (     i     ) * b ) 
            / ( double ) ( n     - 1 );
  }
//
//  Evaluate F at each node, sum and average.
//
//  CHANGE 2) Use the OpenMP timer here to initialize WTIME.
//
  wtime = omp_get_wtime ( );              //added to get initial wtime at 09/06/2016

  total = 0.0;
//
//  CHANGE 3) Fill in these two OpenMP directives here. 
//
# pragma omp parallel private (i) shared (n, x, pi, total)
# pragma omp for reduction (+: total) //added to parallel the for loop at 09/06/2016

  for ( i = 0; i < n; i++ )
  {
    total = total + 50.0 / ( pi * ( 2500.0 * x[i] * x[i] + 1.0 ) );
  }
  total = ( b - a ) * total / ( double ) n;
//
//  CHANGE 4) Call the OpenMP timer here to update WTIME.
//

  wtime = omp_get_wtime ( )-wtime; //added to calculate elapsed wtime at 09/06/2016

//
//  Print quadrature estimate, error, time
//
  error = fabs ( total - exact );

  cout << "\n";
  cout << "  Estimate = " << total << "\n";
  cout << "  Error    = " << error << "\n";
  cout << "  Wallclock time   = " << wtime << "\n";

  delete [] x;

  return 0;
}
