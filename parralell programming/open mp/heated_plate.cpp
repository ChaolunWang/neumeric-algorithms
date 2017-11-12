//----------------heated_plate.cpp-----------------
/*This code is from the lab material of ACSII in fall 2016, directed by Dr. Plewa. It implement a simulation of temperature change of a heated plate
  This program have been modified by Chaolun Wang at 09/06/2016 to let it be able to run parallely. Open_MP was used for the implementation. Also the 
  output of threads number and elapsed wall time has been added for the timing and recording purpose.
*/
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>
# include <omp.h>                    //added by Chaolun Wang at 09/06/2016 to include Open_MP library
using namespace std;

int main ( int argc, char *argv[] );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for HEATED_PLATE.
//
//  Discussion:
//
//    This code solves the steady state heat equation on a rectangular region.
//
//    The sequential version of this program needs approximately
//    18/epsilon iterations to complete. 
//
//
//    The physical region, and the boundary conditions, are suggested
//    by this diagram;
//
//                   W = 0
//             +------------------+
//             |                  |
//    W = 100  |                  | W = 100
//             |                  |
//             +------------------+
//                   W = 100
//
//    The region is covered with a grid of M by N nodes, and an N by N
//    array W is used to record the temperature.  The correspondence between
//    array indices and locations in the region is suggested by giving the
//    indices of the four corners:
//
//                  I = 0
//          [0][0]-------------[0][N-1]
//             |                  |
//      J = 0  |                  |  J = N-1
//             |                  |
//        [M-1][0]-----------[M-1][N-1]
//                  I = M-1
//
//    The steady state solution to the discrete heat equation satisfies the
//    following condition at an interior grid point:
//
//      W[Central] = (1/4) * ( W[North] + W[South] + W[East] + W[West] )
//
//    where "Central" is the index of the grid point, "North" is the index
//    of its immediate neighbor to the "north", and so on.
//   
//    Given an approximate solution of the steady state heat equation, a
//    "better" solution is given by replacing each interior point by the
//    average of its 4 neighbors - in other words, by using the condition
//    as an ASSIGNMENT statement:
//
//      W[Central]  <=  (1/4) * ( W[North] + W[South] + W[East] + W[West] )
//
//    If this process is repeated often enough, the difference between successive 
//    estimates of the solution will go to zero.
//
//    This program carries out such an iteration, using a tolerance specified by
//    the user, and writes the final estimate of the solution to a file that can
//    be used for graphic processing.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 July 2008
//
//  Author:
//
//    Original C version by Michael Quinn.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Michael Quinn,
//    Parallel Programming in C with MPI and OpenMP,
//    McGraw-Hill, 2004,
//    ISBN13: 978-0071232654,
//    LC: QA76.73.C15.Q55.
//
//  Parameters:
//
//    Commandline argument, float EPSILON, the error tolerance.  
//
//  Local parameters:
//
//    Local, double DIFF, the norm of the change in the solution from one iteration
//    to the next.
//
//    Local, double MEAN, the average of the boundary values, used to initialize
//    the values of the solution in the interior.
//
//    Local, double U[M][N], the solution at the previous iteration.
//
//    Local, double W[M][N], the solution computed at the latest iteration.
//
{
# define M 100
# define N 500

  double diff;
  float epsilon;
  int i;
  int iterations;
  int iterations_print;
  int j;
  double mean;
  int success;
  double u[M][N];
  double w[M][N];
  int threads_num;    //added by Chaolun Wang at 09/06/2016 to indicate threads number

  cout << "\n";
  cout << "HEATED_PLATE\n";
  cout << "  C++ version\n";
  cout << "  A program to solve for the steady state temperature distribution\n";
  cout << "  over a rectangular plate.\n";
  cout << "\n";
  cout << "  Spatial grid of " << M << " by " << N << " points.\n";
// 
//  Read EPSILON from the command line or the user.
//
  if ( argc < 2 ) 
  {
    cout << "\n";
    cout << "  Enter EPSILON, the error tolerance:\n";
    cin >> epsilon;
  }
  else
  {
    success = sscanf ( argv[1], "%f", &epsilon );

    if ( success != 1 )
    {
      cout << "\n";
      cout << "HEATED_PLATE\n";
      cout << "  Error reading in the value of EPSILON.\n";
      return 1;
    }
  }

  cout << "\n";
  cout << "  The iteration will be repeated until the change is <= " 
       << epsilon << "\n";
  diff = epsilon;
//start the timing of program
  double wtime = omp_get_wtime ( );            //added by Chaolun Wang at 09/06/2016 to start timing of the program
// 
//  Set the boundary values, which don't change. 
//
  #pragma omp parallel                        //added by Chaolun Wang at 09/06/2016 to parallel the for loop
  #pragma omp for
  for ( i = 1; i < M - 1; i++ )
  {
    threads_num=omp_get_num_threads ( );      //added by Chaolun Wang at 09/06/2016 to get number of threads
    w[i][0] = 100.0;
  }

  #pragma omp parallel                        //added by Chaolun Wang at 09/06/2016 to parallel the for loop
  #pragma omp for
  for ( i = 1; i < M - 1; i++ )
  {
    w[i][N-1] = 100.0;
  }

  #pragma omp parallel                        //added by Chaolun Wang at 09/06/2016 to parallel the for loop
  #pragma omp for
  for ( j = 0; j < N; j++ )
  {
    w[M-1][j] = 100.0;
  }

  #pragma omp parallel                         //added by Chaolun Wang at 09/06/2016 to parallel the for loop
  #pragma omp for
  for ( j = 0; j < N; j++ )
  {
    w[0][j] = 0.0;
  }
//
//  Average the boundary values, to come up with a reasonable
//  initial value for the interior.
//
  mean = 0.0;

  #pragma omp parallel                        //added by Chaolun Wang at 09/06/2016 to parallel the for loop with the reduction on mean
  #pragma omp for reduction(+:mean)
  for ( i = 1; i < M - 1; i++ )
  {
    mean = mean + w[i][0];
  }

  #pragma omp parallel                        //added by Chaolun Wang at 09/06/2016 to parallel the for loop with the reduction on mean
  #pragma omp for reduction(+:mean)
  for ( i = 1; i < M - 1; i++ )
  {
    mean = mean + w[i][N-1];
  }

  #pragma omp parallel                        //added by Chaolun Wang at 09/06/2016 to parallel the for loop with the reduction on mean
  #pragma omp for reduction(+:mean)
  for ( j = 0; j < N; j++ )
  {
    mean = mean + w[M-1][j];
  }

  #pragma omp parallel                        //added by Chaolun Wang at 09/06/2016 to parallel the for loop with the reduction on mean
  #pragma omp for reduction(+:mean)
  for ( j = 0; j < N; j++ )
  {
    mean = mean + w[0][j];
  }
  mean = mean / ( double ) ( 2 * M + 2 * N - 4 );
// 
//  Initialize the interior solution to the mean value.
//
  #pragma omp parallel private (i, j) shared (w, mean) //added by Chaolun Wang at 09/06/2016 to parallel the nested for loop
  #pragma omp for
  for ( i = 1; i < M - 1; i++ )
  {
    for ( j = 1; j < N - 1; j++ )
    {
      w[i][j] = mean;
    }
  }
// 
//  iterate until the  new solution W differs from the old solution U
//  by no more than EPSILON.
//
  iterations = 0;
  iterations_print = 1;
  cout << "\n";
  cout << " Iteration  Change\n";
  cout << "\n";

  while ( epsilon <= diff )                                                //begining of iteration
  {
//
//  Save the old solution in U.
//
    #pragma omp parallel private(i,j) shared(u,w)                          //added by Chaolun Wang at 09/06/2016 to parallel the nested for loop
    #pragma omp for
    for ( i = 0; i < M; i++ ) 
    {
      for ( j = 0; j < N; j++ )
      {
        u[i][j] = w[i][j];
      }
    }
//
//  Determine the new estimate of the solution at the interior points.
//  The new solution W is the average of north, south, east and west neighbors.
//
    #pragma omp parallel private(i,j) shared(u,w)                          //added by Chaolun Wang at 09/06/2016 to parallel the nested for loop
    #pragma omp for
    for ( i = 1; i < M - 1; i++ )
    {
      for ( j = 1; j < N - 1; j++ )
      {
        w[i][j] = ( u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1] ) / 4.0;
      }
    }

    diff = 0.0;
    #pragma omp parallel private(i,j) shared(u,w)                          //added by Chaolun Wang at 09/06/2016 to parallel the nested for loop with the reduction on diff
    #pragma omp for reduction(+:diff)
    for ( i = 1; i < M - 1; i++ )
    {
      for ( j = 1; j < N - 1; j++ )
      {
        diff = diff + fabs ( w[i][j] - u[i][j] );
      }
    }
    diff = diff / ( double ) ( M - 1 ) / ( double ) ( N - 1 );

    iterations++;
    if ( iterations == iterations_print )
    {
      cout << "  " << setw(8) << iterations
           << "  " << diff << "\n";
      iterations_print = 2 * iterations_print;
    }
  }                                                                      //end of iteration

  cout << "\n";
  cout << "  " << setw(8) << iterations
       << "  " << diff << "\n";
  cout << "\n";
  cout << "  Error tolerance achieved.\n";
//end of timing
  wtime = omp_get_wtime ( )-wtime;                                       //added by Chaolun Wang at 09/06/2016 to end the timing and update wtime
// 
//  Terminate.
//
  cout << "\n";
  cout << "HEATED_PLATE:\n";
  cout << "  Normal end of execution.\n";
  cout << "  The number of threads available in parallel:"<<threads_num<<"\n";      //added by Chaolun Wang at 09/06/2016 to output number of threads in parallell
  cout << "  Elapsed wall clock time = " << wtime << "\n";                          //added by Chaolun Wang at 09/06/2016 to output the elapsed wall time
  return 0;

# undef M
# undef N
}
