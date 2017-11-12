//----------------------------------search.cpp-------------------------
/*This is the cpp file which implement the search algorithm. MPI was used
  for the parrallel of the code. This code is modified by Chaolun Wang at 09/15/2016
*/
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <mpi.h> // added by Chaolun Wang at 09/15/2016 to include the MPI library
//IMPORTANT! be sure to distribute the data to the workers! EACH WORKER MAY ONLY GET THE DATA NEEDED
using namespace std;

int main ( int argc, char *argv[]  );
int search ( int data[], int size, int c );  //modified
int f ( int i );

//****************************************************************************80

int main ( int argc, char *argv[]  )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SEARCH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  int b;
  int c;
  int fj;
  int i4_huge =2147483647;
  int j;
  double wtime = 0.0;
  int q;                                            //added for reduce the variable j

//initialize MPI
  MPI::Init ( argc, argv );
//
//  Get the number of processes.
//
  int p = MPI::COMM_WORLD.Get_size (  );
//
//  Determine this processes's rank.
//
  int id = MPI::COMM_WORLD.Get_rank ( );

  MPI::Request request;                              //declear the mpi request, used for check the buffer of Ireceive


  a = 1;
  b = i4_huge/10;
  c = 3081;
  
  int size=(int)((b-a+1)/p);                          //calculate the work amount for each node
  int leftover=b-a+1-size*(p)+size;                   //calculate the work amount of the last node
 
  if(p==1)                                            //exception when only one process exist
    size=leftover;

  int *partialData= new int[leftover];                //dynamically alocate the data chunck


  if(id==0)
  {
    request = MPI::COMM_WORLD.Irecv(partialData, size, MPI::INT, 0, 0);   //the Ireceive was used for commander
  }
  else if(id==p-1)
  {
    MPI::COMM_WORLD.Recv(partialData, leftover, MPI::INT, 0, 0);           //ordinary Receive was used for worker, for the last node
  }  
  else
  {
    MPI::COMM_WORLD.Recv(partialData, size, MPI::INT, 0, 0);               //ordinary Receive was used for worker
  }


  if(id==0){                                       //if the process id is 0 print out the message, start timing
    cout << "\n";
    cout << "SEARCH:\n";
    cout << "  C++ version\n";
    cout << "  Search the integers from A to B\n";
    cout << "  for a value J such that F(J) = C.\n";
    cout << "\n";
    cout << "  A           = " << a << "\n";
    cout << "  B           = " << b << "\n";
    cout << "  C           = " << c << "\n";
    q=-1;


    int *data=new int[b-a+1];
    for(int i=a; i<=b; ++i)
      data[i-a]=i;

    //sent the last part of array to the last worker
    MPI::COMM_WORLD.Send(data+(p-1)*size, leftover, MPI::INT, p-1, 0);

    //sent the part of array to each of the worker except the last one(also the master)
    for(int i=0; i<p-1; ++i)
    {
      MPI::COMM_WORLD.Send(data+i*size, size, MPI::INT, i, 0);
    }
    delete [] data;      //free the memory

    wtime = MPI::Wtime ( );  //start timing
    //generate the array of interger in master process
  }

  if(id==0)
    request.Wait();                                     //wait the Ireceive to get data before calculation

  if(id==p-1)
    j = search ( partialData, leftover, c );            //search for last node
  else
    j = search ( partialData, size, c );                //search for rest of nodes

  delete [] partialData;                                //free the memory of partial data
  MPI::COMM_WORLD.Reduce ( &j, &q, 1, MPI::INT, MPI::MAX, 0 );

  
  if(id==0)
  {
    wtime = MPI::Wtime ( )-wtime;                      //calculate the elapsed wall time
    if ( q == -1 )
    {
      cout << "\n";
      cout << "  No solution was found.\n";
    }
    else
    {
      cout << "\n";
      cout << "  Found     J = " << q << "\n";
      cout << "  Verify F(J) = " << f ( q ) << "\n";
    }

//cout << "  Elapsed time is " << wtime << "\n";
//
//  Terminate.
//
    cout << "\n";
    cout << "SEARCH:\n";
    cout << "  Normal end of execution. time spent: "<<wtime<<"\n";
  }
  MPI::Finalize ( );
  return 0;
}
//****************************************************************************80

int search ( int data[], int size, int c )                                       //changed so that array and array size are used as variable

//****************************************************************************80
//
//  Purpose:
//
//    SEARCH searches integers in [A,B] for a J so that F(J) = C.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, B, the search range.
//
//    Input, int C, the desired function value.
//
//    Output, int SEARCH, the computed solution, or -1
//    if no solution was found.
//
{
  int fi;
  int i;
  int j;

  j = -1;

  for ( i = 0; i <size; i++ )
  {
    fi = f ( data[i] );

    if ( fi == c )
    {
      j = data[i];
      //break;                                                              //IMPORTANT: Break was removed for timing purpose
    }
  }

  return j;
}
//****************************************************************************80

int f ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    F is the function we are analyzing.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the argument.
//
//    Input, int F, the value.
//
{
  int i4_huge = 2147483647;
  int j;
  int k;
  int value;

  value = i;

  for ( j = 1; j <= 5; j++ )
  {
    k = value / 127773;

    value = 16807 * ( value - k * 127773 ) - k * 2836;

    if ( value <= 0 )
    {
      value = value + i4_huge;
    }
  }

  return value;
}

