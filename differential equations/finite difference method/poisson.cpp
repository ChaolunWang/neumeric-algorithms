//--------------------------poisson.cpp--------------------------------
/*This codes implements the poisson equation as a 1D boundary value proplem
  A self made library matrix.h was used to help the calculation 
  Programmed by: Chaolun Wang
	     at: 11/16/2016
*/


#include <vector>
#include "matrix.h"

using namespace std;


int main()
{
	int n=10;
	double a0=0, b0=1;
	double h=(b0-a0)/n;

	//condition 1:
	//initialization
	Matrix<double> A(n+1,n+1), b(n+1, 1), c(n+1, 1);
	A.init(0);
	for(int i=0; i<n+1; ++i)
	{
		A[i][i]=-2;
		if(i-1>=0)
			A[i][i-1]=1;
		if(i+1<=n)
			A[i][i+1]=1;
	}
	A[0][0]=h*h;
	A[0][1]=0;
	A[n][n]=h*h;
	A[n][n-1]=0;
	A=A*(1/h/h);
	cout<<"condition 1: \nThe coefficient matrix is:"<<A;
	b.init(-2);
	b[0][0]=0;
	b[n][0]=3;
	//solve
	c=solve(A,b);
	cout<<"The solution form u(0) to u(1) is:"<<c;

	//condition 2:

	Matrix<double> A2(n+2,n+2), b2(n+2, 1);
	A2.init(0);
	for(int i=0; i<n+2; ++i)
	{
		A2[i][i]=-2;
		if(i-1>=0)
			A2[i][i-1]=1;
		if(i+1<=n+1)
			A2[i][i+1]=1;
	}
	A2[0][0]=h*h;
	A2[0][1]=0;
	A2[n+1][n+1]=h*h;
	A2[n+1][n]=0;
	A2[n+1][n-1]=-1*h*h;
	A2=A2*(1/h/h);
	cout<<"condition 2: \nThe coefficient matrix is:"<<A2;
	b2.init(-2);
	b2[0][0]=0;
	b2[n+1][0]=3*2*h;

	c=solve(A2,b2);
	c.truncate(1,0);
	cout<<"The solution form u(0) to u(1) is:"<<c;

	//condition 3 and 4, the coefficient matrix is singular 
	Matrix<double> A3(n+3,n+3);
	A3.init(0);
	for(int i=0; i<n+3; ++i)
	{
		A3[i][i]=-2;
		if(i-1>=0)
			A3[i][i-1]=1;
		if(i+1<=n+2)
			A3[i][i+1]=1;
	}
	A3[0][0]=h*h;
	A3[0][1]=0;
	A3[0][2]=-1*h*h;
	A3[n+2][n+2]=h*h;
	A3[n+2][n+1]=0;
	A3[n+2][n]=-1*h*h;
	A3=A3*(1/h/h);
	cout<<"condition 3 and 4: \nThe coefficient matrix is:"<<A3;
	cout<<"The matrix is singular\n";

}
