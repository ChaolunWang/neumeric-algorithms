//--------------------------poissonmul.cpp--------------------------------
/*This codes implements the poisson equation as a 1D boundary value proplem
  A self made library matrix.h was used to help the calculation. The equation
  was solved using a weighted jacobi approach with multigrid to remove error
  of different frequency 
  Programmed by: Chaolun Wang
	     at: 12/10/2016
*/


#include <vector>
#include "matrix.h"

using namespace std;

void weightedjac(const Matrix<double>& A, Matrix<double>& x, const Matrix<double>& b, int n, int loops, double w); //the solver using weighted jacobi method
Matrix<double> Iht2h(int n); 								                           //the I2h to h initializer, input is the size of 2h
double l2enorm(const Matrix<double>& a, const Matrix<double>& b, int n); 					   //calculate the l2 error norm

int main()
{
	int n=255, dn=250, N=2000;
	double a0=0, b0=1;
	double h=(b0-a0)/n;
	double u0=1, u1=exp(2), w=0.8;


	//initialize
	Matrix<double> A(n,n), b(n, 1), x(n, 1), I(n, n/2), R(n/2, n), act(n,1), rh(n, 1), r2h(n/2,1), A2h(n/2, n/2), E(n/2, 1);
	A.init(0);
	for(int i=0; i<n; ++i)
	{
		A[i][i]=-2;
		if(i-1>=0)
			A[i][i-1]=1;
		if(i+1<=n-1)
			A[i][i+1]=1;
	}
	A[0][0]=-1*h*h;
	A[0][1]=0;
	A=A*(-1.0/h/h);
	for(int i=1; i<n; ++i)
	{
		b[i][0]=-4*exp(2*(a0+i*h));
	}
	for(int i=0; i<n; ++i)
	{
		act[i][0]=exp(2*(a0+i*h));
	}

	b[0][0]=u0;
	b[n-1][0]=b[n-1][0]+u1/h/h;
	x.init(1);
	I=Iht2h(n);
	R=0.5*I.t();
	//Oridnary weighted Jacobi
	for(int i=0; i<N/dn; ++i)
	{
		weightedjac(A, x, b, n, dn, w);
		cout<<"l2 residual norm for n="<<(i+1)*dn<<": "<<l2enorm(b, A*x, n)<<'\n';
	}
	cout<<"\nUsing weighted Jacobi method, The solution form u(0) to u(1) is:"<<x<<u1<<'\n';
	cout<<"compair to actual result, the l2 error norm is: "<<l2enorm(act, x, n)<<"\n\n";
	//multigrid method
	x.init(1);
	for(int i=0; i<N/dn/2; ++i)
	{
		weightedjac(A, x, b, n, dn, w);
		rh=b-A*x;
		cout<<"l2 residual norm for n="<<(i*2+1)*dn<<"(fine grid): "<<l2enorm(b, A*x, n)<<'\n';
		r2h=R*rh;
		A2h=R*A*I;
		E.init(0);
		weightedjac(A2h, E, r2h, n/2, dn, w);
		x=x+I*E;
		cout<<"l2 residual norm for n="<<(i*2+2)*dn<<"(coarse grid): "<<l2enorm(b, A*x, n)<<'\n';
	}
	cout<<"\nUsing multigrid weighted Jacobi method, The solution form u(0) to u(1) is:"<<x<<u1<<'\n';
	cout<<"compair to actual result, the l2 error norm is: "<<l2enorm(act, x, n)<<'\n';
}

void weightedjac(const Matrix<double>& A, Matrix<double>& x, const Matrix<double>& b, int n, int loops, double w)
{
	Matrix<double> D(n,n), R(n,n);
	D.init(0);
	for(int i=0; i<n; ++i)
		D[i][i]=A[i][i];
	R=A-D;
	for(int i=0; i<n; ++i)
		D[i][i]=1.0/D[i][i];
	for(int i=0; i<loops; ++i)
		x=w*D*(b-R*x)+(1-w)*x;

}

Matrix<double> Iht2h(int n)
{
	Matrix<double> result(n, n/2);
	result.init(0);
	for(int i=0; i<n/2; ++i)
	{
		result[i*2][i]=1;
		result[i*2+1][i]=2;
		result[i*2+2][i]=1;
	}
	return 0.5*result;
}

double l2enorm(const Matrix<double>& a, const Matrix<double>& b, int n)
{
	double result=0;
	for(int i=0; i<n; ++i)
	{
		result+=pow(a[i][0]-b[i][0],2.0);
	}
	return sqrt(result);
}
