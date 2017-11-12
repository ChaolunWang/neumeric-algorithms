//------------------------fs1.cpp------------------------------
/*This is the cpp file which evaluate the first function f(x)=x*x+x
  In this file, the fourier series approximation method was implemented
  The out put will contain the fourier coefficients which will then used
  for the plotting*/
#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

const double TOL=8;                                                       //the max laier of romberg method



double trapezoid(int n, double a, double b, int k, int index, double L);   //trapezoid rule
double romberg(vector<vector<double> > & input, double a, double b, int k, int index, double L);  //romberg method 
double integral(double a, double b, int k, int index);                     //using romberg method to find integration

vector<double> fourier(double L, int n);                                //fourier series for any function any region any n
double f(double x, int n, int index, double L);                         //function to be evaluated                         


int main()
{
	cout<<"Function 1: \n";
	for (int k=1; k<5; ++k)
	{
		vector<double> result=fourier(M_PI, pow(2,k));
		cout<<"\nThe fourier coefficients for n= "<<pow(2,k)<<" are(in order of a0 a1-an b1-bn):\n";
		for(int i=0; i<result.size(); ++i)
			cout<<result[i]<<' ';
		cout<<'\n';
	}


}

double f(double x, int n, int index, double L)
{
	double k=x/M_PI*L;
	if(index==0)
		return k*k+k;                                 //function to be evaluated
	if(index==1)
		return (k*k+k)*cos(n*x);                       //intergrand to calculate a
	if(index==2)
		return (k*k+k)*sin(n*x);                       //intergrand to calculate b
}


double trapezoid(int n, double a, double b, int k, int index, double L)
{
	double result=0, h=(b-a)/n;
	for(int i=1; i<n; ++i)
	    result+=f(a+i*h, k, index, L);
	result+=(f(a, k, index, L)+f(b, k, index, L))/2;
	return result*h;
}

double romberg(vector<vector<double> > & input, double a, double b, int k, int index, double L)
{
	double result;
	int m=input.size(), n;
	if(m==0)
        n=0;
	else
	    n=input.back().size();
	if(m==n)
	{
		result=trapezoid(pow(2,m), a, b, k, index, L);
		vector<double> newline;
		newline.push_back(result);
		input.push_back(newline);
	}
	else
	{
		result=(pow(4,n)*input[m-1][n-1]-input[m-2][n-1])/(pow(4,n)-1);
		input.back().push_back(result);
	}
	return result;
}

double integral(double a, double b, int k, int index, double L)
{
	vector<vector<double> > seed;
	double result;
	double temp;

	do
	{
		temp=romberg(seed, a, b, k, index, L);
		if(seed.size()>TOL)                                       //stopping criterial
			break;
		result=temp;
	} while(true);
	seed.clear();
	return result;
}

vector<double> fourier(double L, int n)
{
	
	vector<double> res;

	res.push_back(integral(-1*M_PI, M_PI, 0, 0, L)/M_PI);
	
	for(int i=1; i<=n; ++i)
		res.push_back(integral(-1*M_PI, M_PI, i, 1, L)/M_PI);
	for(int i=1; i<=n; ++i)
		res.push_back(integral(-1*M_PI, M_PI, i, 2, L)/M_PI);
	return res;
}



