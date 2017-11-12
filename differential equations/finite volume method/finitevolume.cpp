//----------------------finitevolume.cpp---------------------------------
/*This documents involves the implementation of 1D finite voule method, and 
  It has been used on solving the advection problem
  Programmed by: Chaolun Wang
	     at: 11/28/2016
*/

#include<fstream>
#include<iostream>
#include<cmath>

using namespace std;

void finitevol(double data[], double dx, double dt, double a, double b, double v, int n, const char* solution, const char* TV); //finite volume solver
void initialize(double data[], int n, double dx);                      //initialize the initial value
double everagel2norm(double data1[], double data2[], int n);           //calculate average l2 norm

int main()
{
	int n;
	double dt, dx, a, b, v, C;
	a=0;
	b=10;
	n=40;
	v=2;

	//tv=0.4
	C=0.4;
	dx=10.0/n;
	dt=C*dx/v;
	double data[n];
	initialize(data, n, dx);
	finitevol(data, dx, dt, a, b, v, n, "solution0.4", "tv0.4");

	//tv=0.8
	C=0.8;
	dx=10.0/n;
	dt=C*dx/v;
	initialize(data, n, dx);
	finitevol(data, dx, dt, a, b, v, n, "solution0.8", "tv0.8");

	//tv=1
	C=1.0;
	dx=10.0/n;
	dt=C*dx/v;
	initialize(data, n, dx);
	finitevol(data, dx, dt, a, b, v, n, "solution1", "tv1");

	//tv=1.2
	C=1.2;
	dx=10.0/n;
	dt=C*dx/v;
	initialize(data, n, dx);
	finitevol(data, dx, dt, a, b, v, n, "solution1.2", "tv1.2");

	//calculate the order
	cout<<"The following cases, dt=0.005:\n";
	//n=80
	dt=0.005;
	n=80;
	dx=10.0/n;
	double data080[n],data80[n];
	initialize(data080, n, dx);
	initialize(data80, n, dx);
	finitevol(data80, dx, dt, a, b, v, n, "solution80", "tv80");
	cout<<"When n=80, the everage error norm is: "<<everagel2norm(data080, data80, n)<<'\n';

	//n=160
	n=160;
	dx=10.0/n;
	double data016[n],data16[n];
	initialize(data016, n, dx);
	initialize(data16, n, dx);
	finitevol(data16, dx, dt, a, b, v, n, "solution16", "tv16");
	cout<<"When n=160, the everage error norm is: "<<everagel2norm(data016, data16, n)<<'\n';

	//n=320
	n=320;
	dx=10.0/n;
	double data032[n],data32[n];
	initialize(data032, n, dx);
	initialize(data32, n, dx);
	finitevol(data32, dx, dt, a, b, v, n, "solution32", "tv32");
	cout<<"When n=320, the everage error norm is: "<<everagel2norm(data032, data32, n)<<'\n';

	//n=640
	n=640;
	dx=10.0/n;
	double data064[n],data64[n];
	initialize(data064, n, dx);
	initialize(data64, n, dx);
	finitevol(data64, dx, dt, a, b, v, n, "solution64", "tv64");
	cout<<"When n=640, the everage error norm is: "<<everagel2norm(data064, data64, n)<<'\n';
	
}

void initialize(double data[], int n, double dx)
{
	for(int i=0; i<n; ++i)
	{
		if(abs(-5+i*dx)>1)
			data[i]=0;
		else
			data[i]=1;
	}
}

void finitevol(double data[], double dx, double dt, double a, double b, double v, int n, const char* solution, const char* TV)
{
	ofstream sol, tv;
	sol.open (solution);
	tv.open(TV);

	int x=(b-a)/dt;
	for(int i=0; i<n; ++i)
		sol<<data[i]<<' ';
	sol<<'\n';
	double temp=0;
	for(int i=0; i<n; ++i)
	{
		if(i==n-1)
			temp+=abs(data[0]-data[i]);
		else
			temp+=abs(data[i+1]-data[i]);
	}
	tv<<temp<<'\n';
	for(int j=0; j<x; ++j)
	{
		double f[n];
		for(int i=0; i<n; ++i)
			f[i]=v*data[i];
		for(int i=0; i<n; ++i)
		{
			if(i==0)
				data[i]-=dt/dx*(f[0]-f[n-1]);
			else
				data[i]-=dt/dx*(f[i]-f[i-1]);
			sol<<data[i]<<' ';
		}
		sol<<'\n';
		temp=0;
		for(int i=0; i<n; ++i)
		{
			if(i==n-1)
				temp+=abs(data[0]-data[i]);
			else
				temp+=abs(data[i+1]-data[i]);
		}
		tv<<temp<<'\n';
	}


	sol.close();
	tv.close();
}

double everagel2norm(double data1[], double data2[], int n)
{
	double result=0;
	for(int i=0; i<n; ++i)
		result+=pow(data1[i]-data2[i],2);
	result=sqrt(result);
	result/=n;
	return result;
}
