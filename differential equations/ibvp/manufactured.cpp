//-------------------------------manufactured.cpp------------------------
/* This file implements the 1D heat equaiton and the way to use the manufactured
   solution to study the convergence and accuracy
   Programmed by Chaolun Wang
	      at 11/30/2016
*/

#include <iostream>
#include<fstream>
#include<cmath>

using namespace std;

void heat1d(double data[], double dt, double a, double b, double s, double xl, double xr, int n, const char* solution, bool save); //1D heat equation with forward euler used for time marching
void initialize(double data[], int n, double xl, double xr);  //Initialize the algorithm with intial condition
double everagel2norm(double data1[], double data2[], int n);   //calculate the average l2 error norm

int main()
{
	double s=1, xl=0, xr=2*M_PI, a=1, b=1;

	//test space convergence

	cout<<"space convergence:\n";
	int nt=10000;
	double dt=s/nt;
	//time grid fix, spatial grid number =40
	int n=40;
	
	double *data=new double[n+1];
	double datatrue[n+1];

	initialize(data, n, xl, xr);
	initialize(datatrue, n, xl, xr);

	heat1d( data, dt, a, b, s, xl, xr, n, "solutions40.txt", 1);
	cout<<"n=40, The average l2 error norm is: "<<everagel2norm(data, datatrue, n+1)<<'\n';
	delete [] data;

	//time grid fix, spatial grid number =80
	n=80;
	
	double *data2=new double[n+1];
	double datatrue2[n+1];

	initialize(data2, n, xl, xr);
	initialize(datatrue2, n, xl, xr);

	heat1d( data2, dt, a, b, s, xl, xr, n, "solutions80.txt", 1);
	cout<<"n=80, The average l2 error norm is: "<<everagel2norm(data2, datatrue2, n+1)<<'\n';
	delete [] data2;

	//time grid fix, spatial grid number =160
	n=160;
	
	double *data3=new double[n+1];
	double datatrue3[n+1];

	initialize(data3, n, xl, xr);
	initialize(datatrue3, n, xl, xr);

	heat1d( data3, dt, a, b, s, xl, xr, n, "solutions160.txt", 1);
	cout<<"n=160, The average l2 error norm is: "<<everagel2norm(data3, datatrue3, n+1)<<'\n';
	delete [] data3;

	//time grid fix, spatial grid number =320
	n=320;
	
	double *data4=new double[n+1];
	double datatrue4[n+1];

	initialize(data4, n, xl, xr);
	initialize(datatrue4, n, xl, xr);

	heat1d( data4, dt, a, b, s, xl, xr, n, "solutions320.txt", 1);
	cout<<"n=320, The average l2 error norm is: "<<everagel2norm(data4, datatrue4, n+1)<<'\n';
	delete [] data4;

	//test time convergence
	cout<<"time convergence:\n";
	n=100;
	s=1;
	double *datat=new double[n+1];
	double datatruet[n+1];
	initialize(datatruet, n, xl, xr);
	
	//time grid num=5000
	initialize(datat, n, xl, xr);
	nt=500;
	dt=s/nt;
	heat1d( datat, dt, a, b, s, xl, xr, n, "solutiont500.txt", 1);
	cout<<"number of time grids=500, The average l2 error norm is: "<<everagel2norm(datat, datatruet, n+1)<<'\n';

	//time grid num=10000
	initialize(datat, n, xl, xr);
	nt=1000;
	dt=s/nt;
	heat1d( datat, dt, a, b, s, xl, xr, n, "solutiont1000.txt", 1);
	cout<<"number of time grids=1000, The average l2 error norm is: "<<everagel2norm(datat, datatruet, n+1)<<'\n';

	//time grid num=20000
	initialize(datat, n, xl, xr);
	nt=2000;
	dt=s/nt;
	heat1d( datat, dt, a, b, s, xl, xr, n, "solutiont2000.txt", 1);
	cout<<"number of time grids=2000, The average l2 error norm is: "<<everagel2norm(datat, datatruet, n+1)<<'\n';

	//time grid num=40000
	initialize(datat, n, xl, xr);
	nt=4000;
	dt=s/nt;
	heat1d( datat, dt, a, b, s, xl, xr, n, "solutiont4000.txt", 1);
	cout<<"number of time grids=4000, The average l2 error norm is: "<<everagel2norm(datat, datatruet, n+1)<<'\n';

	delete [] datat;
}

void heat1d(double data[], double dt, double a, double b, double s, double xl, double xr, int n, const char* solution, bool save)  //heat equation implementation
{
	ofstream sol;
	if(save)
	{	
		sol.open (solution);
	}

	int x=s/dt;
	double dx=(xr-xl)/n;



	double *temp;
	for(int j=0; j<x; ++j)
	{
		temp=new double[n+1];
		for(int i=1; i<=n-1; ++i)
		{
			temp[i]=dt/dx/dx*(data[i-1]-2*data[i]+data[i+1])+dt*sin(xl+dx*i)+data[i];
		}
		temp[0]=a;
		temp[n]=b;
		delete [] data;
		data=temp;


	}
	if(save)
	{
		for(int i=0; i<=n; ++i)
			sol<<data[i]<<' ';
		sol<<'\n';
	}
	if(save)
		sol.close();

}


void initialize(double data[], int n, double xl, double xr)
{
	double dx=(xr-xl)/n;
	for(int i=0; i<=n; ++i)
		data[i]=sin(xl+dx*i)+1;	
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
