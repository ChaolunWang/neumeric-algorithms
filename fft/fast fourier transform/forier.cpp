//------------------------forier.cpp------------------------
/*this file is for the study of fast forier transform and its use
  in signal processing.
  Programmed by: Chaolun Wang
             at: 11/10/2016
rfftf, rfftb, rffti was used
Data format:
r       r(1) = the sum from i=1 to i=n of r(i)         

        if n is even set l =n/2   , if n is odd set l = (n+1)/2

          then for k = 2,...,l

             r(2*k-2) = the sum from i = 1 to i = n of

                  r(i)*cos((k-1)*(i-1)*2*pi/n)

             r(2*k-1) = the sum from i = 1 to i = n of

                 -r(i)*sin((k-1)*(i-1)*2*pi/n)

        if n is even

             r(n) = the sum from i = 1 to i = n of

                  (-1)**(i-1)*r(i)


*/

#include "fftpack4.h"
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <random>
using namespace std;

void sampling (double (*f)(double), int n, double min, double max, double a[]);
double f1(double t);
double f2(double t);

int main()
{

	int *ifac;
	int n;
	double *r;
	double *rt;
	double *wsave;
	random_device rd;
	mt19937 generator(rd());
	normal_distribution<double> distribution(0, 2);
	//part 1
	for(int k=0; k<3; ++k )
	{
		n=16*pow(2, k);
		ifac = ( int * ) malloc ( 8 * sizeof ( int ) );
		wsave = ( double * ) malloc ( 2 * ( n + 1 ) * sizeof ( double ) );

		rffti ( &n, wsave, ifac );
		r=new double[n];
		sampling(f1, n, -1*M_PI, M_PI, r);

		//Compute the FFT coefficients.
		rfftf ( &n, r, wsave, ifac );
		//Compute inverse FFT of coefficients.  Should get back the original data.

		rfftb ( &n, r, wsave, ifac );
		//normalization
		for ( int i = 0; i < n; i++ )
		{
			r[i] = r[i] /n;
		}
		cout<<"\nThe discrete Fourier approximations for n = "<<n<<" is:\n";
		for ( int i = 0; i < n; i++ )
		{
			cout<<r[i]<<' ';
		}
		cout<<'\n';
 		free ( wsave );
		free ( ifac );
		delete [] r;	
	}
	//n==64 calculate the spectural
	ifac = ( int * ) malloc ( 8 * sizeof ( int ) );
	wsave = ( double * ) malloc ( 2 * ( n + 1 ) * sizeof ( double ) );
	rffti ( &n, wsave, ifac );
	r=new double[n];
	sampling(f1, n, -1*M_PI, M_PI, r);

	//Compute the FFT coefficients.
	rfftf ( &n, r, wsave, ifac );


	cout<<"\nn==64 calculate the spectural: \n";
	rt=new double[n/2];
	rt[0]=abs(r[0]);
	for(int i=1; i<n/2; ++i)
	{
		rt[i]=sqrt(r[i*2-1]*r[i*2-1]+r[i*2]*r[i*2]);
	}

	for(int i=0; i<n/2; ++i)
	{
		cout<<rt[i]<<' ';
	}
	cout<<'\n';
 	free ( wsave );
	free ( ifac );
	delete [] r;	
	delete [] rt;
	




	//part 2
	n=210;
	//clear signal
	ifac = ( int * ) malloc ( 8 * sizeof ( int ) );
	wsave = ( double * ) malloc ( 2 * ( n + 1 ) * sizeof ( double ) );

	rffti ( &n, wsave, ifac );
	r=new double[n];
	sampling(f2, n, 0, 1, r);
	cout<<"\n\nPart 2:\nThe function sample points of clear signal:\n";
	for ( int i = 0; i < n; i++ )
	{
		cout<<r[i]<<' ';
	}
	cout<<'\n';
	//Compute the FFT coefficients.
	rfftf ( &n, r, wsave, ifac );
	free ( wsave );
	free ( ifac );
	rt=new double[n/2];
	rt[0]=abs(r[0]);
	for(int i=1; i<n/2; ++i)
	{
		rt[i]=sqrt(r[i*2-1]*r[i*2-1]+r[i*2]*r[i*2]);
	}
	cout<<"\ncalculate the spectural: \n";
	for(int i=0; i<n/2; ++i)
	{
		cout<<rt[i]<<' ';
	}
	cout<<'\n';
 	delete [] r;
	delete [] rt;
	
	

	//noize signal
	ifac = ( int * ) malloc ( 8 * sizeof ( int ) );
	wsave = ( double * ) malloc ( 2 * ( n + 1 ) * sizeof ( double ) );

	rffti ( &n, wsave, ifac );
	r=new double[n];
	sampling(f2, n, 0, 1, r);
	//add noize with standard deviation to be 2
	for ( int i = 0; i < n; i++ )
	{
		r[i]+=distribution(generator);//*sqrt(3);
	}

	cout<<"\n\nThe function sample points of noize signal:\n";
	for ( int i = 0; i < n; i++ )
	{
		cout<<r[i]<<' ';
	}
	cout<<'\n';
	//Compute the FFT coefficients.
	rfftf ( &n, r, wsave, ifac );
	free ( wsave );
	free ( ifac );
	rt=new double[n/2];
	rt[0]=abs(r[0]);
	for(int i=1; i<n/2; ++i)
	{
		rt[i]=sqrt(r[i*2-1]*r[i*2-1]+r[i*2]*r[i*2]);
	}
	cout<<"\ncalculate the spectural: \n";
	for(int i=0; i<n/2; ++i)
	{
		cout<<rt[i]<<' ';
	}
	cout<<'\n';
 	delete [] r;
	delete [] rt;
}



void sampling (double (*f)(double), int n, double min, double max, double a[])
{
	double h=(max-min)/(n);
	for (int i=0; i<n; ++i)
		a[i]=(*f)(min+h*i);
}

double f1(double t)
{
	return t*t*cos(t);
}

double f2(double t)
{
	return sin(2*M_PI*25*t)+sin(2*M_PI*80*t)+sin(2*M_PI*125*t)+sin(2*M_PI*240*t)+sin(2*M_PI*315*t);
}




