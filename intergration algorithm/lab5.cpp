//----------------------------lab5.cpp-------------------------------------
/*this is the implementation of using Gauss-Kronrod method to solve numerical
  intergral problem and do some analysis related to that.
  Programmed by: Chaolun Wang
             at: 03/19/2016
  (c++ file kronrod.cpp/kronrod.hpp is from John Burkardt:
    http://people.sc.fsu.edu/jburkardt/c_src/kronrod/kronrod.html)
 */

#include "kronrod.hpp"                                              //include the code form John Burkardt
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

const double TOL=0.000001;                                          //set tolerance to 10^-6
const double PI=3.141592653;

double f(double x, double the);                                      //the intergrand function
double k_x(double b, double a, double w, double the, double x);      // function K(x)
double integrateKx(int &n, double x);                                //doing an intergration on K(x) using Gauss-Kronrod

int main()
{
    //for question 1 and 2
	int n;
	double result, step=1/100.0;
	double * w1;                                                     //set up pointers to store weight and nodes
	double * w2;
	double * the;
	ofstream output;
    output.open ("q12result.txt");                  
    
    for(int i=0; i<=100; ++i)
    {
		n=0;                                                         
		if(i<100)
		    result=integrateKx(n, i*step);
		else                                                           //because if i=100 the value of the K(x) will be infinity, the x was setted to be 1-TOL
		    result=integrateKx(n, i*step-TOL);
		output<<i*step<<' '<<result<<' '<<n<<'\n';                    //output data to "q12result.txt"

	}
	output.close();
	n=0;
	//for question 3
	output.open ("q3_1.txt");
	result=integrateKx(n, 0.5);

    step=PI/2/100;
    for(int i=0; i<=100; ++i)                              //for plotting intergrand
        output<<i*step<<' '<<f(0.5, i*step)<<'\n';
    output.close();
	output.open ("q3_2.txt");  
	w1 = new double[n+1];
	w2 = new double[n+1];
	the = new double[n+1];
	kronrod ( n, TOL, the, w1, w2 );
	for(int i=1; i<n+1; i+=2)
	{
		if(i==n)                                              //export gauss node/weight
			output<<PI/2/2*(1+the[i])<<' '<<w2[i]<<'\n';                 
		else
			output<<PI/2/2*(1+the[i])<<' '<<w2[i]<<'\n'<<PI/2/2*(1-1*the[i])<<' '<<w2[i]<<'\n';
	}
    output.close();	
	output.open ("q3_3.txt");  	
	for(int i=0; i<n+1; ++i)                                 //export kernoid node/weight
	{
		if(i==n)
            output<<PI/2/2*(1+the[i])<<' '<<w1[i]<<'\n';
		else
			output<<PI/2/2*(1+the[i])<<' '<<w1[i]<<'\n'<<PI/2/2*(1-1*the[i])<<' '<<w1[i]<<'\n';
	}
    output.close();	
  	delete [] w1;
	delete [] w2;
	delete [] the;	  
	cout<<"claculaion complete, result exported to txt files, execute the lab5plot.m file to plot result\n";
}



double f(double x, double the)                                  //the intergrand function
{
	return 1.0/sqrt(1-x*x*sin(the)*sin(the));
}

double k_x(double b, double a, double w, double the, double x)     // function K(x)
{
	return (b-a)/2*w*f(x,(b-a)/2*the+(a+b)/2);
}

double integrateKx(int &n, double x)                              //doing an intergration on K(x) using Gauss-Kronrod
{
	double err, result;
	double * w1;
	double * w2;
	double * the;
	do
	{
		err=result=0;
		++n;
		w1 = new double[n+1];
		w2 = new double[n+1];
		the = new double[n+1];
		kronrod ( n, TOL, the, w1, w2 );
		for(int i=1; i<n+1; i+=2)
		{
			if(i==n)
			    err+=k_x(PI/2, 0, w2[i], the[i], x);
			else
			    err+=k_x(PI/2, 0, w2[i], -1*the[i], x)+k_x(3.141592653/2, 0, w2[i], the[i], x);      //integrate over 0 and pi/2
		}
		for(int i=0; i<n+1; ++i)
		{
			if(i==n)
			    result+=k_x(PI/2, 0, w1[i], the[i], x);
			else
			    result+=k_x(PI/2, 0, w1[i], the[i], x)+k_x(3.141592653/2, 0, w1[i], -1*the[i], x);    //integrate over 0 and pi/2
		}
		err=abs((err-result)/result);
		delete [] w1;
		delete [] w2;
		delete [] the;
	}while(err>TOL);
	return result;
}









