//--------------------------randomnumber.cpp-----------------------
/*In this code the randomnumber generator using linear congruential method was implemented. Then
  a SND randon number generator was made and the chi-square was calculated. for lab3 of acsII
  Programmed by: Chaolun Wang
	     at: 09/20/2016
*/

#include <cmath>
#include <omp.h>
#include <iostream>
#include <iomanip>
using namespace std;

const int A=1664525;                          //H.W.Lewis rand
const int B=1013904223;
const int E=32;                               //interger is 4 bytes which is 32 bits
const int N=100000000;
const int T=6;
const double pmmu=0.68269;                    //from online z table, the probability in mu+-theta
const double pmmu05=0.19146;                  //from online z table, the probability for m-0.5*theta to 0

double nextrand(int m, int &seed);                                        //calculate the next random number from 0 to 1
double normalDist(double x, double mu, double va);                        //the normal distribution function
double test(double x);                                                    //test function for accept and reject method
double normalRand(int m, int &seed, int precise, double mu, double va);   //random number generator which can generate normal random distributed numbers
double chiSq(int size, int o[], double e[]);                              //function to perform chi square test

int main()
{
	int m=pow(2,E);                       //the number of integer the machine can represent
	int seed=omp_get_wtime()*100;        //generate the seed
	double chi=0;
	double precision=5;
	int i=0;
	double va=1, mu=0;
	double the=sqrt(va);
	double expect[T];
        int observed[T];
        //part1 testing random number generator
	cout<<"The SND random number generator was used to generate 20 numbers:\n";
	for(int i=0; i<20; ++i)
	{		
		cout<<normalRand(m, seed, precision, mu, va)<<'\n';
	}
        //part2 chi-square analysis
	for(int i=0; i<T; ++i)
		observed[i]=0;
	expect[1]=expect[4]=pmmu*N/2.0-pmmu05*N;
	expect[0]=expect[5]=(1.0-pmmu)*N/2.0;
	expect[2]=expect[3]=pmmu05*N;
	for(int i=0; i<N; ++i)
	{		
		double temp=normalRand(m, seed, precision, mu, va);
		if(temp<=mu-the)
			++observed[0];
		else if(temp>mu-the && temp<=mu-0.5*the)
			++observed[1];
		else if(temp<=mu && temp>mu-0.5*the)
			++observed[2];
		else if(temp>mu && temp<mu+0.5*the)
			++observed[3];
		else if(temp<mu+the && temp>=mu+0.5*the)
			++observed[4];
		else
			++observed[5];
	}
        cout<<N<<" random numbers with normal distribution have been generated\n";
	cout<<"For each of the "<<T<<" bins:\n";
	cout<<"region id           observed            expected            \n";
	for(int i=0; i<T; ++i)
		cout<<setw(20)<<left<<i+1<<setw(20)<<observed[i]<<setw(20)<<expect[i]<<'\n';
	cout<<"Performing Chi-square test, result is:\n"<<chiSq(T, observed, expect)<<'\n';
}

//calculate the next random number from 0 to 1
double nextrand(int m, int &seed)
{
	//do
	seed=(A*seed+B)%(m);
	//while (seed==0 ||seed==m);
	return abs((double)seed/m);
}

//the normal distribution function
double normalDist(double x, double mu, double va)
{
	return exp(-1.0*pow(x-mu,2.0)/2.0/va)/(sqrt(2.0*va*M_PI));
}

//test function for accept and reject method
double test(double x)
{
	return normalDist(0,0,1);
}

//random number generator which can generate normal random distributed numbers
double normalRand(int m, int &seed, int precise, double mu, double va)
{
	double y1,y2;
	while(true)
	{
		y2=(nextrand(m, seed)-0.5)*2*precise*sqrt(va)+mu;
		y1=nextrand(m, seed);
		if(y1<normalDist(y2, mu, va)/test(y2))
			return y2;
	}
}

//function to perform chi square test
double chiSq(int size, int o[], double e[])
{
	double result=0;
	for(int i=0; i<size; ++i)
		result+=pow((double)o[i]-e[i],2)/e[i];
	return result;
}
