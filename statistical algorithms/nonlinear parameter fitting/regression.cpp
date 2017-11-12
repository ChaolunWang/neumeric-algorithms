//--------------------------regression.cpp--------------------------------
/*This codes implements the linear-nonlinear model inference algorithm for the two given dataset. The linear algibra was handled with outside lapack functions.
  A self made library matrix.h was used to help the calculation 
  Programmed by: Chaolun Wang
	     at: 10/10/2016
*/


#include <vector>
#include "matrix.h"

using namespace std;
const double TOL=0.0001; 

vector<double> storeDouble(const char* filename);                                                       //can read the file character by character and then converted into double format. 


int main()
{

        //question 1
	cout<<"question1:\n";
  	
	int info;
	char uplo = 'u';
	double at=0;
	vector<double> data=storeDouble("sample1.dat");                                                 //read data form data file
	vector<double> y, x;
	for(int i=0; i<data.size(); i+=2)                                                                //splitting 
	{
		x.push_back(data[i]);
		y.push_back(data[i+1]);
	}

	int n = y.size();
	Matrix<double> A(n, 2),b(n,1), X(2,1);

	for(int i=0; i<n; ++i)                                                                           //splitting 
	{
		A[i][0]=1;
		A[i][1]=x[i];
		b[i][0]=y[i];
	}
	X=(A.t()*A).inv()*A.t()*b;
	cout<<"a0 and a1 was calculated as follow:\n"<<X;




	double sumx=0, sumr=(b-A*X).len();
	for(int i=0; i<x.size(); ++i)
		sumx+=x[i];
	double mx=sumx/x.size();
	sumx=0;
	for(int i=0; i<x.size(); ++i)
		sumx+=(mx-x[i])*(mx-x[i]);

	double t=(X[1][0]-at)/sqrt(sumr/((x.size()-2)*sumx));
	//cout<<"sumx="<<sumx<<" sumr="<<sumr<<'\n';
	cout<<"the T value is "<<t<<" with degrees of freedom "<<x.size()-2<<'\n';



        //question 2
	cout<<"\n\nquestion2:\n";

	x.clear();
	y.clear();
	data=storeDouble("sample2.dat");                                                                //read data form data file

	for(int i=0; i<data.size(); i+=2)                                                                //splitting 
	{
		x.push_back(data[i]);
		y.push_back(data[i+1]);
	}
	n=x.size();

	Matrix<double> J(n,4), Y(n,1), R(n,1);
	for(int i=0; i<n; ++i)
	{
		J[i][0]=1;
		J[i][1]=x[i];
		J[i][2]=x[i]*x[i];
		J[i][3]=x[i]*x[i]*x[i];
		Y[i][1]=y[i];
	}
	//cout<<"jt is: \n"<<J.t();
	Matrix<double> B(4,1);
	B[0][0]=B[1][0]=B[2][0]=B[3][0]=1;                //Guess 1
	double sse, ssep=0;

	for(int ms=0; ms<10000; ++ms)
	{
		R=Y-J*B;
		sse=R.len();
		if(abs(sse-ssep)<TOL)
			break;

		B=B+(J.t()*J).inv()*J.t()*R;

		ssep=sse;
	}
	cout<<"Guess 1: initial values of a0-a3 are all 1: when sse value reach "<<sse<<" the iteration stop, which gives the value of a0-a3 to be:\n"<<B;

	B[0][0]=B[1][0]=B[2][0]=B[3][0]=1000;           //Guess 2
	sse, ssep=0;

	for(int ms=0; ms<10000; ++ms)
	{
		R=Y-J*B;
		sse=R.len();
		if(abs(sse-ssep)<TOL)
			break;

		B=B+(J.t()*J).inv()*J.t()*R;

		ssep=sse;
	}
	cout<<"Guess 2: initial values of a0-a3 are all 1000: when sse value reach "<<sse<<" the iteration stop, which gives the value of a0-a3 to be:\n"<<B;

	return 0;
}


vector<double> storeDouble(const char* filename) //can read the file character by character and then converted into double format.
{   //this method read in the file char by char and then rebuild the double vector.
	char c;
	bool  desimal=false, neg=false;
	int accumulator;
	double temp=0;
	vector<double> result;
	fstream file(filename);
	while(file.get(c)!='\0')
	{
		if(c=='-')
			neg=true;
		else if(c==',' || c=='\n')              //only need to change here if the separator change
		{
			desimal=false;
			if(neg)
				result.push_back(temp*(-1));
			else
				result.push_back(temp);
			temp=0.0;
			neg=false;
		}
		else if(c=='.')
		{
			desimal=true;
			accumulator=-1;
		}
		else if(desimal==false)
			temp=temp*10+c-'0';
		else
		{
			temp+=(c-'0')*pow(10, accumulator);
			--accumulator;
		}
	}
	/*if(neg)
		result.push_back(temp*(-1));
	else
		result.push_back(temp);*/
	file.close();
	return result;
}

