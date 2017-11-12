//-------------------crs.cpp----------------------
/* This is the implementation of the controlled random search algorithm
   Replicated from the Price's paper from 1977. The validation of this 
   method was judjed by comparing with the BFGS mehtod, and also the 
   result from the paper
	 Programmed by: Chaolun Wang
		    at: 10/26/2016
*/


#include <iostream>
#include <cmath>
#include <vector>
#include <time.h>
#include <stdlib.h>
#include <algorithm>
#include <cstring>
#include <fstream>
#include <limits>

using namespace std;

const int N=50;       //cloud size
const int n=2;        //variable number
const int neginf=numeric_limits<int>::min();
const int inf=numeric_limits<int>::max();

double f(double x, double y);   //test function

vector<double> generatingp(const vector<vector<double> > cloud, int n, int N); //generate p value randomly

bool comparithm(const vector<double>& i, const vector<double>& j)       //comparithm used in the partial sort
{
  return i.back()>j.back();
}


int main()
{
	ofstream myfile1, myfile2;
	srand (time(NULL));
	double x, y;
	vector<vector<double> > cloud;
	myfile1.open ("functionvalue.txt");                                  //results stored into output files
	myfile2.open ("cloud.txt");

	for(int i=0; i<N; ++i)
	{
		vector<double> temp(n+1);
		x=(((double)rand())/RAND_MAX-0.5)*10;   //random number form -5 to 5
		y=(((double)rand())/RAND_MAX-0.5)*10;
		temp[0]=x;
		temp[1]=y;
		temp[2]=f(x,y);
		cloud.push_back(temp);
	}

	//recursively execute following step
	for(int i=0; i<4000; ++i)
	{
		if(i==0 || i==4000/3 ||i==(4000*2)/3)
		{
			for(int i=0; i<N; ++i)
			{
				for(int j=0; j<n+1; ++j)
				{
					myfile2<<cloud[i][j]<<' ';
				}
				myfile2<<'\n';
			}
		}
		
		partial_sort (cloud.begin(), cloud.begin()+1, cloud.end(),comparithm);    // a good way to find and move the max point to the index 0, since data was choosen randomly, the changing of order doesn't metter.
		do
		{
			vector<double> temp=generatingp(cloud, n, N);
			if(cloud[0][2]>=temp[2])
			{
				cloud[0]=temp;
				break;
			}
		}
		while(true);
		double min=inf, max=neginf, average=0, value;
		for(int j=0; j<N; ++j)
		{
			value=cloud[j][2];
			if(value<min)
				min=value;
			if(value>max)
				max=value;
			average+=value;
		}
		myfile1<<min<<' '<<max<<' '<<average/N<<'\n';

	}

	for(int i=0; i<N; ++i)
	{
		for(int j=0; j<n+1; ++j)
		{
			myfile2<<cloud[i][j]<<' ';
		}
		myfile2<<'\n';
	}
	myfile1.close();
	myfile2.close();

}



double f(double x, double y)
{
	return 100*pow(y-x*x,2)+pow(6.4*pow(y-0.5,2)-x-0.6,2);
}


vector<double> generatingp(const vector<vector<double> > cloud, int n, int N)
{
	int index, count=0;
	int a[N], k[n+1];
	vector<double> result(3);
	do{
		count=0;
		memset(a, 0, sizeof(a));
	
		while(count<n+1)
		{
			index=rand()%N;
			if(a[index]==0)
			{
				k[count]=index;
				a[index]==1;
				count++;
			}	
		}

		double gx, gy, px, py;
		gx=gy=0;
		for(int i=0; i<n; ++i)
		{
			gx+=cloud[k[i]][0];
			gy+=cloud[k[i]][1];
		}
		gx/=n;
		gy/=n;
		px=gx*2-cloud[k[n]][0];
		py=gy*2-cloud[k[n]][1];
		if(px<5 && px>-5 && py<5 && py>-5)
		{
			result[0]=px;
			result[1]=py;
			result[2]=f(px,py);
			break;
		}
	}while(true);
	return result;
}






