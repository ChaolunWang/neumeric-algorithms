//--------------------------power.cpp--------------------------------
/*This codes tests the statistical power using centain number of sampling
  and compaire the experimental result with theoritical result
  Programmed by: Chaolun Wang
	     at: 10/5/2016
output exampe:
For n=100, the acceptance(1) and rejectance(0) are stored in an array for each of the 100 samples: 
1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 0 1 0 1 0 1 1 0 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 0 0 1 0 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 0 1 1 1 0 1 0 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 
for each n=50 100 250 500 1000 calculate the poportion of rejection:
0.11 0.17 0.39 0.61 0.86 

*/

#include<vector>
#include<fstream>
#include<iostream>
#include<cmath>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */


using namespace std;

vector<bool> storeBool(const char* filename);                                                   //can read the file character by character and then converted into bool format.
double randomSamp(const vector<bool>& data, int n );                                             //get random sample
bool accept(double p, int n);                                                                    //determine whether or not the sample is rejected

int main()
{
	srand (time(NULL));
	vector<bool> data=storeBool("survey1.dat");                                                 //read data form data file
	//cout<<data.size()<<' '<<randomSamp(data,10000)<<'\n';                                     //test read
	cout<<"For n=100, the acceptance(1) and rejectance(0) are stored in an array for each of the 100 samples: \n";	
	bool array[100];
	for (int i=0; i<100; ++i)
		array[i]=accept(randomSamp(data, 100), 100);                                        //calculate the acceptance and rejection and store in array

	for (int i=0; i<100; ++i)
		cout<<array[i]<<' ';
	cout<<'\n';

	cout<<"for each n=50 100 250 500 1000 calculate the poportion of rejection:\n";
	int number[5]={50,100,250,500,1000};
	double result[5]={0,0,0,0,0};
	for(int i=0; i<5; ++i)
	{
		for(int j=0; j<100; ++j)
			result[i]+=(double)accept(randomSamp(data, number[i]), number[i]);    
		result[i]=1.0-result[i]/100;
	}

	for (int i=0; i<5; ++i)
		cout<<result[i]<<' ';
	cout<<'\n';
}






vector<bool> storeBool(const char* filename) //can read the file character by character and then converted into boolean format.
{   //this method read in the file char by char and then rebuild the bool vector.
	char c;
	vector<bool> result;
	fstream file(filename);
	while(file.get(c)!='\0')
	{
		if(c=='0')
			result.push_back(false);
		else if(c=='1')
			result.push_back(true);
	}
	file.close();
	return result;
}

double randomSamp(const vector<bool>& data, int n )
{
	int size=data.size();
	double result=0;
	for (int i=0; i<n; ++i)
	{
		result+=(double)data[rand()%size];
	}
	return result/n;
}

bool accept(double p, int n)                                                                    //determine whether or not the sample is rejected
{
	if (p>0.5+0.98/sqrt(n)||p<0.5-0.98/sqrt(n))
		return false;
	return true;
}
