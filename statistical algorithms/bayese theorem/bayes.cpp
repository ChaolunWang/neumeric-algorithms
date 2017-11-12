//--------------------------bayes.cpp--------------------------------
/*This codes implements the bayes reference for distribution parameters
  Programmed by: Chaolun Wang
	     at: 09/28/2016
*/

#include<vector>
#include<fstream>
#include<iostream>
#include<cmath>

const double ZA=1.96;                                                                            //the value of confidential coefficient is from online table for 95%
using namespace std;

vector<double> storeDouble(const char* filename);                                               //can read the file character by character and then converted into double format.
double postMu(double mu0, double theata0, double theata, double& prev_y, double y, int n);      //to caculate the posterior mean using the parameter estimation technique prev_y: previous mean, y: new data
double confidential(double theata0, double theata, int n);                                      //to calculate the confidential reagion using the prior variance and the data variance

int main()
{
	
	vector<double> data=storeDouble("norm.dat");                                                 //read data form data file
	double sum=0;
	for(int i=0; i<data.size(); ++i)
		sum+=data[i];
	cout<<"The size of data is: "<<data.size()<<" with mean value to be: "<<sum/data.size()<<'\n';   //output the information about data

	//for question 1	
	ofstream myfile;
	myfile.open ("question1.txt");
	double prevousy=0;

	for(int i=0; i<data.size(); ++i)
		myfile<<postMu(0, 0.1, sqrt(0.1), prevousy, data[i], i+1)<<' ';                              //calculate and output the posterior means for the growing size of data. Output file: question1.txt
	int k=1;
	while(confidential(0.1, sqrt(0.1), k)>0.05){                                                     //search the value of n for proper confidence interval for question1
		++k;
	}
	cout<<"Question 1: "<<k-1<<" samples needed before reaching +-0.05\n";
	myfile.close();
	
	//for question 2
	myfile.open ("question2.txt");
	prevousy=0;
	for(int i=0; i<data.size(); ++i)
		myfile<<postMu(0, 1, sqrt(0.1), prevousy, data[i], i+1)<<' ';                                //calculate and output the posterior means for the growing size of data. Output file: question1.txt
	myfile.close();

	k=1;
	while(confidential(1, sqrt(0.1), k)>0.05){                                                       //search the value of n for proper confidence interval for question2
		++k;
	}
	cout<<"Question 2: "<<k-1<<" samples needed before reaching +-0.05\n";
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
		else if(c==',')
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
	if(neg)
		result.push_back(temp*(-1));
	else
		result.push_back(temp);
	file.close();
	return result;
}

double postMu(double mu0, double theata0, double theata, double& prev_y, double y, int n)  //to caculate the posterior mean using the parameter estimation technique prev_y: previous mean, y: new data
{
	prev_y=(prev_y*(n-1.0)+y)/n;
	return (n/(theata*theata)/(1.0/(theata0*theata0)+n/(theata*theata)))*prev_y+(mu0/(theata0*theata0)/(1.0/(theata0*theata0)+n/(theata*theata)));   //the formula in lab report
}

double confidential(double theata0, double theata, int n)                                  //to calculate the confidential reagion using the prior variance and the data variance
{
	return ZA/sqrt(1.0/(theata0*theata0)+n/(theata*theata));                       //the formular in lab report
}
