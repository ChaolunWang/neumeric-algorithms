//---------------------------------pagerank5.cpp---------------------------
/*
    This is the cpp file of the Google Page rank program(for question 5).
    Power method methods was implemented using different termination criteria. 
    The aggressive criteria I proposed to find the order of top ten site is like this:
    If the rank 10 haven't change for more than two iteration, termination will occur.
    Validation of this method was checked by running code multiple times for different
    size of n.
    Programmed by: Chaolun Wang
               at: 02/03/2016
*/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <omp.h>

using namespace std;

const double D=0.85;
const double TOL=0.0000001;

void randomInitializer(double *A, int size);                                        //randomly initialized the matrix
void powerMethod(int n, double * A, double * x, double tolerance);                  //using power method to find the eigenvector in method 2
void aggressivePowerMethod(int n, double * A, double * x, int stableTime);          // power method with the aggressive criteria
double lOneNorm(double * x, int size);                                               //calculate the l1 norm
bool findTop(int * topList, double * x, int size);           //find the top ten, store in topList, if the value in topList is different than original value, return false, else return true;


int main()
{
	//Question 5
	int input[8]={50, 100, 200, 500, 1000, 2000, 3000, 5000};
    double wtime;
    
    cout<<"time spent:\n";
    cout.setf(ios::fixed);
	cout.precision(10); 
    cout<<"n=       Tol=10^(-7)     Tol=10^(-6)     Tol=10^(-5)     Tol=10^(-4)   Aggressive Power Method     \n";
	for(int n=0; n<8; ++n)
	{
        int size=input[n];
	    bool valid;
	    double *matrixA=new double[size*size];
	    double *matrixM=new double[size*size];
	    double *coefficientM=new double[size*size]; 
	    double *rhs=new double[size];
	    double *result=new double[size];	
	    srand(time(NULL));
	    //initializing the matrix
        do
        {   
		    valid=true;
	        randomInitializer(matrixA, size);
	        for(int j=0; j<size; ++j)
	        {
			    int sum=0;
	            for(int i=0; i<size; ++i)
	            {
				    sum+=matrixA[i*size+j];
			    }
			    if(sum==0)
			    {
				    valid=false;
				    break;
			    }
	        }
	    }while(!valid);
        
        for(int k=0; k<4; ++k)
        {
            //Power method at different tolerance
            wtime = omp_get_wtime( ); 
            //calculate M
            for(int j=0; j<size; ++j)
	        {
		        int sum=0;
	            for(int i=0; i<size; ++i)
	            {
			        sum+=matrixA[i*size+j];
		        }
	            for(int i=0; i<size; ++i)
	            {
			        matrixM[i*size+j]=matrixA[i*size+j]/sum;
		        }	
		        result[j]=1.0/size;                   //initial guess	
	        }
	
	        //calculate M^	
            for(int i=0; i<size; ++i)
	        {
	            for(int j=0; j<size; ++j)
	            {
	                coefficientM[i*size+j]=matrixM[i*size+j]*D+(1-D)/size;
		        }
	        }	
	
            powerMethod(size, coefficientM, result, TOL*pow(10, k));   // call aggressive power method

            wtime=omp_get_wtime ( )-wtime;
            if(k==0)
                cout<<setw(5)<<left<<size;
            cout<<"    "<<setw(12)<<wtime;
	    }
	    
	    //Aggressive power method
        wtime = omp_get_wtime ( );
        //calculate M
        for(int j=0; j<size; ++j)
	    {
		    int sum=0;
	        for(int i=0; i<size; ++i)
	        {
			    sum+=matrixA[i*size+j];
		    }
	        for(int i=0; i<size; ++i)
	        {
			    matrixM[i*size+j]=matrixA[i*size+j]/sum;
		    }	
		    result[j]=1.0/size;                   //initial guess set the value to be 1/size	
	    }
	
	    //calculate M^	
        for(int i=0; i<size; ++i)
	    {
	        for(int j=0; j<size; ++j)
	        {
	            coefficientM[i*size+j]=matrixM[i*size+j]*D+(1-D)/size;
		    }
	    }	
	
        aggressivePowerMethod(size, coefficientM, result, 1);   // call aggressive power method
       
        wtime=omp_get_wtime ( )-wtime;
        cout<<"    "<<setw(12)<<wtime<<'\n';
        
        delete [] matrixA;
	    delete [] matrixM;
	    delete [] coefficientM;
	    delete [] rhs;
	    delete [] result;
   }

  //validation of aggressive power method
  // the code was runned again to see if the method can get the same top ten rank for each value of n  
    cout<<"Top ten ID:\n";

	for(int n=0; n<8; ++n)
	{
        int size=input[n], topTen[10];
	    bool valid;
	    double *matrixA=new double[size*size];
	    double *matrixM=new double[size*size];
	    double *coefficientM=new double[size*size]; 
	    double *rhs=new double[size];
	    double *result=new double[size];
	    srand(time(NULL));	
	    //initializing the matrix
        do
        {   
		    valid=true;
	        randomInitializer(matrixA, size);
	        for(int j=0; j<size; ++j)
	        {
			    int sum=0;
	            for(int i=0; i<size; ++i)
	            {
				    sum+=matrixA[i*size+j];
			    }
			    if(sum==0)
			    {
				    valid=false;
				    break;
			    }
	        }
	    }while(!valid);
        
        for(int k=0; k<4; ++k)
        {
            //Power method at different tolerance
            //calculate M
            for(int j=0; j<size; ++j)
	        {
		        int sum=0;
	            for(int i=0; i<size; ++i)
	            {
			        sum+=matrixA[i*size+j];
		        }
	            for(int i=0; i<size; ++i)
	            {
			        matrixM[i*size+j]=matrixA[i*size+j]/sum;
		        }	
		        result[j]=1.0/size;                   //initial guess set the value to be 1/size	
	        }
	
	        //calculate M^	
            for(int i=0; i<size; ++i)
	        {
	            for(int j=0; j<size; ++j)
	            {
	                coefficientM[i*size+j]=matrixM[i*size+j]*D+(1-D)/size;
		        }
	        }	
	
            powerMethod(size, coefficientM, result, TOL*pow(10, k));  // call power method

            findTop(topTen, result, size);                            //find top ten, print out result
            if(k==0)
                cout<<"when n= "<<size<<":\n"<<"ID of top ten webpage was printed out:\n";
            cout<<"Tol=10^(-"<<7-k<<"):       ";
            for(int i=0; i<10; ++i)
                cout<<topTen[i]<<" ";
            cout<<'\n';
	    }
	    
	    //Aggressive power method
        wtime = omp_get_wtime ( );
        //calculate M
        for(int j=0; j<size; ++j)
	    {
		    int sum=0;
	        for(int i=0; i<size; ++i)
	        {
			    sum+=matrixA[i*size+j];
		    }
	        for(int i=0; i<size; ++i)
	        {
			    matrixM[i*size+j]=matrixA[i*size+j]/sum;
		    }	
		    result[j]=1.0/size;                   //initial guess set the value to be 1/size	
	    }
	
	    //calculate M^	
        for(int i=0; i<size; ++i)
	    {
	        for(int j=0; j<size; ++j)
	        {
	            coefficientM[i*size+j]=matrixM[i*size+j]*D+(1-D)/size;
		    }
	    }	
	
        aggressivePowerMethod(size, coefficientM, result, 1);  // call aggressive power method
       
        findTop(topTen, result, size);                          //find top ten, print out result
        cout<<"Aggressive method: ";
        for(int i=0; i<10; ++i)
            cout<<topTen[i]<<" ";
        cout<<'\n';
        
        delete [] matrixA;
	    delete [] matrixM;
	    delete [] coefficientM;
	    delete [] rhs;
	    delete [] result;
   }
}



void randomInitializer(double *A, int size)                        
{
	for(int i=0; i<size; ++i)
	    for(int j=0; j<size; ++j)
	    {
			if(i!=j)
			{	            
	            A[i*size+j]=rand()%2;                                  //initialized by random number 
			}
			else
			    A[i*size+j]=0;
        }
}


void powerMethod(int n, double * A, double * x, double tolerance)
{
	double *tempResult=new double[n];
	double *change=new double[n];	
	do
	{
		for(int i=0; i<n; ++i)
		    change[i]=x[i];
		    
        for(int i=0; i<n; ++i)
	    {
			tempResult[i]=0;
	        for(int j=0; j<n; ++j)
	        {
	            tempResult[i]+=A[i*n+j]*x[j];
		    }
	    }
		double delta=lOneNorm(tempResult, n);
		for(int i=0; i<n; ++i)
		{
			x[i]=tempResult[i]/delta;
		}
		
	   	for(int i=0; i<n; ++i)
		    change[i]=x[i]-change[i];			
	} while(lOneNorm(change,n)>=tolerance);

	delete [] tempResult;
	delete [] change;
}

void aggressivePowerMethod(int n, double * A, double * x, int stableTime)
{
	double *tempResult=new double[n];
	int topList[10], counter=0, number=0;	
	do
	{	    
        for(int i=0; i<n; ++i)
	    {
			tempResult[i]=0;
	        for(int j=0; j<n; ++j)
	        {
	            tempResult[i]+=A[i*n+j]*x[j];
		    }
	    }
		double delta=lOneNorm(tempResult, n);
		for(int i=0; i<n; ++i)
		{
			x[i]=tempResult[i]/delta;
		}
       
		if(findTop(topList, x, n))                //findTop will return true if the value in top list haven't change         
		    counter++;                            //counter will then increase
		else
		    counter=0;
			
	} while(counter<stableTime);                            

	delete [] tempResult;
}

double lOneNorm(double * x, int size)
{
	if(x==NULL)
	    return 0;
	double result=0;
	for(int i=0; i<size; ++i)
	{
		result+=abs(x[i]);
	}
	return result;
}

bool findTop(int * topList, double * x, int size)
{
	bool find, result=true;
	for(int i=0; i<10; ++i)
	{
		int maxId=0;
		do                                                  //initialize the maxId
		{
			find=false;
			for(int n=0; n<i; ++n)
	            if(topList[n]==maxId)
	            {
					++maxId;
					find=true;
					break;	
				}
		}while(find);                                        //this loop can make sure that the initial value of maxId is not inside the previous findTop
		
		for(int k=0; k<size; ++k)                            //find the ID of i th largest value
		{
			if(x[maxId]<x[k])
			{
				find=false;
			    for(int n=0; n<i; ++n)                        // if already in top list, ignore that
			        if(topList[n]==k)
			            find=true;
			    if(!find)                                     //else replace the maxId
			        maxId=k;
			}
		}
		if(topList[i]!=maxId)
		{
			result=false;
			topList[i]=maxId;
		}
	}
	return result;
}
