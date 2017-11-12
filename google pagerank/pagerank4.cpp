//---------------------------------pagerank4.cpp---------------------------
/*
    This is the cpp file of the Google Page rank program(for question 4).
    It implements the page rank algorithm using two different methods. For 
    solving the linear system in Method 1, Gauss Seidel method was used since 
    matrix (I-dM) is diagonal dominant matrix. 
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
const double TOL=0.00000001;

void randomInitializer(double *A, int size);                                            //randomly initialized the matrix
void gaussSeidel(int n, double * A, double * b, double * x, double tolerance);           //using gaussSeidel method to solve system in method 1
void powerMethod(int n, double * A, double * x, double tolerance);                       //using powermethod to find the eigenvector in method 2
double lOneNorm(double * x, int size);                                                   //calculate the l1 norm



int main()
{
	//Question 4
	int input[8]={5, 10, 50, 100, 500, 1000, 2000, 5000};
    double wtime;
    
    cout.setf(ios::fixed);
	cout.precision(10); 
    cout<<"n=       Method 1 time   Method 2 time\n";
	for(int n; n<8; ++n)
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
        
        //Methord1
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
	    }
	
	    //calculate I-DM
        for(int i=0; i<size; ++i)
	    {
	        for(int j=0; j<size; ++j)
	        {
			    if(i!=j)
			        coefficientM[i*size+j]=matrixM[i*size+j]*D*(-1);
			    else
			        coefficientM[i*size+j]=1;
		    }
	    }
	
	    //calculate right hand side;
        for(int i=0; i<size; ++i)
	    {
		    rhs[i]=(1-D)/size;
		    result[i]=1.0/size;                              //initial guess each entry will be 1/n
	    }		
	
	    gaussSeidel(size, coefficientM, rhs, result, TOL);    //call gaussSeidel method
        wtime=omp_get_wtime ( )-wtime;                        //print out result
        cout<<setw(5)<<left<<size<<"    "<<setw(12)<<wtime;
        
	    //Method 2
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
		    result[j]=1.0/size;                   //initial guess	each entry will be 1/n
	    }
	
	    //calculate M^	
        for(int i=0; i<size; ++i)
	    {
	        for(int j=0; j<size; ++j)
	        {
	            coefficientM[i*size+j]=matrixM[i*size+j]*D+(1-D)/size;
		    }
	    }	
	
        powerMethod(size, coefficientM, result, TOL);                   //call powerMethod
       
        wtime=omp_get_wtime ( )-wtime;                                  //print out result
        cout<<"    "<<setw(12)<<wtime<<'\n';
        
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
	            A[i*size+j]=rand()%2;
			}
			else
			    A[i*size+j]=0;
        }
}

void gaussSeidel(int n, double * A, double * b, double * x, double tolerance)                                    //using GaussSeidel's methord to solve Ax=b for n*n matrix, parameter: n:size A:matrix b:right hand side x: unknown                                                    
{
	
	double *change=new double[n];
	do
	{
		for(int i=0; i<n; ++i)
		    change[i]=x[i];
	    for(int i=0; i<n; i++)
	    {
		    double a_ijxj=0.0;
		    for(int j=0; j<n; j++)
		    {
			    if(j!=i)
				a_ijxj+=A[i*n+j]*x[j];
	        }
		    x[i]=(b[i]-a_ijxj)/A[i*(n+1)];
	    }
		for(int i=0; i<n; ++i)
		    change[i]=x[i]-change[i];
    }while(lOneNorm(change,n)>=tolerance);                          //no need to normalized since all entries in x added up to one(l-1 norm of x is 1).
   delete [] change;
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
	} while(lOneNorm(change,n)>=tolerance);                         //no need to normalized since all entries in x added up to one(l-1 norm of x is 1).

	delete [] tempResult;
	delete [] change;
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
