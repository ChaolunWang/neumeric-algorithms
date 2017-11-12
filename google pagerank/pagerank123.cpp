//---------------------------------pagerank123.cpp---------------------------
/*
    This is the cpp file of the Google Page rank program(for question 1, 2, 3).
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

using namespace std;

const double D=0.85;
const double TOL=0.00000001;

void randomInitializer(double *A, int size);                                       //randomly initialized the matrix
void gaussSeidel(int n, double * A, double * b, double * x, double tolerance);     //using gaussSeidel method to solve system in method 1
void powerMethod(int n, double * A, double * x, double tolerance);                 //using powermethod to find the eigenvector in method 2
double lOneNorm(double * x, int size);                                             //calculate the l1 norm


int main()
{
	//Question 1	
	int size;
	bool valid;
	cout<<"Question1:\nPlease input the size of matrix: \n";
	cin>>size;
	double *matrixA=new double[size*size];
	
	srand(time(NULL));
	                                                            //initializing the matrix
    do
    {
		cout<<"Initializing the matrix A;\n";
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
				cout<<"One or more coulmn have no non-zero element, matrix discarded;\n";
				break;
			}
	    }
	}while(!valid);                                                    //if the matrix is invalid, re-generate the matrix
	
	cout<<"Matrix generation successful!\n";
	cout<<"The matrix is:\n";                                         //print out the generated matrix
    for(int i=0; i<size; ++i)
	{
	    for(int j=0; j<size; ++j)
	    {
			cout<<matrixA[i*size+j]<<' ';
		}
		cout<<'\n';
	}

	
	//Question 2
	cout<<"\nQuestion2:\n";
	
	double *matrixM=new double[size*size];
	double *coefficientM=new double[size*size]; 
	double *rhs=new double[size];
	double *result=new double[size];

	
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
	cout<<"The matrix M is:\n";	                            //print out matrix M
    ios_base::fmtflags fstate = cout.flags();
    int oprecision = cout.precision();
    cout.setf(ios::fixed);
	cout.precision(3); 
    for(int i=0; i<size; ++i)
	{
	    for(int j=0; j<size; ++j)
	    {
			cout<<matrixM[i*size+j]<<' ';
		}
		cout<<'\n';
	}
	cout.flags(fstate);
	cout.precision(oprecision);	
	
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
		result[i]=1;                   //initial guess set the value to be 1
	}		
	
	gaussSeidel(size, coefficientM, rhs, result, TOL);            //call gauss seidel method


    cout<<"The page rank(calculated by method 1) result is:\n";
    for(int i=0; i<size; ++i)                                     //print out result
	{
		cout<<result[i]<<'\n';                   
	}	



	
	//Question 3
	cout<<"\nQuestion3:\n";

    for(int i=0; i<size; ++i)
	{
		result[i]=1;                   //initial guess set the value to be 1
	}	
	
	//calculate M^	
    for(int i=0; i<size; ++i)
	{
	    for(int j=0; j<size; ++j)
	    {
	        coefficientM[i*size+j]=matrixM[i*size+j]*D+(1-D)/size;
		}
	}	
	
    powerMethod(size, coefficientM, result, TOL);                      //call power method
    
    cout<<"The page rank(calculated by method 2) result is:\n";
    for(int i=0; i<size; ++i)                     //print out result
	{
		cout<<result[i]<<'\n';                   
	}
			     
	delete [] matrixA;
	delete [] matrixM;
	delete [] coefficientM;
	delete [] rhs;
	delete [] result;

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
    }while(lOneNorm(change,n)>=tolerance);                        //no need to normalized since all entries in x added up to one(l-1 norm of x is 1).
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
	} while(lOneNorm(change,n)>=tolerance);                 //no need to normalized since all entries in x added up to one(l-1 norm of x is 1).

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


