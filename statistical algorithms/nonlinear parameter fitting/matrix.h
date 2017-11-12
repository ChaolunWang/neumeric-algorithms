//---------------------matrix.h-------------------------------
/* This is a template class which is a blueprint of 2D matrix.
   right now only support overloaded << ; * ; [][] ;  = and other matrix operations
   Maybe it is also useful for future homeworks if expanded.
   Programmed by: Chaolun Wang  modefied by: Chaolun Wang
              at: 09/11/2015             at: 10/10/2016
*/

#ifndef _MATRIX_H
#define _MATRIX_H

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

extern "C" int dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);
extern "C" int dpotri_(char *uplo, int *n, double *a, int *lda, int *info);

template<class T>
class Matrix
{
//utility function
friend int check(const char* fileName);                                     //return the number of enturies in a matrix in txt file

template<class U>
friend ostream& operator<<(ostream& os, const Matrix<U> &matrix);           //save matrix in ostream 

template<class U>
friend Matrix<U> operator*(double value, const Matrix<U> &m);                 //matrix multiplication 

public:
	Matrix(int m, int n);                                                   //constructor
	Matrix(const Matrix<T> &m);                                             //copy constructor
	~Matrix();                                                              //destructor
	Matrix<T>& operator=(const Matrix<T> &m);                               //assignment operator
	Matrix<T> operator*(const Matrix<T> &m);                                //matrix multiplication
	Matrix<T> operator*(double value);                                      //matrix multiplication
	Matrix<T> inv();                                			//matrix inverse
	Matrix<T> t();                                			        //matrix transpose
	Matrix<T> operator+(const Matrix<T> &m);                                //matrix +
	Matrix<T> operator-(const Matrix<T> &m);                                //matrix -
	double len();                                			        //matrix eucledian length

	//sub class Row
	class Row                                                               //object Row, can overload second[]
	{
	public:
		Row(T *mRow):mRow(mRow){}
		T& operator[](int n){return mRow[n];}                                //second[] operator for write
		const T& operator[](int n) const{return mRow[n];}                    //second[] operator for read
	private:
		T *mRow;
	};
	
	Row operator[](int m){return Row(&mMatrix[m*mColumn]);}                           //first[] operator for write
	const Row& operator[](int m) const{return Row(mMatrix[m]);}              //first[] operator for read
	bool saveMatrix(const char *fileName, int type=0, int mode=0) const;     //save matrix into: txt file(type=0) , binary file(type=1) mode:append(0) rewrite(1)
	bool readMatrix(const char *fileName, int index=0);                      //read the nth matrix from binary file
	double average() const;                                                  //return the average of all the enturies

private:
    	T *mMatrix;                                                              //inside matrix
	int mRow;                                                                 //row number
	int mColumn;                                                              //column number
 
};


#endif

int check(const char* fileName)                                  
//return the number of enturies in a matrix in txt file
{
	ifstream input;
	int word=0, counter=0;
	char ch;
	input.clear();
	input.open(fileName);
	if(!input)
		return -1;
	while(input.get(ch))
	{
		if(ch=='\r' || ch=='\n' || ch==' ')
		{
			if(word)
			{
			    ++counter;
			    word=0;
			}
		}
		else
		    ++word;
	}
	if(word)
	    ++counter;
	input.close();
	return counter;
}

template<class T>  
ostream& operator<<(ostream& os, const Matrix<T> &matrix)
//save matrix in ostream
{
	ios_base::fmtflags fstate = os.flags();
	int oprecision = os.precision();
	os.setf(ios::fixed);
	os.precision(7);
	for(int i=0; i<matrix.mRow; ++i)
	{
		for(int j=0; j<matrix.mColumn; ++j)
			os<<setw(15)<<left<<matrix.mMatrix[i*matrix.mColumn+j]<<' ';
		os<<'\n';
	}
	os<<'\n';
	os.flags(fstate);
	os.precision(oprecision);
	return os;
}

template<class T>
Matrix<T> operator*(double value, const Matrix<T> &m)
{
	Matrix<T> result(m.mRow, m.mColumn);
	for(int i=0; i<m.mRow;++i)
		for(int j=0; j<m.mColumn; ++j)
		{
			result.mMatrix[i*result.mColumn+j]=value*m.mMatrix[i*m.mColumn+j];
		}
	return result;
}
    
template<class T>
Matrix<T>::Matrix(int m, int n)   
//constructor
{
	mRow=m;
	mColumn=n;
	mMatrix=new T[m*n];
}

template<class T>
Matrix<T>::~Matrix()                                                              
//destructor
{
	delete [] mMatrix;
}

template<class T>
Matrix<T>::Matrix(const Matrix<T> &m)
//copy constructor
{
	mRow=m.mRow;
	mColumn=m.mColumn;
	mMatrix=new T[mRow*mColumn];
	for(int i=0; i<mRow*mColumn; ++i)
		mMatrix[i]=m.mMatrix[i];
}
	
template<class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T> &m)                               
//assignment operator
{
	if(this!=&m)
	{
		delete [] mMatrix;
		mRow=m.mRow;
		mColumn=m.mColumn;
		mMatrix=new T[mRow*mColumn];
		for(int i=0; i<mRow; ++i)
			for(int j=0; j<mColumn; ++j)
				mMatrix[i*mColumn+j]=m.mMatrix[i*mColumn+j];
	}
	return *this;
}

template<class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &m)                                  
//matrix multiplication
{
	Matrix<T> result(mRow, m.mColumn);
	for(int i=0; i<mRow;++i)
		for(int j=0; j<m.mColumn; ++j)
		{
			result.mMatrix[i*result.mColumn+j]=0.0;
			for(int n=0; n<mColumn; ++n)
			    result.mMatrix[i*result.mColumn+j]+=mMatrix[i*mColumn+n]*m.mMatrix[n*m.mColumn+j];
		}
	return result;
}

template<class T>
double Matrix<T>::len()
{
	double result=0;
	for(int i=0; i<mRow;++i)
		for(int j=0; j<mColumn; ++j)
			    result+=mMatrix[i*mColumn+j]*mMatrix[i*mColumn+j];
	return result;
}


template<class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &m)                                  
//matrix multiplication
{
	Matrix<T> result(mRow, mColumn);
	for(int i=0; i<mRow;++i)
		for(int j=0; j<mColumn; ++j)
		{
			    result.mMatrix[i*result.mColumn+j]=mMatrix[i*mColumn+j]+m.mMatrix[i*m.mColumn+j];
		}
	return result;
}

template<class T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &m)                                  
//matrix multiplication
{
	Matrix<T> result(mRow, mColumn);
	for(int i=0; i<mRow;++i)
		for(int j=0; j<m.mColumn; ++j)
		{
			    result.mMatrix[i*result.mColumn+j]=mMatrix[i*mColumn+j]-m.mMatrix[i*m.mColumn+j];
		}
	return result;
}

template<class T>
Matrix<T> Matrix<T>::operator*(double value)                                  
//matrix multiplication
{
	Matrix<T> result(mRow, mColumn);
	for(int i=0; i<mRow;++i)
		for(int j=0; j<mColumn; ++j)
		{
			result.mMatrix[i*result.mColumn+j]=value*mMatrix[i*mColumn+j];

		}
	return result;
}
		 
template<class T>
Matrix<T> Matrix<T>::inv()
{
  	int n = mRow;
	int info;
	char uplo = 'u';
	Matrix<T> result(mRow, mColumn);
	for(int i=0; i<mRow*mColumn; ++i)
		result.mMatrix[i]=mMatrix[i];

	// compute Cholesky factorization
	dpotrf_(&uplo,&n,&result.mMatrix[0],&n,&info);
	if (info != 0) cout << "Error in dpotrf_(): Flag is " << info << endl;
  
  	// use Cholesky factorization to compute the inverse
  	dpotri_(&uplo,&n,&result.mMatrix[0],&n,&info);
	//fill in lower triangular
  	if (info != 0) cout << "Error in dpotri_(): Flag is " << info << endl;
	for(int i=0; i<mRow; ++i)
		for(int j=0; j<mColumn; ++j)
			result.mMatrix[i*mColumn+j]=result.mMatrix[j*mColumn+i];
	return result;
}
 
template<class T>
Matrix<T> Matrix<T>::t()
{
	Matrix<T> result(mColumn, mRow);
	for(int i=0; i<mColumn;++i)
		for(int j=0; j<mRow; ++j)
		{
			result.mMatrix[i*result.mColumn+j]=mMatrix[i+j*mColumn];

		}
	return result;
}

template<class T>   
bool Matrix<T>::saveMatrix(const char *fileName, int type, int mode) const                       
//save matrix into: txt file(type=0) , binary file(type=1) mode:append(0) rewrite(1)
{
	if(type)
	{	
		FILE *fpt;
		if(mode)
		    fpt=fopen(fileName, "wb+");
		else
		    fpt=fopen(fileName, "ab+");
		if(!fpt)
			return false;

                fwrite(mMatrix, sizeof(T), mColumn*mRow, fpt);
        fclose(fpt);
	}
	else
	{
	    ofstream output;
		output.clear();
		if(mode)
		    output.open(fileName);
		else
		    output.open(fileName, std::ios_base::app);
		if(!output)
		    return false;
		output<<*this;
		output.close();
		}
	return true;
}

template<class T>
bool Matrix<T>::readMatrix(const char *fileName, int index)                       
//read a certain matrix from binary file
{
	FILE *fpt;
	fpt=fopen(fileName, "r");
	if(!fpt)
		return false;
	T *buffer=new T[mRow*mColumn*(index+1)];
	fread(buffer, sizeof(T), mRow*mColumn*(index+1), fpt);
	for(int i=0; i<mRow; ++i)
		for(int j=0; j<mColumn; ++j)
		    mMatrix[i][j]=buffer[mRow*mColumn*index+i*mColumn+j];
	delete [] buffer;
	fclose(fpt);
	return true;
}

template<class T>
double Matrix<T>::average() const
//return the average of all the enturies
{
	double result=0;
	for(int i=0; i<mRow; ++i)
		for(int j=0; j<mColumn; ++j)
		    result+=mMatrix[i][j];
	result/=mRow*mColumn;
	return result;
}

//end of matrix.h
