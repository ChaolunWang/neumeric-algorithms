//-----------------bfgs.cpp---------------------------
/*this is the implementation of the bfgs mehtod*/

#include <iostream>
#include <cmath>

using namespace std;
const double E=0.000001;
const double A=1;
struct xvector{
	double x1;
	double x2;
};

struct bmatrix{
	double a1;
	double a2;
	double b1;
	double b2;
};

double f(xvector x);
xvector gradf(xvector x);
bmatrix hessian(xvector x);
bmatrix neginverse(bmatrix b);
xvector operator+(xvector x1, xvector x2);
xvector operator*(double alf, xvector x2);

bmatrix operator+(bmatrix x1, bmatrix x2);
bmatrix operator*(double alf, bmatrix x2);
bmatrix operator*(bmatrix x1, bmatrix x2);

bmatrix operator*(xvector x1, xvector x2);
xvector operator*(bmatrix x1, xvector x2);

double times(xvector x1, xvector x2);

int main()
{


	
	double tol, alf, r=0.9, c=0.9;
	xvector xk, deltaf, p, s, y;
	int count;


//steepest Descent
	cout<<"steepest descent method:\n";
	xk.x1=-3;
	xk.x2=-4;
	count=0;
	do
	{
		deltaf=gradf(xk);
		p=(-1)*deltaf;
		alf=A;
		while(f(xk+alf*p)>(f(xk)+c*alf*times(deltaf, p)))
			alf=r*alf;
		//cout<<"alf: "<<alf<<'\n';
		xk=xk+alf*p;	
		tol=sqrt(times(gradf(xk), gradf(xk)))/(1.0+abs(f(xk)));
		count++;
	}while(tol>=E /*&& count<100*/);
	cout<<"x1= "<<xk.x1<<" x2= "<<xk.x2<<"iteration= "<<count<<'\n';
//newton's mehtod

	cout<<"newton's method:\n";
	xk.x1=-3;
	xk.x2=-4;
	count=0;
	do
	{
		deltaf=gradf(xk);
		p=neginverse(hessian(xk))*deltaf;
		alf=A;
		while(f(xk+alf*p)>(f(xk)+c*alf*times(deltaf, p)))
			alf=r*alf;
		//cout<<"alf: "<<alf<<'\n';
		xk=xk+alf*p;	
		tol=sqrt(times(gradf(xk), gradf(xk)))/(1.0+abs(f(xk)));
	}while(tol>=E /*&& count<100*/);
	cout<<"x1= "<<xk.x1<<" x2= "<<xk.x2<<"iteration= "<<count<<'\n';

//bfgs
	cout<<"bfgs method:\n";
	xk.x1=-3;
	xk.x2=-4;
	bmatrix bk1, bk=hessian(xk);
	count=0;
	do
	{
		deltaf=gradf(xk);
		p=neginverse(bk)*deltaf;
		alf=A;
		while(f(xk+alf*p)>(f(xk)+c*alf*times(deltaf, p)))
			alf=r*alf;
		s=alf*p;
		xk=xk+s;
		y=gradf(xk)+(-1)*deltaf;
		bk1=bk+(1.0/times(y,s))*(y*y)+(-1.0/times(s, bk*s))*(bk*s*s*bk);
		bk=bk1;
		tol=sqrt(times(gradf(xk), gradf(xk)))/(1.0+abs(f(xk)));
		++count;
	}while(tol>=E);
	cout<<"x1= "<<xk.x1<<" x2= "<<xk.x2<<"iteration= "<<count<<'\n';
}

double f(xvector x)
{
	return 100*pow(x.x2-pow(x.x1,2), 2)+pow(1-x.x1,2);
}
xvector gradf(xvector x)
{
	xvector result;
	result.x1=-400*(x.x2-pow(x.x1,2))*x.x1-2*(1-x.x1);
	result.x2=200*(x.x2-pow(x.x1,2));
	return result;
}
bmatrix hessian(xvector x)
{
	bmatrix result;
	result.a1=-400*(x.x2-3*pow(x.x1,2))+2;
	result.a2=result.b1=-400*x.x1;
	result.b2=200;
	return result;
}

bmatrix neginverse(bmatrix b)
{
	bmatrix result;
	double c=-1.0/(b.a1*b.b2-b.a2*b.b1);
	result.a1=c*b.b2;
	result.a2=-1*c*b.a2;
	result.b1=-1*c*b.b1;
	result.b2=c*b.a1;
	return result;
}

xvector operator+(xvector x1, xvector x2)
{
	xvector result;
	result.x1=x1.x1+x2.x1;
	result.x2=x1.x2+ x2.x2;
	return result;
}
xvector operator*(double alf, xvector x2)
{
	xvector result;
	result.x1=x2.x1*alf;
	result.x2=x2.x2*alf;
	return result;
}
bmatrix operator+(bmatrix x1, bmatrix x2)
{
	bmatrix result;
	result.a1=x1.a1+x2.a1;
	result.a2=x1.a2+x2.a2;
	result.b1=x1.b1+x2.b1;
	result.b2=x1.b2+x2.b2;
	return result;
}
bmatrix operator*(double alf, bmatrix x2)
{
	bmatrix result;
	result.a1=alf*x2.a1;
	result.a2=alf*x2.a2;
	result.b1=alf*x2.b1;
	result.b2=alf*x2.b2;
	return result;
}

bmatrix operator*(xvector x1, xvector x2)
{
	bmatrix result;
	result.a1=x1.x1*x2.x1;
	result.a2=x1.x1*x2.x2;
	result.b1=x1.x2*x2.x1;
	result.b2=x1.x2*x2.x2;
	return result;	
}

bmatrix operator*(bmatrix x1, bmatrix x2)
{
	bmatrix result;
	result.a1=x1.a1*x2.a1+x1.a2*x2.b1;
	result.a2=x1.a1*x2.a2+x1.a2*x2.b2;
	result.b1=x1.b1*x2.a1+x1.b2*x2.b1;
	result.b2=x1.b1*x2.a2+x1.b2*x2.b2;
	return result;
}

xvector operator*(bmatrix x1, xvector x2)
{
	xvector result;
	result.x1=x2.x1*x1.a1+x2.x2*x1.a2;
	result.x2=x2.x1*x1.b1+x2.x2*x1.b2;
	return result;
}
double times(xvector x1, xvector x2)
{
	return x1.x1*x2.x1+x1.x2*x2.x2;
}
