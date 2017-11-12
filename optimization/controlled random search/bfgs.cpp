//-----------------bfgs.cpp---------------------------
/*this is the implementation of the steepest descent, Newton's method and bfgs mehtod*/

#include <iostream>
#include <cmath>

using namespace std;

const double E=0.000001;           // the tolerance was set to be 10^-6
const double A=1;                  // initial guess of step size is 1


struct xvector{                    //create structures to represent vector and matrix
	double x1;
	double x2;
};

struct bmatrix{
	double a1;
	double a2;
	double b1;
	double b2;
};

double f(xvector x);                              //evaluated function
xvector gradf(xvector x);                         //find gradient
bmatrix hessian(xvector x);                       //generate hessian matirx
bmatrix neginverse(bmatrix b);                    //inverse a 2*2 matrix
xvector operator+(xvector x1, xvector x2);        //overload other operator for 2*2 matrix
xvector operator*(double alf, xvector x2);

bmatrix operator+(bmatrix x1, bmatrix x2);
bmatrix operator*(double alf, bmatrix x2);
bmatrix operator*(bmatrix x1, bmatrix x2);

bmatrix operator*(xvector x1, xvector x2);
xvector operator*(bmatrix x1, xvector x2);

double times(xvector x1, xvector x2);              //times two vectors

int main()
{


	
	double tol, alf, r=0.5, c=0.001;
	xvector xk, deltaf, p, s, y;
	int count;

        //bfgs method
	cout<<"bfgs method:\n";
	xk.x1=0;
	xk.x2=0;
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
	cout<<"x= "<<xk.x1<<" y= "<<xk.x2<<" f(x,y)= "<<f(xk)<<" iteration= "<<count<<'\n';
}





double f(xvector x)
{
	return 100*pow(x.x2-x.x1*x.x1,2)+pow(6.4*pow(x.x2-0.5,2)-x.x1-0.6,2);
}
xvector gradf(xvector x)
{
	xvector result;
	result.x1=-400*x.x1*(x.x2-x.x1*x.x1)-2*(6.4*pow(x.x2-0.5,2)-x.x1-0.6);
	result.x2=200*(x.x2-x.x1*x.x1)+2*(6.4*pow(x.x2-0.5,2)-x.x1-0.6)*(12.8*x.x2-6.4);
	return result;
}
bmatrix hessian(xvector x)
{
	bmatrix result;
	result.a1=800*x.x1*x.x1+2-400*(x.x2-x.x1*x.x1);
	result.a2=result.b1=-400*x.x1-25.6*(x.x2-0.5);
	result.b2=200+2*pow(12.8*(x.x2-0.5),2)+25.6*(6.4*pow(x.x2-0.5,2)-x.x1-0.6);
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
