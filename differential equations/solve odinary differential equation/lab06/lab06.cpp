//----------------------------lab06.cpp-------------------------------------
/*this is the code for lab06 in which for method, FE, BE, IT and LF was implemented 
 *to solve ode. The magnitude of the eigenvalue of matrix G of the four method was then
  calculated. the data will then be plotted using matlab
 * Programmed by: Chaolun Wang
             at: 04/07/2016
 */



#include <iostream>
#include <fstream>
#include <stdlib.h>     
#include <time.h>
#include <cmath>
#include <complex>
 
using namespace std;

const double H=M_PI*2/32;
const double TA=0;
const double TB=10*M_PI;

//update the value of x and y form x_n, y_n to x_{n+1}  y_{n+1}, the a,b,c,d indecate the four element of matrix G
void update(double &x, double &v, double a, double b, double c, double d);      

int main () {
  ofstream file;
  

  int counter=0;
  double x0=1;
  double v0=0;
  
  file.open ("q11result.txt");                                       // the FE method
  file<<0<<' '<<x0<<' '<<v0<<' '<<0.5<<' '<<0<<'\n';                 //output initial stage
  while(counter*H+TA<TB)
  {
      update(x0,v0,1,H,-1*H,1);                                      //update with matrix G
      counter++;
      file<<counter*H+TA<<' '<<x0<<' '<<v0<<' '<<0.5*(x0*x0+v0*v0)<<' '<<0.5*(x0*x0+v0*v0)-0.5<<'\n';  //output the time, x, v, energy, error
  }
  file.close();
  
  
  counter=0;
  x0=1;
  v0=0;
  
  file.open ("q12result.txt");                                       //the BE method
  file<<0<<' '<<x0<<' '<<v0<<' '<<0.5<<' '<<0<<'\n';                 //output initial stage
  while(counter*H+TA<TB)
  {
      update(x0,v0,1/(1+H*H),H/(1+H*H),-1*H/(1+H*H),1/(1+H*H));      //update with matrix G
      counter++;
      file<<counter*H+TA<<' '<<x0<<' '<<v0<<' '<<0.5*(x0*x0+v0*v0)<<' '<<0.5*(x0*x0+v0*v0)-0.5<<'\n';   //output the time, x, v, energy, error
  }
  file.close();
  
  counter=0;
  x0=1;
  v0=0;
  
  file.open ("q13result.txt");                                        //the IT method
  file<<0<<' '<<x0<<' '<<v0<<' '<<0.5<<' '<<0<<'\n';                  //output initial stage
  while(counter*H+TA<TB)
  {
      update(x0,v0,(-1*H*H+4)/(4+H*H),4*H/(4+H*H),-4*H/(4+H*H),(-1*H*H+4)/(4+H*H));  //update with matrix G
      counter++;
      file<<counter*H+TA<<' '<<x0<<' '<<v0<<' '<<0.5*(x0*x0+v0*v0)<<' '<<0.5*(x0*x0+v0*v0)-0.5<<'\n';         //output the time, x, v, energy, error
  }
  file.close();
  
    counter=0;
  x0=1;
  v0=0;
  
  file.open ("q14result.txt");
  file<<0<<' '<<x0<<' '<<v0<<' '<<0.5<<' '<<0<<'\n';                   //output initial stage
  while(counter*H+TA<TB)
  {
      update(x0,v0,1,H,-1*H,1-H*H);                                    //update with matrix G
      counter++;
      file<<counter*H+TA<<' '<<x0<<' '<<v0<<' '<<0.5*(x0*x0+v0*v0)<<' '<<0.5*(x0*x0+v0*v0)-0.5<<'\n';         //output the time, x, v, energy, error
  }
  file.close();

  cout<<"     Lambda1/Lambda2   Magnitude of Lambda:\n";                                               //calculate the eigen value property of the four method
  complex<double> L11(1, H), L12(1,-1*H);
  cout<<"FE:  "<<abs(L12/L11)<<"                 "<<abs(L11)<<'\n';
  
  complex<double> L21(1.0/(1.0+H*H), H/(1.0+H*H)), L22(1/(1.0+H*H), -1.0*H/(1.0+H*H));
  cout<<"BE:  "<<abs(L22/L21)<<"                 "<<abs(L21)<<'\n';
  
  complex<double> L31((4.0-H*H)/(4.0+H*H), 4.0*H/(4.0+H*H)), L32((4.0-H*H)/(4.0+H*H), -4.0*H/(4.0+H*H));
  cout<<"IT:  "<<abs(L32/L31)<<"                 "<<abs(L31)<<'\n';
 
  complex<double> L41((2.0-H*H)/2.0, sqrt(4.0*H*H-H*H*H*H)/2.0), L42((2.0-H*H)/2.0, sqrt(4.0*H*H-H*H*H*H)/-2.0);
  cout<<"LF:  "<<abs(L32/L31)<<"                 "<<abs(L41)<<'\n';   
  return 0;
}



void update(double &x, double &v, double a, double b, double c, double d) //up date the value of x and y form x_n, y_n to x_{n+1}  y_{n+1}
{
	double x0=x;
	double v0=v;
	x=a*x0+b*v0;
	v=c*x0+d*v0;
}





