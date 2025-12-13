#include <bits/stdc++.h>
using namespace std;

  double f(double x, double y)
  {
      return ((x - y)/2);
  }

int main(){
    fstream fin("RKinput.txt");
     ofstream fout("RKoutput.txt");
   if(!fin)
   {
       cout<<"Can't open input file"<<endl;
       return 1;
   }
    if(!fout)
   {
       cout<<"Can't open output file"<<endl;
       return 1;
   }

 double h,xn,x0,y0;
    fin>>h>>xn>>x0>>y0;

 //"Enter step size h, final xn, initial x0, initial y0:\n";
 cout << "Read values: h=" << h << " xn=" << xn << " x0=" << x0 << " y0=" << y0 << endl;

  double x=x0;
  double y=y0;
  int s=(xn-x0)/h;
    for( int i=0;i<s ;i++)  {
      double k1= h*f(x,y);
      double k2= h*f(x+h/2.0,y+k1/2.0);
      double k3= h*f(x+h/2.0,y+k2/2.0);
      double k4= h*f(x+h,y+k3);
      y=y+(k1+2*k2+2*k3+k4)/6.0;
      x=x+h;

  }
  fout<<"Value of y at x is: "<<y<<endl;
  fin.close();
  fout.close();

return 0;

}
