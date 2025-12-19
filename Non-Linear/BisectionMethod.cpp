#include <bits/stdc++.h>
using namespace std;

double f(double x){
return x*x*x-x*x+2 ;
}

int main() {
     fstream fin("BisecIn.txt");
     ofstream fout("BisecOut.txt");
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

    double x1,x2,e;
    fin>>x1>>x2>>e;
    if(f(x1)*f(x2)>=0)
    {
        fout<<"No root in the given interval"<<endl;
        return 1;

    }
        double x0;
    while (true){
    x0=(x1+x2)/2;
    double f1=f(x1);
    double f2=f(x2);
    double f0=f(x0);
   if( f0== 0 || abs((x2-x1)/x2 )<e)
   {
       fout<<"Root : "<<x0<<endl;
       break;
   }
   if(f1*f0 <0)
   {
       x2=x0;
   }
   else
   {
       x1=x0;
   }
    }

  fin.close();
  fout.close();
    return 0;
}
