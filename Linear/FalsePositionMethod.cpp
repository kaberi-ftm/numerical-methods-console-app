#include <bits/stdc++.h>
using namespace std;

double f(double x){
return x*x-9 ;
}

int main() {
     fstream fin("FalsePIn.txt");
     ofstream fout("FalsePOut.txt");
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
        int iter=0;
    while (iter<1000){
            iter++;
    double f1=f(x1);
    double f2=f(x2);
    if(f1==f2)
    {
        fout<<"Division by zero error in False Position method ."<<endl;
        break;
    }
    x0=x2- f2*(x1-x2)/(f1-f2);
    double f0=f(x0);

   if( abs(f0)< e || abs(x2-x1) < e)
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
   if(iter>1000)
   {
       fout<<"Maximum iteration exceeded ."<<endl;
       break;
   }


    }
    fin.close();
  fout.close();
    return 0;
}
