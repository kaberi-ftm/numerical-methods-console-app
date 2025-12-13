#include <iostream>
#include <bits/stdc++.h>
#include<vector>
#include<cmath>
using namespace std;
int n;
vector<double> a;
double f(double x)
{ double sum = 0;
    for (int i = 0; i <= n; i++) {
        sum+=a[i]*pow(x,n-i);
    }
    return sum;
}
double df(double x){
double sum = 0;
    for (int i = 0; i < n; i++) {
        sum+=a[i]*(n-i)*pow(x,n-i-1);
    }
    return sum;

}

int main()
{
    ifstream fin("NewtonRaphson_input.txt");
    ofstream fout("NewtonRaphson_output.txt");
    fout<<"NO of Degree of the equation : ";
    fout<<endl<<endl;
    fin>>n;
    a.resize(n + 1);
    for (int i = 0; i <= n; i++) {
        fin>>a[i];
    }

    fout<<"The Function f = ";
     for (int i=0; i<=n; i++) {
        fout <<a[i]<<"x^"<<(n - i);
        if(i!=n)fout<< " + ";
    }

    fout<<endl<<endl;
    

   double xmax=1;
    for(int i=1;i<=n;i++)
    {
        xmax=max(xmax,1+fabs(a[i]/a[0]));
    }

double step=0.45,a=-xmax,b;
double e=1e-3;
int c=0;
while(a<xmax){
    b=a+step;
    if(f(a)*f(b)<0){

        c++;
        fout<<c<<" :"<<endl;

        double x1=a,x2=b,xr;
        int it=0;
        do{
                 it++;
            x2= x1-(f(x1)/df(x1));
            xr=x1;
            x1=x2;

       

        }while(fabs(f(x2))>e && fabs(x2-xr)>e);

        fout<<"Root"<<c<<"  :"<<x2<<endl;
        fout<<"NUmber of iterations :"<<it<<endl;
        fout<<"Search Interval for root"<<c<<" ["<<a<<","<<b<<"]"<<endl<<endl<<endl;

    }
    a=b;
}




     return 0;
}
