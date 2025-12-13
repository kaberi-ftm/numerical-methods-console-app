#include<bits/stdc++.h>
using namespace std;

const double EPS=1e-3;
const double STEP=0.45;

double f(double x,const vector<double>&a)
{
    double r=0;
    int n=a.size()-1;
    for(int i=0;i<=n;i++)
    {
        r+=a[i]*pow(x,n-i);
    }
    return r;
}

double secant(double x0,double x1,const vector<double>&a,int &it)
{
    double x2;
    it=0;
    while(true)
    {
        double f0=f(x0,a);

        double f1=f(x1,a);
        x2=x1-f1*(x1-x0)/(f1-f0);

        it++;
        if(fabs(x2-x1)<EPS && fabs(f(x2,a))<EPS)
        {
            break;
        }
        x0=x1;
        x1=x2;
    }
    return x2;
}

int main()
{
    ifstream fin("Secant_input.txt");
    ofstream fout("Secant_output.txt");

    int n;
    fin>>n;

    vector<double>a(n+1);


    for(int i=0;i<=n;i++)
    {
        fin>>a[i];
    }

    fout<<"Function: ";

    for(int i=0;i<=n;i++)
    {
        if(a[i]==0){continue;}
        if(i!=0 && a[i]>0){fout<<"+";}
        fout<<a[i];
        if(n-i>0){fout<<"x";}
        if(n-i>1){fout<<"^"<<(n-i);}
    }
    fout<<endl<<endl;

    double xmax=0;
    for(int i=1;i<=n;i++)
    {
        xmax=max(xmax,fabs(a[i]/a[0]));
    }
    xmax=1+xmax;

    int c=0;
    for(double x=-xmax; x<xmax; x+=STEP)
    {
        double l=x;
        double r=x+STEP;
        if(f(l,a)*f(r,a)<=0)
        {
            c++;
            int it;
            double root=secant(l,r,a,it);
            
            fout<<"Root "<<c<<": "<<root<<endl;

            fout<<"Search interval for root "<<c<<" = ["<<l<<","<<r<<"]"<<endl;

            fout<<"Iteration needed for root "<<c<<" = "<<it<<endl;

              fout<<endl<<endl;
        }
    }

    fin.close();
    fout.close();
    return 0;
}
