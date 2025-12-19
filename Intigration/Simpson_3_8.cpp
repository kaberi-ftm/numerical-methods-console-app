#include <bits/stdc++.h>
using namespace std;


double f(double x)
{
    return 1.0/(1.0+x*x);
}

int main()
{
    ifstream fin("Simpson_3_8_input.txt");
    ofstream fout("Simpson_3_8_output.txt");

    double a, b;
    int n;
   fout<<"Upper limit : ";
    fin >> a ;
    fout<<endl;
   fout<<"Lower limit : ";
    fin>>b;
    fout<<endl;
fout<<" Number of intervals( multiple of 3) : ";
    fin>>n;
fout<<endl;
    double h=(b-a)/n;
    double sum = 0.0;

    for(int i=0;i<=n;i++)
    {
         double x=a+i*h;

        if(i == 0||i == n)
            sum+=f(x);
        else if(i%3 == 0)
            sum+=2*f(x);
        else
            sum+=3*f(x);
    }

    double integral=(3*h*sum) / 8;

    fout<<fixed<<setprecision(6);
    fout<<"Integral using Simpson's 3/8 Rule = "<<integral<<endl;

    fin.close();
    fout.close();

    return 0;
}
