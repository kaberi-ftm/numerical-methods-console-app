#include <bits/stdc++.h>
using namespace std;

int main()
{
    ifstream fin("DividedDifference_input.txt");
    ofstream fout("DividedDifference_output.txt");

    int n;
    fin>>n;

    double x[50], y[50][50], X;

    for(int i=0;i<n;i++)
        fin>>x[i];

    for(int i=0;i<n;i++)
        fin>>y[i][0];

    fin>>X;



    for(int j=1;j<n;j++)
        for(int i= 0;i < n-j;i++)
            y[i][j]=(y[i+1][j-1]-y[i][j-1])/(x[i+j]-x[i]);



    fout<<"Newton's Divided Difference interpolation Table :"<<endl;
    for(int i = 0;i< n;i++)
    {
        fout<<setw(7)<<x[i];
        for(int j = 0;j<n-i;j++)
            fout<<setw(7)<<y[i][j];
        fout<<endl;
    }



    double ans =y[0][0];
    double term=1;

    for(int i = 1;i<n;i++)
    {
        term *=(X-x[i-1]);
        ans +=term*y[0][i];
    }


    fout<<ans<<endl;
  

    fin.close();
    fout.close();

    return 0;
}
