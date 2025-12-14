#include <bits/stdc++.h>
using namespace std;

double fact(int n)
{
    double f = 1;
    for(int i = 2; i<= n; i++)
    {
        f *= i;
    }
    return f;
}

int main()
{
    ifstream fin("Backward_Interpolation_input.txt");
    ofstream fout("Backward_Interpolation_output.txt");

    int n;
    fin>>n;

    double x[50],y[50][50],X;

    for(int i= 0; i<n; i++)
        fin>>x[i];

    for(int i= 0; i < n; i++)
        fin>>y[i][0];

    fin>>X;



    for(int j=1;j< n;j++)
        for(int i=0;i<n-j;i++)
            y[i][j]= y[i+1][j-1]- y[i][j-1];



    fout << "Backward Interpolation Table :"<< endl<<endl;
    for(int i = 0; i < n; i++)
    {
        fout << setw(8) << x[i];
        for(int j = 0; j < n - i; j++)
            fout << setw(8) << y[i][j];
        fout << endl;
    }
fout <<endl<<endl;


    double h =x[1]-x[0];

    double u =(X-x[n-1])/h;



    double ans = y[n - 1][0];
    double u_prod = 1;

    for(int i = 1;i<n;i++)
    {
        u_prod *=(u +(i-1));
        ans += (u_prod*y[n-i-1][i])/fact(i);
    }



    fout << ans << endl<<endl;

   

    fin.close();
    fout.close();

    return 0;
}
