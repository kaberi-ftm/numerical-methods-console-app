#include <iostream>
#include <bits/stdc++.h>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

const int m = 100;

void print_matrix(double mat[][m], int n, ofstream &fout)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            fout << setw(10) << mat[i][j];
        fout << endl;
    }
}

bool lu(double mat[][m], double l[][m], double u[][m], int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            l[i][j] = u[i][j] = 0;

    for (int i = 0; i < n; i++)
    {
        for (int k = i; k < n; k++)
        {
            double sum = 0;
            for (int j = 0; j < i; j++)
                sum += l[i][j] * u[j][k];
            u[i][k] = mat[i][k] - sum;
        }

        for (int k = i; k < n; k++)
        {
            if (i == k)
                l[i][i] = 1;
            else
            {
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum += l[k][j] * u[j][i];
                if (fabs(u[i][i]) < 1e-9)
                    return false;
                l[k][i] = (mat[k][i] - sum) / u[i][i];
            }
        }
    }

    if (fabs(u[n - 1][n - 1]) < 1e-9)
        return false;

    return true;
}

void ForwardSubstitution(double l[][m], double b[], double y[], int n)
{
    for (int i = 0; i < n; i++)
    {
        double sum = 0;
        for (int j = 0; j < i; j++)
            sum += l[i][j] * y[j];
        y[i] = b[i] - sum;
    }
}

void BackwardSubstitution(double u[][m], double x[], double y[], int n)
{
    for (int i = n - 1; i >= 0; i--)
    {
        double sum = 0;
        for (int j = i + 1; j < n; j++)
            sum += u[i][j] * x[j];
        x[i] = (y[i] - sum) / u[i][i];
    }
}

int main()
{
    ifstream fin("LU_Decomposition_input.txt");
    ofstream fout("LU_Decomposition_output.txt");

    int n;
    double a[m][m], l[m][m], u[m][m];
    double b[m], x[m], y[m];

    
        fout << "Enter num of eqn: ";
        fin >> n;

        fout << "Enter the augumented matrix :" << endl;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                fin >> a[i][j];
            fin >> b[i];
        }

        bool c = lu(a, l, u, n);

        if (c)
        {
            fout << "Lower Triangular Matrix :" << endl;
            print_matrix(l, n, fout);

        fout<<endl<<endl;

            fout << "Upper Triangular Matrix :" << endl;
            print_matrix(u, n, fout);

        fout<<endl<<endl;

            ForwardSubstitution(l, b, y, n);
            BackwardSubstitution(u, x, y, n);

            fout << "The solution " << endl;
            for (int i = 0; i < n; i++)
                fout << "x" << i + 1 << " = " << x[i] << endl;

            fout<<endl<<endl;

            fout << "The solution is unique" << endl;
        }
        else
        {
            fout << " The solution is infinite/no solution" << endl;
        }

     
    return 0;
}
