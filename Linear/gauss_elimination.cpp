#include <iostream>
#include <vector>
#include <bits/stdc++.h>
#include <fstream>

using namespace std;

int main()
{
    ifstream fin("gauss_elimination_input.txt");
    ofstream fout("gauss_elimination_output.txt");

    
        int n;
        fout<<"Enter Number of equations :" << endl;
        fin>>n;

        vector<vector<double>> a(n, vector<double>(n + 1));
        fout << "Enter the augmented matrix:" << endl;
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= n; j++)
                fin >> a[i][j];

        for (int i = 0; i < n - 1; i++)
        {
            if (fabs(a[i][i]) < 1e-9)
            {
                bool s = false;
                for (int k = i + 1; k < n; k++)
                {
                    if (fabs(a[k][i]) > 1e-9)
                    {
                        swap(a[i], a[k]);
                        s = true;
                        break;
                    }
                }
                if (!s) continue;
            }

            for (int k = i + 1; k < n; k++)
            {
                double f = a[k][i] / a[i][i];
                for (int j = i; j <= n; j++)
                    a[k][j] -= f * a[i][j];
            }
        }

        bool noSol = false, infinite = false;
        for (int i = 0; i < n; i++)
        {
            bool ze = true;
            for (int j = 0; j < n; j++)
                if (fabs(a[i][j]) > 1e-9)
                {
                    ze = false;
                }

            if (ze && fabs(a[i][n]) > 1e-9)
            {
                noSol = true;
            }
            else if (ze && fabs(a[i][n]) < 1e-9)
            {
                infinite = true;
            }
        }

        fout << endl;
        fout << "Echalon form:" << endl;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j <= n; j++)
                fout << fixed << setprecision(3) << setw(8) << a[i][j];
            fout << endl;
        }
        fout << endl << endl;

        if (noSol)
        {
            fout << "The system has no solution" << endl << endl << endl;
        }
        else if (infinite)
        {
            fout << "The system has infinite solution" << endl;
        }
        else
        {
            vector<double> x(n);
            for (int i = n - 1; i >= 0; i--)
            {
                x[i] = a[i][n];
                for (int j = i + 1; j < n; j++)
                    x[i] -= a[i][j] * x[j];
                x[i] /= a[i][i];
            }

            for (int i = 0; i < n; i++)
            {
                fout << "x" << i + 1 << " = "
                     << fixed << setprecision(3) << x[i] << endl;
            }
            fout << "The system has unique solution" << endl;
        }

      

    return 0;
}
