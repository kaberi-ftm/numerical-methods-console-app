#include <bits/stdc++.h>
using namespace std;

int main()
{
    ifstream fin("RegressionPolyIn.txt");
    ofstream fout("RegressionPolyOut.txt");

    if (!fin || !fout)
    {
        cout << "File can't be opened ." << endl;
        return 1;
    }

    fout << fixed << setprecision(3);

    int n;
    fin >> n;
    vector<double> x(n), y(n);
    for (int i = 0; i < n; i++)
       {
            fin >> x[i] >> y[i];
       }

    double sx = 0, sx2 = 0, sx3 = 0, sx4 = 0;
    double sy = 0, sxy = 0, sx2y = 0;

    for (int i = 0; i < n; i++)
    {
        sx += x[i];
        sx2 += x[i] * x[i];
        sx3 += x[i] * x[i] * x[i];
        sx4 += x[i] * x[i] * x[i] * x[i];
        sy += y[i];
        sxy += x[i] * y[i];
        sx2y += x[i] * x[i] * y[i];
    }

    int eq = 3;
    vector<vector<double>> a(eq, vector<double>(eq + 1));

    a[0] = {(double)n, sx,  sx2, sy   };
    a[1] = { sx,  sx2, sx3, sxy  };
    a[2] = { sx2, sx3, sx4, sx2y };

    for (int i = 0; i < eq - 1; i++)
    {
        if (fabs(a[i][i]) < 1e-9)
        {
            for (int k = i + 1; k < eq; k++)
            {
                if (fabs(a[k][i]) > 1e-9)
                {
                    swap(a[i], a[k]);
                    break;
                }
            }
        }

        for (int k = i + 1; k < eq; k++)
        {
            double factor = a[k][i] / a[i][i];
            for (int j = i; j <= eq; j++)
                a[k][j] -= factor * a[i][j];
        }
    }

    vector<double> coef(eq);
    for (int i = eq - 1; i >= 0; i--)
    {
        coef[i] = a[i][eq];
        for (int j = i + 1; j < eq; j++)
            coef[i] -= a[i][j] * coef[j];
        coef[i] /= a[i][i];
    }

    fout << "Polynomial Regression Equation:" << endl;
    fout << "y = " << coef[0] << " + " << coef[1]
         << "x + " << coef[2] << "x^2" << endl;

    return 0;
}
