#include <bits/stdc++.h>
using namespace std;

int main()
{
    ifstream fin("RegressionTranIn.txt");
    ofstream fout("RegressionTranOut.txt");

    if (!fin || !fout)
    {
        cout << "File can't be opened ." << endl;
        return 1;
    }

    fout << fixed << setprecision(3);

    int n;
    fin >> n;

    vector<double> x(n), y(n), Y(n);
    for (int i = 0; i < n; i++)
    {
        fin >> x[i] >> y[i];
        Y[i] = log(y[i]);
    }

    double sx = 0, sY = 0, sxY = 0, sx2 = 0;

    for (int i = 0; i < n; i++)
    {
        sx += x[i];
        sY += Y[i];
        sxY += x[i] * Y[i];
        sx2 += x[i] * x[i];
    }

    double b = (n * sxY - sx * sY) / (n * sx2 - sx * sx);
    double A = (sY - b * sx) / n;
    double a = exp(A);

    fout << "Exponential Regression Equation:" << endl;
    fout << "y = " << a << " * e^(" << b << "x)" << endl;

    return 0;
}
