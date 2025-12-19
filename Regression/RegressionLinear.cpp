#include <bits/stdc++.h>
using namespace std;

int main()
{
    ifstream fin("RegressionlinearIn.txt");
    ofstream fout("RegressionlinearOut.txt");

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

    double sx = 0,sy = 0,sxy = 0,sx2 = 0;

    for (int i = 0; i < n; i++)
    {
        sx+= x[i];
        sy+= y[i];
        sxy+= x[i] * y[i];
        sx2+= x[i]*x[i];
    }

    double b = (n * sxy - sx * sy) / (n * sx2 - sx * sx);
    double a = (sy - b * sx) / n;

    fout << "Linear Regression Equation: " << endl;
    fout << "y = " << a << " + " << b << "x" << endl;

    return 0;
}
