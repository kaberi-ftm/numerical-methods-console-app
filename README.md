- [Solution of Linear Equations](#solution-of-linear-equations)
  - [Gauss Elimination Method](#gauss-elimination-method)
    - [Theory](#gauss-elimination-theory)
    - [Code](#gauss-elimination-code)
    - [Input](#gauss-elimination-input)
    - [Output](#gauss-elimination-output)

  - [Gauss Jordan Elimination Method](#gauss-jordan-elimination-method)
    - [Theory](#gauss-jordan-theory)
    - [Code](#gauss-jordan-code)
    - [Input](#gauss-jordan-input)
    - [Output](#gauss-jordan-output)

  - [LU Decomposition Method](#lu-decomposition-method)
    - [Theory](#lu-decomposition-theory)
    - [Code](#lu-decomposition-code)
    - [Input](#lu-decomposition-input)
    - [Output](#lu-decomposition-output)

  - [Inverse Matrix Method](#inverse-matrix-method)
    - [Theory](#inverse-matrix-theory)
    - [Code](#inverse-matrix-code)
    - [Input](#inverse-matrix-input)
    - [Output](#inverse-matrix-output)

- [Solution of Non-Linear Equations](#solution-of-non-linear-equations)
  - [Bisection Method](#bisection-method)
    - [Theory](#bisection-theory)
    - [Code](#bisection-code)
    - [Input](#bisection-input)
    - [Output](#bisection-output)

  - [False Position Method](#false-position-method)
    - [Theory](#false-position-theory)
    - [Code](#false-position-code)
    - [Input](#false-position-input)
    - [Output](#false-position-output)

  - [Newton Raphson Method](#newton-raphson-method)
    - [Theory](#newton-raphson-theory)
    - [Code](#newton-raphson-code)
    - [Input](#newton-raphson-input)
    - [Output](#newton-raphson-output)

  - [Secant Method](#secant-method)
    - [Theory](#secant-theory)
    - [Code](#secant-code)
    - [Input](#secant-input)
    - [Output](#secant-output)

- [Interpolation](#interpolation)
  - [Forward Interpolation](#forward-interpolation)
    - [Theory](#forward-interpolation-theory)
    - [Code](#forward-interpolation-code)
    - [Input](#forward-interpolation-input)
    - [Output](#forward-interpolation-output)

  - [Backward Interpolation](#backward-interpolation)
    - [Theory](#backward-interpolation-theory)
    - [Code](#backward-interpolation-code)
    - [Input](#backward-interpolation-input)
    - [Output](#backward-interpolation-output)

  - [Divided Difference Method](#divided-difference-method)
    - [Theory](#divided-difference-theory)
    - [Code](#divided-difference-code)
    - [Input](#divided-difference-input)
    - [Output](#divided-difference-output)

- [Numerical Differentiation](#numerical-differentiation)
  - [Differentiation Using Interpolation](#differentiation-using-interpolation)
    - [Theory](#differentiation-theory)
    - [Code](#differentiation-code)
    - [Input](#differentiation-input)
    - [Output](#differentiation-output)

- [Numerical Integration](#numerical-integration)
  - [Simpson’s 1/3 Rule](#simpsons-13-rule)
    - [Theory](#simpson-13-theory)
    - [Code](#simpson-13-code)
    - [Input](#simpson-13-input)
    - [Output](#simpson-13-output)

  - [Simpson’s 3/8 Rule](#simpsons-38-rule)
    - [Theory](#simpson-38-theory)
    - [Code](#simpson-38-code)
    - [Input](#simpson-38-input)
    - [Output](#simpson-38-output)

- [Differential Equations](#differential-equations)
  - [Runge Kutta 4th Order Method](#runge-kutta-4th-order-method)
    - [Theory](#rk4-theory)
    - [Code](#rk4-code)
    - [Input](#rk4-input)
    - [Output](#rk4-output)

- [Regression Analysis](#regression-analysis)
  - [Linear Regression](#linear-regression)
    - [Theory](#linear-regression-theory)
    - [Code](#linear-regression-code)
    - [Input](#linear-regression-input)
    - [Output](#linear-regression-output)

  - [Polynomial Regression](#polynomial-regression)
    - [Theory](#polynomial-regression-theory)
    - [Code](#polynomial-regression-code)
    - [Input](#polynomial-regression-input)
    - [Output](#polynomial-regression-output)

  - [Transcendental Regression](#transcendental-regression)
    - [Theory](#transcendental-regression-theory)
    - [Code](#transcendental-regression-code)
    - [Input](#transcendental-regression-input)
    - [Output](#transcendental-regression-output)

### Solution of Linear Equations

### Gauss Elimination Method

#### Gauss Elimination Theory


#### Gauss Elimination Code
```cpp
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
```



---
