
---

[Solution of Linear Equations  :](#solution-of-linear-equations)
----


| Algorithm | 📄 Concept | 💻 Source | 📥 Test Data | 📤 Result |
| :--- | :---: | :---: | :---: | :---: |
| [**Gauss Elimination**](#gauss-elimination-method) | [Theory](#gauss-elimination-theory) | [Code](#gauss-elimination-code) | [Input](#gauss-elimination-input) | [Output](#gauss-elimination-output) |
| [**Gauss Jordan Elimination**](#gauss-jordan-elimination-method) | [Theory](#gauss-jordan-theory) | [Code](#gauss-jordan-code) | [Input](#gauss-jordan-input) | [Output](#gauss-jordan-output) |
| [**LU Decomposition**](#lu-decomposition-method) | [Theory](#lu-decomposition-theory) | [Code](#lu-decomposition-code) | [Input](#lu-decomposition-input) | [Output](#lu-decomposition-output) |
| [**Inverse Matrix**](#inverse-matrix-method) | [Theory](#inverse-matrix-theory) | [Code](#inverse-matrix-code) | [Input](#inverse-matrix-input) | [Output](#inverse-matrix-output) |

---

[Solution of NonLinear Equations  :](#solution-of-linear-equations)
----

| Algorithm | 📄 Concept | 💻 Source | 📥 Test Data | 📤 Result |
| :--- | :---: | :---: | :---: | :---: |
| [**Bisection Method**](#bisection-method) | [Theory](#bisection-theory) | [Code](#bisection-code) | [Input](#bisection-input) | [Output](#bisection-output) |
| [**False Position Method**](#false-position-method) | [Theory](#false-position-theory) | [Code](#false-position-code) | [Input](#false-position-input) | [Output](#false-position-output) |
| [**Newton Raphson Method**](#newton-raphson-method) | [Theory](#newton-raphson-theory) | [Code](#newton-raphson-code) | [Input](#newton-raphson-input) | [Output](#newton-raphson-output) |
| [**Secant Method**](#secant-method) | [Theory](#secant-theory) | [Code](#secant-code) | [Input](#secant-input) | [Output](#secant-output) |

---

[Solution of Interpolations  :](#solution-of-linear-equations)
----

| Algorithm | 📄 Concept | 💻 Source | 📥 Test Data | 📤 Result |
| :--- | :---: | :---: | :---: | :---: |
| [**Forward Interpolation**](#forward-interpolation-method) | [Theory](#forward-interpolation-theory) | [Code](#forward-interpolation-code) | [Input](#forward-interpolation-input) | [Output](#forward-interpolation-output) |
| [**Backward Interpolation**](#Backward-Interpolation) | [Theory](#backward-interpolation-theory) | [Code](#backward-interpolation-code) | [Input](#backward-interpolation-input) | [Output](#backward-interpolation-output) |
| [**Divided Difference Method**](#divided-difference-method) | [Theory](#divided-difference-theory) | [Code](#divided-difference-code) | [Input](#divided-difference-input) | [Output](#divided-difference-output) |

---

[Solution of Differentiation  :](#solution-of-linear-equations)
----

| Algorithm | 📄 Concept | 💻 Source | 📥 Test Data | 📤 Result |
| :--- | :---: | :---: | :---: | :---: |
| [**Differentiation Using Interpolation**](#differentiation-using-interpolation) | [Theory](#differentiation-theory) | [Code](#differentiation-code) | [Input](#differentiation-input) | [Output](#differentiation-output) |

---

[Solution of Intigration  :](#solution-of-linear-equations)
----

| Algorithm | 📄 Concept | 💻 Source | 📥 Test Data | 📤 Result |
| :--- | :---: | :---: | :---: | :---: |
| [**Simpson’s 1/3 Rule**](#simpsons-13-rule) | [Theory](#simpson-13-theory) | [Code](#simpson-13-code) | [Input](#simpson-13-input) | [Output](#simpson-13-output) |
| [**Simpson’s 3/8 Rule**](#simpsons-38-rule) | [Theory](#simpson-38-theory) | [Code](#simpson-38-code) | [Input](#simpson-38-input) | [Output](#simpson-38-output) |

---

[Solution of Differential Equations  :](#solution-of-differential-equations)
----

| Algorithm | 📄 Concept | 💻 Source | 📥 Test Data | 📤 Result |
| :--- | :---: | :---: | :---: | :---: |
| [**Runge Kutta 4th Order Method**](#runge-kutta-4th-order-method) | [Theory](#rk4-theory) | [Code](#rk4-code) | [Input](#rk4-input) | [Output](#rk4-output) |

---

[Solution of Regression  :](#solution-of-linear-equations)
----

| Algorithm | 📄 Concept | 💻 Source | 📥 Test Data | 📤 Result |
| :--- | :---: | :---: | :---: | :---: |
| [**Linear Regression**](#linear-regression) | [Theory](#linear-regression-theory) | [Code](#linear-regression-code) | [Input](#linear-regression-input) | [Output](#linear-regression-output) |
| [**Polynomial Regression**](#polynomial-regression) | [Theory](#polynomial-regression-theory) | [Code](#polynomial-regression-code) | [Input](#polynomial-regression-input) | [Output](#polynomial-regression-output) |
| [**Transcendental Regression**](#transcendental-regression) | [Theory](#transcendental-regression-theory) | [Code](#transcendental-regression-code) | [Input](#transcendental-regression-input) | [Output](#transcendental-regression-output) |


## Solution of Linear Equations

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
#### Gauss Elimination Input
```
3
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3

```
#### Gauss Elimination Output
```
Enter Number of equations :
Enter the augmented matrix:

Echalon form:
   2.000   1.000  -1.000   8.000
   0.000   0.500   0.500   1.000
   0.000   0.000  -1.000   1.000


x1 = 2.000
x2 = 3.000
x3 = -1.000
The system has unique solution
```

### Gauss Jordan elimination Method
#### Gauss Jordan Theory
#### Gauss Jordan Code
#### Gauss Jordan Input
#### Gauss Jordan Output

### LU Decomposition Method
#### LU Decomposition Theory
#### LU Decomposition Code
#### LU Decomposition Input
#### LU Decomposition Output

### Inverse Matrix Method
#### Inverse Matrix Theory
#### Inverse Matrix Code
#### Inverse Matrix Input
#### Inverse Matrix Output

### bisection-method
#### bisection Theory
#### bisection Code
#### bisection Input
#### bisection Output

### false-position-method
#### false-position Theory
#### false-position Code
#### false-position Input
#### false-position Output

### newton-raphson-method
#### newton-raphson-Theory
#### newton-raphson-Code
#### newton-raphson-Input
#### newton-raphson-Output

### Secant Method
#### Secant Theory
#### Secant Code
#### Secant Input
#### Secant Output

### forward-interpolation-method
#### Forward Interpolation Theory
#### Forward Interpolation Code
#### Forward Interpolation Input
#### Forward Interpolation Output

### Backward-Interpolation 
#### Backward Interpolation Theory
#### Backward Interpolation Code
#### Backward Interpolation Input
#### Backward Interpolation Output

### divided-difference-method
#### Divided Difference Theory
#### Divided Difference Code
#### Divided Difference Input
#### Divided Difference Output

### Differentiation  Method
#### Differentiation Theory
#### Differentiation Code
#### Differentiation Input
#### Differentiation Output

### Simpson 1/3  Method
#### Simpson 1/3 Theory
#### Simpson 1/3 Code
#### Simpson 1/3 Input
#### Simpson 1/3 Output

### Simpson 3/8  Method
#### Simpson 3/8 Theory
#### Simpson 3/8 Code
#### Simpson 3/8 Input
#### Simpson 3/8 Output

### RK4  Method
#### RK4 Theory
#### RK4 Code
#### RK4 Input
#### RK4 Output

### Linear Regression  Method
#### Linear Regression Theory
#### Linear Regression Code
#### Linear Regression Input
#### Linear Regression Output

### Polynomial Regression  Method
#### Polynomial Regression Theory
#### Polynomial Regression Code
#### Polynomial Regression Input
#### Polynomial Regression Output

### Transcendental Regression  Method
#### Transcendental Regression Theory
#### Transcendental Regression Code
#### Transcendental Regression Input
#### Transcendental Regression Output




---
