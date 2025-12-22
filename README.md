
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

[Solution of NonLinear Equations  :](#solution-of-nonlinear-equations)
----

| Algorithm | 📄 Concept | 💻 Source | 📥 Test Data | 📤 Result |
| :--- | :---: | :---: | :---: | :---: |
| [**Bisection Method**](#bisection-method) | [Theory](#bisection-theory) | [Code](#bisection-code) | [Input](#bisection-input) | [Output](#bisection-output) |
| [**False Position Method**](#false-position-method) | [Theory](#false-position-theory) | [Code](#false-position-code) | [Input](#false-position-input) | [Output](#false-position-output) |
| [**Newton Raphson Method**](#newton-raphson-method) | [Theory](#newton-raphson-theory) | [Code](#newton-raphson-code) | [Input](#newton-raphson-input) | [Output](#newton-raphson-output) |
| [**Secant Method**](#secant-method) | [Theory](#secant-theory) | [Code](#secant-code) | [Input](#secant-input) | [Output](#secant-output) |

---

[Solution of Interpolations  :](#solution-of-interpolations)
----

| Algorithm | 📄 Concept | 💻 Source | 📥 Test Data | 📤 Result |
| :--- | :---: | :---: | :---: | :---: |
| [**Forward Interpolation**](#forward-interpolation-method) | [Theory](#forward-interpolation-theory) | [Code](#forward-interpolation-code) | [Input](#forward-interpolation-input) | [Output](#forward-interpolation-output) |
| [**Backward Interpolation**](#backward-interpolation-method) | [Theory](#backward-interpolation-theory) | [Code](#backward-interpolation-code) | [Input](#backward-interpolation-input) | [Output](#backward-interpolation-output) |
| [**Divided Difference Method**](#divided-difference-method) | [Theory](#divided-difference-theory) | [Code](#divided-difference-code) | [Input](#divided-difference-input) | [Output](#divided-difference-output) |

---

[Solution of Differentiation :](#solution-of-differentiation)
----

| Algorithm | 📄 Concept | 💻 Source | 📥 Test Data | 📤 Result |
| :--- | :---: | :---: | :---: | :---: |
| [**Differentiation Using Interpolation**](#differentiation-using-interpolation-method) | [Theory](#differentiation-theory) | [Code](#differentiation-code) | [Input](#differentiation-input) | [Output](#differentiation-output) |

---

[Solution of Integration  :](#solution-of-integration)
----

| Algorithm | 📄 Concept | 💻 Source | 📥 Test Data | 📤 Result |
| :--- | :---: | :---: | :---: | :---: |
| [**Simpson’s 1/3 Rule**](#simpsons-13-rule-method) | [Theory](#simpson-13-theory) | [Code](#simpson-13-code) | [Input](#simpson-13-input) | [Output](#simpson-13-output) |
| [**Simpson’s 3/8 Rule**](#simpsons-38-rule-method) | [Theory](#simpson-38-theory) | [Code](#simpson-38-code) | [Input](#simpson-38-input) | [Output](#simpson-38-output) |

---

[Solution of Differential Equations  :](#solution-of-differential-equations)
----

| Algorithm | 📄 Concept | 💻 Source | 📥 Test Data | 📤 Result |
| :--- | :---: | :---: | :---: | :---: |
| [**Runge Kutta 4th Order Method**](#runge-kutta-4th-order-method) | [Theory](#rk4-theory) | [Code](#rk4-code) | [Input](#rk4-input) | [Output](#rk4-output) |

---

[Solution of Regression  :](#solution-of-regression)
----

| Algorithm | 📄 Concept | 💻 Source | 📥 Test Data | 📤 Result |
| :--- | :---: | :---: | :---: | :---: |
| [**Linear Regression**](#linear-regression-method) | [Theory](#linear-regression-theory) | [Code](#linear-regression-code) | [Input](#linear-regression-input) | [Output](#linear-regression-output) |
| [**Polynomial Regression**](#polynomial-regression-method) | [Theory](#polynomial-regression-theory) | [Code](#polynomial-regression-code) | [Input](#polynomial-regression-input) | [Output](#polynomial-regression-output) |
| [**Transcendental Regression**](#transcendental-regression-method) | [Theory](#transcendental-regression-theory) | [Code](#transcendental-regression-code) | [Input](#transcendental-regression-input) | [Output](#transcendental-regression-output) |


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

### Gauss Jordan Elimination Method
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
```cpp
#include <bits/stdc++.h>
using namespace std;
#define MAX 100
 ofstream fout;
void getCofactor(double mat[MAX][MAX], double temp[MAX][MAX],int p, int q, int n)
{
    int i = 0, j = 0;
    for (int row = 0; row < n; row++)
        {
        for (int col = 0; col< n; col++)
        {
            if (row != p &&col!= q) {
                temp[i][j++] = mat[row][col];

                if (j == n-1)
                    {
                    j=0;
                    i++;
                }
            }
        }
    }
}

double determinant(double mat[MAX][MAX], int n)
{
    if (n==1)
       {
        return mat[0][0];
       }

    double det = 0;
    double temp[MAX][MAX];
    int sign = 1;

    for (int f = 0; f < n; f++) {
        getCofactor(mat, temp, 0, f, n);
        det += sign * mat[0][f] * determinant(temp, n - 1);
        sign = -sign;
    }

    return det;
}

void adjoint(double mat[MAX][MAX], double adj[MAX][MAX], int n)
{
    if (n == 1) {
        adj[0][0] = 1;
        return;
    }

    int sign = 1;
    double temp[MAX][MAX];

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            getCofactor(mat, temp, i, j, n);
            sign = ((i + j) % 2 == 0) ? 1 : -1;
            adj[j][i] = sign * determinant(temp, n - 1);
        }
    }
}

bool inverse(double mat[MAX][MAX], double inv[MAX][MAX], int n)
{
    double det = determinant(mat, n);

    if (det == 0) {
        fout << "Inverse does not exist as determinant(mat)= 0."<<endl;
        return false;
    }

    double adj[MAX][MAX];
    adjoint(mat, adj, n);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inv[i][j] = adj[i][j] / det;

    return true;
}

int main()
{
     fstream fin("MatrixIn.txt");
     ofstream fout("MatrixOut.txt");
   if(!fin)
   {
       cout<<"Can't open input file"<<endl;
       return 1;
   }
    if(!fout)
   {
       cout<<"Can't open output file"<<endl;
       return 1;
   }
    int n;
    double mat[MAX][MAX], B[MAX],inv[MAX][MAX],X[MAX];
    fin >> n;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            fin >> mat[i][j];

        for (int i = 0; i < n; i++)
        fin >> B[i];

    if (inverse(mat, inv, n)) {
        fout << "Inverse Matrix:"<<endl;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                fout << fixed << setprecision(3) << inv[i][j] << " ";
            fout << endl;
        }
    }
    for (int i = 0; i < n; i++) {
        X[i] = 0;
        for (int j = 0; j < n; j++)
            X[i] += inv[i][j] * B[j];
    }
    fout << "Solution Vector X:\n";
    for (int i = 0; i < n; i++)
        fout << "x" << i + 1 << " = "
             << fixed << setprecision(3) << X[i] << "\n";
    fin.close();
    fout.close();

    return 0;
}
```
#### Inverse Matrix Input
```
3
2 1 -1
-3 -1 2
-2 1 2
8 -11 -3
```
#### Inverse Matrix Output
```
Inverse Matrix:
4.000 3.000 -1.000 
-2.000 -2.000 1.000 
5.000 4.000 -1.000 
Solution Vector X:
x1 = 2.000
x2 = 3.000
x3 = -1.000
```
## Solution of NonLinear Equations
### Bisection Method
#### Bisection Theory
#### Bisection Code

```cpp
#include <bits/stdc++.h>
using namespace std;

double f(double x){
return x*x*x-x*x+2 ;
}

int main() {
     fstream fin("BisecIn.txt");
     ofstream fout("BisecOut.txt");
   if(!fin)
   {
       cout<<"Can't open input file"<<endl;
       return 1;
   }
    if(!fout)
   {
       cout<<"Can't open output file"<<endl;
       return 1;
   }

    double x1,x2,e;
    fin>>x1>>x2>>e;
    if(f(x1)*f(x2)>=0)
    {
        fout<<"No root in the given interval"<<endl;
        return 1;

    }
        double x0;
    while (true){
    x0=(x1+x2)/2;
    double f1=f(x1);
    double f2=f(x2);
    double f0=f(x0);
   if( f0== 0 || abs((x2-x1)/x2 )<e)
   {
       fout<<"Root : "<<x0<<endl;
       break;
   }
   if(f1*f0 <0)
   {
       x2=x0;
   }
   else
   {
       x1=x0;
   }
    }

  fin.close();
  fout.close();
    return 0;
}
```
#### Bisection Input
```
-2 2 0.001
```
#### Bisection Output

```
Root : -1
```

### False Position Method
#### False Position Theory
#### False Position Code

```cpp
#include <bits/stdc++.h>
using namespace std;

double f(double x){
return x*x-9 ;
}

int main() {
     fstream fin("FalsePIn.txt");
     ofstream fout("FalsePOut.txt");
   if(!fin)
   {
       cout<<"Can't open input file"<<endl;
       return 1;
   }
    if(!fout)
   {
       cout<<"Can't open output file"<<endl;
       return 1;
   }

    double x1,x2,e;
    fin>>x1>>x2>>e;
    if(f(x1)*f(x2)>=0)
    {
        fout<<"No root in the given interval"<<endl;
        return 1;

    }
        double x0;
        int iter=0;
    while (iter<1000){
            iter++;
    double f1=f(x1);
    double f2=f(x2);
    if(f1==f2)
    {
        fout<<"Division by zero error in False Position method ."<<endl;
        break;
    }
    x0=x2- f2*(x1-x2)/(f1-f2);
    double f0=f(x0);

   if( abs(f0)< e || abs(x2-x1) < e)
   {
       fout<<"Root : "<<x0<<endl;
       break;
   }
   if(f1*f0 <0)
   {
       x2=x0;

   }
   else

   {
       x1=x0;

   }
   if(iter>1000)
   {
       fout<<"Maximum iteration exceeded ."<<endl;
       break;
   }


    }
    fin.close();
  fout.close();
    return 0;
}
```
#### False Position Input
```
0 5 0.001
```
#### False Position Output
```
Root : 2.99991
```

### Newton Raphson Method
#### Newton Raphson Theory
#### Newton Raphson Code
#### Newton Raphson Input
#### Newton Raphson Output

### Secant Method
#### Secant Theory
#### Secant Code
#### Secant Input
#### Secant Output

## Solution of Interpolations 
### Forward Interpolation Method
#### Forward Interpolation Theory
#### Forward Interpolation Code
#### Forward Interpolation Input
#### Forward Interpolation Output

### Backward Interpolation Method 
#### Backward Interpolation Theory
#### Backward Interpolation Code
#### Backward Interpolation Input
#### Backward Interpolation Output

### Divided Difference Method
#### Divided Difference Theory
#### Divided Difference Code
#### Divided Difference Input
#### Divided Difference Output
## Solution of Differentiation
### Differentiation Using Interpolation Method
#### Differentiation Theory
#### Differentiation Code
```cpp
#include <bits/stdc++.h>
using namespace std;

int main() {
    ifstream fin("Differentiation_Input.txt");
    ofstream fout("Differentiation_Output.txt");

    if (!fin || !fout) {
        cout << "File can't be opened." << endl;
        return 1;
    }

    fout << fixed << setprecision(6);

    int n;
    fin >> n;

    double x[50], y[50], X;
    for (int i = 0; i < n; i++) fin >> x[i];
    for (int i = 0; i < n; i++) fin >> y[i];
    fin >> X;

    double h = x[1] - x[0];

    double fd[50][50];
    for (int i = 0; i < n; i++) fd[i][0] = y[i];
    for (int j = 1; j < n; j++)
        for (int i = 0; i < n - j; i++)
            fd[i][j] = fd[i + 1][j - 1] - fd[i][j - 1];

    fout << "Forward Difference Table:\n";
    for (int i = 0; i < n; i++) {
        fout << setw(8) << x[i];
        for (int j = 0; j < n - i; j++)
            fout << setw(12) << fd[i][j];
        fout << endl;
    }

    int idx = -1;
    for (int i = 0; i < n; i++)
        if (abs(x[i] - X) < 1e-6) {
            idx = i;
            break;
        }

    double der1, der2;

    if (idx == 0) {
        der1 = (-3*y[0] + 4*y[1] - y[2]) / (2*h);
        der2 = (y[0] - 2*y[1] + y[2]) / (h*h);
    } else if (idx == n - 1) {
        der1 = (3*y[n-1] - 4*y[n-2] + y[n-3]) / (2*h);
        der2 = (y[n-1] - 2*y[n-2] + y[n-3]) / (h*h);
    } else {
        der1 = (y[idx+1] - y[idx-1]) / (2*h);
        der2 = (y[idx+1] - 2*y[idx] + y[idx-1]) / (h*h);
    }

    fout << "\nFirst derivative at X = " << X << " is: " << der1 << endl;
    fout << "Second derivative at X = " << X << " is: " << der2 << endl;

    return 0;
}
```
#### Differentiation Input
```
5
1 2 3 4 5
1 8 27 64 125
3
```
#### Differentiation Output
```
Forward Difference Table:
1.000000    1.000000    7.000000   12.000000    6.000000    0.000000
2.000000    8.000000   19.000000   18.000000    6.000000
3.000000   27.000000   37.000000   24.000000
4.000000   64.000000   61.000000
5.000000  125.000000

First derivative at X = 3.000000 is: 28.000000
Second derivative at X = 3.000000 is: 18.000000
```
## Solution of Integration
### Simpson’s 1/3 Rule Method
#### Simpson 1/3 Theory
#### Simpson 1/3 Code
#### Simpson 1/3 Input
#### Simpson 1/3 Output

### Simpson’s 3/8 Rule Method
#### Simpson 3/8 Theory
#### Simpson 3/8 Code
#### Simpson 3/8 Input
#### Simpson 3/8 Output
## Solution of Differential Equations
### Runge Kutta 4th Order Method
#### RK4 Theory
#### RK4 Code
```cpp
#include <bits/stdc++.h>
using namespace std;

  double f(double x, double y)
  {
      return ((x - y)/2);
  }

int main(){
    fstream fin("RKinput.txt");
     ofstream fout("RKoutput.txt");
   if(!fin)
   {
       cout<<"Can't open input file"<<endl;
       return 1;
   }
    if(!fout)
   {
       cout<<"Can't open output file"<<endl;
       return 1;
   }

 double h,xn,x0,y0;
    fin>>h>>xn>>x0>>y0;

 //"Enter step size h, final xn, initial x0, initial y0:\n";
 cout << "Read values: h=" << h << " xn=" << xn << " x0=" << x0 << " y0=" << y0 << endl;

  double x=x0;
  double y=y0;
  int s=(xn-x0)/h;
    for( int i=0;i<s ;i++)  {
      double k1= h*f(x,y);
      double k2= h*f(x+h/2.0,y+k1/2.0);
      double k3= h*f(x+h/2.0,y+k2/2.0);
      double k4= h*f(x+h,y+k3);
      y=y+(k1+2*k2+2*k3+k4)/6.0;
      x=x+h;

  }
  fout<<"Value of y at x is: "<<y<<endl;
  fin.close();
  fout.close();

return 0;

}
```
#### RK4 Input
```
0.2 2 0 1
```
#### RK4 Output
```
Value of y at x is: 1.10364
```

## Solution of Regression
### Linear Regression Method
#### Linear Regression Theory
#### Linear Regression Code
```cpp
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
```
#### Linear Regression Input
```
7
1 3
2 4
3 4
4 5
5 8
6 9
7 10
```
#### Linear Regression Output
```
Linear Regression Equation: 
y = 1.143 + 1.250x
```

### Polynomial Regression Method
#### Polynomial Regression Theory
#### Polynomial Regression Code
```cpp
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
```
#### Polynomial Regression Input
```
6
1 4
2 8
3 14
4 22
5 32
6 44
```
#### Polynomial Regression Output
```
Polynomial Regression Equation:
y = 2.000 + 1.000x + 1.000x^2
```

### Transcendental Regression Method
#### Transcendental Regression Theory
#### Transcendental Regression Code
```cpp
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
```
#### Transcendental Regression Input
```
6
1 3.0
2 4.9
3 8.1
4 13.5
5 22.3
6 36.6
```
#### Transcendental Regression Output
```
Exponential Regression Equation:
y = 1.807 * e^(0.502x)
```



---
