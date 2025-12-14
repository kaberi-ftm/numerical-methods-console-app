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
