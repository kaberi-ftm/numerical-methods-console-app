#include <bits/stdc++.h>
using namespace std;

int main() {
    ifstream fin("DiffInput.txt");
    ofstream fout("DiffOutput.txt");

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
