#ifdef _WIN32
#include <windows.h>
#endif
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

const double PI = acos(-1.0);

double f(double x) { return cos(x) / sin(x); }                         // ctg x
double f4(double x) {                                                  // f''''(x)
    double s = sin(x), c = cos(x);
    return 8.0 * c * (2.0 + c * c) / pow(s, 5);
}

// Интерполяционный многочлен Лагранжа
double lagrange(const vector<double>& X, const vector<double>& Y, double x) {
    int n = X.size();
    double res = 0;
    for (int i = 0; i < n; i++) {
        double li = 1;
        for (int j = 0; j < n; j++)
            if (j != i) li *= (x - X[j]) / (X[i] - X[j]);
        res += Y[i] * li;
    }
    return res;
}

// Разделённые разности для многочлена Ньютона
vector<double> dividedDiff(const vector<double>& X, const vector<double>& Y) {
    int n = X.size();
    vector<double> a = Y;
    for (int k = 1; k < n; k++)
        for (int i = n - 1; i >= k; i--)
            a[i] = (a[i] - a[i - 1]) / (X[i] - X[i - k]);
    return a;   // a[k] = f(x0; ...; xk)
}

double newton(const vector<double>& X, const vector<double>& a, double x) {
    int n = X.size();
    double res = a[n - 1];
    for (int i = n - 2; i >= 0; i--)
        res = res * (x - X[i]) + a[i];      // схема Горнера
    return res;
}

void solve(const string& name, const vector<double>& X, double xs) {
    int n = X.size();
    vector<double> Y(n);
    for (int i = 0; i < n; i++) Y[i] = f(X[i]);

    vector<double> a = dividedDiff(X, Y);
    double L = lagrange(X, Y, xs);
    double N = newton(X, a, xs);
    double exact = f(xs);

    // Практическая погрешность
    double errL = fabs(exact - L);
    double errN = fabs(exact - N);

    // Теоретическая оценка (3.9): M_{n+1}/(n+1)! * |omega(x*)|
    double omega = 1;
    for (int i = 0; i < n; i++) omega *= (xs - X[i]);
    double M = 0;
    int steps = 100000;
    double lo = X.front(), hi = X.back();
    for (int i = 0; i <= steps; i++) {
        double x = lo + (hi - lo) * i / steps;
        M = max(M, fabs(f4(x)));
    }
    double est = M / 24.0 * fabs(omega);   // (n+1)! = 4! = 24

    cout << "=== Набор " << name << " ===\n";
    cout << "Узлы X_i:  ";
    for (double x : X) cout << x << "  ";
    cout << "\nY_i:       ";
    for (double y : Y) cout << y << "  ";
    cout << "\n\nРазделённые разности: ";
    for (double v : a) cout << v << "  ";
    cout << "\n\nX* = " << xs << "\n";
    cout << "Точное значение ctg(X*)  = " << exact << "\n";
    cout << "Лагранж   L3(X*)         = " << L << ",  |погр.| = " << errL << "\n";
    cout << "Ньютон    N3(X*)         = " << N << ",  |погр.| = " << errN << "\n";
    cout << "M4 = max|f''''|          = " << M << "\n";
    cout << "|omega(X*)|              = " << fabs(omega) << "\n";
    cout << "Оценка (3.9)             = " << est << "\n";
    cout << "Оценка >= факт. погр.?   : " << (est >= errL ? "да" : "нет") << "\n\n";
}

int main() {
    #ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
#endif
    cout << fixed << setprecision(8);
    cout << fixed << setprecision(8);
    double xs = PI / 3;

    vector<double> Xa = {PI / 8, 2 * PI / 8, 3 * PI / 8, 4 * PI / 8};
    vector<double> Xb = {PI / 8, 5 * PI / 16, 3 * PI / 8, PI / 2};

    solve("а", Xa, xs);
    solve("б", Xb, xs);
    return 0;
}