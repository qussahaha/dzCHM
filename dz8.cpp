#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <streambuf>
#include <string>
#ifdef _WIN32
#include <windows.h>   // обязательно ДО "using namespace std" (иначе конфликт std::byte)
#endif
using namespace std;

#ifdef _WIN32

// Вывод UTF-8 в консоль Windows через WriteConsoleW (не зависит от кодовой страницы)
class ConsoleBuf : public streambuf {
    string buf;
    HANDLE h;
    bool isConsole;
    void flushBuf() {
        if (buf.empty()) return;
        DWORD written;
        if (isConsole) {
            int n = MultiByteToWideChar(CP_UTF8, 0, buf.data(), (int)buf.size(), NULL, 0);
            wstring w(n, L'\0');
            MultiByteToWideChar(CP_UTF8, 0, buf.data(), (int)buf.size(), &w[0], n);
            WriteConsoleW(h, w.c_str(), (DWORD)w.size(), &written, NULL);
        } else {
            WriteFile(h, buf.data(), (DWORD)buf.size(), &written, NULL);  // вывод в файл/конвейер
        }
        buf.clear();
    }
public:
    ConsoleBuf() {
        h = GetStdHandle(STD_OUTPUT_HANDLE);
        DWORD mode;
        isConsole = GetConsoleMode(h, &mode) != 0;
    }
    ~ConsoleBuf() { flushBuf(); }
protected:
    int_type overflow(int_type c) override {
        if (!traits_type::eq_int_type(c, traits_type::eof())) {
            buf.push_back((char)c);
            if (c == '\n') flushBuf();
        }
        return traits_type::not_eof(c);
    }
    streamsize xsputn(const char* s, streamsize n) override {
        buf.append(s, (size_t)n);
        flushBuf();
        return n;
    }
    int sync() override { flushBuf(); return 0; }
};
#endif

typedef vector<double> Vec;
typedef vector<Vec> Mat;

// ---------- Общие функции ----------

// Метод Гаусса с выбором главного элемента
Vec gauss(Mat A, Vec b) {
    int n = b.size();
    for (int i = 0; i < n; i++) {
        int p = i;
        for (int r = i + 1; r < n; r++)
            if (fabs(A[r][i]) > fabs(A[p][i])) p = r;
        swap(A[i], A[p]); swap(b[i], b[p]);
        for (int r = i + 1; r < n; r++) {
            double q = A[r][i] / A[i][i];
            for (int c = i; c < n; c++) A[r][c] -= q * A[i][c];
            b[r] -= q * b[i];
        }
    }
    Vec x(n);
    for (int i = n - 1; i >= 0; i--) {
        double s = b[i];
        for (int c = i + 1; c < n; c++) s -= A[i][c] * x[c];
        x[i] = s / A[i][i];
    }
    return x;
}

// Значение многочлена по схеме Горнера
double poly(const Vec& a, double x) {
    double s = 0;
    for (int i = (int)a.size() - 1; i >= 0; i--) s = s * x + a[i];
    return s;
}

// ---------- Задание 1: кубический сплайн ----------

void task1() {
    cout << "===== Задание 1. Кубический сплайн =====\n";
    Vec x = {0.0, 0.9, 1.8, 2.7, 3.6};
    Vec f = {0.0, 0.36892, 0.85408, 1.7856, 6.3138};
    double X = 1.5;

    int n = x.size() - 1;
    Vec h(n);
    for (int i = 0; i < n; i++) h[i] = x[i + 1] - x[i];

    // Трёхдиагональная система для c_1..c_{n-1}; c_0 = c_n = 0
    int m = n - 1;
    Vec a(m), b(m), c(m), d(m);
    for (int k = 0; k < m; k++) {
        int i = k + 1;
        a[k] = (k > 0) ? h[i - 1] : 0;
        b[k] = 2 * (h[i - 1] + h[i]);
        c[k] = (k < m - 1) ? h[i] : 0;
        d[k] = 3 * ((f[i + 1] - f[i]) / h[i] - (f[i] - f[i - 1]) / h[i - 1]);
    }
    // Прогонка
    Vec P(m), Q(m);
    P[0] = -c[0] / b[0]; Q[0] = d[0] / b[0];
    for (int k = 1; k < m; k++) {
        double den = b[k] + a[k] * P[k - 1];
        P[k] = -c[k] / den;
        Q[k] = (d[k] - a[k] * Q[k - 1]) / den;
    }
    Vec cc(n + 1, 0.0);
    cc[m] = Q[m - 1];
    for (int k = m - 2; k >= 0; k--) cc[k + 1] = P[k] * cc[k + 2] + Q[k];

    Vec A(n), B(n), C(n), D(n);
    for (int i = 0; i < n; i++) {
        A[i] = f[i];
        C[i] = cc[i];
        D[i] = (cc[i + 1] - cc[i]) / (3 * h[i]);
        B[i] = (f[i + 1] - f[i]) / h[i] - h[i] * (cc[i + 1] + 2 * cc[i]) / 3;
    }

    cout << "S_i(x) = a + b(x-x_i) + c(x-x_i)^2 + d(x-x_i)^3\n";
    for (int i = 0; i < n; i++)
        cout << "[" << x[i] << ", " << x[i + 1] << "]: a=" << A[i] << " b=" << B[i]
             << " c=" << C[i] << " d=" << D[i] << "\n";

    int k = 0;
    while (k < n - 1 && X > x[k + 1]) k++;
    double t = X - x[k];
    double S = A[k] + B[k] * t + C[k] * t * t + D[k] * t * t * t;
    cout << "X* = " << X << " лежит на [" << x[k] << ", " << x[k + 1] << "]\n";
    cout << "S(X*) = " << S << "\n\n";
}

// ---------- Задание 2: МНК ----------

// Нормальная система МНК для многочлена степени deg
Vec fit(const Vec& x, const Vec& y, int deg) {
    int m = deg + 1, N = x.size();
    Mat A(m, Vec(m));
    Vec b(m);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            double s = 0;
            for (int k = 0; k < N; k++) s += pow(x[k], i + j);
            A[i][j] = s;
        }
        double s = 0;
        for (int k = 0; k < N; k++) s += y[k] * pow(x[k], i);
        b[i] = s;
    }
    return gauss(A, b);
}

void task2() {
    cout << "===== Задание 2. МНК =====\n";
    Vec x = {-0.9, 0.0, 0.9, 1.8, 2.7, 3.6};
    Vec y = {-0.36892, 0.0, 0.36892, 0.85408, 1.7856, 6.3138};

    Vec a1 = fit(x, y, 1), a2 = fit(x, y, 2);
    const Vec* res[2] = {&a1, &a2};

    for (int deg = 1; deg <= 2; deg++) {
        const Vec& a = *res[deg - 1];
        double Phi = 0;
        for (size_t k = 0; k < x.size(); k++) {
            double e = poly(a, x[k]) - y[k];
            Phi += e * e;
        }
        cout << "Степень " << deg << ": P(x) = ";
        for (int i = 0; i <= deg; i++)
            cout << (i ? " + (" : "") << a[i] << (i ? ")" : "") << (i ? "*x^" + to_string(i) : "");
        cout << "\n  Phi = " << Phi << "\n";
    }

    // Данные для графиков
    ofstream fp("points.dat"), fc("curves.dat");
    for (size_t k = 0; k < x.size(); k++) fp << x[k] << " " << y[k] << "\n";
    for (double t = -1.0; t <= 3.7; t += 0.05)
        fc << t << " " << poly(a1, t) << " " << poly(a2, t) << "\n";
    fp.close(); fc.close();

    ofstream g("plot.gp");
    g << "set terminal pngcairo size 900,600\nset output 'plot.png'\n"
         "set grid; set xlabel 'x'; set ylabel 'y'; set key left top\n"
         "plot 'points.dat' with points pt 7 title 'y_i', "
         "'curves.dat' using 1:2 with lines lw 2 title 'P1(x)', "
         "'curves.dat' using 1:3 with lines lw 2 title 'P2(x)'\n";
    g.close();

    // Попытка построить график автоматически (если установлен gnuplot)
#ifdef _WIN32
    const char* cmd = "gnuplot plot.gp >nul 2>nul";
#else
    const char* cmd = "gnuplot plot.gp >/dev/null 2>&1";
#endif
    if (system(cmd) == 0)
        cout << "График сохранён в plot.png\n";
    else
        cout << "gnuplot не найден: данные в points.dat, curves.dat, скрипт plot.gp\n";
}

int main() {
#ifdef _WIN32
    ConsoleBuf consoleBuf;
    streambuf* oldBuf = cout.rdbuf(&consoleBuf);
#endif
    cout << fixed << setprecision(6);
    task1();
    task2();
#ifdef _WIN32
    cout.flush();
    cout.rdbuf(oldBuf);
#endif
    return 0;
}
