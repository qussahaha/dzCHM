#include <iostream>
#include <cmath>
#include <iomanip>
#ifdef _WIN32
#include <windows.h>
#endif

using namespace std;

// Исходная функция: f(x) = sqrt(1 - x^2) - e^x + 0.1
double f(double x) {
    return sqrt(1.0 - x * x) - exp(x) + 0.1;
}

// Первая производная: f'(x) = -x / sqrt(1 - x^2) - e^x
double df(double x) {
    return -x / sqrt(1.0 - x * x) - exp(x);
}

// Вторая производная: f''(x) = -1 / (1 - x^2)^(3/2) - e^x
double d2f(double x) {
    double term = 1.0 - x * x;
    return -1.0 / (term * sqrt(term)) - exp(x);
}

// Функция для метода простой итерации:
// Из уравнения e^x = sqrt(1 - x^2) + 0.1 выражаем x:
// x = phi(x) = ln(sqrt(1 - x^2) + 0.1)
double phi(double x) {
    return log(sqrt(1.0 - x * x) + 0.1);
}

// Производная phi'(x) для проверки условия сходимости |phi'(x)| < 1
double dphi(double x) {
    double root_term = sqrt(1.0 - x * x);
    return (-x / root_term) / (root_term + 0.1);
}

// 1. Метод простой итерации
double simpleIterationMethod(double x0, double eps, int& iterations) {
    double x_prev = x0;
    double x_next = phi(x_prev);
    iterations = 1;

    cout << "\n--- Метод простой итерации ---\n";
    cout << left << setw(10) << "Iter"
         << setw(18) << "x_n"
         << setw(18) << "|x_n - x_{n-1}|"
         << setw(18) << "|f(x_n)|" << "\n";
    cout << string(64, '-') << "\n";

    cout << left << setw(10) << iterations
         << fixed << setprecision(10)
         << setw(18) << x_next
         << scientific << setprecision(4)
         << setw(18) << fabs(x_next - x_prev)
         << setw(18) << fabs(f(x_next)) << "\n";

    while (fabs(x_next - x_prev) > eps) {
        x_prev = x_next;
        x_next = phi(x_prev);
        iterations++;

        cout << left << setw(10) << iterations
             << fixed << setprecision(10)
             << setw(18) << x_next
             << scientific << setprecision(4)
             << setw(18) << fabs(x_next - x_prev)
             << setw(18) << fabs(f(x_next)) << "\n";
    }

    return x_next;
}

// 2. Метод Ньютона (метод касательных)
double newtonMethod(double x0, double eps, int& iterations) {
    double x_prev = x0;
    double x_next = x_prev - f(x_prev) / df(x_prev);
    iterations = 1;

    cout << "\n--- Метод Ньютона (касательных) ---\n";
    cout << left << setw(10) << "Iter"
         << setw(18) << "x_n"
         << setw(18) << "|x_n - x_{n-1}|"
         << setw(18) << "|f(x_n)|" << "\n";
    cout << string(64, '-') << "\n";

    cout << left << setw(10) << iterations
         << fixed << setprecision(10)
         << setw(18) << x_next
         << scientific << setprecision(4)
         << setw(18) << fabs(x_next - x_prev)
         << setw(18) << fabs(f(x_next)) << "\n";

    while (fabs(x_next - x_prev) > eps) {
        x_prev = x_next;
        x_next = x_prev - f(x_prev) / df(x_prev);
        iterations++;

        cout << left << setw(10) << iterations
             << fixed << setprecision(10)
             << setw(18) << x_next
             << scientific << setprecision(4)
             << setw(18) << fabs(x_next - x_prev)
             << setw(18) << fabs(f(x_next)) << "\n";
    }

    return x_next;
}

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
#endif

    // Заданная точность вычислений
    const double eps = 1e-6;

    // Графически определённый отрезок локализации положительного корня: [0.0, 0.2]
    // Начальное приближение x0 = 0.1
    const double x0 = 0.1;

    cout << "Решение нелинейного уравнения (Вариант 3):\n";
    cout << "sqrt(1 - x^2) - e^x + 0.1 = 0\n\n";
    cout << fixed << setprecision(6);
    cout << "Начальное приближение x0 = " << x0 << "\n";
    cout << "Заданная точность eps    = " << scientific << setprecision(1) << eps << "\n";
    cout << fixed << setprecision(6);
    cout << "Проверка условий в точке x0:\n";
    cout << "  |phi'(x0)|      = " << fabs(dphi(x0)) << " < 1 (условие сходимости МПИ выполнено)\n";
    cout << "  f(x0) * f''(x0) = " << f(x0) * d2f(x0) << " > 0 (условие Фурье для метода Ньютона выполнено)\n";

    int iter_mpi = 0;
    double root_mpi = simpleIterationMethod(x0, eps, iter_mpi);

    int iter_newton = 0;
    double root_newton = newtonMethod(x0, eps, iter_newton);

    cout << "\n================ ИТОГОВЫЕ РЕЗУЛЬТАТЫ ================\n";
    cout << fixed << setprecision(8);
    cout << "1. Метод простой итерации:\n";
    cout << "   Корень x = " << root_mpi << " (найден за " << iter_mpi << " итер., |f(x)| = "
         << scientific << setprecision(2) << fabs(f(root_mpi)) << ")\n";

    cout << fixed << setprecision(8);
    cout << "2. Метод Ньютона:\n";
    cout << "   Корень x = " << root_newton << " (найден за " << iter_newton << " итер., |f(x)| = "
         << scientific << setprecision(2) << fabs(f(root_newton)) << ")\n";

    return 0;
}