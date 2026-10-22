#include <iostream>
#include <cmath>
#include <iomanip>
#include <clocale>

#ifdef _WIN32
#include <windows.h>
#endif

using namespace std;

// Настройка консоли для корректного вывода кириллицы
void setupRussianLocale() {
#ifdef _WIN32
    // Устанавливаем кодировку консоли UTF-8 для современных компиляторов Windows
    SetConsoleCP(65001);
    SetConsoleOutputCP(65001);
#endif
    setlocale(LC_ALL, "ru_RU.UTF-8");
}

// Первое уравнение системы: (x1^2 + a^2)*x2 - a^3 = 0
double f1(double x1, double x2, double a) {
    return (x1 * x1 + a * a) * x2 - a * a * a;
}

// Второе уравнение системы: (x1 - a/2)^2 + (x2 - a/2)^2 - a^2 = 0
double f2(double x1, double x2, double a) {
    double half_a = a / 2.0;
    return pow(x1 - half_a, 2) + pow(x2 - half_a, 2) - a * a;
}

int main() {
    setupRussianLocale();

    // Параметры для Варианта 3
    const double a = 4.0;
    const double eps = 1e-6; // Точность вычислений
    const int max_iter = 100; // Максимальное число итераций

    // Начальное приближение (определено графически для x1 > 0, x2 > 0)
    double x1 = 5.5;
    double x2 = 1.5;

    cout << "==========================================================" << endl;
    cout << " Решение системы нелинейных уравнений (Вариант 3, a = " << a << ")" << endl;
    cout << " Метод: Метод Ньютона" << endl;
    cout << " Начальное приближение: x1 = " << x1 << ", x2 = " << x2 << endl;
    cout << " Заданная точность: " << eps << endl;
    cout << "==========================================================" << endl << endl;

    cout << left << setw(10) << "Итерация" 
         << setw(16) << "x1" 
         << setw(16) << "x2" 
         << setw(16) << "Погрешность" << endl;
    cout << "----------------------------------------------------------" << endl;

    int iter = 0;
    double error = 0.0;

    do {
        iter++;

        // Значения функций в текущей точке
        double F1 = f1(x1, x2, a);
        double F2 = f2(x1, x2, a);

        // Элементы матрицы Якоби (частные производные)
        double J11 = 2.0 * x1 * x2;             // df1 / dx1
        double J12 = x1 * x1 + a * a;           // df1 / dx2
        double J21 = 2.0 * (x1 - a / 2.0);      // df2 / dx1
        double J22 = 2.0 * (x2 - a / 2.0);      // df2 / dx2

        // Определитель матрицы Якоби
        double det = J11 * J22 - J12 * J21;

        if (fabs(det) < 1e-12) {
            cout << "Ошибка: Определитель матрицы Якоби близок к нулю!" << endl;
            return 1;
        }

        // Находим приращения dx1 и dx2 по формулам Крамера для системы J * dx = -F
        double dx1 = (-F1 * J22 - (-F2) * J12) / det;
        double dx2 = (J11 * (-F2) - J21 * (-F1)) / det;

        // Обновляем значения неизвестных
        x1 += dx1;
        x2 += dx2;

        // Оценка текущей погрешности (максимальное по модулю приращение)
        error = max(fabs(dx1), fabs(dx2));

        cout << left << setw(10) << iter 
             << fixed << setprecision(7) << setw(16) << x1 
             << setw(16) << x2 
             << scientific << setprecision(3) << setw(16) << error << endl;

        if (iter >= max_iter) {
            cout << "\nПревышено максимальное число итераций!" << endl;
            break;
        }

    } while (error > eps);

    cout << "----------------------------------------------------------" << endl;
    cout << fixed << setprecision(6);
    cout << "\nРешение успешно найдено за " << iter << " итерации(й):" << endl;
    cout << "x1 = " << x1 << endl;
    cout << "x2 = " << x2 << endl;

    cout << "\nПроверка подстановкой в исходную систему (невязка):" << endl;
    cout << scientific << setprecision(3);
    cout << "f1(x1, x2) = " << f1(x1, x2, a) << " (должно быть ~ 0)" << endl;
    cout << "f2(x1, x2) = " << f2(x1, x2, a) << " (должно быть ~ 0)" << endl;

    return 0;
}