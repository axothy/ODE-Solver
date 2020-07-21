#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>

using namespace std;

double Func(double F, double m, double k, double v) {
    return (F - k * (v * v)) / m;
}

int main()
{
    setlocale(LC_ALL, "rus");

    //Объявлем переменные 
    int i, methodNumber;
    const int n = 100;
    double F, m, k, v0;
    double h;
    double K1, K2, K3, K4;

    //Откроем файл
    ifstream in("Data.txt");

    // Если мы не можем открыть этот файл для чтения его содержимого
    if (!in)
    {
        // То выводим следующее сообщение об ошибке и выполняем exit()
        cerr << "Невозможно открыть файл Data.txt..." << endl;
        exit(1);
    }

    // Прочитаем переменные из открытого файла
    in >> F >> m >> k >> v0 >> h >> methodNumber;
    in.close(); // Закрываем файл

    std::cout << " Cила F:                " << F << "\n";
    std::cout << " Масса m:               " << m << std::endl;
    std::cout << " Коэффициент k:         " << k << std::endl;
    std::cout << " Начальная скорость v0: " << v0 << std::endl;
    std::cout << " Шаг h:                 " << h << std::endl;
    std::cout << " Итерации n:            " << n << std::endl;
    std::cout << " Метод решения:         " << methodNumber << std::endl;

    std::cout << "" << std::endl;
    std::cout << " ================================================" << std::endl;
    std::cout << "" << std::endl;

    //Объявляем массивы

    double v[n + 1];
    double v_analytical[n];
    double t[n + 1];
    double percent[n];

    v[1] = 0;
    t[1] = 0;
    v_analytical[1] = 0;

    for (i = 2;i <= n;i++) {
        t[i] = t[i - 1] + h;
        v_analytical[i] = sqrt(F / k) * abs(tanh(((1 / 2) * log(abs(((1 + (sqrt(k / F) * v0))) / (1 - (sqrt(k / F) * v0))))) - (sqrt(k * F) / m) * t[i]));
    }

    //Конструкция ветвления относительно переменной method
    if (methodNumber == 1)
    {
        std::cout << "Решаем методом Эйлера..." << std::endl;
        for (i = 2; i <= n; i++)
        {
            t[i] = t[i - 1] + h;
            v[i] = v[i - 1] + h * Func(F, m, k, v[i - 1]);
        }
    }

    if (methodNumber == 2)
    {
        std::cout << "Решаем методом Рунге-Кутты 2 порядка..." << std::endl;
        for (i = 2;i <= n;i++) {
            K1 = h * Func(F, m, k, v[i - 1]);
            K2 = h * Func(F, m, k, (v[i - 1] + K1 / 2.0));
            t[i] = t[i - 1] + h;
            v[i] = v[i - 1] + K2;
        }
    }

    if (methodNumber == 3)
    {
        std::cout << "Решаем методом Рунге-Кутты 4 порядка..." << std::endl;
        for (i = 2;i <= n;i++) {
            K1 = h * Func(F, m, k, v[i - 1]);
            K2 = h * Func(F, m, k, (v[i - 1] + K1 / 2.0));
            K3 = h * Func(F, m, k, (v[i - 1] + K2 / 2.0));
            K4 = h * Func(F, m, k, (v[i - 1] + K3));
            t[i] = t[i - 1] + h;
            v[i] = v[i - 1] + (K1 + 2 * K2 + 2 * K3 + K4) / 6.0;
        }
    }

    if (methodNumber == 4)
    {
        std::cout << "Решаем аналитически..." << std::endl;
        for (i = 2;i <= n;i++) {
            v[i] = v_analytical[i];
        }
    }

    //Открываем файл для записи решения

    ofstream fout("Solution.txt");
    fout << "Решение методом: " << methodNumber << endl;
    for (int i = 1; i <= n; i++) {
        fout << "t=" << t[i] << "     " << "v=" << v[i] << endl;
    }
    std::cout << "Файл Solution.txt был успешно создан." << std::endl;
    fout.close(); // Закрываем файл

    for (i = 1; i <= n; i++)
    {
        percent[i] = (100*abs(v_analytical[i] - v[i]))/v[i];
    }

    //Открываем файл для записи погрешностей

    ofstream crout("Calculation Errors.txt");
    crout << "Подсчет отклонения в процентах" << endl;
    for (int i = 1; i <= n; i++) {
        crout << percent[i] << "%" << endl;
    }
    std::cout << "Файл Calculcation Errors.txt был успешно создан." << std::endl;
    crout.close(); // Закрываем файл

    system("pause");
    return 0;
}
