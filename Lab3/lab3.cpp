// g++ ./lab-3.cpp -fopenmp=libomp && a.exe

//------------------------------------------------------------
// Программа решения уравнений Пуассона методом Гаусса-Зейделя
//------------------------------------------------------------
#include <iostream>
#include <cstring>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include "locale.h"

using namespace std;

// Функции решения уравнения (результаты всех версий должны быть идентичны!)
int Calc_ser(double **u, double **f, int N, double eps); // последовательная
int Calc_blk(double **u, double **f, int N, double eps); // блочная последовательная
int Calc_par(double **u, double **f, int N, double eps); // параллельная (блочная)

// Инициализация массивов
void Init(double **u, double **f, int N);
void clear_arr(double** arr, int N);
double **new_arr(int N);
void delete_arr(double **arr, int N);

// Вывод части массива для контроля
void Output(double **u, int N);

int main(int argc, char **argv)
{
    double **u = NULL, **f = NULL;

    const int N = 600;       // Количество точек сетки по каждой размерности
    const double eps = 0.0001; // Точность вычислений
    int icnt;                // Количество итераций
    double stime = -1;       // Время решения

    f = new_arr(N);     // Выделение памяти под правую часть значений уравнения
    u = new_arr(N + 2); // Выделение памяти под неизвестные и краевые условия

    //	 Последовательная реализация

    cout << "\n\t*** Serial version ***\n";
    Init(u, f, N);                 // Инициализация краевых условий и правой части уравнения


    stime = omp_get_wtime();
    icnt = Calc_ser(u, f, N, eps); // Вызов функции расчета по методу Гаусса-Зейделя
    cout << "Solution time = " << omp_get_wtime() - stime << endl;
    cout << "Iterations =    " << icnt << endl;
    cout << "Results:\n";
    Output(u, N); // Вывод результатов на экран

    //	 Последовательная блочная реализация
    cout << "\n\t*** Block serial version ***\n";
    Init(u, f, N);                 // Инициализация краевых условий и правой части

    stime = omp_get_wtime();
    icnt = Calc_blk(u, f, N, eps); // Вызов блочной функции расчета
    cout << "Solution time = " << omp_get_wtime() - stime << endl;
    cout << "Iterations =    " << icnt << endl;
    cout << "Results:\n";
    Output(u, N);

    //   Параллельная реализация (блочная)
    cout << "\n\t*** Parallel version ***\n";
    Init(u, f, N);                 // Инициализация краевых условий и правой части

    stime = omp_get_wtime();
    icnt = Calc_par(u, f, N, eps); // Вызов параллельной функции расчета
    cout << "Solution time = " << omp_get_wtime() - stime << endl;
    cout << "Iterations =    " << icnt << endl;
    cout << "Results:\n";
    Output(u, N);

    // Освобождение памяти массивов
    delete_arr(f, N);
    delete_arr(u, N + 2);

    return 0;
}

// Последовательная функция, реализующая алгоритм Гаусса-Зейделя
// Входные параметры: массив неизвестных и краевых значений, массив правых частей, количество точек сетки по каждому направлению, точность вычислений
int Calc_ser(double **u, double **f, int N, double eps)
{
    double max;               // Максимальная ошибка на итерации
    double h = 1.0 / (N + 1); // Величина шага
    int icnt = 0;             // Количество итераций

    do
    {
        icnt++;
        max = 0;
        for (int i = 1; i <= N; i++)
            for (int j = 1; j <= N; j++)
            {
                double u0 = u[i][j];
                u[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j] // Расчет по формуле Гаусса-Зейделя
                                  + u[i][j - 1] + u[i][j + 1] - h * h * f[i - 1][j - 1]);
                double d = fabs(u[i][j] - u0); // Разность нового значения неизвестной и значения с предыдущей итерации
                if (d > max)                   // Поиск максимальной ошибки
                    max = d;
            }
    } while (max > eps);

    return icnt;
}

// Последовательная функция, реализующая блочный алгоритм Гаусса-Зейделя
int Calc_blk(double **u, double **f, int N, double eps)
{
    double max;
    double h = 1.0 / (N + 1);
    int icnt = 0;

    const int BlockSize = 20; // Размер блока
    int bcnt;                 // Количество блоков в ряд

    if (N % BlockSize == 0) // Если количество точек по каждому из направлений сетки делится нацело на размер блока, то проводятся вычисления
    {
        bcnt = N / BlockSize;
        do
        {
            icnt++;
            max = 0;

            int bi, bj, i, j, ii, jj;
            double temp, diff;

            for (ii = 0; ii < bcnt; ii++) {
                for (jj = 0; jj <= ii; jj++) {
                    bi = ii - jj;
                    bj = jj;

                    // Циклы по элементам внутри блока
                    for (i = bi * BlockSize + 1; i <= (bi + 1) * BlockSize; i++) {
                        for (j = bj * BlockSize + 1; j <= (bj + 1) * BlockSize; j++) {
                            temp = u[i][j];
                            u[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1] - h * h * f[i - 1][j - 1]);
                            diff = fabs(u[i][j] - temp);

                            if (diff > max)
                                max = diff;
                        }
                    }
                }
            }

            for (ii = bcnt-2; ii >= 0; ii--) {
                for (jj = 0; jj <= ii; jj++) {
                    bi = jj - ii - 1 + bcnt;
                    bj = bcnt - 1 - jj;

                    // Циклы по элементам внутри блока
                    for (i = bi * BlockSize + 1; i <= (bi + 1) * BlockSize; i++) {
                        for (j = bj * BlockSize + 1; j <= (bj + 1) * BlockSize; j++) {
                            temp = u[i][j];
                            u[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1] - h * h * f[i - 1][j - 1]);
                            diff = fabs(u[i][j] - temp);

                            if (diff > max)
                                max = diff;
                        }
                    }
                }
            }

            // здесь реализаия: циклы по блокам, в них по элементам

        } while (max > eps);
    }
    else
    {
        cout << "Error!!!" << endl;
    }

    return icnt;
}

// Параллельная реализия блочного алгоритма Гаусса-Зейделя
int Calc_par(double **u, double **f, int N, double eps)
{
    double max;
    double h = 1.0 / (N + 1);
    int icnt = 0;

    const int BlockSize = 20; // Размер блока
    int bcnt;                 // Количество блоков в ряд

    if (N % BlockSize == 0) // Если количество точек по каждому из направлений сетки делится нацело на размер блока, то проводятся вычисления
    {
        bcnt = N / BlockSize;
        // а сколько блоков в волне?
        do
        {
            icnt++;
            max = 0;

            int bi, bj, i, j, ii, jj;
            double temp, diff;

            for (ii = 0; ii < bcnt; ii++) {
                #pragma omp parallel for schedule(auto) \
                                         private(i, j, bi, bj, temp, diff)
                    for (jj = 0; jj <= ii; jj++) {
                        bi = ii - jj;
                        bj = jj;

                        // Циклы по элементам внутри блока
                        for (i = bi * BlockSize + 1; i <= (bi + 1) * BlockSize; i++) {
                            for (j = bj * BlockSize + 1; j <= (bj + 1) * BlockSize; j++) {
                                temp = u[i][j];
                                u[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1] 
                                                  - h * h * f[i - 1][j - 1]);
                                diff = fabs(u[i][j] - temp);

                                if (diff > max)
                                    max = diff;
                            }
                        }
                    }
            }

            for (ii = bcnt-2; ii >= 0; ii--) {
                #pragma omp parallel for schedule(auto) \
                                         private(i, j, bi, bj, temp, diff)
                    for (jj = 0; jj <= ii; jj++) {
                        bi = jj - ii - 1 + bcnt;
                        bj = bcnt - 1 - jj;

                        // Циклы по элементам внутри блока
                        for (i = bi * BlockSize + 1; i <= (bi + 1) * BlockSize; i++) {
                            for (j = bj * BlockSize + 1; j <= (bj + 1) * BlockSize; j++) {
                                temp = u[i][j];
                                u[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1] 
                                                  - h * h * f[i - 1][j - 1]);
                                diff = fabs(u[i][j] - temp);

                                if (diff > max)
                                    max = diff;
                            }
                        }
                    }
                }   
        } while (max > eps);
    }
    else
    {
        cout << "Error!!!" << endl;
    }
    return icnt;
}


// Функция выделения памяти под 2D массив
double **new_arr(int N)
{
    double **f = new double *[N];
    for (int i = 0; i < N; i++)
    {
        f[i] = new double[N];
    }
    return f;
}

void clear_arr(double** arr, int N) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            arr[i][j] = 0.0;
        }
}

// Функция освобождения памяти 2D массива
void delete_arr(double **arr, int N)
{
    for (int i = 0; i < N; i++)
    {
        delete[] arr[i];
    }
    delete[] arr;
}

// Задание граничных значений
double G(double x, double y)
{
    if (x == 0)
        return 1 - 2 * y;
    if (x == 1)
        return -1 + 2 * y;
    if (y == 0)
        return 1 - 2 * x;
    if (y == 1)
        return -1 + 2 * x;

    return 0;
}

// Задание правой части
double F(double x, double y)
{
    return 2.2;
}

// Инициализация массивов правой части и краевых условий
void Init(double **u, double **f, int N)
{
    memset(&u[0][0], 0, sizeof(double) * N);
    memset(&f[0][0], 0, sizeof(double) * N);

    double h = 1.0 / (N + 1);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
            f[i][j] = F((i + 1) * h, (j + 1) * h);
    }
    for (int i = 1; i < N + 1; i++)
    {
        for (int j = 1; j < N + 1; j++)
            u[i][j] = 0.2;
        u[i][0] = G(i * h, 0);
        u[i][N + 1] = G(i * h, (N + 1) * h);
    }
    for (int j = 0; j < N + 2; j++)
    {
        u[0][j] = G(0, j * h);
        u[N + 1][j] = G((N + 1) * h, j * h);
    }
}

// Функция вывода прореженной матрицы решения
void Output(double **u, int N)
{
    const int K = 5;
    cout << fixed << setprecision(8);
    for (int i = 0; i <= K; i++)
    {
        for (int j = 0; j <= K; j++)
            cout << setw(12) << u[i * (N + 1) / K][j * (N + 1) / K];
        cout << endl;
    }
}


