#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <windows.h>

int** a;
int N;

double get_time() {
    LARGE_INTEGER frequency, counter;
    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&counter);
    return (double)counter.QuadPart / frequency.QuadPart;
}

void fill_row_wise_zeros() {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            a[i][j] = 0;
}

void fill_row_wise_sequential() {
    int k = 1;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            a[i][j] = k++;
}

void fill_row_wise_random() {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            a[i][j] = rand();
}

void fill_column_wise_zeros() {
    for (int j = 0; j < N; j++)
        for (int i = 0; i < N; i++)
            a[i][j] = 0;
}

void fill_column_wise_sequential() {
    int k = 1;
    for (int j = 0; j < N; j++)
        for (int i = 0; i < N; i++)
            a[i][j] = k++;
}

void fill_column_wise_random() {
    for (int j = 0; j < N; j++)
        for (int i = 0; i < N; i++)
            a[i][j] = rand();
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        printf("Usage: %s <matrix_size>\n", argv[0]);
        return 1;
    }
    N = atoi(argv[1]);
    if (N <= 0) {
        printf("Invalid matrix size\n");
        return 1;
    }

    a = (int**)malloc(N * sizeof(int*));
    if (!a) {
        perror("malloc failed");
        return 1;
    }
    for (int i = 0; i < N; i++) {
        a[i] = (int*)malloc(N * sizeof(int));
        if (!a[i]) {
            perror("malloc failed");
            for (int j = 0; j < i; j++) free(a[j]);
            free(a);
            return 1;
        }
    }

    srand(time(NULL));

    double start, end;
    const char* labels[] = {
        "Row-wise zeros",
        "Row-wise sequential",
        "Row-wise random",
        "Column-wise zeros",
        "Column-wise sequential",
        "Column-wise random"
    };
    void (*functions[])() = {
        fill_row_wise_zeros,
        fill_row_wise_sequential,
        fill_row_wise_random,
        fill_column_wise_zeros,
        fill_column_wise_sequential,
        fill_column_wise_random
    };

    for (int i = 0; i < 6; i++) {
        start = get_time();
        functions[i]();
        end = get_time();
        printf("%s: %.6f s\n", labels[i], end - start);
    }

    for (int i = 0; i < N; i++) free(a[i]);
    free(a);
    return 0;
}