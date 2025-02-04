#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <time.h>
#include <windows.h>
#include <stdlib.h>
#include <string.h>
#define ARR_SIZE 20000 
static int k = 1;
enum ChooseFillValue
{
    ZERO,
    LINEAR,
    RANDOM
};
enum RuleFill {
    COLUMN,
    LINE
};

int a[ARR_SIZE][ARR_SIZE];
int zero_fill() {
    return 0;
}

int linear_fill() {
   
    return k++;
}

int random_fill() {
    return rand();
}
int main() {
    system("chcp 1251 > null");
    int i, j;
    enum ChooseFillValue chooseFillVal = LINEAR;
    enum RuleFill rule = COLUMN;
    int (*value)() = NULL;
    char msg[50];

    switch (chooseFillVal) {
    case ZERO:
        value = zero_fill;
        strcpy(msg, "Заполнение нулями");
        break;
    case LINEAR:
        value = linear_fill;
        strcpy(msg, "Заполнение линейно");
        break;
    case RANDOM:
        value = random_fill;
        strcpy(msg, "Заполнение рандомно");
        break;
    }
    const char *ruleStr = rule == COLUMN ? "Заполнение по столбцам" : "Заполнение построчно";
    printf("Выбрано:\n%s\n%s", ruleStr,msg);
    ULONGLONG start,end;
    //clock_t start, end;
    switch (rule) {
    case COLUMN:
        start= GetTickCount64();
        //start = clock();
        for (i = 0; i < ARR_SIZE; i++)
            for (j = 0; j < ARR_SIZE; j++)
                a[i][j] = value();
        end = GetTickCount64();
        //end=clock();
        // printf("%lf s.", (float)(end - start) / CLOCKS_PER_SEC);
        printf("%lld ms.", end - start);
        break;
    case LINE:
        start = GetTickCount64();
        //start = clock();
        for (i = 0; i < ARR_SIZE; i++)
            for (j = 0; j < ARR_SIZE; j++)
                a[j][i] = value();
        end = GetTickCount64();
        //end=clock();
       // printf("%lf s.", (float)(end - start) / CLOCKS_PER_SEC);
        printf("%lld ms.", end - start);
        break;
    }
   
    return 0;
}
