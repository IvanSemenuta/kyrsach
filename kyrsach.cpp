#include <iostream>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <locale.h>

using namespace std;

#define input_filename (char*)"InputMatrix.txt"


#pragma warning(disable:4996)


void freeMatrix(double** arr, int N) {
    int i;

    for (i = 0; i < N; i++) {
        free(arr[i]);
        arr[i] = NULL;
    }
    arr = NULL;
    free(arr);
}


void printExpandedMatrix(double** arr, double* freeRow, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%7.2f ", arr[i][j]);

        }
        printf("| %-7.2f\n", freeRow[i]);
    }
    printf("\n");
}


void copyArra(double** copy, double** arr, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            copy[i][j] = arr[i][j];
        }
    }
}


double determinant(double** arr, int N) {
    double** subMatrix = (double**)malloc((N - 1) * sizeof(double*));;
    int i, p, j;

    for (i = 0; i < N - 1; i++) {
        subMatrix[i] = (double*)malloc((N - 1) * sizeof(double));
    }

    if (N == 1) return arr[0][0];
    if (N == 2) return arr[0][0] * arr[1][1] - arr[0][1] * arr[1][0];

    double det = 0;
    for (p = 0; p < N; p++) {
        for (i = 1; i < N; i++) {
            int subCol = 0;
            for (j = 0; j < N; j++) {
                if (j == p)
                    continue;
                subMatrix[i - 1][subCol++] = arr[i][j];
            }
        }
        det += (p % 2 == 0 ? 1 : -1) * arr[0][p] * determinant(subMatrix, N - 1);
    }

    for (i = 0; i < N - 1; i++)
        free(subMatrix[i]);
    free(subMatrix);

    return det;
}


void cramer(double** arr, double* freeRow, double* result, int N) {
    double** copy = (double**)malloc(N * sizeof(double*));
    int i, j;
    double det_i,
        det = determinant(arr, N);

    for (i = 0; i < N; i++) {
        copy[i] = (double*)malloc(N * sizeof(double));
    }

    printf("\n----------------------------- Метод крамера -----------------------------\n");

    if (det == 0) {
        printf("\nВизначник: %-.3lf (розв'язок цим методом неможливий)\n", det);
        for (i = 0; i < N; i++) {
            result[i] = -1;
        }
    }
    else {

        for (i = 0; i < N; i++) {
            printf("      det(%d)\t", i + 1);
        }
        printf("\n");
        for (i = 0; i < N; i++) {
            printf("x%d = -------\t", i + 1);
        }
        printf("\n");
        for (i = 0; i < N; i++) {
            printf("       det\t");
        }

        printf("\n");
        printf("\ndet = %.2lf != 0\n", det);
        for (i = 0; i < N; i++) {
            copyArra(copy, arr, N);
            for (j = 0; j < N; j++) {
                copy[j][i] = freeRow[j];
            }

            det_i = determinant(copy, N);
            printf("det(%d) = %-10.2lf\t", i + 1, det_i);
            result[i] = det_i / det;
        }

        printf("\n");
        for (i = 0; i < N; i++) {
            printf("x%d = %-10.4lf\t\t", i + 1, result[i]);
        }
    }
    printf("\n-------------------------------------------------------------------------\n\n");

    freeMatrix(copy, N);
}


void gaussElimination(double** arr, double* freeRow, double* result, int N) {
    double* copy_freeRow = (double*)malloc(N * sizeof(double));
    double** copy_arr = (double**)malloc(N * sizeof(double));
    int i, j, rank = N, g;

    printf("\n------------------------------ Метод Гауса ------------------------------\n\n");

    for (i = 0; i < N; i++) {
        copy_arr[i] = (double*)malloc(N * sizeof(double));
        copy_freeRow[i] = freeRow[i];
    }
    copyArra(copy_arr, arr, N);

    for (i = 0; i < N; i++) {
        int pivotRow = i;
        for (j = i + 1; j < N; j++) {
            if (fabs(copy_arr[j][i]) > fabs(copy_arr[pivotRow][i])) {
                pivotRow = j;
            }
        }

        if (fabs(copy_arr[pivotRow][i]) < 0.000001) {
            rank--;
            continue;
        }

        if (pivotRow != i) {
            double* tempRow = copy_arr[i];
            copy_arr[i] = copy_arr[pivotRow];
            copy_arr[pivotRow] = tempRow;

            double tempConst = copy_freeRow[i];
            copy_freeRow[i] = copy_freeRow[pivotRow];
            copy_freeRow[pivotRow] = tempConst;
        }

        for (g = i + 1; g < N; g++) {
            double factor = copy_arr[g][i] / copy_arr[i][i];
            for (int k = i; k < N; k++) {
                copy_arr[g][k] -= factor * copy_arr[i][k];
            }
            copy_freeRow[g] -= factor * copy_freeRow[i];
        }

        if (rank == N) printExpandedMatrix(copy_arr, copy_freeRow, N);
    }

    printf("Ранг матриці: %d ", rank);

    if (rank < N) {
        printf("(розв'язок за цим методом неможливий)\n");
        for (i = 0; i < N; i++) {
            result[i] = -1;
        }
    }
    else {
        // Виконуємо зворотний хід для знаходження невідомих
        for (i = N - 1; i >= 0; i--) {
            result[i] = copy_freeRow[i];
            for (j = i + 1; j < N; j++) {
                result[i] -= copy_arr[i][j] * result[j];
            }
            result[i] /= copy_arr[i][i];
        }
    }
    printf("\n-------------------------------------------------------------------------\n\n");

    freeMatrix(copy_arr, N);
    free(copy_freeRow);
}


void readMatrixFromFile(double*** arr, double** freeRow, int* N) {
    int i, j;

    FILE* file = fopen(input_filename, "r");
    if (!file) {
        printf("Помилка відкриття файлу!\n");
        return;
    }

    fscanf(file, "%d", N);

    *freeRow = (double*)malloc(*N * sizeof(double));
    *arr = (double**)malloc(*N * sizeof(double*));
    for (i = 0; i < *N; i++) {
        (*arr)[i] = (double*)malloc(*N * sizeof(double));
    }


    for (i = 0; i < *N; i++) {
        for (j = 0; j < *N; j++) {
            fscanf(file, "%lf", &(*arr)[i][j]);
        }
    }

    for (i = 0; i < *N; i++) {
        fscanf(file, "%lf", &(*freeRow)[i]);
    }

    fclose(file);
}


void readMartixFromConsole(double*** arr, double** freeRow, int* N) {
    int userInputType = 0,
        userTryAgain = 0,
        i, j;

    // ------------------------------------------- Entering size -------------------------------------------
    printf("\nВведіть розмір квадратної матриці: ");
    scanf("%d", N);

    *freeRow = (double*)malloc(*N * sizeof(double));
    *arr = (double**)malloc(*N * sizeof(double*));
    for (i = 0; i < *N; i++) {
        (*arr)[i] = (double*)malloc(*N * sizeof(double));
    }


    // ----------------------------------------- Input type (entering) -----------------------------------------
    printf("\nВкажіть як ви хочете ввести дані:\n");
    printf("1 - Власні дані\n");
    printf("2 - Випадкові \n");

    do {
        scanf("%d", &userInputType);

        if (userInputType != 1 && userInputType != 2) {
            printf("Такого вводу немає!\n");
        }
    } while (userInputType != 1 && userInputType != 2);
    while (getchar() != '\n');

    // ----------------------------------------- Input type (calculating) -----------------------------------------
    switch (userInputType) {
    case 1: {
        printf("Введіть елементи матриці по рядках (0 1 2 3):\n");
        for (i = 0; i < *N; i++) {
            for (j = 0; j < *N; j++) {
                scanf("%lf", &(*arr)[i][j]);
                //printf("%-5.2f ", (*arr)[i][j]);
            }
        }
        printf("Введіть рядок вільних членів (0 1 2 3):\n");
        for (i = 0; i < *N; i++) {
            scanf("%lf", &(*freeRow)[i]);
        }
        break;
    }
    case 2: {
        srand(time(0));
        for (i = 0; i < *N; i++) {
            for (j = 0; j < *N; j++) {
                (*arr)[i][j] = 10 + rand() % 10 - 10;
                printf("%-5.2f ", (*arr)[i][j]);
            }
            printf("\n");
        }

        for (i = 0; i < *N; i++) {
            (*freeRow)[i] = 10 + rand() % 10 - 10;
        }
        break;
    }
    }
}


void writeResultToFile(double** arr, double** results, double* freeRow, int N) {
    int i, j;

    FILE* file = fopen(input_filename, "w");
    if (!file) {
        printf("Помилка відкриття файлу!\n");
        return;
    }
    fprintf(file, "%d\n", N);
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            fprintf(file, "%0.2lf ", arr[i][j]);
        }
        fprintf(file, "\n");
    }

    fprintf(file, "\n");
    for (i = 0; i < N; i++) {
        fprintf(file, "%0.2lf ", freeRow[i]);
    }

    fprintf(file, "\nКрамер: ");
    for (i = 0; i < N; i++) {
        fprintf(file, "%0.3lf ", results[0][i]);
    }


    fprintf(file, "\nГаусс: ");
    for (i = 0; i < N; i++) {
        fprintf(file, "%0.3lf ", results[1][i]);
    }

    fclose(file);
}


int main() {


    int userInputFrom = 0,
        index,
        userMethod = 0,
        N, i, j;

    double** arr = nullptr,
        * freeRow = nullptr,
        ** results = nullptr;

    system("chcp 1251 & cls");

    do {

        do {
            printf("З чим ви хочете працювати:\n");
            printf("1 - Консоль\n");
            printf("2 - Файл\n");
            printf("3 - Вийти\n");

            scanf("%d", &userInputFrom);

            switch (userInputFrom) {
            case 1: { readMartixFromConsole(&arr, &freeRow, &N); break; }
            case 2: { readMatrixFromFile(&arr, &freeRow, &N); break; }
            case 3: { return 0; }
            }
        } while (userInputFrom < 1 || userInputFrom > 3);


        system("cls");
        printf("Вхідна матриці:\n");
        printExpandedMatrix(arr, freeRow, N);

        // ------------------------------------------- Введення даних -------------------------------------------
        printf("\nВведіть метод розвязку СЛАР:\n");
        printf("0 - Вийти\n");
        printf("1 - Метод крамера\n");
        printf("2 - Метод гауса\n");
        printf("3 - Табличний вигляд усіх методів\n");

        do {
            scanf("%d", &userMethod);
            if (userMethod < 0 || userMethod > 3) {
                printf("Такого методу не існує\n");
                while (getchar() != 0);
                system("cls");
            }
        } while (userMethod < 0 || userMethod > 3);

        // ----------------------------------------- Розподілення пам'яті -----------------------------------------

        index = 0;
        if (userMethod == 3) {
            results = (double**)malloc(2 * sizeof(double**));
            for (i = 0; i < 2; i++) {
                results[i] = (double*)malloc(N * sizeof(double*));
            }
        }
        else {
            results = (double**)malloc(1 * sizeof(double**));
            results[0] = (double*)malloc(N * sizeof(double*));
        }

        // ----------------------------------------- Методи обчислення -----------------------------------------

        if (userMethod == 1 || userMethod == 3) {
            cramer(arr, freeRow, results[index], N);

            if (userMethod == 3)
                index++;
        }
        if (userMethod == 2 || userMethod == 3) {
            gaussElimination(arr, freeRow, results[index], N);

            if (userMethod == 3)
                index++;
        }

        // ---------------------------------------------- Вивід ---------------------------------------------
        if (userMethod != 3) {
            printf("Розвязок СЛАР:\n");
            for (i = 0; i < N; i++) {
                printf("x%d: %7.4f\n", i, results[0][i]);
            }
        }
        else {
            printf("\t  Крамер\t  Гаус\n");
            for (i = 0; i < N; i++) {
                printf("x%d:\t", i + 1);
                for (j = 0; j < index; j++) {
                    printf("%7.3f\t\t", results[j][i]);
                }
                printf("\n");
            }
            writeResultToFile(arr, results, freeRow, N);
        }


        system("pause");
        system("cls");

        // ---------------------------------------------- Очислення пам'яті ---------------------------------------------
        freeMatrix(results, index);
        freeMatrix(arr, N);
        free(freeRow);

    } while (userMethod != 0);

    return 0;
}
