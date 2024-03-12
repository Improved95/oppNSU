#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define PI 3.14159265358979323846
#define N 2200

const double epsilon = 0.00001;
const double tao = 0.0005;

void printMatrix(double *matrix) {
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            printf("%f ", matrix[i * N + j]);
        }
        printf("\n");
    }
}

void printVector(double *vector) {
    for (size_t i = 0; i < N; ++i) {
        printf("%f ", vector[i]);
    }
    printf("\n");
}

void setZeroVector(double *vector) {
    memset(vector, 0, N * sizeof(double));
}

void subVector(double *vector1, double *vector2) {
    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i) {
        vector1[i] -= vector2[i];
    } 
}

void mulMatrixVector(double *matrix, double *inputVector, double *outputVector) {
    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            outputVector[i] += matrix[i * N + j] * inputVector[j];
        }
    }
}

double getNorm(double *vector) {
    double res = 0;
    #pragma omp parallel for reduction(+:res)
    for (size_t i = 0; i < N; ++i) {
        double a = vector[i];
        res += (a * a);
    }
    return res;
}

int main() {
    double *matrixA = malloc(N * N * sizeof(double));
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            matrixA[i * N + j] = 1.0;
            if (i == j) {
                matrixA[i * N + j] = 2.0;
            }
        }
    }

    double *vectorU = calloc(sizeof(double), N);
    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i) {
        vectorU[i] = sin(2 * PI * (i + 1) / N);
    }
    // printVector(vectorU);
    // printf("\n");

    double *vectorX = calloc(sizeof(double), N);
    double *vectorB = calloc(sizeof(double), N);
    double *vectorAxn_b = calloc(sizeof(double), N);

    setZeroVector(vectorB);
    mulMatrixVector(matrixA, vectorU, vectorB);

    double normB = getNorm(vectorB);
    
    for(size_t k = 0; 1; ++k) {
        setZeroVector(vectorAxn_b);
        mulMatrixVector(matrixA, vectorX, vectorAxn_b);
        subVector(vectorAxn_b, vectorB);

        double normAx_b = getNorm(vectorAxn_b);
        if (normAx_b / normB < epsilon * epsilon) {
            // printf("%ld\n", k);
            break;
        }

        #pragma omp parallel for
        for (size_t i = 0; i < N; ++i) {
            vectorX[i] = vectorX[i] - (tao * vectorAxn_b[i]);
        }
    }

    // printVector(vectorX);

    free(matrixA);
    free(vectorU);
    free(vectorX);
    free(vectorB);
    free(vectorAxn_b);
    return 0;
}
