#include <iostream>
#include <cmath>
#include <omp.h>
#include <cstdlib>
#include <cstddef>
#include <stdio.h>

const long int N = 500;
double pi = 3.1415926535;
const double epsilon = 0.00001;
double tau = 0.001;
double norm = 0;

double EuclideanNorm(const double* u) {
    double norm = 0;
    for (int i = 0; i < N; i++) {
        norm += u[i] * u[i];
    }
    return sqrt(norm);
}

void sub(double* a, double* b, double* c) {
    for (int i = 0; i < N; i++) {
        c[i] = a[i] - b[i];
    }
}

void Mul(double* A, double* b, double* result, int n) {
    unsigned int i, j;
    for (i = 0; i < n; i++) {
        result[i] = 0;
        for (j = 0; j < n; j++) {
            result[i] += A[i * n + j] * b[j];
        }
    }
}


void ScalMul(double* A, double tau) {
    int i;
    for (i = 0; i < N; ++i) {
        A[i] = A[i] * tau;
    }
}


double drand(double low, double high) {
    double f = (double)rand() / RAND_MAX;
    return low + f * (high - low);
}

double rand_double() {
    return (double)rand() / RAND_MAX * 4.0 - 2.0;
}

double* Vector_U_builder() {
    double* U = new double[N];
    for (int j = 0; j < N; j++) {
        U[j] = (double)sin((2 * j * pi) / N);
    }
    return U;
}

double* Vector_X_builder() {
    double* X = new double[N];
    for (int j = 0; j < N; j++) {
        X[j] = 0.0;
    }
    return X;
}

int main(int argc, char* argv[]) {
    //srand(100);
    //omp_set_num_threads(4);

    double* X = Vector_X_builder();
    double* U = Vector_U_builder();
    double* Ax = new double[N];

    double normAxb = 0; // ||A*xn - b||
    double normb = 0;
    double saveRes = 10000;
    double res = 0;

    //double end;

    
    //matrix A
    double* A = new double[N * N];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++)
        {
            if (j == i) {
                A[i * N + j] = 2.0;
            }
            else {
                A[i * N + j] = 1.0;
            }
        }
    }
    //

    //vector B
    double* B = new double[N];
    for (int i = 0; i < N; i++) {
        B[i] = 0;
        for (int j = 0; j < N; j++) {
            B[i] += A[i * N + j] * U[j];
        }
    }
    //

    Mul(A, X, Ax, N); // A*xn
    sub(Ax, B, Ax); // A*xn - b
    normAxb = EuclideanNorm(Ax); // ||A*xn - b||
    normb = EuclideanNorm(B);
    ScalMul(Ax, tau); // TAU*(A*xn - b)

    double* nextX = new double[N];

    sub(X, Ax, nextX); // xn - TAU * (A*xn - b)
    //saveRes = normAxb / normb;
    res = normAxb / normb;

    int countIt = 1;

    double start = clock();
    //unsigned int i, j;

    #pragma omp parallel num_threads(4)
    {
        int i = omp_get_thread_num();
        printf_s("Hello from thread %d\n", i);
        //printf("asf");
        //double norm;num_threads(MAX)
        while (res > epsilon) {
            #pragma omp master
            norm = 0;

            #pragma omp for
            for (int i = 0; i < N; i++) {
                X[i] = nextX[i];
            }

            #pragma omp for
            for (int i = 0; i < N; i++) {
                Ax[i] = 0;
                for (int j = 0; j < N; j++) {
                    Ax[i] += A[i * N + j] * X[j];
                }
            }

            #pragma omp for
            for (int i = 0; i < N; i++) {
                Ax[i] = Ax[i] - B[i];
            }

            #pragma omp for reduction(+:norm)
            for (int i = 0; i < N; i++) {
                norm += Ax[i] * Ax[i];
            }

            #pragma omp single
            norm = sqrt(norm);

            #pragma omp for
            for (int i = 0; i < N; ++i) {
                Ax[i] = Ax[i] * tau;
            }

            #pragma omp for
            for (int i = 0; i < N; i++) {
                nextX[i] = X[i] - Ax[i];
            }

            #pragma omp single
            {
                res = norm / normb;
                countIt++;
                if (countIt == 10) {
                    if (saveRes < res) {
                        tau = (-1) * tau;
                    }
                    countIt = 0;
                    saveRes = res;
                }
            }
        }
    }

    bool Check = true;
    for (int i = 0; i < N; i++) {
        if (X[i] - U[i] > 0.001 || X[i] - U[i] < -0.001) {
            Check = false;
        }
        //printf("%f ", X[i] - U[i]);
    }

    if (Check) {
        printf("ItIsTRUE ");
        printf("Time: %f \n", (double)(clock() - start) / 1000);
    }
    else
    {
        printf("FALSE");
    }

    delete []Ax;
    delete[] nextX;
    delete[] X;
    delete[] B;
    delete[] A;
    return 0;
}