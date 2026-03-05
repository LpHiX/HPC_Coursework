/**
 * @file serial_solver.cpp
 * @author Martin Leung
 */
#include <iostream>
#include "serial_solver.h"

#define F77NAME(x) x##_
extern "C"
{
    double F77NAME(dnrm2)(
        const int &n, 
        double *x, 
        const int &incx);
    double F77NAME(idamax)(
        const int &n, 
        double *x, 
        const int &incx);
}

SerialSolver::SerialSolver(int Nx, int Ny, int Nz, int test, double epsilon): 
    Nx(Nx),
    Ny(Ny),
    Nz(Nz),
    nx(Nx - 2),
    ny(Ny - 2),
    nz(Nz - 2),
    test(test),
    epsilon(epsilon),
    hx(1.0 / (Nx - 1)),
    hy(1.0 / (Ny - 1)),
    hz(1.0 / (Nz - 1)),
    hx2(1/(hx*hx)),
    hy2(1/(hy*hy)),
    hz2(1/(hz*hz))
    {
        u   = new double[Nx * Ny * Nz];
        ddu = new double[nx * ny * nz];
        f   = new double[nx * ny * nz];
        r   = new double[nx * ny * nz];
        
        if (test == 1){
            initialize_test_1();
        }
    }

int SerialSolver::uIndex(int i, int j, int k){
    return (i * Ny + j) * Nz + k;
}
int SerialSolver::u2rIndex(int i, int j, int k){
    return ((i-1) * ny + (j-1)) * nz + (k-1);
}

double SerialSolver::run_solver(){
    for (int i = 1; i < Nx-1; i++){
        for (int j = 1; j < Ny-1; j++){
            for (int k = 1; k < Nz-1; k++){
                int reduced_index = u2rIndex(i, j, k);
                f[reduced_index] = 6;
                ddu[reduced_index] =(
                    +     u[uIndex(i + 1, j    , k    )]
                    - 2 * u[uIndex(i    , j    , k    )]
                    +     u[uIndex(i - 1, j    , k    )]) * hx2 + (
                    +     u[uIndex(i    , j + 1, k    )]
                    - 2 * u[uIndex(i    , j    , k    )]
                    +     u[uIndex(i    , j - 1, k    )]) * hy2 + (
                    +     u[uIndex(i    , j    , k + 1)]
                    - 2 * u[uIndex(i    , j    , k    )]
                    +     u[uIndex(i    , j    , k - 1)]) * hz2;
                r[reduced_index] = f[reduced_index] - ddu[reduced_index];
            }
        }
    }
    residual = F77NAME(dnrm2)(nx * ny * nz, r, 1);

    return residual;
}

void SerialSolver::initialize_test_1(){
    for (int i = 0; i < Nx; i++){
        double x = i * hx;
        for (int j = 0; j < Ny; j++){
            double y = j * hy;
            for (int k = 0; k < Nz; k++){
                double z = k * hz;
                u[uIndex(i, j, k)] = x*x + y*y + z*z;
            }
        }
    }
}

SerialSolver::~SerialSolver(){
    delete[] u;
    delete[] ddu;
    delete[] f;
    delete[] r;
}