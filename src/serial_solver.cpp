/**
 * @file serial_solver.cpp
 * @author Martin Leung
 */
#include <iostream>
#include "serial_solver.h"
#include <cmath>

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
    max_iter(1e6),
    epsilon(epsilon),
    hx(1.0 / (Nx - 1)),
    hy(1.0 / (Ny - 1)),
    hz(1.0 / (Nz - 1)),
    hx2(1/(hx*hx)),
    hy2(1/(hy*hy)),
    hz2(1/(hz*hz)),
    j_coeff(0.5 / (hx2 + hy2 + hz2))
    {
        u   = new double[Nx * Ny * Nz]();
        u2  = new double[Nx * Ny * Nz]();
        ddu = new double[nx * ny * nz];
        f   = new double[nx * ny * nz];
        r   = new double[nx * ny * nz];
        
        initialize_test(test);
    }

int SerialSolver::uIndex(int i, int j, int k){
    return (i * Ny + j) * Nz + k;
}
int SerialSolver::u2rIndex(int i, int j, int k){
    return ((i-1) * ny + (j-1)) * nz + (k-1);
}

int SerialSolver::run_solver(){
    int iterations = 0;
    double residual = get_residual();

    while (residual > epsilon && iterations < max_iter){
        // std::cout << "Iter: " << iterations << ", residual: " << residual << std::endl;
        for (int i = 1; i < Nx-1; i++){
            for (int j = 1; j < Ny-1; j++){
                for (int k = 1; k < Nz-1; k++){
                    u2[uIndex(i,j,k)] = j_coeff * (
                        + (u[uIndex(i + 1, j    , k    )] + u[uIndex(i - 1, j    , k    )]) * hx2 
                        + (u[uIndex(i    , j + 1, k    )] + u[uIndex(i    , j - 1, k    )]) * hy2
                        + (u[uIndex(i    , j    , k + 1)] + u[uIndex(i    , j    , k - 1)]) * hz2 - f[u2rIndex(i,j,k)]);
                }
            }
        }
        double* temp = u;
        u = u2;
        u2 = temp;

        residual = get_residual();
        iterations++;
    }
    if (iterations >= max_iter){
        std::cout << "Failed to converge within " << max_iter << " iterations, residual: " << residual << std::endl;
        return 1;
    }
    std::cout << "Converged in " << iterations << " iterations, residual: " << residual << std::endl;
    return 0;
}

double SerialSolver::get_residual(){
    for (int i = 1; i < Nx-1; i++){
        for (int j = 1; j < Ny-1; j++){
            for (int k = 1; k < Nz-1; k++){
                int reduced_index = u2rIndex(i, j, k);
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
    return F77NAME(dnrm2)(nx * ny * nz, r, 1);;
}

void SerialSolver::initialize_test(int test){
    switch (test ) {
        case 1:
            for (int i = 0; i < Nx; i++){
                double x = i * hx;
                for (int j = 0; j < Ny; j++){
                    double y = j * hy;
                    for (int k = 0; k < Nz; k++){
                        double z = k * hz;
                        double scalar = x*x + y*y + z*z;
                        u[uIndex(i, j, k)] = scalar;
                        u2[uIndex(i, j, k)] = scalar;
                    }
                }
            }
            for (int i = 1; i < Nx-1; i++){
                for (int j = 1; j < Ny-1; j++){
                    for (int k = 1; k < Nz-1; k++){
                        f[u2rIndex(i, j, k)] = 6;
                        u[uIndex(i,j,k)] = 0;
                    }
                }
            }
            break;
        case 2:
            for (int i = 0; i < Nx; i++){
                double x = i * hx;
                for (int j = 0; j < Ny; j++){
                    double y = j * hy;
                    for (int k = 0; k < Nz; k++){
                        double z = k * hz;
                        double scalar = sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z);
                        u[uIndex(i, j, k)] = scalar;
                        u2[uIndex(i, j, k)] = scalar;
                    }
                }
            }
            for (int i = 1; i < Nx-1; i++){
                double x = i * hx;
                for (int j = 1; j < Ny-1; j++){
                    double y = j * hy;
                    for (int k = 1; k < Nz-1; k++){
                        double z = k * hz;
                        f[u2rIndex(i, j, k)] = - 3 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z);
                        u[uIndex(i,j,k)] = 0;
                    }
                }
            }
            break;
        case 3:
            for (int i = 0; i < Nx; i++){
                double x = i * hx;
                for (int j = 0; j < Ny; j++){
                    double y = j * hy;
                    for (int k = 0; k < Nz; k++){
                        double z = k * hz;
                        double scalar = sin(M_PI * x) * sin(4 * M_PI * y) * sin(8 * M_PI * z);
                        u[uIndex(i, j, k)] = scalar;
                        u2[uIndex(i, j, k)] = scalar;
                    }
                }
            }
            for (int i = 1; i < Nx-1; i++){
                double x = i * hx;
                for (int j = 1; j < Ny-1; j++){
                    double y = j * hy;
                    for (int k = 1; k < Nz-1; k++){
                        double z = k * hz;
                        f[u2rIndex(i,j,k)] = - 81 * M_PI * M_PI * sin(M_PI * x) * sin(4 * M_PI * y) * sin(8 * M_PI * z);
                        u[uIndex(i,j,k)] = 0;
                    }
                }
            }
            break;
        case 4:
            for (int i = 1; i < Nx-1; i++){
                double x = i * hx;
                for (int j = 1; j < Ny-1; j++){
                    double y = j * hy;
                    for (int k = 1; k < Nz-1; k++){
                        double z = k * hz;
                        f[u2rIndex(i,j,k)] = 100 * exp(-100 * ((x - 0.5)*(x - 0.5)+(y - 0.5)*(y - 0.5)+(z - 0.5)*(z - 0.5)));
                    }
                }
            }
            break;
        case 5:
            for (int i = 1; i < Nx-1; i++){
                double x = i * hx;
                for (int j = 1; j < Ny-1; j++){
                    for (int k = 1; k < Nz-1; k++){
                        if(x < 0.5) {f[u2rIndex(i,j,k)] = -1;}
                        else {f[u2rIndex(i,j,k)] = 1;}
                    }
                }
            }
            break;
    }
}

SerialSolver::~SerialSolver(){
    delete[] u;
    delete[] u2;
    delete[] ddu;
    delete[] f;
    delete[] r;
}