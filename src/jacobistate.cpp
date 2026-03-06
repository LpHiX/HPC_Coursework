/**
 * @file jacobi_localstate.cpp
 * @author Martin Leung
 */
#include <jacobi_localstate.h>
#include <cmath>
/**
 * Math Explanation:
 * 
 * N=11,
 * Serial :
 * 
 * 0 1 2 3 4 5 6 7 8 9 10
 * 
 * MPI :
 * N=7         N=6
 * 0 1 2 3 4 5 6
 *           5 6 7 8 7 19
 * Shared boundary: 5,6
 * Aim: Pass as little information as possible. Meaning: for loop might get a little complicated.
 */



JacobiLocalState::JacobiLocalState(int Nx, int Ny, int Nz, double* f): 
    Nx(Nx),
    Ny(Ny),
    Nz(Nz),
    nx(Nx - 2),
    ny(Ny - 2),
    nz(Nz - 2),
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
        r   = new double[nx * ny * nz];        
    }

void JacobiLocalState::jacobi_step(){
    for (int i = 1; i < Nx-1; i++){
        for (int j = 1; j < Ny-1; j++){
            for (int k = 1; k < Nz-1; k++){
                int redind = u2rIndex(i,j,k);
                u2[uIndex(i,j,k)] = j_coeff * (
                    + (u[uIndex(i + 1, j    , k    )] + u[uIndex(i - 1, j    , k    )]) * hx2 
                    + (u[uIndex(i    , j + 1, k    )] + u[uIndex(i    , j - 1, k    )]) * hy2
                    + (u[uIndex(i    , j    , k + 1)] + u[uIndex(i    , j    , k - 1)]) * hz2 - f[redind]);
            }
        }
    }
    double* temp = u;
    u = u2;
    u2 = temp;
}

double JacobiLocalState::get_residual() const{
    double sum_residual = 0;
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
                double r_step = f[reduced_index] - ddu[reduced_index];
                sum_residual += r_step*r_step;
            }
        }
    }
    return sqrt(sum_residual);
}

void JacobiLocalState::get_boundary_constants(Direction plane_dir, bool plane_sign, int &offset, int &M, int &N, int &iMul, int &jMul) const {
    offset = 0;
    switch (plane_dir) {
        case DIR_X:
            M = Ny;
            N = Nz;
            iMul = Nz;
            jMul = 1;
            if (plane_sign){offset = (Nx - 1)* Ny * Nz;}
            break;
        case DIR_Y:
            M = Nx;
            N = Nz;
            iMul = Ny * Nz;
            jMul = 1;
            if (plane_sign){offset = (Ny - 1) * Nz;}
        case DIR_Z:
            M = Nx;
            N = Ny;
            iMul = Ny * Nz;
            jMul = Ny;
            if (plane_sign){offset = Nz - 1;}
            break;
        default:
            throw "Not a possible direction";
    }
}

void JacobiLocalState::set_u_boundary(double *plane, Direction plane_dir, bool plane_sign){
    int offset, M, N, iMul, jMul;
    get_boundary_constants(plane_dir, plane_sign, offset, M, N, iMul, jMul);
    

    for (int im = 0; im < M; im++){
        for (int jn = 0; jn < N; jn++){
            u[offset + im * iMul + jn * jMul] = plane[im * N + jn];
        }
    }
}

void JacobiLocalState::get_u_boundary(double *&plane, Direction plane_dir, bool plane_sign) const{
    int offset, M, N, iMul, jMul;
    get_boundary_constants(plane_dir, plane_sign, offset, M, N, iMul, jMul);
    
    for (int im = 0; im < M; im++){
        for (int jn = 0; jn < N; jn++){
            plane[im * N + jn] = u[offset + im * iMul + jn * jMul];
        }
    }
}

JacobiLocalState::~JacobiLocalState(){
    delete[] u;
    delete[] u2;
    delete[] ddu;
    delete[] f;
    delete[] r;
}