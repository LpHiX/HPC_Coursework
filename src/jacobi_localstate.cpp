/**
 * @file jacobi_localstate.cpp
 * @author Martin Leung
 */
#include "jacobi_localstate.h"
#include <cmath>




JacobiLocalState::JacobiLocalState(int Nx, int Ny, int Nz, int lNx, int lNy, int lNz, double* lf): 
Nx(Nx),
Ny(Ny),
Nz(Nz),
lNx(lNx),
lNy(lNy),
lNz(lNz),
lnx(lNx - 2),
lny(lNy - 2),
lnz(lNz - 2),
hx(1.0 / (Nx - 1)),
hy(1.0 / (Ny - 1)),
hz(1.0 / (Nz - 1)),
hx2(1/(hx*hx)),
hy2(1/(hy*hy)),
hz2(1/(hz*hz)),
j_coeff(0.5 / (hx2 + hy2 + hz2)),
lf(lf)
{
    lu   = new double[lNx * lNy * lNz]();
    lu2  = new double[lNx * lNy * lNz]();
    lddu = new double[lnx * lny * lnz];
    lr   = new double[lnx * lny * lnz];        
}

void JacobiLocalState::jacobi_step(){
    for (int k = 1; k < lNz-1; k++){
        for (int j = 1; j < lNy-1; j++){
            for (int i = 1; i < lNx-1; i++){
                int redind = u2rIndex(i,j,k);
                lu2[uIndex(i,j,k)] = j_coeff * (
                    + (lu[uIndex(i + 1, j    , k    )] + lu[uIndex(i - 1, j    , k    )]) * hx2 
                    + (lu[uIndex(i    , j + 1, k    )] + lu[uIndex(i    , j - 1, k    )]) * hy2
                    + (lu[uIndex(i    , j    , k + 1)] + lu[uIndex(i    , j    , k - 1)]) * hz2 - lf[redind]);
            }
        }
    }
    double* temp = lu;
    lu = lu2;
    lu2 = temp;
}

double JacobiLocalState::get_residualsquared() const{
    double sum_residual = 0;
    for (int k = 1; k < lNz-1; k++){
        for (int j = 1; j < lNy-1; j++){
            for (int i = 1; i < lNx-1; i++){
                int reduced_index = u2rIndex(i, j, k);
                lddu[reduced_index] =(
                    +     lu[uIndex(i + 1, j    , k    )]
                    - 2 * lu[uIndex(i    , j    , k    )]
                    +     lu[uIndex(i - 1, j    , k    )]) * hx2 + (
                    +     lu[uIndex(i    , j + 1, k    )]
                    - 2 * lu[uIndex(i    , j    , k    )]
                    +     lu[uIndex(i    , j - 1, k    )]) * hy2 + (
                    +     lu[uIndex(i    , j    , k + 1)]
                    - 2 * lu[uIndex(i    , j    , k    )]
                    +     lu[uIndex(i    , j    , k - 1)]) * hz2;
                double r_step = lf[reduced_index] - lddu[reduced_index];
                sum_residual += r_step*r_step;
            }
        }
    }
    return sum_residual;
}

int JacobiLocalState::get_planesize(Direction plane_dir) const{
    int M, N, dummy_offset, dummy_i, dummy_j;
    get_boundary_constants(plane_dir, false, dummy_offset, M, N, dummy_i, dummy_j);
    return M * N;
}

void JacobiLocalState::get_boundary_constants(Direction plane_dir, bool is_sending, int &offset, int &M, int &N, int &iMul, int &jMul) const {

    // // Writing it all back out to figure it out one by one, this is easier for me to understand.
    // // X positive -> z,y
    // int offset = lNx - 1;
    // int M = lNz;
    // int N = lNy;
    // int iMul = lNy * lNx;
    // int jMul = lNx;
    // for (int im = 0; im < M; im++){
    //     for (int jn = 0; jn < N; jn++){
    //         u[offset + im * lNy * lNx  + jn * lNx] = plane[im * N + jn];
    //     }
    // }
    // // Y positive -> z,x
    // int offset = (lNy - 1) * lNx;
    // int M = lNz;
    // int N = lNx;
    // int iMul = lNy * lNx;
    // int jMul = 1;
    // for (int im = 0; im < M; im++){
    //     for (int jn = 0; jn < N; jn++){
    //         u[offset + im * lNy * lNx  + jn] = plane[im * N + jn];
    //     }
    // }
    // // Z positive -> x,y
    // int offset = (lNz - 1) * lNy * lNx;
    // int M = lNy;
    // int N = lNx;
    // int iMul = lNx
    // int jMul = 1;
    // for (int im = 0; im < M; im++){
    //     for (int jn = 0; jn < N; jn++){
    //         u[offset + im * lNy * lNx  + jn * lNx] = plane[im * N + jn];
    //     }
    // }
    offset = 0;
    int shift = 0;
    switch (plane_dir) {
        case POS_X:
            offset = lNx - 1;
            shift = - 1;
        case NEG_X:
            if (plane_dir == NEG_X) shift = 1;
            M = lNz;
            N = lNy;
            iMul = lNy * lNx;
            jMul = lNx;
            break;
        case POS_Y:
            offset = (lNy - 1) * lNx;
            shift = - lNx;
        case NEG_Y:
            if (plane_dir == NEG_Y) shift = lNx;
            M = lNz;
            N = lNx;
            iMul = lNy * lNx;
            jMul = 1;
            break;
        case POS_Z:
            offset = (lNz - 1) * lNy * lNx;
            shift = - lNy * lNx;
        case NEG_Z:
            if (plane_dir == NEG_Z) shift = lNy * lNx;
            M = lNy;
            N = lNx;
            iMul = lNx;
            jMul = 1;
            break;
        default:
            throw "Not a possible direction";
    }
    if (is_sending) offset += shift;
}

void JacobiLocalState::set_u_boundary(double *plane, Direction plane_dir){
    int offset, M, N, iMul, jMul;
    get_boundary_constants(plane_dir, false, offset, M, N, iMul, jMul);
    

    for (int im = 0; im < M; im++){
        for (int jn = 0; jn < N; jn++){
            lu[offset + im * iMul + jn * jMul] = plane[im * N + jn];
        }
    }
}

void JacobiLocalState::get_u_boundary(double *&plane, Direction plane_dir) const{
    int offset, M, N, iMul, jMul;
    get_boundary_constants(plane_dir, true, offset, M, N, iMul, jMul);
    
    for (int im = 0; im < M; im++){
        for (int jn = 0; jn < N; jn++){
            plane[im * N + jn] = lu[offset + im * iMul + jn * jMul];
        }
    }
}

JacobiLocalState::~JacobiLocalState(){
    delete[] lu;
    delete[] lu2;
    delete[] lddu;
    delete[] lf;
    delete[] lr;
}