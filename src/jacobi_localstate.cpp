/**
 * @file jacobi_localstate.cpp
 * @author Martin Leung
 */
#include "jacobi_localstate.h"
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
    for (int i = 1; i < lNx-1; i++){
        for (int j = 1; j < lNy-1; j++){
            for (int k = 1; k < lNz-1; k++){
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
    for (int i = 1; i < lNx-1; i++){
        for (int j = 1; j < lNy-1; j++){
            for (int k = 1; k < lNz-1; k++){
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
    offset = 0;
    int shift = 0;
    switch (plane_dir) {
        case POS_X:
            offset = (lNx - 1)* lNy * lNz;
            shift = - lNy * lNz;
        case NEG_X:
            if (plane_dir == NEG_X) shift = lNy * lNz;
            M = lNy;
            N = lNz;
            iMul = lNz;
            jMul = 1;
            break;
        case POS_Y:
            offset = (lNy - 1) * lNz;
            shift = - lNz;
        case NEG_Y:
            if (plane_dir == NEG_Y) shift = lNz;
            M = lNx;
            N = lNz;
            iMul = lNy * lNz;
            jMul = 1;
            break;
        case POS_Z:
            offset = lNz - 1;
            shift = 1;
        case NEG_Z:
            if (plane_dir == NEG_Z) shift = -1;
            M = lNx;
            N = lNy;
            iMul = lNy * lNz;
            jMul = lNz;
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