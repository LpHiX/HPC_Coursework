/**
 * @file jacobi_localstate.h
 * @author Martin Leung
 */
#pragma once

class JacobiLocalState{
public:
    const int Nx, Ny, Nz, nx, ny, nz;
    const double hx, hy, hz, hx2, hy2, hz2, j_coeff;
    
    JacobiLocalState(int Nx, int Ny, int Nz, double* f);
    JacobiLocalState() = delete;
    ~JacobiLocalState();

    inline int uIndex(int i, int j, int k) const { return (i * Ny + j) * Nz + k; }
    inline int u2rIndex(int i, int j, int k) const { return ((i-1) * ny + (j-1)) * nz + (k-1); }

    double get_residual() const;
    void set_u_boundary(double *plane, Direction plane_dir, bool plane_sign);
private:
    double *u, *u2, *ddu, *f, *r;
    void jacobi_step();
    void get_u_boundary(double *&plane, Direction plane_dir, bool plane_sign) const;
    void get_boundary_constants(Direction plane_dir, bool plane_sign, int &offset, int &M, int &N, int &iMul, int &jMul) const;
};

enum Direction {
    DIR_X,
    DIR_Y,
    DIR_Z
};