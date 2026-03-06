#pragma once

class JacobiState{
public:
    const int Nx, Ny, Nz, nx, ny, nz;
    const double hx, hy, hz, hx2, hy2, hz2, j_coeff;
    
    JacobiState(int Nx, int Ny, int Nz, double* f);
    JacobiState() = delete;
    ~JacobiState();

    inline int uIndex(int i, int j, int k) const { return (i * Ny + j) * Nz + k; }
    inline int u2rIndex(int i, int j, int k) const { return ((i-1) * ny + (j-1)) * nz + (k-1); }

    double get_residual() const;
    void set_u_boundary(double *plane, Direction plane_dir, bool plane_sign);
private:
    double *u, *u2, *ddu, *f, *r;
    void jacobi_step();
};

enum Direction {
    DIR_X,
    DIR_Y,
    DIR_Z
};