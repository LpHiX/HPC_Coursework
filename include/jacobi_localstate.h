/**
 * @file jacobi_localstate.h
 * @author Martin Leung
 */
#pragma once

enum Direction {
    POS_X,
    NEG_X,
    POS_Y,
    NEG_Y,
    POS_Z,
    NEG_Z
};

class JacobiLocalState{
public:
    const int Nx, Ny, Nz, lNx, lNy, lNz, lnx, lny, lnz;
    const double hx, hy, hz, hx2, hy2, hz2, j_coeff;
    
    JacobiLocalState(int Nx, int Ny, int Nz, int lNx, int lNy, int lNz, double* lf);
    JacobiLocalState() = delete;
    ~JacobiLocalState();

    inline int uIndex(int i, int j, int k) const { return (i * lNy + j) * lNz + k; }
    inline int u2rIndex(int i, int j, int k) const { return ((i-1) * lny + (j-1)) * lnz + (k-1); }

    double get_residualsquared() const;
    void set_u_boundary(double *plane, Direction plane_dir);
    void jacobi_step();
    int get_planesize(Direction plane_dir) const;
    void get_u_boundary(double *&plane, Direction plane_dir) const;
    void get_boundary_constants(Direction plane_dir, bool is_sending, int &offset, int &M, int &N, int &iMul, int &jMul) const;
private:
    double *lu, *lu2, *lddu, *lf, *lr;

};

