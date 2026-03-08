/**
 * @file jacobi_localstate.h
 * @author Martin Leung
 */
#pragma once
#include <string>

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
    double *lf, *lu, *lu2, *lddu, *lr;
    
    JacobiLocalState(int Nx, int Ny, int Nz, int lNx, int lNy, int lNz);
    JacobiLocalState() = delete;
    ~JacobiLocalState();

    inline int uIndex(int i, int j, int k) const { return (k * lNy + j) * lNx + i; }
    inline int u2rIndex(int i, int j, int k) const { return ((k-1) * lny + (j-1)) * lnx + (i-1); }

    double get_residualsquared() const;
    void set_u_boundary(double *plane, Direction plane_dir);
    void jacobi_step();
    int get_planesize(Direction plane_dir) const;
    void get_u_boundary(double *&plane, Direction plane_dir) const;
    void get_boundary_constants(Direction plane_dir, bool is_sending, int &offset, int &M, int &N, int &iMul, int &jMul) const;
    void apply_test_conditions(int test, int start_i, int start_j, int start_k);
    void read_set_forcing(std::string filename, int start_i, int start_j, int start_k);
    void pack_solution_u2();
    // void write_solution(std::string filename) const;
};

