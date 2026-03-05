/**
 * @file serial_solver.h
 * @author Martin Leung
 */
#pragma once

class SerialSolver{
public:
    const int Nx, Ny, Nz, nx, ny, nz, test, max_iter;
    const double epsilon, hx, hy, hz, hx2, hy2, hz2, j_coeff;
    
    SerialSolver(int Nx, int Ny, int Nz, int test, double epsilon);
    SerialSolver() = delete;
    ~SerialSolver();

    const double* get_u() const { return u; }
    const double* get_ddu() const { return ddu; }
    const double* get_f() const { return f; }
    const double* get_r() const { return r; }
    double get_residual();
    int run_solver();
private:
    double *u, *u2, *ddu, *f, *r;
    void initialize_test(int test);
    int uIndex(int i,int j, int k);
    int u2rIndex(int i,int j, int k);

};