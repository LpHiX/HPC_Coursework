/**
 * @file serial_solver.h
 * @author Martin Leung
 */
#pragma once

class SerialSolver{
public:
    SerialSolver(int Nx, int Ny, int Nz, int test, double epsilon);
    SerialSolver() = delete;
    double run_solver();
private:
    int Nx, Ny, Nz, nx, ny, nz, test;
    double epsilon, hx, hy, hz, hx2, hy2, hz2;
    double *u, *ddu, *f, *r;
    void initialize_test_1();
    int uIndex(int i,int j, int k);
    int u2rIndex(int i,int j, int k);
};