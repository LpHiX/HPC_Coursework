/**
 * @file serial_solver.h
 * @author Martin Leung
 */
#pragma once
#include <memory>
#include "jacobi_localstate.h"
#include <string>

class SerialSolver{
public:
    const int Nx, Ny, Nz;
    const double epsilon;
    
    SerialSolver(int Nx, int Ny, int Nz, double epsilon);
    SerialSolver() = delete;
    // ~SerialSolver();

    // const double* get_u() const { return u; }
    // const double* get_u2() const { return u2; }
    // const double* get_ddu() const { return ddu; }
    // const double* get_f() const { return f; }
    // const double* get_r() const { return r; }

    int solve();
    double get_residual();
    void initialize(int test);
    void initialize(std::string filename);
    void write_solution(std::string filename);
private:
    std::unique_ptr<JacobiLocalState> localstate;
};