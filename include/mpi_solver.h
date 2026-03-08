/**
 * @file mpi_solver.h
 * @author Martin Leung
 */

#pragma once

#include "jacobi_localstate.h"
#include <string>
#include <memory>
#include <mpi.h>

class MPISolver{
public:
    const int Nx, Ny, Nz, nx, ny, nz, Px, Py, Pz;
    const double epsilon;

    MPISolver(int Nx, int Ny, int Nz, int Px, int Py, int Pz, double epsilon);
    MPISolver() = delete;
    ~MPISolver();

    void solve();
    void initialize(int test);
    void initialize(std::string filename);
    void write_solution(std::string filename, int test);
    double get_residual() const;
private:
    int lNx, lNy, lNz, lnx, lny, lnz;
    MPI_Comm gridcomm;
    int gridrank;
    int neighbours_ranks[6];
    double* send_buffers[6] = {nullptr};
    double* recv_buffers[6] = {nullptr};
    std::unique_ptr<JacobiLocalState> localstate;

    // void exchange();
    int start_i, start_j, start_k;
};

inline std::string DirectionNames[6] = {"POS_X","NEG_X","POS_Y","NEG_Y","POS_Z","NEG_Z"};