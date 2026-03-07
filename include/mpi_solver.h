/**
 * @file mpi_solver.h
 * @author Martin Leung
 */

#pragma once

#include "jacobi_localstate.h"
#include <string>
#include <memory>

class MPISolver{
public:
    const int Nx, Ny, Nz, Px, Py, Pz;
    const double epsilon;

    MPISolver(int Nx, int Ny, int Nz, int Px, int Py, int Pz, double epsilon);
    MPISolver() = delete;
    ~MPISolver();

    void solve();
private:
    int neighbours_ranks[6];
    double* transfer_plane_x;
    double* transfer_plane_y;
    double* transfer_plane_z;
    std::unique_ptr<JacobiLocalState> localstate;
};

inline std::string DirectionNames[6] = {"POS_X","NEG_X","POS_Y","NEG_Y","POS_Z","NEG_Z"};