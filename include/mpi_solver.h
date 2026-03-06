/**
 * @file mpi_solver.h
 * @author Martin Leung
 */

#pragma once

class MPISolver{
public:
    MPISolver(int Nx, int Ny, int Nz, Px, Py, Pz, rank);
    const int Nx, Ny, Nz, Px, Py, Pz, rank;
}