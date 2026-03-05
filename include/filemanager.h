/**
 * @file filemanager.h
 * @author Martin Leung
 */
#pragma once

#include <string>
#include "serial_solver.h"

// class FileManager{
// public:
    void write_solution(SerialSolver &ss, std::string filename);
    void write_sample_forcing(const int Nx, const int Ny, const int Nz, std::string filename);
    void read_forcing(std::string filename, int &Nx, int &Ny, int &Nz, double *&f);
// };