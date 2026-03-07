/**
 * @file serial_solver.cpp
 * @author Martin Leung
 */
#include <iostream>
#include "serial_solver.h"
#include <cmath>
#include <iomanip>
#include <boost/timer/timer.hpp>

SerialSolver::SerialSolver(int Nx, int Ny, int Nz, double epsilon): 
    Nx(Nx),
    Ny(Ny),
    Nz(Nz),
    epsilon(epsilon)
    {
        localstate = std::make_unique<JacobiLocalState>(Nx, Ny, Nz, Nx, Ny, Nz);
    }
    


int SerialSolver::solve(){
    boost::timer::cpu_timer wall_timer;
    double last_report_time = -2;

    int iterations = 0;
    double residual = get_residual();

    while (residual > epsilon){
        localstate->jacobi_step();
        
        if (iterations % 100 == 0){
            residual = get_residual();
            double current_time = wall_timer.elapsed().wall / 1e9; // Convert nanoseconds to seconds        
            if (current_time - last_report_time >= 2.0) {
                std::cout 
                << std::setw(15) << "Step: " << std::setw(10) << iterations 
                << std::setw(15) << "Time: " << std::setw(10) << current_time << "s"
                << std::setw(15) << "Residual: " << std::setw(15) << residual << std::endl;
                last_report_time = current_time;
            }
        }
        iterations++;
    }
    std::cout << "Converged in " << iterations << " iterations, residual: " << residual << std::endl;
    return 0;
}

double SerialSolver::get_residual(){ return sqrt(localstate->get_residualsquared());}

void SerialSolver::initialize_test(int test){
    localstate->apply_test_conditions(test, 0, 0, 0);
}