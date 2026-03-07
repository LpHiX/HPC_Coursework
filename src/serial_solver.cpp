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

        double* lf = new double[Nx*Ny*Nz];   
        for (int i = 0; i < Nx*Ny*Nz; i++){lf[i]=6;}

        localstate = std::make_unique<JacobiLocalState>(Nx, Ny, Nz, Nx, Ny, Nz, lf);
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

// void SerialSolver::initialize_test(int test){
//     switch (test ) {
//         case 0:
//             break;
//         case 1:
//             // for (int i = 0; i < Nx; i++){
//             //     double x = i * hx;
//             //     for (int j = 0; j < Ny; j++){
//             //         double y = j * hy;
//             //         for (int k = 0; k < Nz; k++){
//             //             double z = k * hz;
//             //             double scalar = x*x + y*y + z*z;
//             //             u[uIndex(i, j, k)] = scalar;
//             //             u2[uIndex(i, j, k)] = scalar;
//             //         }
//             //     }
//             // }
//             for (int i = 1; i < Nx-1; i++){
//                 for (int j = 1; j < Ny-1; j++){
//                     for (int k = 1; k < Nz-1; k++){
//                         f[u2rIndex(i, j, k)] = 6;
//                         // u[uIndex(i,j,k)] = 0;
//                     }
//                 }
//             }
//             break;
//         case 2:
//             for (int i = 0; i < Nx; i++){
//                 double x = i * hx;
//                 for (int j = 0; j < Ny; j++){
//                     double y = j * hy;
//                     for (int k = 0; k < Nz; k++){
//                         double z = k * hz;
//                         double scalar = sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z);
//                         u[uIndex(i, j, k)] = scalar;
//                         u2[uIndex(i, j, k)] = scalar;
//                     }
//                 }
//             }
//             for (int i = 1; i < Nx-1; i++){
//                 double x = i * hx;
//                 for (int j = 1; j < Ny-1; j++){
//                     double y = j * hy;
//                     for (int k = 1; k < Nz-1; k++){
//                         double z = k * hz;
//                         f[u2rIndex(i, j, k)] = - 3 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z);
//                         u[uIndex(i,j,k)] = 0;
//                     }
//                 }
//             }
//             break;
//         case 3:
//             for (int i = 0; i < Nx; i++){
//                 double x = i * hx;
//                 for (int j = 0; j < Ny; j++){
//                     double y = j * hy;
//                     for (int k = 0; k < Nz; k++){
//                         double z = k * hz;
//                         double scalar = sin(M_PI * x) * sin(4 * M_PI * y) * sin(8 * M_PI * z);
//                         u[uIndex(i, j, k)] = scalar;
//                         u2[uIndex(i, j, k)] = scalar;
//                     }
//                 }
//             }
//             for (int i = 1; i < Nx-1; i++){
//                 double x = i * hx;
//                 for (int j = 1; j < Ny-1; j++){
//                     double y = j * hy;
//                     for (int k = 1; k < Nz-1; k++){
//                         double z = k * hz;
//                         f[u2rIndex(i,j,k)] = - 81 * M_PI * M_PI * sin(M_PI * x) * sin(4 * M_PI * y) * sin(8 * M_PI * z);
//                         u[uIndex(i,j,k)] = 0;
//                     }
//                 }
//             }
//             break;
//         case 4:
//             for (int i = 1; i < Nx-1; i++){
//                 double x = i * hx;
//                 for (int j = 1; j < Ny-1; j++){
//                     double y = j * hy;
//                     for (int k = 1; k < Nz-1; k++){
//                         double z = k * hz;
//                         f[u2rIndex(i,j,k)] = 100 * exp(-100 * ((x - 0.5)*(x - 0.5)+(y - 0.5)*(y - 0.5)+(z - 0.5)*(z - 0.5)));
//                     }
//                 }
//             }
//             break;
//         case 5:
//             for (int i = 1; i < Nx-1; i++){
//                 double x = i * hx;
//                 for (int j = 1; j < Ny-1; j++){
//                     for (int k = 1; k < Nz-1; k++){
//                         if(x < 0.5) {f[u2rIndex(i,j,k)] = -1;}
//                         else {f[u2rIndex(i,j,k)] = 1;}
//                     }
//                 }
//             }
//             break;
//     }
// }

// SerialSolver::~SerialSolver(){
//     delete[] u;
//     delete[] u2;
//     delete[] ddu;
//     delete[] f;
//     delete[] r;
// }