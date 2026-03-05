/**
 * @file filemanager.cpp
 * @author Martin Leung
 */
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "filemanager.h"
#include <cmath>
#include <limits>

void write_solution(SerialSolver &ss, std::string filename){
    std::ofstream fileOutput(filename, std::ios::out | std::ios::trunc);
    fileOutput.precision(std::numeric_limits<double>::max_digits10);
    fileOutput << ss.Nx << " " << ss.Ny << " " << ss.Nz << std::endl;

    const double *u = ss.get_u();

    for (int i = 0; i < ss.Nx; i++){
        double x = i * ss.hx;
        for (int j = 0; j < ss.Ny; j++){
            double y = j * ss.hy;
            for (int k = 0; k < ss.Nz; k++){
                double z = k * ss.hz;
                fileOutput << x << " " << y << " " << z << " " << u[ss.uIndex(i,j,k)] << std::endl;
            }
        }
    }
    fileOutput.close();
}

void write_sample_forcing(const int Nx, const int Ny, const int Nz, std::string filename){
    std::ofstream fileOutput(filename, std::ios::out | std::ios::trunc);
    // fileOutput.precision(7);
    fileOutput << Nx << " " << Ny << " " << Nz << std::endl;
    fileOutput.precision(std::numeric_limits<double>::max_digits10);

    const double hx = 1.0 / (Nx - 1);
    const double hy = 1.0 / (Ny - 1);
    const double hz = 1.0 / (Nz - 1);

    for (int i = 1; i < Nx-1; i++){
        double x = i * hx;
        for (int j = 1; j < Ny-1; j++){
            double y = j * hy;
            for (int k = 1; k < Nz-1; k++){
                double z = k * hz;
                fileOutput << x << " " << y << " " << z << " " << - 3 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z) << std::endl;
            }
        }
    }
    fileOutput.close();
}

void read_forcing(std::string filename, int &Nx, int &Ny, int &Nz, double *&f){
    if (f != nullptr){
        throw "Tried to read forcing array into a initialized array";
    }

    std::ifstream fileInput(filename);
    if (fileInput.good()){
        std::string line;

        std::getline(fileInput, line);
        
        std::string Nxstr, Nystr, Nzstr;
        std::stringstream s(line);
        std::getline(s, Nxstr, ' ');
        std::getline(s, Nystr, ' ');
        std::getline(s, Nzstr, ' ');
        Nx = std::stoi(Nxstr);
        Ny = std::stoi(Nystr);
        Nz = std::stoi(Nzstr);
        const int nx = Nx - 2;
        const int ny = Ny - 2;
        const int nz = Nz - 2;
        const double hx = 1.0 / (Nx - 1);
        const double hy = 1.0 / (Ny - 1);
        const double hz = 1.0 / (Nz - 1);

        std::cout << "Read file with Nx:" << Nx << ", Ny:" << Ny << ", Nz:" << Nz << std::endl;
        f = new double[nx*ny*nz];

        while (std::getline(fileInput, line)){
            std::stringstream s(line);
            std::string xstr, ystr, zstr, fstr;
            std::getline(s, xstr, ' ');
            std::getline(s, ystr, ' ');
            std::getline(s, zstr, ' ');
            std::getline(s, fstr, ' ');

            double x,y,z,f_value;
            x = std::stod(xstr);
            y = std::stod(ystr);
            z = std::stod(zstr);
            f_value = std::stod(fstr);

            int i = std::round(x / hx);
            int j = std::round(y / hy);
            int k = std::round(z / hz);

            int index = ((i-1) * ny + (j-1)) * nz + (k-1);
            f[index] = f_value;

            // if (index % 10000){std::cout << x << " " << y << " " << z << " " << f_value << std::endl;}
        }
    }
    fileInput.close();
}