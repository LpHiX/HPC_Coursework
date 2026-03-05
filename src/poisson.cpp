/**
 * @file poisson.cpp
 * @author Martin Leung
 */
#include <iostream>
#include <boost/program_options.hpp>
using namespace std;
namespace po = boost::program_options;

#include "sayhi.h"

#define F77NAME(x) x##_
extern "C"
{
    double F77NAME(dnrm2)(
        const int &n, 
        double *x, 
        const int &incx);
    double F77NAME(idamax)(
        const int &n, 
        double *x, 
        const int &incx);
}

int main(int argc, char* argv[]){
    po::options_description opts("Availble options");
    opts.add_options()
        ("help", "Print available options")
        ("forcing", po::value<string>(), "Input forcing file")
        ("test", po::value<int>()->default_value(1), "Test case to use (1-5)")
        ("Nx", po::value<int>()->default_value(32), "Number of grid points (x)")
        ("Ny", po::value<int>()->default_value(32), "Number of grid points (y)")
        ("Nz", po::value<int>()->default_value(32), "Number of grid points (z)")
        ("epsilon", po::value<double>()->default_value(1e-8), "Residual threshold");
        

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << opts << endl;
        return 1;
    }

    string forcing;

    // Not implemented yet:
    [[maybe_unused]] const int test          = vm["test"].as<int>();
    [[maybe_unused]] const double epsilon    = vm["epsilon"].as<double>();

    const int Nx            = vm["Nx"].as<int>();
    const int Ny            = vm["Ny"].as<int>();
    const int Nz            = vm["Nz"].as<int>();
    const int nx = Nx - 2;
    const int ny = Ny - 2;
    const int nz = Nz - 2;

    if (vm.count("forcing")){
        forcing = vm["forcing"].as<string>();
    }

    
    double hx = 1.0 / (Nx - 1);
    double hy = 1.0 / (Ny - 1);
    double hz = 1.0 / (Nz - 1);

    double* u   = new double[nx * ny * nz];
    double* ddu = new double[nx * ny * nz];
    double* f   = new double[nx * ny * nz];
    double* r   = new double[nx * ny * nz];


    for (int i = 0; i < Nx; i++){
        double x = i * hx;
        for (int j = 0; j < Ny; j++){
            double y = j * hy;
            for (int k = 0; k < Nz; k++){
                double z = k * hz;
                u[(i * Nx + j) * Ny + k] = x*x + y*y + z*z;
            }
        }
    }
    cout << "Initialization Complete." << endl;

    for (int i = 1; i < Nx-1; i++){
        for (int j = 1; j < Ny-1; j++){
            for (int k = 1; k < Nz-1; k++){
                int reduced_index = ((i - 1) * nx + (j - 1)) * ny + (k - 1);
                f[reduced_index] = 6;
                ddu[reduced_index] =(
                    +     (u[((i + 1) * Nx + (j    )) * Ny + (k    )])
                    - 2 * (u[((i    ) * Nx + (j    )) * Ny + (k    )])
                    +     (u[((i - 1) * Nx + (j    )) * Ny + (k    )])) / (hx*hx) + (
                    +     (u[((i    ) * Nx + (j + 1)) * Ny + (k    )])
                    - 2 * (u[((i    ) * Nx + (j    )) * Ny + (k    )])
                    +     (u[((i    ) * Nx + (j - 1)) * Ny + (k    )])) / (hy*hy) + (
                    +     (u[((i    ) * Nx + (j    )) * Ny + (k + 1)])
                    - 2 * (u[((i    ) * Nx + (j    )) * Ny + (k    )])
                    +     (u[((i    ) * Nx + (j    )) * Ny + (k - 1)])) / (hz*hz);
                r[reduced_index] = f[reduced_index] - ddu[reduced_index];
            }
        }
    }

    // cout << u[nx * ny * nz - 1] << endl;
    // cout << ddu[0] << endl;
    // cout << ddu[5000] << endl;
    // cout << f[0] << endl;
    // cout << f[5000] << endl;
    // cout << r[0] << endl;
    // cout << r[5000] << endl;
    // int maxarg = F77NAME(idamax)((nx-2) * (ny-2) * (nz-2), r, 1);
    // cout << "Max arg: " << maxarg << endl;
    // cout << "Max r: " << r[maxarg] << endl;

    cout << "ddu and r Calculated." << endl;

    cout << "Residual: " << F77NAME(dnrm2)(nx * ny * nz - 1, r, 1) << endl;

    cout << "End of program." << endl;

    delete[] u;
    delete[] ddu;
    delete[] f;
    delete[] r;

    sayHi();
    return 0;
}


