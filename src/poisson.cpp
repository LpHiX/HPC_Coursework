/**
 * @file poisson.cpp
 * @author Martin Leung
 */
#include <iostream>
#include <string>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "serial_solver.h"



int main(int argc, char* argv[]){
    po::options_description opts("Availble options");
    opts.add_options()
        ("help", "Print available options")
        ("forcing", po::value<std::string>(), "Input forcing file")
        ("test", po::value<int>()->default_value(1), "Test case to use (1-5)")
        ("Nx", po::value<int>()->default_value(32), "Number of grid points (x)")
        ("Ny", po::value<int>()->default_value(32), "Number of grid points (y)")
        ("Nz", po::value<int>()->default_value(32), "Number of grid points (z)")
        ("epsilon", po::value<double>()->default_value(1e-8), "Residual threshold");
        

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << opts << std::endl;
        return 1;
    }

    std::string forcing;

    // Not implemented yet:
    const int test          = vm["test"].as<int>();
    const double epsilon    = vm["epsilon"].as<double>();

    const int Nx            = vm["Nx"].as<int>();
    const int Ny            = vm["Ny"].as<int>();
    const int Nz            = vm["Nz"].as<int>();

    if (vm.count("forcing")){
        forcing = vm["forcing"].as<std::string>();
    }

    
    SerialSolver ss = SerialSolver(Nx, Ny, Nz, test, epsilon);
    std::cout << "Residual: " << ss.run_solver() << std::endl;
    return 0;
}


