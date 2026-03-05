/**
 * @file poisson.cpp
 * @author Martin Leung
 */
#include <iostream>
#include <string>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "serial_solver.h"
#include "filemanager.h"
#include <boost/timer/timer.hpp>


int main(int argc, char* argv[]){
    boost::timer::auto_cpu_timer t;
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
    int test          = vm["test"].as<int>();
    const double epsilon    = vm["epsilon"].as<double>();

    int Nx            = vm["Nx"].as<int>();
    int Ny            = vm["Ny"].as<int>();
    int Nz            = vm["Nz"].as<int>();

    double* f = nullptr;

    if (vm.count("forcing")){
        // if (vm.count("test")){
        //     std::cout << "Can't have both --forcing and --test at the same time" << std::endl;
        //     return 1;
        // }


        forcing = vm["forcing"].as<std::string>();
        read_forcing(forcing, Nx, Ny, Nz, f);
        test = 0;
    }
    write_sample_forcing(32, 32, 32, "testcase2forcing.txt");

    SerialSolver ss = SerialSolver(Nx, Ny, Nz, test, epsilon, f);
    ss.run_solver();
    write_solution(ss, "solution.txt");

    // std::cout << "Residual: " << ss.get_residual() << std::endl;
    return 0;
}


