/**
 * @file poisson.cpp
 * @author Martin Leung
 */
#include <iostream>
#include <string>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "serial_solver.h"
#include <boost/timer/timer.hpp>
#include <fstream>


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

    if (vm.count("forcing")){
        if (!vm["test"].defaulted()){
            std::cout << "Can't have both --forcing and --test at the same time" << std::endl;
            return 1;
        }
        forcing = vm["forcing"].as<std::string>();
        test = 0;

        std::ifstream fileInput(forcing);
        if (!fileInput.good()) {
            std::cout << "Error opening forcing file." << std::endl;
            return 1;
        }
        fileInput >> Nx >> Ny >> Nz;
        fileInput.close();
    }
    if (test < 0 || test > 5){
        throw std::runtime_error("(--test " + std::to_string(test) + ") not defined");
    }
    // write_sample_forcing(32, 32, 32, "testcase2forcing.txt");

    SerialSolver ss = SerialSolver(Nx, Ny, Nz, epsilon);
    if (test != 0){
        ss.initialize(test);
    } else {
        ss.initialize(forcing);
    }
    ss.solve();
    ss.write_solution("solution.txt");

    // write_solution(ss, "solution.txt");

    // std::cout << "Residual: " << ss.get_residual() << std::endl;
    return 0;
}


