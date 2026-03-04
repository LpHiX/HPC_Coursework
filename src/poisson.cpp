/**
 * @file poisson.cpp
 * @author Martin Leung
 */
#include <iostream>
#include <boost/program_options.hpp>
using namespace std;
namespace po = boost::program_options;

#include "sayhi.h"


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
        

    sayHi();
    return 0;
}


