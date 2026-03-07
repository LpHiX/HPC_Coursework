/**
 * @file poisson-mpi.cpp
 * @author Martin Leung
 */

#include <iostream>
#include <string>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "mpi_solver.h"
#include <boost/timer/timer.hpp>
#include <mpi.h>

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int rank, size, retval_rank, retval_size;
    retval_rank = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    retval_size = MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (retval_rank == MPI_ERR_COMM || retval_size == MPI_ERR_COMM)
    {
        std::cout << "Invalid Communicator" << std::endl;
        return 1;
    }

    po::options_description opts("Availble options");
    opts.add_options()
        ("help", "Print available options")
        ("forcing", po::value<std::string>(), "Input forcing file")
        ("test", po::value<int>()->default_value(1), "Test case to use (1-5)")
        ("Nx", po::value<int>()->default_value(32), "Number of grid points (x)")
        ("Ny", po::value<int>()->default_value(32), "Number of grid points (y)")
        ("Nz", po::value<int>()->default_value(32), "Number of grid points (z)")
        ("epsilon", po::value<double>()->default_value(1e-8), "Residual threshold")   
        ("Px", po::value<int>()->default_value(1), "Number of processes (x)")
        ("Py", po::value<int>()->default_value(1), "Number of processes (y)")
        ("Pz", po::value<int>()->default_value(1), "Number of processes (z)");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
        if (rank == 0)
        {
            std::cout << opts << std::endl;
        }
        MPI_Finalize();
        return 0;
    }

    std::string forcing;

    // Not implemented yet:
    // int test          = vm["test"].as<int>();
    const double epsilon    = vm["epsilon"].as<double>();

    int Nx = vm["Nx"].as<int>();
    int Ny = vm["Ny"].as<int>();
    int Nz = vm["Nz"].as<int>();
    int Px = vm["Px"].as<int>();
    int Py = vm["Py"].as<int>();
    int Pz = vm["Pz"].as<int>();

    if (Px * Py * Pz != size){
        if (rank == 0) {
            std::cout << "Error: Px * Py * Pz =" << Px*Py*Pz 
                      << " does not equal total MPI processes (" << size << ")." << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    if (rank == 0){
        boost::timer::auto_cpu_timer t;
    }

    // double *f = nullptr;

    if (vm.count("forcing"))
    {
        // if (vm.count("test")){
        //     std::cout << "Can't have both --forcing and --test at the same time" << std::endl;
        //     return 1;
        // }

        forcing = vm["forcing"].as<std::string>();
        // read_forcing(forcing, Nx, Ny, Nz, f);
        // test = 0;
    }

    MPISolver solver(Nx, Ny, Nz, Px, Py, Pz, epsilon);
    solver.solve();

    MPI_Finalize();
    return 0;
}
