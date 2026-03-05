#include "serial_solver.h"
#include <iostream>
#define BOOST_TEST_MODULE SerialSolverTest
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_SUITE(SerialSolverSuite)

    BOOST_AUTO_TEST_CASE( TestCase1 )
{
    SerialSolver ss(32, 32, 32, 1, 1e-7);
    ss.run_solver();

    const double* u = ss.get_u();
    const double* ddu = ss.get_ddu();
    const double* f = ss.get_f();
    const double* r = ss.get_r();

    BOOST_CHECK_SMALL(ss.get_residual(), ss.epsilon*10);
    BOOST_CHECK_SMALL(u[0], ss.epsilon * 10);
    BOOST_CHECK_CLOSE(u[ss.Nx * ss.Ny * ss.Nz - 1], 3, ss.epsilon*10);
    BOOST_CHECK_CLOSE(ddu[0], 6, ss.epsilon*10);
    BOOST_CHECK_CLOSE(ddu[ss.nx * ss.ny * ss.nz - 1], 6, ss.epsilon*10);
    BOOST_CHECK_CLOSE(f[0], 6, ss.epsilon*10);
    BOOST_CHECK_CLOSE(f[ss.nx * ss.ny * ss.nz - 1], 6, ss.epsilon*10);
    BOOST_CHECK_SMALL(r[0], 1e-7);
    BOOST_CHECK_SMALL(r[ss.nx * ss.ny * ss.nz - 1], ss.epsilon*10);
    
}

BOOST_AUTO_TEST_CASE(TestCase1_DifferentIndex)
{
    SerialSolver ss(64, 32, 16, 1, 1e-7);
    ss.run_solver();

    const double* u = ss.get_u();
    const double* ddu = ss.get_ddu();
    const double* f = ss.get_f();
    const double* r = ss.get_r();

    BOOST_CHECK_SMALL(ss.get_residual(), ss.epsilon*10);
    BOOST_CHECK_SMALL(u[0], ss.epsilon * 10);
    BOOST_CHECK_CLOSE(u[ss.Nx * ss.Ny * ss.Nz - 1], 3, ss.epsilon*10);
    BOOST_CHECK_CLOSE(ddu[0], 6, ss.epsilon*10);
    BOOST_CHECK_CLOSE(ddu[ss.nx * ss.ny * ss.nz - 1], 6, ss.epsilon*10);
    BOOST_CHECK_CLOSE(f[0], 6, ss.epsilon*10);
    BOOST_CHECK_CLOSE(f[ss.nx * ss.ny * ss.nz - 1], 6, ss.epsilon*10);
    BOOST_CHECK_SMALL(r[0], 1e-7);
    BOOST_CHECK_SMALL(r[ss.nx * ss.ny * ss.nz - 1], ss.epsilon*10);
}

BOOST_AUTO_TEST_SUITE_END()