/**
 * @file serial_solver_test.cpp
 * @author Martin Leung
 */
#include "serial_solver.h"
#include "filemanager.h"
#include <iostream>
#define BOOST_TEST_MODULE SerialSolverTest
#include <boost/test/included/unit_test.hpp>

#define F77NAME(x) x##_
extern "C"
{
    int F77NAME(idamax)(
        const int &n, 
        const double *x, 
        const int &incx);
    double F77NAME(dnrm2)(
        const int &n, 
        const double *x, 
        const int &incx);
    double F77NAME(dcopy)(
        const int &n, 
        const double *x, 
        const int &incx, 
        const double *y, 
        const int &incy); 
    double F77NAME(daxpy)(
        const int &n,
        const double &alpha,
        const double *x, 
        const int &incx, 
        const double *y, 
        const int &incy); 
}

void verify_test1(SerialSolver& ss){
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

    int maxarg = F77NAME(idamax)(ss.nx * ss.ny * ss.nz, r, 1);
    BOOST_CHECK_SMALL(r[maxarg-1], ss.epsilon*10);
    std::cout << "Max error: " << r[maxarg-1] << std::endl;
}

void max_error(SerialSolver& ss){
    double* error = new double[ss.Nx*ss.Ny*ss.Nz];

    const double* u2 = ss.get_u2();

    F77NAME(dcopy)(ss.Nx*ss.Ny*ss.Nz, u2, 1, error, 1);

    ss.run_solver();
    BOOST_CHECK_SMALL(ss.get_residual(), ss.epsilon*10);

    const double* u = ss.get_u();

    F77NAME(daxpy)(ss.Nx*ss.Ny*ss.Nz, -1.0, u, 1, error, 1);
    int maxarg = F77NAME(idamax)(ss.Nx * ss.Ny * ss.Nz, error, 1);
    // BOOST_CHECK_SMALL(error[maxarg- 1], ss.epsilon*10);
    std::cout << "Max error: " << error[maxarg- 1] << std::endl;

    delete[] error;
}

BOOST_AUTO_TEST_SUITE(SerialSolverSuite)

    BOOST_AUTO_TEST_CASE( TestCase1 )
{
    SerialSolver ss(32, 32, 32, 1, 1e-8);
    verify_test1(ss);
}

BOOST_AUTO_TEST_CASE(TestCase1_DifferentIndex)
{
    SerialSolver ss(64, 32, 16, 1, 1e-8);
    verify_test1(ss);
}
BOOST_AUTO_TEST_CASE( TestCase2 )
{
    SerialSolver ss(64, 64, 64, 2, 1e-8);
    max_error(ss);
}
BOOST_AUTO_TEST_CASE( TestCase3 )
{
    SerialSolver ss(32, 64, 128, 3, 1e-8);
    max_error(ss);

}
BOOST_AUTO_TEST_CASE( TestCase4 )
{
    SerialSolver ss(64, 64, 64, 4, 1e-8);
    ss.run_solver();
    BOOST_CHECK_SMALL(ss.get_residual(), ss.epsilon*10);

}
BOOST_AUTO_TEST_CASE( TestCase5 )
{
    SerialSolver ss(64, 64, 64, 5, 1e-8);
    ss.run_solver();
    BOOST_CHECK_SMALL(ss.get_residual(), ss.epsilon*10);

}

BOOST_AUTO_TEST_CASE( FileRead )
{
    SerialSolver ss = SerialSolver(32, 32, 32, 2, 1e-8);
    ss.run_solver();
    
    double* error = new double[ss.Nx*ss.Ny*ss.Nz];
    const double* u = ss.get_u();
    F77NAME(dcopy)(ss.Nx*ss.Ny*ss.Nz, u, 1, error, 1);

    double* f = nullptr;
    int a,b,c;
    read_forcing("testcase2forcing.txt", a, b, c, f);

    SerialSolver ss2 = SerialSolver(32, 32, 32, 0, 1e-8, f);
    ss2.run_solver();

    const double* u2 = ss2.get_u();

    F77NAME(daxpy)(ss.Nx*ss.Ny*ss.Nz, -1.0, u2, 1, error, 1);
    double rsdl = F77NAME(dnrm2)(32*32*32, error, 1);
    BOOST_CHECK_SMALL(rsdl, ss.epsilon*10);
    std::cout << "Difference between the two solutions: " << rsdl << std::endl;

    // int maxarg = F77NAME(idamax)(ss.Nx * ss.Ny * ss.Nz, error, 1);
    // BOOST_CHECK_SMALL(error[maxarg- 1], ss.epsilon*10);
    // std::cout << "Max error: " << error[maxarg- 1] << std::endl;

    delete[] error;
}

BOOST_AUTO_TEST_SUITE_END()