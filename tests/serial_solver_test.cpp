#include "serial_solver.h"

#define BOOST_TEST_MODULE SerialSolverTest
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE ( SerialSolverSuite )

BOOST_AUTO_TEST_CASE ( ResidualTestCase1 ){
    SerialSolver ss(32, 32, 32, 1, 1e-7);
    BOOST_CHECK_SMALL ( ss.run_solver(), 1e-7);
}

BOOST_AUTO_TEST_CASE ( ResidualTestCase1_DifferentIndex ){
    SerialSolver ss(64, 32, 16, 1, 1e-7);
    BOOST_CHECK_SMALL ( ss.run_solver(), 1e-7);
}

BOOST_AUTO_TEST_SUITE_END()