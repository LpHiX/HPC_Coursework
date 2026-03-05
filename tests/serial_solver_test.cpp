#include "serial_solver.h"

#define BOOST_TEST_MODULE SerialSolverTest
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE ( ResidualTest ){
    SerialSolver ss(32, 32, 32, 1, 1e-7);
    BOOST_CHECK_SMALL ( ss.run_solver(), 1e-7);
}