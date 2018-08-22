#include <cmath>
#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "geometry.h"

//! \brief Check geometry class for 2D case
TEST_CASE("Geometry is checked for 2D case", "[geometry][2D]") {

  // Dimension
  const unsigned Dim = 2;
  // Tolerance
  const double Tolerance = 1.E-7;

  SECTION("Check inverse rotation matrix") {

    auto geometry = std::make_shared<mpm::Geometry<Dim>>();

    Eigen::Matrix<double, 2, 1> angles;
    // clang-format off
    angles << 45 * M_PI / 180,   // alpha
              30 * M_PI / 180;   // beta

    Eigen::Matrix<double, 2, 2> inverse_rotation_matrix;
    // clang-format off
    inverse_rotation_matrix <<  0.258819045102521, 0.965925826289068,
                               -0.965925826289068, 0.258819045102521;
    // clang-format on
    auto check_inverse_rotation_matrix =
        geometry->compute_inverse_rotation_matrix(angles);
    REQUIRE(check_inverse_rotation_matrix.cols() == 2);
    REQUIRE(check_inverse_rotation_matrix.rows() == 2);
    REQUIRE(check_inverse_rotation_matrix(0, 0) ==
            Approx(inverse_rotation_matrix(0, 0)).epsilon(Tolerance));
    REQUIRE(check_inverse_rotation_matrix(1, 0) ==
            Approx(inverse_rotation_matrix(1, 0)).epsilon(Tolerance));
    REQUIRE(check_inverse_rotation_matrix(0, 1) ==
            Approx(inverse_rotation_matrix(0, 1)).epsilon(Tolerance));
    REQUIRE(check_inverse_rotation_matrix(1, 1) ==
            Approx(inverse_rotation_matrix(1, 1)).epsilon(Tolerance));
  }
}

//! \brief Check cell class for 3D case
TEST_CASE("Geometry is checked for 3D case", "[geometry][3D]") {

  // Dimension
  const unsigned Dim = 3;
  // Tolerance
  const double Tolerance = 1.E-7;

  SECTION("Check inverse rotation matrix") {

    auto geometry = std::make_shared<mpm::Geometry<Dim>>();

    Eigen::Matrix<double, 3, 1> angles;
    // clang-format off
    angles << 45 * M_PI / 180,   // alpha 
              30 * M_PI / 180,   // beta
              60 * M_PI / 180;   // gamma
    // clang-format on

    Eigen::Matrix<double, 3, 3> inverse_rotation_matrix;
    // clang-format off
    inverse_rotation_matrix <<  0.435595740399158,  0.789149130992431,  0.433012701892219,
                               -0.659739608441171, -0.047367172745376,  0.75,
                                0.612372435695794, -0.612372435695795,  0.5;
    // clang-format on
    auto check_inverse_rotation_matrix =
        geometry->compute_inverse_rotation_matrix(angles);
    REQUIRE(check_inverse_rotation_matrix.cols() == 3);
    REQUIRE(check_inverse_rotation_matrix.rows() == 3);
    REQUIRE(check_inverse_rotation_matrix(0, 0) ==
            Approx(inverse_rotation_matrix(0, 0)).epsilon(Tolerance));
    REQUIRE(check_inverse_rotation_matrix(0, 1) ==
            Approx(inverse_rotation_matrix(0, 1)).epsilon(Tolerance));
    REQUIRE(check_inverse_rotation_matrix(0, 2) ==
            Approx(inverse_rotation_matrix(0, 2)).epsilon(Tolerance));
    REQUIRE(check_inverse_rotation_matrix(1, 0) ==
            Approx(inverse_rotation_matrix(1, 0)).epsilon(Tolerance));
    REQUIRE(check_inverse_rotation_matrix(1, 1) ==
            Approx(inverse_rotation_matrix(1, 1)).epsilon(Tolerance));
    REQUIRE(check_inverse_rotation_matrix(1, 2) ==
            Approx(inverse_rotation_matrix(1, 2)).epsilon(Tolerance));
    REQUIRE(check_inverse_rotation_matrix(2, 0) ==
            Approx(inverse_rotation_matrix(2, 0)).epsilon(Tolerance));
    REQUIRE(check_inverse_rotation_matrix(2, 1) ==
            Approx(inverse_rotation_matrix(2, 1)).epsilon(Tolerance));
    REQUIRE(check_inverse_rotation_matrix(2, 2) ==
            Approx(inverse_rotation_matrix(2, 2)).epsilon(Tolerance));
  }
}
