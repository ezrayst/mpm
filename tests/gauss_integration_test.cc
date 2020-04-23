#include <cmath>
#include <limits>
#include <memory>
#include <vector>

#include "Eigen/Dense"
#include "catch.hpp"

#include "cell.h"
#include "element.h"
#include "factory.h"
#include "geometry.h"
#include "hexahedron_element.h"
#include "hexahedron_quadrature.h"
#include "node.h"
#include "quadrilateral_element.h"
#include "quadrilateral_quadrature.h"

//! \brief Check Gauss quadrature rule for 2d Quadrilateral
TEST_CASE("Random gauss integration is checked for 2D Quadrilateral",
          "[RandomGaussQuad][2D]") {

  // Dimension
  const unsigned Dim = 2;
  // Degrees of freedom
  const unsigned Dof = 2;
  // Number of phases
  const unsigned Nphases = 1;
  // Number of nodes per cell
  const unsigned Nnodes = 4;
  // Tolerance
  const double Tolerance = 1.E-9;

  // Shape function
  // 4-noded quadrilateral shape functions
  std::shared_ptr<mpm::Element<Dim>> element =
      Factory<mpm::Element<Dim>>::instance()->create("ED2Q4");

  Eigen::Vector2d coords;
  coords.setZero();

  coords << -1.0, -1.0;
  std::shared_ptr<mpm::NodeBase<Dim>> node0 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

  coords << 1.0, -1.0;
  std::shared_ptr<mpm::NodeBase<Dim>> node1 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

  coords << 1.0, 1.0;
  std::shared_ptr<mpm::NodeBase<Dim>> node2 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

  coords << -1.0, 1.0;
  std::shared_ptr<mpm::NodeBase<Dim>> node3 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

  // Create cell
  mpm::Index id = 0;
  auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element, false);
  // Add nodes
  cell->add_node(0, node0);
  cell->add_node(1, node1);
  cell->add_node(2, node2);
  cell->add_node(3, node3);

  SECTION("polynomial order of one") {
    // For polynomail nrder of one, the minimum number of quadrature points is
    // four.
    const unsigned porder = 1;
    const unsigned nterms = 4;

    const double val_1_by_sqrt3 = 1. / std::sqrt(3.);
    Eigen::Vector2d qpoint0;
    qpoint0 << -val_1_by_sqrt3, -val_1_by_sqrt3;
    Eigen::Vector2d qpoint1;
    qpoint1 << val_1_by_sqrt3, -val_1_by_sqrt3;
    Eigen::Vector2d qpoint2;
    qpoint2 << val_1_by_sqrt3, val_1_by_sqrt3;
    Eigen::Vector2d qpoint3;
    qpoint3 << -val_1_by_sqrt3, val_1_by_sqrt3;

    std::vector<Eigen::Vector2d> qpoints;
    qpoints.emplace_back(qpoint0);
    qpoints.emplace_back(qpoint1);
    qpoints.emplace_back(qpoint2);
    qpoints.emplace_back(qpoint3);

    std::shared_ptr<mpm::Quadrature<2>> quadrature =
        std::make_shared<mpm::QuadrilateralQuadrature<2, 4>>();
    auto gauss_weights = quadrature->compute_gauss_weights(qpoints, porder);

    // Check the size of gauss weights vector
    REQUIRE(gauss_weights.rows() == 4);
    // Check gauss weights
    REQUIRE(gauss_weights(0) == Approx(1.0).epsilon(Tolerance));
    REQUIRE(gauss_weights(1) == Approx(1.0).epsilon(Tolerance));
    REQUIRE(gauss_weights(2) == Approx(1.0).epsilon(Tolerance));
    REQUIRE(gauss_weights(3) == Approx(1.0).epsilon(Tolerance));
  }
}

//! \brief Check Gauss quadrature rule for 2d Triangle
TEST_CASE("Random gauss integration is checked for 2D Triangle",
          "[RandomGaussTri][2D]") {

  // Dimension
  const unsigned Dim = 2;
  // Degrees of freedom
  const unsigned Dof = 2;
  // Number of phases
  const unsigned Nphases = 1;
  // Number of nodes per cell
  const unsigned Nnodes = 4;
  // Tolerance
  const double Tolerance = 1.E-9;

  // Create a triangle
  Eigen::Matrix<double, 2, 1> v1, v2, v3;
  v1 << 0.5, 0.5;
  v2 << 1.2, 1.1;
  v3 << 0.9, 0.8;
  std::vector<Eigen::Matrix<double, 2, 1>> tri_vertices{v1, v2, v3};

  // Generate random points in the triangle
  unsigned numpoints = 12;
  std::vector<Eigen::Matrix<double, 2, 1>> points;
  mpm::geometry::shape3verticesface::generate_random_points(tri_vertices,
                                                            numpoints, &points);

  // Compute reference coordinates of points
  std::vector<Eigen::Matrix<double, 2, 1>> xi;
  std::for_each(
      points.begin(), points.end(),
      [=, &xi](const Eigen::Matrix<double, 2, 1>& point) {
        xi.emplace_back(
            mpm::geometry::shape3verticesface::reference_coordinates<Dim>(
                tri_vertices, point));
      });

  // Compute gauss gauss
  const unsigned porder = 2;
  auto gauss_weights =
      mpm::geometry::shape3verticesface::compute_gauss_weights(xi, porder);

  const unsigned N_nodes = 3;
  Eigen::VectorXd jacobian_determinant(xi.size());
  int p = 0;
  for (auto& local_coord : xi) {
    jacobian_determinant(p) = std::fabs(
        (mpm::geometry::shape3verticesface::jacobian_matrix<Dim, N_nodes>(
             tri_vertices, local_coord))
            .determinant());
    ++p;
  }

  // Compute triangle area (0th order polynomial)
  double area = (jacobian_determinant.cwiseProduct(gauss_weights)).sum();

  // Get the triangle area
  double tri_area = mpm::geometry::shape3verticesface::volume(tri_vertices);
  REQUIRE(area == Approx(tri_area).epsilon(Tolerance));
}