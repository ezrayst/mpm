#ifndef MPM_TRIANGLE_QUADRATURE_H_
#define MPM_TRIANGLE_QUADRATURE_H_

#include <exception>

#include <Eigen/Dense>

#include "quadrature.h"

//! MPM namespace
namespace mpm {

// Triangle quadrature class derived from Quadrature base class
//! \brief Quadrature (gauss points) for a triangle element
//! \tparam Tdim Dimension
//! \tparam Tnquadratures number of quadratures
template <unsigned Tdim, unsigned Tnquadratures>
class TriangleQuadrature : public Quadrature<Tdim> {

 public:
  TriangleQuadrature() : Quadrature<Tdim>() {
    static_assert(Tdim == 2, "Invalid dimension for a triangle element");
    static_assert(((Tnquadratures == 1) || (Tnquadratures == 3)),
                  "Invalid number of quadratures");
  }

  //! Return quadrature points
  //! \param[out] qpoints Quadrature points in local coordinates
  Eigen::MatrixXd quadratures() const override;

  //! Return weights
  //! \param[out] weights Weights for quadrature points
  Eigen::VectorXd weights() const override;

  //! Compute gauss weights
  //! \param[in] xi Reference coordinates of quadrature points
  //! \param[in] porder Order of integral polynomial
  Eigen::VectorXd compute_gauss_weights(
      const std::vector<Eigen::Matrix<double, Tdim, 1>>& xi,
      const unsigned porder) const override;
};

}  // namespace mpm

#include "triangle_quadrature.tcc"

#endif  // MPM_TRIANGLE_QUADRATURE_H_
