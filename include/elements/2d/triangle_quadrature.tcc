// Getting the quadratures for Tnquadratures = 1
template <>
inline Eigen::MatrixXd mpm::TriangleQuadrature<2, 1>::quadratures() const {
  Eigen::Matrix<double, 2, 1> quadratures;
  quadratures(0, 0) = 1. / 3.;
  quadratures(1, 0) = 1. / 3.;

  return quadratures;
}

// Getting the weights for Tnquadratures = 1
template <>
inline Eigen::VectorXd mpm::TriangleQuadrature<2, 1>::weights() const {
  Eigen::VectorXd weights(1);
  weights(0) = 0.5;

  return weights;
}

// Getting the quadratures for Tnquadratures = 3
template <>
inline Eigen::MatrixXd mpm::TriangleQuadrature<2, 3>::quadratures() const {
  Eigen::Matrix<double, 2, 3> quadratures;
  const double val_1_by_sqrt3 = 1. / std::sqrt(3.);

  quadratures(0, 0) = 1. / 6.;
  quadratures(1, 0) = 1. / 6.;
  quadratures(0, 1) = 2. / 3.;
  quadratures(1, 1) = 1. / 6.;
  quadratures(0, 2) = 1. / 6.;
  quadratures(1, 2) = 2. / 3.;

  return quadratures;
}

// Getting the weights for Tnquadratures = 3
template <>
inline Eigen::VectorXd mpm::TriangleQuadrature<2, 3>::weights() const {
  Eigen::VectorXd weights(3);
  weights(0) = 1. / 3.;
  weights(1) = 1. / 3.;
  weights(2) = 1. / 3.;

  return weights;
}

// Compute Gauss weights of a given set of points
template <unsigned Tdim, unsigned Tnquadratures>
Eigen::VectorXd
    mpm::TriangleQuadrature<Tdim, Tnquadratures>::compute_gauss_weights(
        const std::vector<Eigen::Matrix<double, Tdim, 1>>& xi,
        const unsigned porder) const {

  // number of monomials
  unsigned nterms = pow((porder + 1), Tdim);
  if (nterms != Tnquadratures)
    throw std::runtime_error(
        "Invalid polynomial order, cannot compute gauss weights");

  // Check minimum number of quadrature points
  if (xi.size() < nterms)
    throw std::runtime_error(
        "Number of quadrature points is less than minimum required");

  // Fill 'A' matrix with monomials
  Eigen::MatrixXd polyA(nterms, xi.size());
  int p = 0;
  for (auto& local_coord : xi) {
    polyA.col(p) =
        mpm::Polynomial::evaluate_monomials<Tdim>(porder, local_coord);
    ++p;
  }
  // transpose of polynomial 'A' matrix
  Eigen::MatrixXd polyA_trans = polyA.transpose();

  // return gauss_weights;
  return (
      polyA_trans * ((polyA * polyA_trans).inverse()) *
      mpm::Polynomial::DefiniteIntegral<Tdim,
                                        Tnquadratures>::Tri_Definite_Integrals);
}