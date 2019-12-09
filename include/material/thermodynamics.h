#ifndef MPM_MATERIAL_THERMODYNAMICS_H_
#define MPM_MATERIAL_THERMODYNAMICS_H_

#include <iostream>
#include <limits>

#include <cmath>

#include "Eigen/Dense"

#include "material.h"

namespace mpm {

//! Thermodynamics class
//! \brief Mohr Coulomb material model
//! \details Mohr Coulomb material model with softening
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Thermodynamics : public Material<Tdim> {
 public:
  //! Define a vector of 6 dof
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  //! Define a Matrix of 6 x 6
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

  //! Failure state
  enum FailureState { Elastic = 0, Yield = 1 };

  //! Constructor with id and material properties
  //! \param[in] material_properties Material properties
  Thermodynamics(unsigned id, const Json& material_properties);

  //! Destructor
  ~Thermodynamics() override{};

  //! Delete copy constructor
  Thermodynamics(const Thermodynamics&) = delete;

  //! Delete assignement operator
  Thermodynamics& operator=(const Thermodynamics&) = delete;

  //! Initialise history variables
  //! \retval state_vars State variables with history
  mpm::dense_map initialise_state_variables() override;

  //! Thermodynamic pressure
  //! \param[in] volumetric_strain dVolumetric_strain
  //! \retval pressure Pressure for volumetric strain
  double thermodynamic_pressure(double volumetric_strain) const override {
    return 0;
  };

  //! Compute stress
  //! \param[in] stress Stress
  //! \param[in] strain Strain
  //! \param[in] particle Constant point to particle base
  //! \param[in] state_vars History-dependent state variables
  //! \retval updated_stress Updated value of stress
  Vector6d compute_stress(const Vector6d& stress, const Vector6d& strain,
                          const ParticleBase<Tdim>* ptr,
                          mpm::dense_map* state_vars) override;

  //! Compute stress invariants (p, q, etc)
  //! \param[in] stress Stress
  //! \param[in] state_vars History-dependent state variables
  //! \retval status of computation of stress invariants
  bool compute_strain_invariants(const Vector6d& strain,
                                 mpm::dense_map* state_vars);

  //! Compute state variables (void ratio, p_image, e_image, etc)
  //! \param[in] stress Stress
  //! \param[in] state_vars History-dependent state variables
  //! \retval status of computation of stress invariants
  bool compute_state_variables(const Vector6d& stress, const Vector6d& dstrain,
                               mpm::dense_map* state_vars);

 protected:
  //! material id
  using Material<Tdim>::id_;
  //! Material properties
  using Material<Tdim>::properties_;
  //! Logger
  using Material<Tdim>::console_;

 private:
  //! Elastic stiffness matrix
  Matrix6x6 de_;
  //! Plastic stiffness matrix
  Matrix6x6 dp_;
  //! Density
  double density_{std::numeric_limits<double>::max()};
  //! Density dependency index B1
  double B1_{std::numeric_limits<double>::max()};
  //! Elastic bulk modulus B0
  double elastic_bulk_modulus_{std::numeric_limits<double>::max()};
  //! Xi, twice the ratio of shear to bulk modulus
  double xi_modulus_ratio_{std::numeric_limits<double>::max()};
  //! Cohesion related parameter
  double cohesion_{std::numeric_limits<double>::max()};
  //! Non-linear index m
  double m_{std::numeric_limits<double>::max()};
  //! Stress induced anisotropy index
  double eta_{std::numeric_limits<double>::max()};
  //! Liquefaction control b
  double b_{std::numeric_limits<double>::max()};
  //! Coefficient for granular kinetic fluctuation induced by volumetric
  //! deformation
  double mv_{std::numeric_limits<double>::max()};
  //! Coefficient for granular kinetic fluctuation induced by shear deformation
  double ms_{std::numeric_limits<double>::max()};
  //! Scale coefficient for volumetric dissipation
  double lambdav_{std::numeric_limits<double>::max()};
  //! Slope of critical state line parameter
  double alpha_{std::numeric_limits<double>::max()};
  //! Maximum value parameter for residual strain
  double h_{std::numeric_limits<double>::max()};

};  // Thermodynamics class
}  // namespace mpm

#include "thermodynamics.tcc"

#endif  // MPM_MATERIAL_THERMODYNAMICS_H_
