//! Constructor with id and material properties
template <unsigned Tdim>
mpm::Thermodynamics<Tdim>::Thermodynamics(unsigned id,
                                          const Json& material_properties)
    : Material<Tdim>(id, material_properties) {
  try {
    // Density
    density_ = material_properties["density"].template get<double>();
    // Density dependency index B1
    B1_ = material_properties["B1"].template get<double>();
    // Elastic bulk modulus B0
    elastic_bulk_modulus_ =
        material_properties["elastic_bulk_modulus"].template get<double>();
    // Xi, twice the ratio of shear to bulk modulus
    xi_modulus_ratio_ =
        material_properties["xi_modulus_ratio"].template get<double>();
    // Cohesion related parameter
    cohesion_ = material_properties["cohesion"].template get<double>();
    // Non-linear index m
    m_ = material_properties["m"].template get<double>();
    // Stress induced anisotropy index
    eta_ = material_properties["eta"].template get<double>();
    // Liquefaction control b
    b_ = material_properties["b"].template get<double>();
    // Coefficient for granular kinetic fluctuation induced by volumetric
    // deformation
    mv_ = material_properties["mv"].template get<double>();
    // Coefficient for granular kinetic fluctuation induced by shear deformation
    ms_ = material_properties["ms"].template get<double>();
    // Scale coefficient for volumetric dissipation
    lambdav_ = material_properties["lambdav"].template get<double>();
    // Slope of critical state line parameter
    alpha_ = material_properties["alpha"].template get<double>();
    // Maximum value parameter for residual strain
    h_ = material_properties["h"].template get<double>();

    // Properties
    properties_ = material_properties;

  } catch (std::exception& except) {
    console_->error("Material parameter not set: {}\n", except.what());
  }
}

//! Initialise state variables
template <unsigned Tdim>
mpm::dense_map mpm::Thermodynamics<Tdim>::initialise_state_variables() {
  mpm::dense_map state_vars = {// Elastic potential density
                               {"omega_e", 0.},
                               // Liquefaction state variable
                               {"sc", 0.},
                               // Liquefaction parameter
                               {"k", 0.},
                               // Parameter r
                               {"r", 0.5},
                               // Parameter ts
                               {"ts", 0.35},
                               // Elastic volumetric strain
                               {"epse_v", 1.E-15},
                               // Elastic strain second invariant
                               {"epse2", 1.E-15},
                               // Elastic strain third invariant
                               {"epse3", 1.E-15},
                               // Equivalent plastic shear strain
                               {"epss", 1.E-15},
                               // Plastic strain components
                               {"plastic_strain0", 0.},
                               {"plastic_strain1", 0.},
                               {"plastic_strain2", 0.},
                               {"plastic_strain3", 0.},
                               {"plastic_strain4", 0.},
                               {"plastic_strain5", 0.}};

  return state_vars;
}

//! Compute strain invariants
template <unsigned Tdim>
bool mpm::Thermodynamics<Tdim>::compute_strain_invariants(
    const Vector6d& strain, mpm::dense_map* state_vars) {

  // Compute volumetric strain (first invariant)
  double epse_v = strain(0) + strain(1) + strain(2);

  if (abs(epse_v) > 1.E-15)
    (*state_vars)["epse_v"] = epse_v;
  else
    (*state_vars)["epse_v"] = 1.E-15;

  // Compute deviatoric strain (note in our code we use Engineering Strain)
  Vector6d deviatoric_strain = Vector6d::Zero();
  deviatoric_strain(0) = strain(0) - epse_v / 3.;
  deviatoric_strain(1) = strain(1) - epse_v / 3.;
  deviatoric_strain(2) = strain(2) - epse_v / 3.;
  deviatoric_strain(3) = strain(3) / 2.;
  deviatoric_strain(4) = strain(4) / 2.;
  deviatoric_strain(5) = strain(5) / 2.;

  // Compute second invariant (check convention 6 or 3)
  double epse2 =
      sqrt((pow((deviatoric_strain(0) - deviatoric_strain(1)), 2.) +
            pow((deviatoric_strain(1) - deviatoric_strain(2)), 2.) +
            pow((deviatoric_strain(2) - deviatoric_strain(0)), 2.)) /
               6.0 +
           pow(deviatoric_strain(3), 2) + pow(deviatoric_strain(4), 2) +
           pow(deviatoric_strain(5), 2));

  if (abs(epse2) > 1.E-15)
    (*state_vars)["epse2"] = epse2;
  else
    (*state_vars)["epse2"] = 1.E-15;

  // Compute third invariant
  double epse3 =
      (deviatoric_strain(0) * deviatoric_strain(1) * deviatoric_strain(2)) -
      (deviatoric_strain(2) * pow(deviatoric_strain(3), 2)) -
      (deviatoric_strain(0) * pow(deviatoric_strain(4), 2)) -
      (deviatoric_strain(1) * pow(deviatoric_strain(5), 2)) +
      ((2 * deviatoric_strain(3) * deviatoric_strain(4) *
        deviatoric_strain(5)));

  if (abs(epse3) > 1.E-15)
    (*state_vars)["epse3"] = epse3;
  else
    (*state_vars)["epse3"] = 1.E-15;

  return true;
}

//! Compute state parameters
template <unsigned Tdim>
bool mpm::Thermodynamics<Tdim>::compute_state_variables(
    const Vector6d& stress, const Vector6d& strain,
    mpm::dense_map* state_vars) {

  // Get state variables
  const double epse_v = (*state_vars)["epse_v"];
  const double epse2 = (*state_vars)["epse2"];
  const double epse3 = (*state_vars)["epse3"];

  // Compute elastic potential density
  double omega_e;

  if (epse_v <= 1.E-15 || epse2 <= 1.E-15)
    omega_e = elastic_bulk_modulus_ * exp(B1_ * density_) *
              pow((epse_v + cohesion_), m_) *
              (pow(epse_v, 2.) / (m_ + 2.) +
               xi_modulus_ratio_ * pow(epse2, 2.) + eta_);
  else
    omega_e =
        elastic_bulk_modulus_ * exp(B1_ * density_) *
        pow((epse_v + cohesion_), m_) *
        (pow(epse_v, 2.) / (m_ + 2.) + xi_modulus_ratio_ * pow(epse2, 2.) +
         eta_ * pow(epse3, 5.) / epse_v / pow(epse2, 2.));

  // Update state variables
  (*state_vars)["omega_e"] = omega_e;

  // Compute k (need to double check)
  double k = 6. / pow(omega_e, 1. / 3.);

  // Update state variables
  (*state_vars)["k"] = k;

  // Get state variables
  const double epss = (*state_vars)["epss"];

  // Compute sc
  // double sc =
  //     epss / (epss + b_) * (1 - tanh(k * pow(omega_e, 1. / 3.) - 3.)) / 2.;
  double sc = 0.0;
  
  // Update state variables
  (*state_vars)["sc"] = sc;

  return true;
}

//! Compute plastic tensor
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::Thermodynamics<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& strain,
    const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars) {

  // Compute state and strain invariants
  this->compute_strain_invariants(strain, state_vars);
  this->compute_state_variables(stress, strain, state_vars);

  // Get state variables
  const double epse_v = (*state_vars)["epse_v"];
  const double epse2 = (*state_vars)["epse2"];
  const double epse3 = (*state_vars)["epse3"];
  const double epss = (*state_vars)["epss"];
  const double sc = (*state_vars)["sc"];
  double zeta = 2. / 3.;
  // Need to check zeta definition

  // Compute deviatoric strain (note in our code we use Engineering Strain)
  Vector6d deviatoric_strain = Vector6d::Zero();
  deviatoric_strain(0) = strain(0) - epse_v / 3.;
  deviatoric_strain(1) = strain(1) - epse_v / 3.;
  deviatoric_strain(2) = strain(2) - epse_v / 3.;
  deviatoric_strain(3) = strain(3) / 2.;
  deviatoric_strain(4) = strain(4) / 2.;
  deviatoric_strain(5) = strain(5) / 2.;

  // Compute derivative in terms of p (p0 in matlab)
  double p0;
  double B = elastic_bulk_modulus_ * exp(B1_ * density_);
  if (epse2 <= 1.E-15)
    p0 = B * pow((epse_v + cohesion_), m_) * epse_v +
         B * m_ * pow((epse_v + cohesion_), (m_ - 1.)) * xi_modulus_ratio_ *
             pow(epse2, 2.);
  else
    p0 = B * pow((epse_v + cohesion_), m_) * epse_v +
         B * m_ * pow((epse_v + cohesion_), (m_ - 1.)) * xi_modulus_ratio_ *
             pow(epse2, 2.) +
         B * zeta * (m_ - 1) * pow((epse_v + cohesion_), (m_ - 2.)) *
             pow(epse3, 5.) / pow(epse2, 2.);

  // Compute cij
  Vector6d cij = Vector6d::Zero();
  cij(0) = pow(deviatoric_strain(0), 2.) + pow(deviatoric_strain(3), 2.) +
           pow(deviatoric_strain(5), 2.) - pow(epse2, 2. / 3.);
  cij(1) = pow(deviatoric_strain(3), 2.) + pow(deviatoric_strain(1), 2.) +
           pow(deviatoric_strain(4), 2.) - pow(epse2, 2. / 3.);
  cij(2) = pow(deviatoric_strain(5), 2.) + pow(deviatoric_strain(4), 2.) +
           pow(deviatoric_strain(2), 2.) - pow(epse2, 2. / 3.);
  cij(3) = deviatoric_strain(0) * deviatoric_strain(3) +
           deviatoric_strain(3) * deviatoric_strain(1) +
           deviatoric_strain(5) * deviatoric_strain(4);
  cij(4) = deviatoric_strain(3) * deviatoric_strain(5) +
           deviatoric_strain(1) * deviatoric_strain(4) +
           deviatoric_strain(4) * deviatoric_strain(2);
  cij(5) = deviatoric_strain(5) * deviatoric_strain(0) +
           deviatoric_strain(4) * deviatoric_strain(3) +
           deviatoric_strain(2) * deviatoric_strain(5);

  // Compute stress
  Vector6d updated_stress = Vector6d::Zero();
  Vector6d one_volumetric;
  one_volumetric << 1., 1., 1., 0., 0., 0.;
  Vector6d Eij = Vector6d::Zero();

  if (epse2 <= 1.E-15) {

    Eij = 2 * xi_modulus_ratio_ * deviatoric_strain;
    updated_stress = (1. - sc) * ((p0 * one_volumetric) +
                                  B * pow((epse_v + cohesion_), m_) * Eij);
  } else {
    Eij = (2 * xi_modulus_ratio_ -
           2 * zeta * pow(epse3, 5.) / pow(epse2, 4.) / (epse_v + cohesion_)) *
              deviatoric_strain +
          (5 * zeta * pow(epse3, 2.) / pow(epse2, 2.) / (epse_v + cohesion_)) *
              cij;
    updated_stress = (1. - sc) * ((p0 * one_volumetric) +
                                  B * pow((epse_v + cohesion_), m_) * Eij);
  }

  return updated_stress;
}