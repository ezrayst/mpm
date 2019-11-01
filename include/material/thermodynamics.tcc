//! Constructor with id and material properties
template <unsigned Tdim>
mpm::Thermodynamics<Tdim>::Thermodynamics(unsigned id,
                                          const Json& material_properties)
    : Material<Tdim>(id, material_properties) {
  try {
    // Density B1
    density_ = material_properties["density"].template get<double>();
    // Elastic bulk modulus B0
    elastic_bulk_modulus_ =
        material_properties["elastic_bulk_modulus"].template get<double>();
    // Xi, twice the ratio of shear to bulk modulus
    xi_modulus_ratio_ =
        material_properties["xi_modulus_ratio"].template get<double>();
    // Cohesion
    cohesion_ = material_properties["cohesion"].template get<double>();
    // Non-linear index m
    m_ = material_properties["m"].template get<double>();
    // Stress induced anisotropy index
    eta_ = material_properties["eta"].template get<double>();
    // Liquefaction control b
    b_ = material_properties["b"].template get<double>() * M_PI / 180.;
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
                               // Mean stress
                               {"p", 0.},
                               // Deviatoric stress
                               {"q", 0.},
                               // Lode angle
                               {"lode_angle", 0.},
                               // J2
                               {"j2", 0.},
                               // J3
                               {"j3", 0.},
                               // Elastic volumetric strain
                               {"epse_v", 0.},
                               // Elastic strain second invariant
                               {"epse2", 0.},
                               // Elastic strain third invariant
                               {"epse3", 0.},
                               // Equivalent plastic shear strain
                               {"epss", 0.},
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
    const Vector6d& stress, mpm::dense_map* state_vars) {

  return true;
}

//! Compute state parameters
template <unsigned Tdim>
bool mpm::Thermodynamics<Tdim>::compute_state_variables(
    const Vector6d& stress, const Vector6d& dstrain,
    mpm::dense_map* state_vars) {

  return true;
}

//! Compute plastic tensor
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::Thermodynamics<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars) {

  Vector6d updated_stress;

  return updated_stress;
}