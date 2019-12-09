#include <fstream>
#include <iostream>
#include <limits>

#include <cmath>

#include "Eigen/Dense"
#include "catch.hpp"
#include "json.hpp"

#include "cell.h"
#include "material/material.h"
#include "node.h"
#include "particle.h"

//! Check Thermodynamics class in 3D
TEST_CASE("Thermodynamics is checked in 3D", "[material][Thermodynamics][3D]") {
  // Tolerance
  const double Tolerance = 1.E-7;

  const unsigned Dim = 3;

  // Add particle
  mpm::Index pid = 0;
  Eigen::Matrix<double, Dim, 1> coords;
  coords.setZero();
  auto particle = std::make_shared<mpm::Particle<Dim>>(pid, coords);

  // Initialise material
  Json jmaterial;
  jmaterial["density"] = 1800.;
  jmaterial["B1"] = 0.0098;
  jmaterial["elastic_bulk_modulus"] = 5500.;
  jmaterial["xi_modulus_ratio"] = 2. / 3.;
  jmaterial["cohesion"] = 0.0;
  jmaterial["m"] = 0.5;
  jmaterial["eta"] = 1000.;
  jmaterial["b"] = 0.01;
  jmaterial["mv"] = 8.5E+3;
  jmaterial["ms"] = 1.7E+4;
  jmaterial["lambdav"] = 1.000;
  jmaterial["alpha"] = 7.0;
  jmaterial["h"] = 3.5;

  //! Check for id = 0
  SECTION("Thermodynamics id is zero") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "Thermodynamics3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);
  }

  SECTION("Thermodynamics id is positive") {
    //! Check for id is a positive value
    unsigned id = std::numeric_limits<unsigned>::max();
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "Thermodynamics3D", std::move(id), jmaterial);
    REQUIRE(material->id() == std::numeric_limits<unsigned>::max());
  }

  // ! Check failed initialisation
  SECTION("Thermodynamics failed initialisation") {
    unsigned id = 0;
    // Initialise material
    Json jmaterial;
    jmaterial["density"] = 1000.;
    jmaterial["poisson_ratio"] = 0.3;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "Thermodynamics3D", std::move(id), jmaterial);
  }

  //! Check material properties
  SECTION("Thermodynamics check material properties") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "Thermodynamics3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Get material properties
    REQUIRE(material->template property<double>("density") ==
            Approx(jmaterial["density"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("B1") ==
            Approx(jmaterial["B1"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("elastic_bulk_modulus") ==
            Approx(jmaterial["elastic_bulk_modulus"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("xi_modulus_ratio") ==
            Approx(jmaterial["xi_modulus_ratio"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("cohesion") ==
            Approx(jmaterial["cohesion"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("m") ==
            Approx(jmaterial["m"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("eta") ==
            Approx(jmaterial["eta"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("b") ==
            Approx(jmaterial["b"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("mv") ==
            Approx(jmaterial["mv"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("ms") ==
            Approx(jmaterial["ms"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("lambdav") ==
            Approx(jmaterial["lambdav"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("alpha") ==
            Approx(jmaterial["alpha"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("h") ==
            Approx(jmaterial["h"]).epsilon(Tolerance));

    // Check if state variable is initialised
    SECTION("State variable is initialised") {
      mpm::dense_map state_variables = material->initialise_state_variables();
      REQUIRE(state_variables.empty() == false);
      REQUIRE(state_variables.at("omega_e") == Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("sc") == Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("k") == Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("r") == Approx(0.5).epsilon(Tolerance));
      REQUIRE(state_variables.at("ts") == Approx(0.35).epsilon(Tolerance));
      REQUIRE(state_variables.at("epse_v") ==
              Approx(1.E-15).epsilon(Tolerance));
      REQUIRE(state_variables.at("epse2") == Approx(1.E-15).epsilon(Tolerance));
      REQUIRE(state_variables.at("epse3") == Approx(1.E-15).epsilon(Tolerance));
      REQUIRE(state_variables.at("epss") == Approx(1.E-15).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain0") ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain1") ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain2") ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain3") ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain4") ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain5") ==
              Approx(0.).epsilon(Tolerance));
    }
  }

  SECTION("Thermodynamics check stresses") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "Thermodynamics3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    REQUIRE(stress(0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.).epsilon(Tolerance));

    stress(0) = 0.0;
    stress(1) = 0.0;
    stress(2) = 0.0;
    stress(3) = 0.0;
    stress(4) = 0.0;
    stress(5) = 0.0;

    // Initialise strain
    mpm::Material<Dim>::Vector6d strain;
    strain.setZero();
    strain(0) = 0.0010000;
    strain(1) = 0.0000000;
    strain(2) = 0.0000000;
    strain(3) = 0.0000000;
    strain(4) = 0.0000000;
    strain(5) = 0.0000000;

    // Compute updated stress
    mpm::dense_map state_vars = material->initialise_state_variables();

    // stress =
    //     material->compute_stress(stress, strain, particle.get(), &state_vars);
    // // Check stressees
    // REQUIRE(stress(0) == Approx(1.63461538461538e+04).epsilon(Tolerance));
    // REQUIRE(stress(1) == Approx(1.25000000000000e+04).epsilon(Tolerance));
    // REQUIRE(stress(2) == Approx(0.86538461538462e+04).epsilon(Tolerance));
    // REQUIRE(stress(3) == Approx(0.000000e+00).epsilon(Tolerance));
    // REQUIRE(stress(4) == Approx(0.000000e+00).epsilon(Tolerance));
    // REQUIRE(stress(5) == Approx(0.000000e+00).epsilon(Tolerance));

    // Initialise strain
    // strain(0) = 0.0020000;
    // strain(1) = 0.0010000;
    // strain(2) = 0.0000000;
    // strain(3) = 0.0000100;
    // strain(4) = 0.0000000;
    // strain(5) = 0.0000000;

    // Reset stress
    stress.setZero();

    // Compute updated stress
    stress =
        material->compute_stress(stress, strain, particle.get(), &state_vars);

    // Check stressees
    REQUIRE(stress(0) == Approx(1.63461538461538e+04).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(1.25000000000000e+04).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(0.86538461538462e+04).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(3.84615384615385e+01).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.00000000000000e+00).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.00000000000000e+00).epsilon(Tolerance));
  }
}
