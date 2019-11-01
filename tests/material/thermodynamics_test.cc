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
  auto particle = std::make_shared<mpm::Particle<Dim, 1>>(pid, coords);

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

  //! Check failed initialisation
  // SECTION("Thermodynamics failed initialisation") {
  //   unsigned id = 0;
  //   // Initialise material
  //   Json jmaterial;
  //   jmaterial["density"] = 1000.;
  //   jmaterial["poisson_ratio"] = 0.3;
  //   auto material =
  //       Factory<mpm::Material<Dim>, unsigned, const
  //       Json&>::instance()->create(
  //           "Thermodynamics3D", std::move(id), jmaterial);
  // }

  //! Check material properties
  // SECTION("Thermodynamics check material properties") {
  //   unsigned id = 0;
  //   auto material =
  //       Factory<mpm::Material<Dim>, unsigned, const
  //       Json&>::instance()->create(
  //           "Thermodynamics3D", std::move(id), jmaterial);
  //   REQUIRE(material->id() == 0);

  //   // Get material properties
  //   REQUIRE(material->template property<double>("density") ==
  //           Approx(jmaterial["density"]).epsilon(Tolerance));
  //   REQUIRE(material->template property<double>("youngs_modulus") ==
  //           Approx(jmaterial["youngs_modulus"]).epsilon(Tolerance));
  //   REQUIRE(material->template property<double>("poisson_ratio") ==
  //           Approx(jmaterial["poisson_ratio"]).epsilon(Tolerance));
  //   REQUIRE(material->template property<double>("shear_modulus_constant") ==
  //           Approx(jmaterial["shear_modulus_constant"]).epsilon(Tolerance));
  //   REQUIRE(material->template property<double>("shear_modulus_exponent") ==
  //           Approx(jmaterial["shear_modulus_exponent"]).epsilon(Tolerance));
  //   REQUIRE(material->template property<double>("reference_pressure") ==
  //           Approx(jmaterial["reference_pressure"]).epsilon(Tolerance));
  //   // REQUIRE(material->template property<double>("friction_cs") ==
  //   //         Approx(jmaterial["friction_cs"] * pi_constant /
  //   //         180).epsilon(Tolerance));
  //   REQUIRE(material->template property<double>("N") ==
  //           Approx(jmaterial["N"]).epsilon(Tolerance));
  //   REQUIRE(material->template property<double>("e_min") ==
  //           Approx(jmaterial["e_min"]).epsilon(Tolerance));
  //   REQUIRE(material->template property<double>("e_max") ==
  //           Approx(jmaterial["e_max"]).epsilon(Tolerance));
  //   REQUIRE(material->template property<double>("crushing_pressure") ==
  //           Approx(jmaterial["crushing_pressure"]).epsilon(Tolerance));
  //   REQUIRE(material->template property<double>("chi") ==
  //           Approx(jmaterial["chi"]).epsilon(Tolerance));
  //   REQUIRE(material->template property<double>("hardening_modulus") ==
  //           Approx(jmaterial["hardening_modulus"]).epsilon(Tolerance));
  //   REQUIRE(material->template property<double>("void_ratio_initial") ==
  //           Approx(jmaterial["void_ratio_initial"]).epsilon(Tolerance));
  //   REQUIRE(material->template property<double>("p_image_initial") ==
  //           Approx(jmaterial["p_image_initial"]).epsilon(Tolerance));
  //   REQUIRE(material->template property<double>("p_cohesion_initial") ==
  //           Approx(jmaterial["p_cohesion_initial"]).epsilon(Tolerance));
  //   REQUIRE(material->template property<double>("m") ==
  //           Approx(jmaterial["m"]).epsilon(Tolerance));
  //   REQUIRE(material->template property<bool>("bond_model") ==
  //           jmaterial["bond_model"]);

  //   // // Check if state variable is initialised
  //   // SECTION("State variable is initialised") {
  //   //   mpm::dense_map state_variables =
  //   //   material->initialise_state_variables();
  //   REQUIRE(state_variables.empty()
  //   //   == false); REQUIRE(state_variables.at("p") ==
  //   //   Approx(0.).epsilon(Tolerance)); REQUIRE(state_variables.at("q") ==
  //   //   Approx(0.).epsilon(Tolerance));
  //   //   REQUIRE(state_variables.at("void_ratio") ==
  //   //   Approx(jmaterial["void_ratio_initial"]).epsilon(Tolerance));
  //   //   REQUIRE(state_variables.at("e_image") ==
  //   //   Approx(0.9350064171).epsilon(Tolerance));
  //   //   REQUIRE(state_variables.at("p_image") ==
  //   //   Approx(8701.46).epsilon(Tolerance));
  //   //   REQUIRE(state_variables.at("psi_image") ==
  //   //   Approx(-0.0850064171).epsilon(Tolerance));

  //   // }
  // }

  // SECTION("Thermodynamics check stresses") {
  //   unsigned id = 0;
  //   auto material =
  //       Factory<mpm::Material<Dim>, unsigned, const
  //       Json&>::instance()->create(
  //           "Thermodynamics3D", std::move(id), jmaterial);
  //   REQUIRE(material->id() == 0);

  //   // Initialise stress
  //   mpm::Material<Dim>::Vector6d stress;
  //   stress.setZero();
  //   REQUIRE(stress(0) == Approx(0.).epsilon(Tolerance));
  //   REQUIRE(stress(1) == Approx(0.).epsilon(Tolerance));
  //   REQUIRE(stress(2) == Approx(0.).epsilon(Tolerance));
  //   REQUIRE(stress(3) == Approx(0.).epsilon(Tolerance));
  //   REQUIRE(stress(4) == Approx(0.).epsilon(Tolerance));
  //   REQUIRE(stress(5) == Approx(0.).epsilon(Tolerance));

  //   stress(0) = -50000;
  //   stress(1) = -50000;
  //   stress(2) = -50000;
  //   stress(3) = 0;
  //   stress(4) = 0;
  //   stress(5) = 0;

  //   // Initialise dstrain
  //   unsigned multiplier = 10;
  //   mpm::Material<Dim>::Vector6d dstrain;
  //   dstrain.setZero();
  //   dstrain(0) = 0.00025000 / multiplier;
  //   dstrain(1) = -0.0005000 / multiplier;
  //   dstrain(2) = 0.0002500 / multiplier;
  //   dstrain(3) = 0.0000000 / multiplier;
  //   dstrain(4) = 0.0000000 / multiplier;
  //   dstrain(5) = 0.0000000 / multiplier;

  //   // Compute updated stress
  //   mpm::dense_map state_vars = material->initialise_state_variables();
  //   std::ofstream myfile;
  //   std::ofstream myfile2;
  //   myfile.open("Thermodynamics_stress.txt");
  //   myfile2.open("Thermodynamics_state.txt");

  //   // Write initial states
  //   myfile << stress(0) << '\t' << stress(1) << '\t' << stress(2) << '\t'
  //          << stress(3) << '\t' << stress(4) << '\t' << stress(5) << '\n';
  //   myfile2 << (state_vars).at("p_image") << '\t' <<
  //   (state_vars).at("e_image")
  //           << '\t' << (state_vars).at("void_ratio") << '\t'
  //           << (state_vars).at("lode_angle") << '\t'
  //           << (state_vars).at("M_theta") << '\t'
  //           << (state_vars).at("p_cohesion") << '\t' <<
  //           (state_vars).at("zeta")
  //           << '\n';
  //   // Loop
  //   for (unsigned i = 0; i < multiplier * 1000 - 1; ++i) {
  //     stress = material->compute_stress(stress, dstrain, particle.get(),
  //                                       &state_vars);

  //     if (i % multiplier == 0) {
  //       myfile << stress(0) << '\t' << stress(1) << '\t' << stress(2) << '\t'
  //              << stress(3) << '\t' << stress(4) << '\t' << stress(5) <<
  //              '\n';
  //       myfile2 << (state_vars).at("p_image") << '\t'
  //               << (state_vars).at("e_image") << '\t'
  //               << (state_vars).at("void_ratio") << '\t'
  //               << (state_vars).at("lode_angle") << '\t'
  //               << (state_vars).at("M_theta") << '\t'
  //               << (state_vars).at("p_cohesion") << '\t'
  //               << (state_vars).at("zeta") << '\n';
  //     }
  //   }
  //   myfile.close();

  //   // // Check stressees
  //   // REQUIRE(stress(0) == Approx(1.63461538461538e+04).epsilon(Tolerance));
  //   // REQUIRE(stress(1) == Approx(1.25000000000000e+04).epsilon(Tolerance));
  //   // REQUIRE(stress(2) == Approx(0.86538461538462e+04).epsilon(Tolerance));
  //   // REQUIRE(stress(3) == Approx(0.000000e+00).epsilon(Tolerance));
  //   // REQUIRE(stress(4) == Approx(0.000000e+00).epsilon(Tolerance));
  //   // REQUIRE(stress(5) == Approx(0.000000e+00).epsilon(Tolerance));

  //   // // Initialise strain
  //   // strain(0) = 0.0010000;
  //   // strain(1) = 0.0005000;
  //   // strain(2) = 0.0000000;
  //   // strain(3) = 0.0000100;
  //   // strain(4) = 0.0000000;
  //   // strain(5) = 0.0000000;

  //   // // Reset stress
  //   // stress.setZero();

  //   // // Compute updated stress
  //   // stress =
  //   //     material->compute_stress(stress, strain, particle.get(),
  //   //     &state_vars);

  //   // // Check stressees
  //   // REQUIRE(stress(0) == Approx(1.63461538461538e+04).epsilon(Tolerance));
  //   // REQUIRE(stress(1) == Approx(1.25000000000000e+04).epsilon(Tolerance));
  //   // REQUIRE(stress(2) == Approx(0.86538461538462e+04).epsilon(Tolerance));
  //   // REQUIRE(stress(3) == Approx(3.84615384615385e+01).epsilon(Tolerance));
  //   // REQUIRE(stress(4) == Approx(0.00000000000000e+00).epsilon(Tolerance));
  //   // REQUIRE(stress(5) == Approx(0.00000000000000e+00).epsilon(Tolerance));
  // }
}
