#include <fstream>
#include <iostream>
#include <limits>
#include <string>

#include <cmath>

#include "Eigen/Dense"
#include "catch.hpp"
#include "json.hpp"

#include "cell.h"
#include "materials/material.h"
#include "node.h"
#include "particle.h"

//! Check NorSand class in 3D
TEST_CASE("NorSand is checked in 3D", "[material][NorSand][3D]") {
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

  // Testing
  jmaterial["density"] = 2000.;
  jmaterial["poisson_ratio"] = 0.2;
  jmaterial["reference_pressure"] = 1.0E+5;
  jmaterial["friction_cs"] = 30.;
  jmaterial["N"] = 0.3;
  jmaterial["lambda"] = 0.09;
  jmaterial["kappa"] = 0.03;
  jmaterial["gamma"] = 0.8;
  jmaterial["chi"] = 3.5;
  jmaterial["hardening_modulus"] = 200.0;
  jmaterial["void_ratio_initial"] = 0.75;
  // jmaterial["void_ratio_initial"] = 1.20;
  jmaterial["p_image_initial"] = 1.0E+4;
  jmaterial["bond_model"] = false;
  jmaterial["p_cohesion_initial"] = 2.0E+3;
  jmaterial["p_dilation_initial"] = 5.0E+3;
  jmaterial["m_cohesion"] = 2.0;
  jmaterial["m_dilation"] = 2.0;
  jmaterial["m_modulus"] = 1000;

  // Castlegate
  // jmaterial["density"] = 2000.;
  // jmaterial["poisson_ratio"] = 0.17;
  // jmaterial["reference_pressure"] = 1.0E+5;
  // jmaterial["friction_cs"] = 27.;
  // jmaterial["N"] = 0.3;
  // jmaterial["lambda"] = 0.11;
  // jmaterial["kappa"] = 0.008;
  // jmaterial["gamma"] = 1.00;
  // jmaterial["chi"] = 3.5;
  // jmaterial["hardening_modulus"] = 600.0;
  // jmaterial["void_ratio_initial"] = 0.38;
  // // jmaterial["p_image_initial"] = 1.5E+6;
  // jmaterial["p_image_initial"] = 3.0E+6;
  // // jmaterial["p_image_initial"] = 7.5E+6;
  // // jmaterial["p_image_initial"] = 1.5E+7;
  // // jmaterial["p_image_initial"] = 2.4E+7;
  // jmaterial["bond_model"] = true;
  // jmaterial["p_cohesion_initial"] = 1.2E+7;
  // jmaterial["p_dilation_initial"] = 1.2E+7;
  // jmaterial["m_cohesion"] = 10;
  // jmaterial["m_dilation"] = 1;
  // jmaterial["m_modulus"] = 100;

  double pi_constant = M_PI;

  // Check triaxial drain test
  SECTION("NorSand check drained stresses") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "NorSand3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();

    // // Testing with 20 kPa
    stress(0) = -20000;
    stress(1) = -20000;
    stress(2) = -20000 - 1;
    stress(3) = 0;
    stress(4) = 0;
    stress(5) = 0;

    // 1000 psi (pore pressure 500 psi)
    // stress(0) = -3447000;
    // stress(1) = -3447000;
    // stress(2) = -3447000 - 1;
    // stress(3) = 0;
    // stress(4) = 0;
    // stress(5) = 0;

    // // 1500 psi (pore pressure 500 psi)
    // stress(0) = -6895000;
    // stress(1) = -6895000;
    // stress(2) = -6895000 - 1;
    // stress(3) = 0;
    // stress(4) = 0;
    // stress(5) = 0;

    // // 3000 psi (pore pressure 500 psi)
    // stress(0) = -17237000;
    // stress(1) = -17237000;
    // stress(2) = -17237000 - 1;
    // stress(3) = 0;
    // stress(4) = 0;
    // stress(5) = 0;

    // // 5500 psi (pore pressure 500 psi)
    // stress(0) = -34474000;
    // stress(1) = -34474000;
    // stress(2) = -34474000 - 1;
    // stress(3) = 0;
    // stress(4) = 0;
    // stress(5) = 0;

    // // 8500 psi (pore pressure 500 psi)
    // stress(0) = -55158000;
    // stress(1) = -55158000;
    // stress(2) = -55158000 - 1;
    // stress(3) = 0;
    // stress(4) = 0;
    // stress(5) = 0;

    // Check dmatrix
    const double poisson_ratio_ = jmaterial["poisson_ratio"];
    const double kappa = jmaterial["kappa"];
    const double e_init = jmaterial["void_ratio_initial"];
    const double p_dilation = jmaterial["m_dilation"];
    const double p_cohesion = jmaterial["m_cohesion"];
    const double m_modulus = jmaterial["m_modulus"];
    const double p = (stress(0) + stress(1) + stress(2)) / 3.0;
    const double bulk_modulus_ =
        (1. + e_init) / kappa * p + m_modulus * (p_cohesion + p_dilation);
    const double G = 3. * bulk_modulus_ * (1. - 2. * poisson_ratio_) /
                     (2.0 * (1. + poisson_ratio_));
    const double a1 = bulk_modulus_ + (4.0 / 3.0) * G;
    const double a2 = bulk_modulus_ - (2.0 / 3.0) * G;

    // Initialise dstrain (DRAINED)
    // unsigned multiplier = 50;
    // mpm::Material<Dim>::Vector6d dstrain;
    // dstrain.setZero();
    // dstrain(2) = -0.0005 / multiplier;
    // dstrain(0) = -1 * dstrain(2) * a2 / (a2 + a1);
    // dstrain(1) = -1 * dstrain(2) * a2 / (a2 + a1);
    // dstrain(3) = 0.0000000 / multiplier;
    // dstrain(4) = 0.0000000 / multiplier;
    // dstrain(5) = 0.0000000 / multiplier;

    // Initialise dstrain (UNDRAINED)
    unsigned multiplier = 50;
    mpm::Material<Dim>::Vector6d dstrain;
    dstrain.setZero();
    dstrain(2) = -0.0005 / multiplier;
    dstrain(0) = -0.5 * dstrain(2);
    dstrain(1) = -0.5 * dstrain(2);
    dstrain(3) = 0.0000000 / multiplier;
    dstrain(4) = 0.0000000 / multiplier;
    dstrain(5) = 0.0000000 / multiplier;

    // // Castlegate
    // unsigned multiplier = 10;
    // mpm::Material<Dim>::Vector6d dstrain;
    // dstrain.setZero();
    // dstrain(2) = -0.0001 / multiplier;
    // dstrain(0) = -1 * dstrain(2) * a2 / (a2 + a1);
    // dstrain(1) = -1 * dstrain(2) * a2 / (a2 + a1);
    // dstrain(3) = 0.0000000 / multiplier;
    // dstrain(4) = 0.0000000 / multiplier;
    // dstrain(5) = 0.0000000 / multiplier;

    // Compute updated stress
    mpm::dense_map state_vars = material->initialise_state_variables();
    std::ofstream myfile;
    std::ofstream myfile2;
    myfile.open("norsand_stress.txt");
    myfile2.open("norsand_state.txt");

    // Write initial states
    myfile << stress(0) << '\t' << stress(1) << '\t' << stress(2) << '\t'
           << stress(3) << '\t' << stress(4) << '\t' << stress(5) << '\n';
    myfile2 << (state_vars).at("p_image") << '\t' << (state_vars).at("e_image")
            << '\t' << (state_vars).at("void_ratio") << '\t'
            << (state_vars).at("p_cohesion") << '\t'
            << (state_vars).at("p_dilation") << '\n';
    // Loop
    for (unsigned i = 0; i < multiplier * 1000 - 1; ++i) {
      stress = material->compute_stress(stress, dstrain, particle.get(),
                                        &state_vars);

      if (i % multiplier == 0) {
        myfile << stress(0) << '\t' << stress(1) << '\t' << stress(2) << '\t'
               << stress(3) << '\t' << stress(4) << '\t' << stress(5) << '\n';
        myfile2 << (state_vars).at("p_image") << '\t'
                << (state_vars).at("e_image") << '\t'
                << (state_vars).at("void_ratio") << '\t'
                << (state_vars).at("p_cohesion") << '\t'
                << (state_vars).at("p_dilation") << '\n';
      }
    }
    myfile.close();
  }
}