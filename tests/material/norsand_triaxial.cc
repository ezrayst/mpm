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

  // Castlegate
  jmaterial["density"] = 2000.;
  jmaterial["poisson_ratio"] = 0.17;
  jmaterial["reference_pressure"] = 1.0E+5;
  jmaterial["friction_cs"] = 27.;
  jmaterial["N"] = 0.3;
  jmaterial["lambda"] = 0.11;
  jmaterial["kappa"] = 0.008;
  jmaterial["gamma"] = 1.00;
  jmaterial["chi"] = 3.5;
  jmaterial["hardening_modulus"] = 600.0;
  jmaterial["void_ratio_initial"] = 0.38;
  // jmaterial["p_image_initial"] = 1.5E+6;
  jmaterial["p_image_initial"] = 3.0E+6;
  // jmaterial["p_image_initial"] = 7.5E+6;
  // jmaterial["p_image_initial"] = 1.5E+7;
  // jmaterial["p_image_initial"] = 2.4E+7;
  jmaterial["bond_model"] = true;
  jmaterial["p_cohesion_initial"] = 1.2E+7;
  jmaterial["p_dilation_initial"] = 1.2E+7;
  jmaterial["m_cohesion"] = 10;
  jmaterial["m_dilation"] = 1;
  jmaterial["m_modulus"] = 100;

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

    // 1000 psi (pore pressure 500 psi)
    // stress(0) = -3447000;
    // stress(1) = -3447000;
    // stress(2) = -3447000 - 1;
    // stress(3) = 0;
    // stress(4) = 0;
    // stress(5) = 0;

    // // 1500 psi (pore pressure 500 psi)
    stress(0) = -6895000;
    stress(1) = -6895000;
    stress(2) = -6895000 - 1;
    stress(3) = 0;
    stress(4) = 0;
    stress(5) = 0;

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

    // Initialise dstrain
    unsigned multiplier = 10;
    mpm::Material<Dim>::Vector6d dstrain;
    dstrain.setZero();
    dstrain(2) = -0.0001 / multiplier;
    dstrain(0) = -1 * dstrain(2) * a2 / (a2 + a1);
    dstrain(1) = -1 * dstrain(2) * a2 / (a2 + a1);
    dstrain(3) = 0.0000000 / multiplier;
    dstrain(4) = 0.0000000 / multiplier;
    dstrain(5) = 0.0000000 / multiplier;

    // Compute updated stress
    mpm::dense_map state_vars = material->initialise_state_variables();
    std::ofstream myfile;
    std::ofstream myfile2;
    myfile.open("norsand_stress.txt");
    myfile2.open("norsand_state.txt");

    // // Compute stress invariants
    // stress_neg = -stress;

    // // Compute mean stress_neg p
    // double mean_p = (stress_neg(0) + stress_neg(1) + stress_neg(2)) / 3.;

    // // Compute J2
    // double j2 = (std::pow((stress_neg(0) - stress_neg(1)), 2) +
    //              std::pow((stress_neg(1) - stress_neg(2)), 2) +
    //              std::pow((stress_neg(0) - stress_neg(2)), 2)) /
    //                 6.0 +
    //             std::pow(stress_neg(3), 2) + std::pow(stress_neg(4), 2) +
    //             std::pow(stress_neg(5), 2);

    // // Compute q
    // double deviatoric_q = std::sqrt(3 * j2);

    // // Compute the deviatoric stress_neg
    // Vector6d dev_stress_neg = stress_neg;
    // for (unsigned i = 0; i < 3; ++i) dev_stress_neg(i) -= mean_p;

    // // Compute J3
    // double j3 = (dev_stress_neg(0) * dev_stress_neg(1) * dev_stress_neg(2)) -
    //             (dev_stress_neg(2) * std::pow(dev_stress_neg(3), 2)) +
    //             ((2 * dev_stress_neg(3) * dev_stress_neg(4) *
    //             dev_stress_neg(5)) -
    //              (dev_stress_neg(0) * std::pow(dev_stress_neg(4), 2)) -
    //              (dev_stress_neg(1) * std::pow(dev_stress_neg(5), 2)));

    // // Compute Lode angle value
    // double lode_angle_val = (3. * std::sqrt(3.) / 2.) * (j3 /
    // std::pow(j2, 1.5)); if (lode_angle_val > 1.0) lode_angle_val = 1.0; if
    // (lode_angle_val < -1.0) lode_angle_val = -1.0;

    // // Compute Lode angle (sin convention)
    // double lode_angle = (1. / 3.) * asin(lode_angle_val);
    // if (lode_angle > M_PI / 6.) lode_angle = M_PI / 6.;
    // if (lode_angle < -M_PI / 6.) lode_angle = -M_PI / 6.;

    // // Compute M_theta (Jefferies and Shuttle, 2011)
    // double cos_lode_angle = cos(3. / 2. * lode_angle + M_PI / 4.);
    // double M_theta = Mtc_ - std::pow(Mtc_, 2) / (3. + Mtc_) * cos_lode_angle;

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

      // // Compute stress invariants
      // stress_neg = -stress;

      // // Compute mean stress_neg p
      // mean_p = (stress_neg(0) + stress_neg(1) + stress_neg(2)) / 3.;

      // // Compute J2
      // j2 = (std::pow((stress_neg(0) - stress_neg(1)), 2) +
      //              std::pow((stress_neg(1) - stress_neg(2)), 2) +
      //              std::pow((stress_neg(0) - stress_neg(2)), 2)) /
      //                 6.0 +
      //             std::pow(stress_neg(3), 2) + std::pow(stress_neg(4), 2) +
      //             std::pow(stress_neg(5), 2);

      // // Compute q
      // deviatoric_q = std::sqrt(3 * j2);

      // // Compute the deviatoric stress_neg
      // Vector6d dev_stress_neg = stress_neg;
      // for (unsigned i = 0; i < 3; ++i) dev_stress_neg(i) -= mean_p;

      // // Compute J3
      // j3 = (dev_stress_neg(0) * dev_stress_neg(1) * dev_stress_neg(2)) -
      //             (dev_stress_neg(2) * std::pow(dev_stress_neg(3), 2)) +
      //             ((2 * dev_stress_neg(3) * dev_stress_neg(4) *
      //             dev_stress_neg(5)) -
      //              (dev_stress_neg(0) * std::pow(dev_stress_neg(4), 2)) -
      //              (dev_stress_neg(1) * std::pow(dev_stress_neg(5), 2)));

      // // Compute Lode angle value
      // lode_angle_val = (3. * std::sqrt(3.) / 2.) * (j3 / std::pow(j2, 1.5));
      // if (lode_angle_val > 1.0) lode_angle_val = 1.0;
      // if (lode_angle_val < -1.0) lode_angle_val = -1.0;

      // // Compute Lode angle (sin convention)
      // lode_angle = (1. / 3.) * asin(lode_angle_val);
      // if (lode_angle > M_PI / 6.) lode_angle = M_PI / 6.;
      // if (lode_angle < -M_PI / 6.) lode_angle = -M_PI / 6.;

      // // Compute M_theta (Jefferies and Shuttle, 2011)
      // cos_lode_angle = cos(3. / 2. * lode_angle + M_PI / 4.);
      // M_theta = Mtc_ - std::pow(Mtc_, 2) / (3. + Mtc_) * cos_lode_angle;

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