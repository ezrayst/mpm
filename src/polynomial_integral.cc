#include "polynomial.h"

//! 2D Quadrilateral
//! Polynomial order 0 and number of terms 1
//! 4 monomials: 1
//! Integration of monomials over [-1,-1] x [1,1] unit square
template <>
const Eigen::Matrix<double, 1, 1>
    mpm::Polynomial::DefiniteIntegral<2, 1>::Square_Definite_Integrals =
        (Eigen::MatrixXd(1, 1) << 4.0).finished();

//! 2D Quadrilateral
//! Polynomial order 1 and number of terms 4
//! 4 monomials: 1, y, x, xy
//! Integration of monomials over [-1,-1] x [1,1] unit square
template <>
const Eigen::Matrix<double, 4, 1>
    mpm::Polynomial::DefiniteIntegral<2, 4>::Square_Definite_Integrals =
        (Eigen::MatrixXd(4, 1) << 4.0, 0., 0., 0.).finished();

//! 2D Quadrilateral
//! Polynomial order 2 and number of terms 9
//! 9 monomials: 1, y, y^2, x, xy, xy^2, x^2, x^2 y, x^2 y^2
//! Integration of monomials over [-1,-1] x [1,1] unit square
template <>
const Eigen::Matrix<double, 9, 1>
    mpm::Polynomial::DefiniteIntegral<2, 9>::Square_Definite_Integrals =
        (Eigen::MatrixXd(9, 1) << 4., 0., 4. / 3, 0., 0., 0., 4. / 3, 0.,
         4. / 9)
            .finished();

//! 2D Quadrilateral
//! Polynomial order 3 and number of terms 16
//! 16 monomials: 1, y, y^2, y^3, x, xy, xy^2, xy^3, x^2, x^2 y, x^2 y^2, x^2
//! y^3, x^3, x^3 y, x^3 y^2, x^3 y^3
//! Integration of monomials over [-1,-1] x [1,1] unit square
template <>
const Eigen::Matrix<double, 16, 1>
    mpm::Polynomial::DefiniteIntegral<2, 16>::Square_Definite_Integrals =
        (Eigen::MatrixXd(16, 1) << 4.0, 0., 4. / 3, 0., 0., 0., 0., 0., 4. / 3,
         0., 4. / 9, 0., 0., 0., 0., 0.)
            .finished();

//! 2D Quadrilateral
//! Polynomial order 4 and number of terms 25
//! 25 monomials:
//! Integration of monomials over [-1,-1] x [1,1] unit square
template <>
const Eigen::Matrix<double, 25, 1>
    mpm::Polynomial::DefiniteIntegral<2, 25>::Square_Definite_Integrals =
        Eigen::Matrix<double, 25, 1>::Zero();

//! 3D Hexahedron
//! Polynomial order 0 and number of terms 1
//! 4 monomials: 1
//! Integration of monomials over [-1,-1, -1] x [1,1,1] unit cube
template <>
const Eigen::Matrix<double, 1, 1>
    mpm::Polynomial::DefiniteIntegral<3, 1>::Square_Definite_Integrals =
        (Eigen::MatrixXd(1, 1) << 8.0).finished();

//! 3D Hexahedron
//! Polynomial order 1 and number of terms 8
//! 4 monomials: 1, z, y, yz, x, xz, xy, xyz
//! Integration of monomials over [-1,-1, -1] x [1,1,1] unit cube
template <>
const Eigen::Matrix<double, 8, 1>
    mpm::Polynomial::DefiniteIntegral<3, 8>::Square_Definite_Integrals =
        (Eigen::MatrixXd(8, 1) << 8.0, 0., 0., 0., 0., 0., 0., 0.).finished();

//! 3D Hexahedron
//! Polynomial order 2 and number of terms 27
//! 27 monomials: 1, z, z^2, y, yz, yz^2, y^2, y^2 z, y^2 z^2, x, xz, x z^2, x
//! y, xyz, xyz^2, xy^2, xy^2 z, xy^2 z^2, x^2, x^2 z, x^2 z^2, x^2 y, x^2 yz,
//! x^2 yz^2, x^2 y^2, y^2 z, x^2 y^2 z^2
//! Integration of monomials over [-1,-1, -1] x [1,1, 1] unit cube
template <>
const Eigen::Matrix<double, 27, 1>
    mpm::Polynomial::DefiniteIntegral<3, 27>::Square_Definite_Integrals =
        (Eigen::MatrixXd(27, 1) << 8., 0., 8. / 3., 0., 0., 0., 8. / 3., 0.,
         8. / 9., 0., 0., 0., 0., 0., 0., 0., 0., 0., 8. / 3., 0., 8. / 9., 0.,
         0., 0., 8. / 9., 0., 8. / 27.)
            .finished();

//! 3D Hexahedron
//! Polynomial order 3 and number of terms 64
//! 64 monomials: 1, z, z^2, z^3, y, yz, yz^2, yz^3, y^2, y^2 z, y^2 z^2,
//! y^2 z^3, y^3, y^3 z, y^3 z^2, y^3 z^3,
//! x, xz, xz^2, xz^3, xy, xyz, xyz^2, xyz^3, xy^2, xy^2 z, xy^2 z^2,
//! xy^2 z^3, xy^3, xy^3 z, xy^3 z^2, xy^3 z^3,
//! x^2, x^2 z, x^2 z^2, x^2 z^3, x^2 y, x^2 yz, x^2 yz^2, x^2 yz^3, x^2 y^2,
//! x^2 y^2 z, x^2 y^2 z^2, x^2 y^2 z^3, x^2 y^3, x^2 y^3 z, x^2 y^3 z^2,
//! x^2 y^3 z^3,
//! x^3, x^2 z, x^3 z^2, x^3 z^3, x^3 y, x^2 yz, x^3 yz^2, x^3 yz^3, x^3 y^2,
//! x^3 y^2 z, x^3 y^2 z^2, x^3 y^2 z^3, x^3 y^3, x^3 y^3 z, x^3 y^3 z^2,
//! x^3 y^3 z^3,
//! Integration of monomials over [-1,-1, -1] x [1,1, 1] unit cube
template <>
const Eigen::Matrix<double, 64, 1>
    mpm::Polynomial::DefiniteIntegral<3, 64>::Square_Definite_Integrals =
        (Eigen::MatrixXd(64, 1) << 8., 0., 8. / 3., 0., 0., 0., 0., 0., 8. / 3.,
         0., 8. / 9., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
         0., 0., 0., 0., 0., 0., 0., 8. / 3., 0., 8. / 9., 0., 0., 0., 0., 0.,
         8. / 9., 0., 8. / 27., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
         0., 0., 0., 0., 0., 0., 0., 0., 0.)
            .finished();

//! 3D Hexahedron
//! Polynomial order 4 and number of terms 125
//! 125 monomials:
//! Integration of monomials over [-1,-1, -1] x [1,1, 1] unit cube
template <>
const Eigen::Matrix<double, 125, 1>
    mpm::Polynomial::DefiniteIntegral<3, 125>::Square_Definite_Integrals =
        Eigen::Matrix<double, 125, 1>::Zero();

//! 2D Triangle
//! Polynomial order 0 and number of terms 1
//! 4 monomials: 1
//! Integration of monomials over [0,0] x [1,1-x] unit triangle
template <>
const Eigen::Matrix<double, 1, 1>
    mpm::Polynomial::DefiniteIntegral<2, 1>::Tri_Definite_Integrals =
        (Eigen::MatrixXd(1, 1) << 1. / 2.).finished();

//! 2D Triangle
//! Polynomial order 1 and number of terms 3
//! 4 monomials: 1, x, y
//! Integration of monomials over [0,0] x [1,1-x] unit triangle
template <>
const Eigen::Matrix<double, 3, 1>
    mpm::Polynomial::DefiniteIntegral<2, 3>::Tri_Definite_Integrals =
        (Eigen::MatrixXd(3, 1) << 1. / 2., 1. / 6., 1. / 6.).finished();

//! 2D Triangle
//! Polynomial order 1 and number of terms 4
//! 4 monomials: 1, y, x, xy
//! Integration of monomials over [0,0] x [1,1-x] unit triangle
template <>
const Eigen::Matrix<double, 4, 1> mpm::Polynomial::DefiniteIntegral<
    2, 4>::Tri_Definite_Integrals =
    (Eigen::MatrixXd(4, 1) << 1. / 2., 1. / 6., 1. / 6., 1. / 24.).finished();

//! 2D Triangle
//! Polynomial order 2 and number of terms 9
//! 9 monomials: 1, y, y^2, x, xy, xy^2, x^2, x^2 y, x^2 y^2
//! Integration of monomials over [0,0] x [1,1-x] unit triangle
template <>
const Eigen::Matrix<double, 9, 1>
    mpm::Polynomial::DefiniteIntegral<2, 9>::Tri_Definite_Integrals =
        (Eigen::MatrixXd(9, 1) << 1. / 2., 1. / 6., 1. / 12., 1. / 6., 1. / 24.,
         1. / 60., 1. / 12., 1. / 60., 1. / 180.)
            .finished();

//! 2D Triangle
//! Polynomial order 3 and number of terms 16
//! 16 monomials: 1, y, y^2, y^3, x, xy, xy^2, xy^3, x^2, x^2 y, x^2 y^2, x^2
//! y^3, x^3, x^3 y, x^3 y^2, x^3 y^3
//! Integration of monomials over [0,0] x [1,1-x] unit triangle
template <>
const Eigen::Matrix<double, 16, 1>
    mpm::Polynomial::DefiniteIntegral<2, 16>::Tri_Definite_Integrals =
        (Eigen::MatrixXd(16, 1) << 1. / 2., 1. / 6., 1. / 12., 1. / 20.,
         1. / 6., 1. / 24., 1. / 60., 1. / 120., 1. / 12., 1. / 60., 1. / 180.,
         1. / 420., 1. / 20., 1. / 120., 1. / 420., 1. / 1120.)
            .finished();

//! 2D Triangle
//! Polynomial order 4 and number of terms 25
//! 9 monomials:
//! Integration of monomials over [0,0] x [1,1-x] unit triangle
template <>
const Eigen::Matrix<double, 25, 1>
    mpm::Polynomial::DefiniteIntegral<2, 25>::Tri_Definite_Integrals =
        Eigen::Matrix<double, 25, 1>::Zero();

//! 3D Tetrahedron
//! Polynomial order 1 and number of terms 8
//! 4 monomials: 1, z, y, yz, x, xz, xy, xyz
//! Integration of monomials over [0,0,0] x [0, 1-x, 1-x-y] unit tetrahedron
template <>
const Eigen::Matrix<double, 8, 1>
    mpm::Polynomial::DefiniteIntegral<3, 8>::Tri_Definite_Integrals =
        (Eigen::MatrixXd(8, 1) << 1. / 6., 1. / 24., 1. / 24., 1. / 120.,
         1. / 24., 1. / 120., 1. / 120., 1. / 720.)
            .finished();

//! 3D Tetrahedron
//! Polynomial order 2 and number of terms 27
//! 27 monomials: 1, z, z^2, y, yz, yz^2, y^2, y^2 z, y^2 z^2, x, xz, x z^2, x
//! y, xyz, xyz^2, xy^2, xy^2 z, xy^2 z^2, x^2, x^2 z, x^2 z^2, x^2 y, x^2 yz,
//! x^2 yz^2, x^2 y^2, x^2 y^2 z, x^2 y^2 z^2
//! Integration of monomials over [0,0,0] x [0, 1-x, 1-x-y] unit tetrahedron
template <>
const Eigen::Matrix<double, 27, 1>
    mpm::Polynomial::DefiniteIntegral<3, 27>::Tri_Definite_Integrals =
        (Eigen::MatrixXd(27, 1) << 1. / 6., 1. / 24., 1. / 60., 1. / 24.,
         1. / 120., 1. / 360., 1. / 60., 1. / 360., 1. / 1260., 1. / 24.,
         1. / 120., 1. / 360., 1. / 120., 1. / 720., 1. / 2520., 1. / 360.,
         1. / 2520., 1. / 10080., 1. / 60., 1. / 360., 1. / 1260., 1. / 360.,
         1. / 2520., 1. / 10080., 1. / 1260., 1. / 10080., 1. / 45360.)
            .finished();

//! 3D Triangle
//! Polynomial order 3 and number of terms 64
template <>
const Eigen::Matrix<double, 64, 1>
    mpm::Polynomial::DefiniteIntegral<3, 64>::Tri_Definite_Integrals =
        Eigen::Matrix<double, 64, 1>::Zero();

//! 3D Triangle
//! Polynomial order 4 and number of terms 125
template <>
const Eigen::Matrix<double, 125, 1>
    mpm::Polynomial::DefiniteIntegral<3, 125>::Tri_Definite_Integrals =
        Eigen::Matrix<double, 125, 1>::Zero();
