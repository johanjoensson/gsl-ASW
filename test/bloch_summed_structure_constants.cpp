#include "../src/structure_const.h"
#include "../src/crystal.h"
#include "../src/atom.h"
#include "../src/log_mesh.h"
#include "GSLpp/vector.h"
#include <gtest/gtest.h>
#include <array>

namespace {
std::array<GSL::Vector, 5> R_vecs = {GSL::Vector{0, 0, 0}, {0.5, 0, 0}, {0, 0.5, 0},
                                     {0, 0, 0.5}, {0.5, 0.5, 0.5}};
std::array<GSL::Vector, 5> K_vecs = {GSL::Vector{0, 0, 0}, {0.5, 0, 0}, {0, 0.5, 0},
                                     {0, 0, 0.5}, {0.5, 0.5, 0.5}};

Crystal_t<3, Atom> set_up_crystal()
{
    double kappa = std::sqrt(0.015);
    GSL::Vector a = {2.0, 0.0, 0.0}, b = {0.0, 2.0, 0.0}, c = {0.0, 0.0, 2.0};
    Crystal_t<3, Atom> cr(Lattice_t<3>{a, b, c});
    double Rmax = calc_Rmax(cr.volume(), kappa, lm {5, 0}, 5e-14);
	double Kmax = calc_Kmax(cr.volume(), kappa, lm {5, 0}, 5e-14);
    cr.set_Rn(Rmax);
    cr.set_Kn(Kmax);

	Atom at{Logarithmic_mesh(), {0,0,0}};
    at.set_Z(1);
    cr.set_size({1, 1, 1});
	cr.add_sites({{0, 0, 0}});
	cr.add_basis({at});

    return cr;
}
}

TEST(BlochStructureConstants, LExchange1)
{
    Crystal_t<3, Atom> cr = set_up_crystal();
    Bloch_summed_structure_constant B1(cr, {0,0}, {1,0});
    Bloch_summed_structure_constant B2(cr, {1,0}, {0,0});

    for(GSL::Vector tau : R_vecs){
        for(GSL::Vector k : K_vecs){
            ASSERT_DOUBLE_EQ(B1(tau, k).re(), -B2(tau, k).re());
            ASSERT_DOUBLE_EQ(B1(tau, k).im(), -B2(tau, k).im());
        }
    }
}

TEST(BlochStructureConstants, LExchange2)
{
    Crystal_t<3, Atom> cr = set_up_crystal();
    Bloch_summed_structure_constant B1(cr, {2,0}, {1,0});
    Bloch_summed_structure_constant B2(cr, {1,0}, {2,0});

    for(GSL::Vector tau : R_vecs){
        for(GSL::Vector k : K_vecs){
            ASSERT_DOUBLE_EQ(B1(tau, k).re(), -B2(tau, k).re());
            ASSERT_DOUBLE_EQ(B1(tau, k).im(), -B2(tau, k).im());
        }
    }
}

TEST(BlochStructureConstants, LExchange3)
{
    Crystal_t<3, Atom> cr = set_up_crystal();
    Bloch_summed_structure_constant B1(cr, {2,0}, {0,0});
    Bloch_summed_structure_constant B2(cr, {0,0}, {2,0});

    for(GSL::Vector tau : R_vecs){
        for(GSL::Vector k : K_vecs){
            ASSERT_DOUBLE_EQ(B1(tau, k).re(), B2(tau, k).re());
            ASSERT_DOUBLE_EQ(B1(tau, k).im(), B2(tau, k).im());
        }
    }
}

TEST(BlochStructureConstants, Negate1)
{
    Crystal_t<3, Atom> cr = set_up_crystal();
    Bloch_summed_structure_constant B(cr, {1,0}, {0,0});

    for(GSL::Vector tau : R_vecs){
        for(GSL::Vector k : K_vecs){
            ASSERT_NEAR(B(-tau, k).re(), -B(tau, -k).re(), 5e-14);
            ASSERT_NEAR(B(-tau, k).im(), -B(tau, -k).im(), 5e-14);
        }
    }
}

TEST(BlochStructureConstants, Negate2)
{
    Crystal_t<3, Atom> cr = set_up_crystal();
    Bloch_summed_structure_constant B(cr, {1,0}, {2,0});

    for(GSL::Vector tau : R_vecs){
        for(GSL::Vector k : K_vecs){
            ASSERT_NEAR(B(-tau, k).re(), -B(tau, -k).re(), 5e-14);
            ASSERT_NEAR(B(-tau, k).im(), -B(tau, -k).im(), 5e-14);
        }
    }
}

TEST(BlochStructureConstants, Negate3)
{
    Crystal_t<3, Atom> cr = set_up_crystal();
    Bloch_summed_structure_constant B(cr, {2,0}, {0,0});

    for(GSL::Vector tau : R_vecs){
        for(GSL::Vector k : K_vecs){
            ASSERT_NEAR(B(-tau, k).re(), B(tau, -k).re(), 5e-14);
            ASSERT_NEAR(B(-tau, k).im(), B(tau, -k).im(), 5e-14);
        }
    }
}

TEST(BlochStructureConstants, Conjugate1)
{
    Crystal_t<3, Atom> cr = set_up_crystal();
    Bloch_summed_structure_constant B(cr, {1,0}, {0,0});

    for(GSL::Vector tau : R_vecs){
        for(GSL::Vector k : K_vecs){
            ASSERT_NEAR(B(tau, k).re(), B(tau, -k).re(), 5e-14);
            ASSERT_NEAR(B(tau, k).im(), -B(tau, -k).im(), 5e-14);
        }
    }
}

TEST(BlochStructureConstants, Conjugate2)
{
    Crystal_t<3, Atom> cr = set_up_crystal();
    Bloch_summed_structure_constant B(cr, {2,0}, {1,0});

    for(GSL::Vector tau : R_vecs){
        for(GSL::Vector k : K_vecs){
            ASSERT_NEAR(B(tau, k).re(), B(tau, -k).re(), 5e-14);
            ASSERT_NEAR(B(tau, k).im(), -B(tau, -k).im(), 5e-14);
        }
    }
}

TEST(BlochStructureConstants, Conjugate3)
{
    Crystal_t<3, Atom> cr = set_up_crystal();
    Bloch_summed_structure_constant B(cr, {2,0}, {0,0});

    for(GSL::Vector tau : R_vecs){
        for(GSL::Vector k : K_vecs){
            ASSERT_NEAR(B(tau, k).re(), B(tau, -k).re(), 5e-14);
            ASSERT_NEAR(B(tau, k).im(), -B(tau, -k).im(), 5e-14);
        }
    }
}
