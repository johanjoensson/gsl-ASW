#include "../src/bloch_sum.h"
#include "../src/crystal.h"
#include "../src/atom.h"
#include "../src/log_mesh.h"
#include "GSLpp/vector.h"
#include <gtest/gtest.h>
#include <array>
#include <iostream>

namespace {
const std::array<GSL::Vector, 5> R_vecs = {GSL::Vector{0, 0, 0}, {0.5, 0, 0}, {0, 0.5, 0},
                                     {0, 0, 0.5}, {0.5, 0.5, 0.5}};
const std::array<GSL::Vector, 5> K_vecs = {GSL::Vector{0, 0, 0}, {0.5, 0, 0}, {0, 0.5, 0},
                                     {0, 0, 0.5}, {0.5, 0.5, 0.5}};
Crystal_t<3, Atom> set_up_crystal()
{
    double kappa = std::sqrt(0.015);
    GSL::Vector a = {3.0, 0.0, 0.0}, b = {0.0, 3.0, 0.0}, c = {0.0, 0.0, 3.0};
    Crystal_t<3, Atom> cr(Lattice_t<3>{a, b, c});
    double Rmax = calc_Rmax(cr.volume(), kappa, lm {5, 0}, 1e-14);
    double Kmax = calc_Kmax(cr.volume(), kappa, lm {5, 0}, 1e-14);
    cr.set_Rn(Rmax);
    cr.set_Kn(Kmax);

    Atom at{{0,0,0}};
    at.set_Z(1);
    cr.set_size({1, 1, 1});
	cr.add_sites({{0, 0, 0}});
	cr.add_basis({at});

    return cr;
}
}

TEST(BlochSum, Negate0)
{
    Crystal_t<3, Atom> cr{set_up_crystal()};
    Bloch_sum D;

    for(auto tau : R_vecs){
        for(auto k : K_vecs){
            // EXPECT_DOUBLE_EQ(D(-tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), D(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re());
            // EXPECT_DOUBLE_EQ(D(-tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), D(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im());
            // ASSERT_DOUBLE_EQ(D(-tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), D(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re());
            // ASSERT_DOUBLE_EQ(D(-tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), D(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im());
            ASSERT_NEAR(D({0, 0}, std::sqrt(0.015), cr, -tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), D({0, 0}, std::sqrt(0.015), cr, tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re(), 5e-12);
            ASSERT_NEAR(D({0, 0}, std::sqrt(0.015), cr, -tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), D({0, 0}, std::sqrt(0.015), cr, tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im(), 5e-12);
        }
    }
}

TEST(BlochSum, Negate1)
{
    Crystal_t<3, Atom> cr = set_up_crystal();
    Bloch_sum D;

    for(GSL::Vector tau : R_vecs){
        for(GSL::Vector k : K_vecs){
            // ASSERT_DOUBLE_EQ(D(-tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), -D(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re());
            // ASSERT_DOUBLE_EQ(D(-tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), -D(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im());
            ASSERT_NEAR(D({1,0}, std::sqrt(0.015), cr, -tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), -D({1,0}, std::sqrt(0.015), cr, tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re(), 5e-12);
            ASSERT_NEAR(D({1,0}, std::sqrt(0.015), cr, -tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), -D({1,0}, std::sqrt(0.015), cr, tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im(), 5e-12);
        }
    }
}

TEST(BlochSum, Negate2)
{
    Crystal_t<3, Atom> cr = set_up_crystal();
    Bloch_sum D;

    for(GSL::Vector tau : R_vecs){
        for(GSL::Vector k : K_vecs){
           // ASSERT_DOUBLE_EQ(D(-tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), -D(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re());
           // ASSERT_DOUBLE_EQ(D(-tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), -D(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im());
            ASSERT_NEAR(D({1,1}, std::sqrt(0.015), cr, -tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), -D({1,1}, std::sqrt(0.015), cr, tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re(), 5e-12);
            ASSERT_NEAR(D({1,1}, std::sqrt(0.015), cr, -tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), -D({1,1}, std::sqrt(0.015), cr, tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im(), 5e-12);
        }
    }
}

TEST(BlochSum, Negate3)
{
    Crystal_t<3, Atom> cr = set_up_crystal();
    Bloch_sum D;

    for(GSL::Vector tau : R_vecs){
        for(GSL::Vector k : K_vecs){
           // ASSERT_DOUBLE_EQ(D(-tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), D(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re());
           // ASSERT_DOUBLE_EQ(D(-tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), D(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im());
            ASSERT_NEAR(D({2,0}, std::sqrt(0.015), cr, -tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), D({2,0}, std::sqrt(0.015), cr, tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re(), 5e-12);
            ASSERT_NEAR(D({2,0}, std::sqrt(0.015), cr, -tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), D({2,0}, std::sqrt(0.015), cr, tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im(), 5e-12);
        }
    }
}

TEST(BlochSum, Conjugate0)
{
    Crystal_t<3, Atom> cr = set_up_crystal();
    Bloch_sum D;

    for(GSL::Vector tau : R_vecs){
        for(GSL::Vector k : K_vecs){
           // ASSERT_DOUBLE_EQ(D(tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), D(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re());
           // ASSERT_DOUBLE_EQ(D(tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), -D(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im());
            ASSERT_NEAR(D({0,0}, std::sqrt(0.015), cr, tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), D({0,0}, std::sqrt(0.015), cr, tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re(), 5e-12);
            ASSERT_NEAR(D({0,0}, std::sqrt(0.015), cr, tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), -D({0,0}, std::sqrt(0.015), cr, tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im(), 5e-12);
        }
    }
}

TEST(BlochSum, Conjugate1)
{
    Crystal_t<3, Atom> cr = set_up_crystal();
    Bloch_sum D;

    for(GSL::Vector tau : R_vecs){
        for(GSL::Vector k : K_vecs){
           // ASSERT_DOUBLE_EQ(D(tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), D(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re());
           // ASSERT_DOUBLE_EQ(D(tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), -D(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im());
            ASSERT_NEAR(D({1,0}, std::sqrt(0.015), cr, tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), D({1,0}, std::sqrt(0.015), cr, tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re(), 5e-12);
            ASSERT_NEAR(D({1,0}, std::sqrt(0.015), cr, tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), -D({1,0}, std::sqrt(0.015), cr, tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im(), 5e-12);
        }
    }
}

TEST(BlochSum, Conjugate2)
{
    Crystal_t<3, Atom> cr = set_up_crystal();
    Bloch_sum D;

    for(GSL::Vector tau : R_vecs){
        for(GSL::Vector k : K_vecs){
           // ASSERT_DOUBLE_EQ(D(tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), D(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re());
           // ASSERT_DOUBLE_EQ(D(tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), -D(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im());
            ASSERT_NEAR(D({1,1}, std::sqrt(0.015), cr, tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), D({1,1}, std::sqrt(0.015), cr, tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re(), 5e-12);
            ASSERT_NEAR(D({1,1}, std::sqrt(0.015), cr, tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), -D({1,1}, std::sqrt(0.015), cr, tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im(), 5e-12);
        }
    }
}

TEST(BlochSum, Conjugate3)
{
    Crystal_t<3, Atom> cr = set_up_crystal();
    Bloch_sum D;

    for(GSL::Vector tau : R_vecs){
        for(GSL::Vector k : K_vecs){
           // ASSERT_DOUBLE_EQ(D(tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), D(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re());
           // ASSERT_DOUBLE_EQ(D(tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), -D(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im());
            ASSERT_NEAR(D({2,0}, std::sqrt(0.015), cr, tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), D({2,0}, std::sqrt(0.015), cr, tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re(), 5e-12);
            ASSERT_NEAR(D({2,0}, std::sqrt(0.015), cr, tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), -D({2,0}, std::sqrt(0.015), cr, tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im(), 5e-12);
        }
    }
}
