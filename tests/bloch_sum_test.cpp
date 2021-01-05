#include "../src/bloch_sum.h"
#include "../src/crystal.h"
#include "../src/atom.h"
#include "../src/log_mesh.h"
#include "GSLpp/vector.h"
#include <gtest/gtest.h>
#include <array>
#include <iostream>

#define TOL 1e-10
#define KAPPA std::sqrt(0.015)

namespace {
const std::array<GSL::Vector, 7> R_vecs = {GSL::Vector{0.001, 0, 0}, {0, 0.001, 0},
                                     {0, 0, 0.001}, {0.001, 0.001, 0}, {0.001, 0, 0.001}, {0, 0.001, 0.001}, {0.001, 0.001, 0.001}};
const std::array<GSL::Vector, 7> K_vecs = {GSL::Vector{0.001, 0, 0}, {0, 0.001, 0},
                                     {0, 0, 0.001}, {0.001, 0.001, 0}, {0.001, 0, 0.001}, {0, 0.001, 0.001}, {0.001, 0.001, 0.001}};
Crystal_t<3, Atom> set_up_crystal()
{
    double kappa = KAPPA;
    GSL::Vector a = {3.0, 0.0, 0.0}, b = {0.0, 3.0, 0.0}, c = {0.0, 0.0, 3.0};
    Crystal_t<3, Atom> cr(Lattice_t<3>{a, b, c});
    double Rmax = calc_Rmax(cr.volume(), kappa, lm {6, 0}, 1e-14);
    double Kmax = calc_Kmax(cr.volume(), kappa, lm {6, 0}, 1e-14);
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

TEST(BlochSum, Negate)
{
    Crystal_t<3, Atom> cr{set_up_crystal()};
    Bloch_sum D1, D2;
    for(lm l = {4, 0, 0}; l != lm {5, 0, 0}; l++){
        for(int i = -1000; i <= 1000; i += 50){
            for(auto tau : R_vecs){
                for(auto k : K_vecs){
                            // EXPECT_DOUBLE_EQ(D(-tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), D(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re());
                            // EXPECT_DOUBLE_EQ(D(-tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), D(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im());
                            // ASSERT_DOUBLE_EQ(D(-tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), D(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re());
                            // ASSERT_DOUBLE_EQ(D(-tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), D(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im());
                        auto val1 = D1(l, KAPPA, cr, -i*(tau*cr.lat().lat()), i*(k*cr.lat().recip_lat()));
                        auto val2 = GSL::pow_int(-1, l.l)*D2(l, KAPPA, cr, i*(tau*cr.lat().lat()), -i*(k*cr.lat().recip_lat()));
                        ASSERT_NEAR(val1.re(), val2.re(), TOL);
                        ASSERT_NEAR(val1.im(), val2.im(), TOL);
                }
            }
        }
    }

}

TEST(BlochSum, Conjugate)
{
    Crystal_t<3, Atom> cr = set_up_crystal();
    Bloch_sum D;
    for(lm l = {4, 0, 0}; l != lm {5, 0, 0}; l++){
        for(int i = -1000; i <= 1000; i += 50){
            for(GSL::Vector tau : R_vecs){
                for(GSL::Vector k : K_vecs){
                            // ASSERT_DOUBLE_EQ(D(tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), D(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re());
                            // ASSERT_DOUBLE_EQ(D(tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), -D(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im());
                    ASSERT_NEAR(D(l, KAPPA, cr, i*tau*cr.lat().lat(), i*k*cr.lat().recip_lat()).re(), D(l, KAPPA, cr, i*tau*cr.lat().lat(), -i*k*cr.lat().recip_lat()).re(), TOL);
                    ASSERT_NEAR(D(l, KAPPA, cr, i*tau*cr.lat().lat(), i*k*cr.lat().recip_lat()).im(), -D(l, KAPPA, cr, i*tau*cr.lat().lat(), -i*k*cr.lat().recip_lat()).im(), TOL);
                }
            }
        }
    }
}

TEST(BlochSum, NegateDot)
{
    Crystal_t<3, Atom> cr{set_up_crystal()};
    Bloch_sum D1, D2;

    for(lm l = {4, 0, 0}; l != lm {5, 0, 0}; l += 50){
        for(int i = -1000; i <= 1000; i++){
            for(auto tau : R_vecs){
                    for(auto k : K_vecs){
                            // EXPECT_DOUBLE_EQ(D.dot(-tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), D.dot(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re());
                            // EXPECT_DOUBLE_EQ(D.dot(-tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), D.dot(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im());
                            // ASSERT_DOUBLE_EQ(D.dot(-tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), D.dot(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re());
                            // ASSERT_DOUBLE_EQ(D.dot(-tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), D.dot(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im());
                        auto val1 = D1.dot(l, KAPPA, cr, -i*(tau*cr.lat().lat()), i*(k*cr.lat().recip_lat()));
                        auto val2 = GSL::pow_int(-1, l.l)*D2.dot(l, KAPPA, cr, i*(tau*cr.lat().lat()), -i*(k*cr.lat().recip_lat()));
                        ASSERT_NEAR(val1.re(), val2.re(), TOL);
                        ASSERT_NEAR(val1.im(), val2.im(), TOL);
                    }
            }
        }
    }

}

TEST(BlochSum, ConjugateDot)
{
    Crystal_t<3, Atom> cr = set_up_crystal();
    Bloch_sum D;

    for(lm l = {4, 0, 0}; l != lm {5, 0, 0}; l++){
        for(int i = -1000; i <= 1000; i += 50){
            for(GSL::Vector tau : R_vecs){
                    for(GSL::Vector k : K_vecs){
                            // ASSERT_DOUBLE_EQ(D.dot(tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), D.dot(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re());
                            // ASSERT_DOUBLE_EQ(D.dot(tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), -D.dot(tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im());
                            ASSERT_NEAR(D.dot(l, KAPPA, cr, i*tau*cr.lat().lat(), i*k*cr.lat().recip_lat()).re(), D.dot(l, KAPPA, cr, i*tau*cr.lat().lat(), -i*k*cr.lat().recip_lat()).re(), TOL);
                            ASSERT_NEAR(D.dot(l, KAPPA, cr, i*tau*cr.lat().lat(), i*k*cr.lat().recip_lat()).im(), -D.dot(l, KAPPA, cr, i*tau*cr.lat().lat(), -i*k*cr.lat().recip_lat()).im(), TOL);
                    }
        }    }
    }
}
