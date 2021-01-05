#include "../src/structure_const.h"
#include "../src/crystal.h"
#include "../src/atom.h"
#include "../src/log_mesh.h"
#include "GSLpp/vector.h"
#include <gtest/gtest.h>
#include <array>

#define TOL 1e-12
#define KAPPA std::sqrt(0.015)
namespace {
std::array<GSL::Vector, 5> R_vecs = {GSL::Vector{0, 0, 0}, {0.5, 0, 0}, {0, 0.5, 0},
                                     {0, 0, 0.5}, {0.5, 0.5, 0.5}};
std::array<GSL::Vector, 5> K_vecs = {GSL::Vector{0, 0, 0}, {0.5, 0, 0}, {0, 0.5, 0},
                                     {0, 0, 0.5}, {0.5, 0.5, 0.5}};

Crystal_t<3, Atom> set_up_crystal()
{
//    double kappa = KAPPA;
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

TEST(BlochStructureConstants, LExchange)
{
    Crystal_t<3, Atom> cr = set_up_crystal();
    Bloch_summed_structure_constant B1;
    Bloch_summed_structure_constant B2;

    for(lm l1 = {4, 0, 0}; l1 != lm {5, 0, 0}; l1++){
        for(lm l2 = {4, 0, 0}; l2 != lm {5, 0, 0}; l2++){
            for(const auto& tau : R_vecs){
                for(const auto& k : K_vecs){
                    EXPECT_NEAR(B1(cr, {4, 4}, l1, l2, KAPPA, tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), GSL::pow_int(-1, l1.l - l2.l)*B2(cr, {4, 4}, l2, l1, KAPPA, tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), TOL)
                        << "l1  = " << l1 << ", l2 = " << l2 << "\n tau = " << tau << ", k = " << k << "\n";
                    EXPECT_NEAR(B1(cr, {4, 4}, l1, l2, KAPPA, tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), GSL::pow_int(-1, l1.l - l2.l)*B2(cr, {4, 4}, l2, l1, KAPPA, tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), TOL);
                }
            }
        }
    }
}

TEST(BlochStructureConstants, Negate)
{
    Crystal_t<3, Atom> cr = set_up_crystal();
    Bloch_summed_structure_constant B;
    for(lm l1 = {4, 0, 0}; l1 != lm {5, 0, 0}; l1++){
        for(lm l2 = {4, 0, 0}; l2 != lm {5, 0, 0}; l2++){
            for(const auto& tau : R_vecs){
                for(const auto& k : K_vecs){
                    ASSERT_NEAR(B(cr, {4, 0}, l1, l2, KAPPA, -tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), GSL::pow_int(-1, l1.l - l2.l)*B(cr, {4, 0}, l1, l2, KAPPA, tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re(), TOL);
                    ASSERT_NEAR(B(cr, {4, 0}, l1, l2, KAPPA, -tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), GSL::pow_int(-1, l1.l - l2.l)*B(cr, {4, 0}, l1, l2, KAPPA, tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im(), TOL);
                }
            }
        }
    }
}

TEST(BlochStructureConstants, Conjugate)
{
    Crystal_t<3, Atom> cr = set_up_crystal();
    Bloch_summed_structure_constant B;
    for(lm l1 = {4, 0, 0}; l1 != lm {5, 0, 0}; l1++){
        for(lm l2 = {4, 0, 0}; l2 != lm {5, 0, 0}; l2++){
            for(const auto& tau : R_vecs){
                for(const auto& k : K_vecs){
                    ASSERT_NEAR(B(cr, {4, 0}, l1, l2, KAPPA, tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), B(cr, {4, 0}, l1, l2, KAPPA, tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re(), TOL);
                    ASSERT_NEAR(B(cr, {4, 0}, l1, l2, KAPPA, tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), -B(cr, {4, 0}, l1, l2, KAPPA, tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im(), TOL);
                }
            }
        }
    }
}

TEST(BlochStructureConstants, LExchangeDot)
{
    Crystal_t<3, Atom> cr = set_up_crystal();
    Bloch_summed_structure_constant B1;
    Bloch_summed_structure_constant B2;
    for(lm l1 = {4, 0, 0}; l1 != lm {5, 0, 0}; l1++){
        for(lm l2 = {4, 0, 0}; l2 != lm {5, 0, 0}; l2++){
            for(const auto& tau : R_vecs){
                for(const auto& k : K_vecs){
                    ASSERT_NEAR(B1(cr, {4, 0}, l1, l2, KAPPA, tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), GSL::pow_int(-1, l1.l - l2.l)*B2(cr, {4, 0}, l2, l1, KAPPA, tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), TOL);
                    ASSERT_NEAR(B1(cr, {4, 0}, l1, l2, KAPPA, tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), GSL::pow_int(-1, l1.l - l2.l)*B2(cr, {4, 0}, l2, l1, KAPPA, tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), TOL);
                }
            }
        }
    }
}

TEST(BlochStructureConstants, NegateDot)
{
    Crystal_t<3, Atom> cr = set_up_crystal();
    Bloch_summed_structure_constant B;
    for(lm l1 = {4, 0, 0}; l1 != lm {5, 0, 0}; l1++){
        for(lm l2 = {4, 0, 0}; l2 != lm {5, 0, 0}; l2++){
            for(const auto& tau : R_vecs){
                for(const auto& k : K_vecs){
                    ASSERT_NEAR(B.dot(cr, {4, 0}, l1, l2, KAPPA, -tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), GSL::pow_int(-1, l1.l - l2.l)*B.dot(cr, {4, 0}, l1, l2, KAPPA, tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re(), TOL);
                    ASSERT_NEAR(B.dot(cr, {4, 0}, l1, l2, KAPPA, -tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), GSL::pow_int(-1, l1.l - l2.l)*B.dot(cr, {4, 0}, l1, l2, KAPPA, tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im(), TOL);
                }
            }
        }
    }
}

TEST(BlochStructureConstants, ConjugateDot)
{
    Crystal_t<3, Atom> cr = set_up_crystal();
    Bloch_summed_structure_constant B;
    for(lm l1 = {4, 0, 0}; l1 != lm {5, 0, 0}; l1++){
        for(lm l2 = {4, 0, 0}; l2 != lm {5, 0, 0}; l2++){
            for(const auto& tau : R_vecs){
                for(const auto& k : K_vecs){
                    ASSERT_NEAR(B.dot(cr, {4, 0}, l1, l2, KAPPA, tau*cr.lat().lat(), k*cr.lat().recip_lat()).re(), B.dot(cr, {4, 0}, l1, l2, KAPPA, tau*cr.lat().lat(), -k*cr.lat().recip_lat()).re(), TOL);
                    ASSERT_NEAR(B.dot(cr, {4, 0}, l1, l2, KAPPA, tau*cr.lat().lat(), k*cr.lat().recip_lat()).im(), -B.dot(cr, {4, 0}, l1, l2, KAPPA, tau*cr.lat().lat(), -k*cr.lat().recip_lat()).im(), TOL);
                }
            }
        }
    }
}
