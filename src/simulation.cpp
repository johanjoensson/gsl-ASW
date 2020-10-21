#include "simulation.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <set>
#include <algorithm>
#include "atom.h"
#include "utils.h"
#include "envelope_fun.h"
#include "structure_const.h"
#include "GSLpp/special_functions.h"
#include "GSLpp/eigen.h"
#include "GSLpp/linalg.h"

void Simulation::set_up_crystal()
{
    double at_vol = 0;
    std::vector<std::vector<Site_t<3>>> nn = cryst.calc_nearest_neighbours();
    // Calculate MT radii
    for(const auto& site : cryst.sites()){
        std::cout << "Nearest neighbor " << nn[site.index()][0].pos() << "\n";
        cryst.atom(site).set_MT(nn[site.index()][0].pos().norm<double>()/2);
        at_vol += GSL::pow_int(cryst.atom(site).get_MT(), 3);
        std::cout << "RMT = " << cryst.atom(site).get_MT() << "\n";
    }
    at_vol *= 4*M_PI/3;
    // Set AS radii, atomic sphere volumes should add up to the crystal volume
    for(const auto& site : cryst.sites()){
        cryst.atom(site).set_AS( std::cbrt(cryst.volume()/at_vol) *
        cryst.atom(site).get_MT());
    }
}

void Simulation::add_states(const Site_t<3>& center, const double kappa)
{
    std::set<lm> p;
    lm l = {1, 0, 0};
    size_t nel = 0;
    size_t Z = cryst.atom(center).get_Z();
    for(; nel < Z || /*l.l > 0*/ l.m > -l.l; l++, nel += 2){
        p.insert(l);
    }
    l = lm {l.n - 1, l.n -2, l.n - 2};
/*    while(l.l <= 3) {
        p.insert(l);
        if(l.m >= l.l){
            l.l++;
            l.m = -l.l;
        }else{
            l.m++;
        }
    }
*/
    nel = 0;
    for(lm tmp : p){
        if(nel + static_cast<size_t>(/*2*tmp.n*tmp.n*/ 2*(2*tmp.l + 1)) >= Z){
            std::cout << "valence : " << tmp << std::endl;
            // basis_valence.push_back(Augmented_spherical_wave(Hs_m, Bs_m, center, kappa, tmp, UP));
            // basis_valence.push_back(Augmented_spherical_wave(Hs_m, Bs_m, center, kappa, tmp, DOWN));
            basis_valence.push_back(Augmented_spherical_wave(center, kappa, tmp, UP, false));
            basis_valence.push_back(Augmented_spherical_wave(center, kappa, tmp, DOWN, false));
        }else{
            std::cout << "core : " << tmp << std::endl;
            basis_core.push_back(Augmented_spherical_wave(center, kappa, tmp, UP, true));
            basis_core.push_back(Augmented_spherical_wave(center, kappa, tmp, DOWN, true));
        }
        nel += 2;
    }
    std::stable_sort(basis_valence.begin(), basis_valence.end(),
        [](const Augmented_spherical_wave& w1, const Augmented_spherical_wave& w2){
            return w1.s() < w2.s();
        });
    std::stable_sort(basis_core.begin(), basis_core.end(),
        [](const Augmented_spherical_wave& w1, const Augmented_spherical_wave& w2){
            return w1.s() < w2.s();
        });
}

void Simulation::set_up_basis()
{
    std::cout << "Setting up logarithmic meshes" << std::endl;
    // Set logarithmic meshes for the atoms and count number of electrons
    for(size_t i = 0; i < cryst.atoms().size(); i++){
        // auto index = std::distance(cryst.atoms().begin(), std::find(cryst.atoms().begin(), cryst.atoms().end(), at));
        // at_meshes[static_cast<size_t>(index)] = Logarithmic_mesh(at.get_AS(),
        at_meshes[i] = Logarithmic_mesh(cryst.atoms(i).get_AS(),
            static_cast<size_t>(1000*(0.5 + static_cast<double>(cryst.atoms(i).get_Z())/118)));
    }
    set_up_augmented_functions();
    // Divide electrons into core and valence states
    for(const auto& site: cryst.sites()){
        for(auto kappa : kappas){
            add_states(site, kappa);
        }
    }
    std::cout << "# valence states " << basis_valence.size() << std::endl;
}

void Simulation::set_up_augmented_functions()
{
    std::cout << "Setting up augmented functions" << std::endl;
    lm lH{1, 0, 0};
    size_t nel = 0;
    for(size_t i = 0; i < cryst.atoms().size(); i++){;
        size_t z = cryst.atoms(i).get_Z();
        Hs_m.push_back(Hankel_container());
        Bs_m.push_back(Bessel_container());
        while(nel < z){
            int l_max = 3;
            if(nel + static_cast<size_t>(2*lH.n*lH.n) < z){
                l_max = lH.n - 1;
            }

            for(int n = lH.n, l = lH.l; l <= std::min(l_max, lH.n - 1); l++){
            for(auto kappa : kappas){
                for(auto s : {spin{UP}, spin{DOWN}}){
                    Hs_m.back().add_function(Augmented_Hankel({n, l, -l}, kappa, s,
                        at_meshes[i]));
                }
            }
            }
            nel += 2*static_cast<size_t>(lH.n*lH.n);
            lH += lH.n*lH.n;
        }
        // lm lB{1, 0, 0};
        nel = 0;
        // while(nel < z){
            // if(nel + static_cast<size_t>(2*lB.n*lB.n) > z){
            // for(int n = lB.n, l = lB.l; l < lB.n; l++){
            for(lm l = {lH.n - 1, 0,  0}; l.n == lH.n - 1; l += 2*l.l + 1){
                for(auto kappa : kappas){
                    for(auto s : {spin{UP}, spin{DOWN}}){
                        Bs_m.back().add_function(Augmented_Bessel(l, kappa, s,
                            at_meshes[i]));
                    }
                }
            // }
            }

            // nel += static_cast<size_t>(2*lB.n*lB.n);
            // lB += lB.n*lB.n;
        // }
        nel = 0;
        lH = {1, 0, 0};
        std::cout << "Hs_m contains " << Hs_m.back().size() << " elements\n";
        std::cout << "Bs_m contains " << Bs_m.back().size() << " elements\n";
    }
}

void Simulation::set_up_potential(const XC_FUN func)
{
    std::cout << "Setting up potential" << std::endl;
    // Set up initial potential
    pot.set_xc_fun(func);
    pot.initial_pot(cryst.volume());
}

void Simulation::init_augmented_functions()
{
    std::cout << "Initialising augmented Hankel functions" << std::endl;
    size_t N = basis_valence.size();
    for(size_t i = 0; i < cryst.atoms().size(); i++){
        for(auto& H : Hs_m[i]){
            std::cout << "H " << H.l() << std::flush;
            bool core = false;
            if(H.s() == UP){
                core = H.l() < basis_valence[0].l();
            }else{
                core = H.l() < basis_valence[N/2].l();
            }
            H.update(pot.sphere(i), -GSL::pow_int(static_cast<double>(cryst.atoms(i).get_Z())/H.l().n, 2), core);
            H.EH() -= pot.MT0();
            std::cout << ", EH = " << H.EH() << "\n";
        }
    }

    std::cout << "Initialising augmented Bessel functions" << std::endl;
    for(size_t i = 0; i < cryst.atoms().size(); i++){
        for(auto& J : Bs_m[i]){
            if(J.l().m != -J.l().l){
                std::cout << "Unexpected\n";
                continue;
            }
            J.update(pot.sphere(i), -GSL::pow_int(static_cast<double>(cryst.atoms(i).get_Z())/J.l().n, 2), false);
            J.EJ() -= pot.MT0();
            std::cout << "J " << J.l() << ", JH = " << J.EJ() << "\n";
        }
    }
}

Simulation::Simulation(const Crystal_t<3, Atom>& crystal, const XC_FUN func, const std::vector<double> kappas_n,
 std::function<double(const size_t, const double)> at_pot)
 : kappas(kappas_n), Hs_m(), Bs_m(), cryst(crystal), at_meshes(cryst.atoms().size()), pot(cryst, at_meshes, at_pot), n_m(cryst, at_meshes), basis_valence(), basis_core(), k_eigenvals(), XH1(),
   XS1(), XH2(), XS2(), XH3(), XS3(), B_m()
{

    set_up_crystal();
    set_up_basis();
    set_up_potential(func);
    init_augmented_functions();

    size_t nH = 0;
    size_t nB = 0;

    auto H_i = std::max_element(Hs_m.begin(), Hs_m.end(),
    [](const Hankel_container& Hs1, const Hankel_container& Hs2){return Hs1.size() < Hs2.size();});
    auto B_i = std::max_element(Bs_m.begin(), Bs_m.end(),
    [](const Bessel_container& Bs1, const Bessel_container& Bs2){return Bs1.size() < Bs2.size();});
    if(H_i != Hs_m.end()){
        nH = H_i->size();
    }
    if(B_i != Bs_m.end()){
        nB = B_i->size();
    }
    XH1 = std::vector<GSL::Matrix>(cryst.sites().size(), GSL::Matrix(nH, nH));
    XS1 = std::vector<GSL::Matrix>(cryst.sites().size(), GSL::Matrix(nH, nH));

    XH2 = std::vector<GSL::Matrix>(cryst.sites().size(), GSL::Matrix(nH, nB));
    XS2 = std::vector<GSL::Matrix>(cryst.sites().size(), GSL::Matrix(nH, nB));

    XH3 = std::vector<GSL::Matrix>(cryst.sites().size(), GSL::Matrix(nB, nB));
    XS3 = std::vector<GSL::Matrix>(cryst.sites().size(), GSL::Matrix(nB, nB));
}

double Simulation::X_H1(const Augmented_Hankel& Ht1, const Augmented_Hankel& Ht2, const Site_t<3>& site)
{
    double res = 0;
    const double k1 = Ht1.kappa(), k2 = Ht2.kappa();
    const Envelope_Hankel H1(Ht1.l(), k1);
    const Envelope_Hankel H2(Ht2.l(), k2);
    res += Ht2.EH()*augmented_integral(Ht1, Ht2);
    res += -k2*k2*off_atomic_integral(H1, H2, at_meshes[cryst.atom_index(site)].r_back());
    return res;
}

double Simulation::X_H2(const Augmented_Hankel& Ht1, const Augmented_Bessel& Jt2, const Site_t<3>& site)
{
    double res = 0;
    const double k1 = Ht1.kappa(), k2 = Jt2.kappa();
    const Envelope_Hankel H1(Ht1.l(), k1);
    const Envelope_Bessel J2(Jt2.l(), k2);
    if(std::abs(Ht1.EH() - Jt2.EJ()) < 5e-4){
        res += Jt2.EJ()*augmented_integral(Ht1, Ht1);
    }else{
//        std::cout << "Error in interal = " << ( 1./(Ht1.EH() - Jt2.EJ()) - augmented_integral(Ht1, Jt2) ) * ( Ht1.EH() - Jt2.EJ() ) << "\n";
        // res += Jt2.EJ()/(Ht1.EH() - Jt2.EJ());
        res += Jt2.EJ()*augmented_integral(Ht1, Jt2);

    }
    res -= -k2*k2*atomic_integral(H1, J2, at_meshes[cryst.atom_index(site)].r_back());
    if(std::abs(k1 - k2) > 1e-10){
        res += -k2*k2/(-k2*k2 + k1*k1);
    }

    return res;
}

double Simulation::X_H3(const Augmented_Bessel& Jt1, const Augmented_Bessel& Jt2, const Site_t<3>& site)
{
    double res = 0;
    const double k1 = Jt1.kappa(), k2 = Jt2.kappa();
    const Envelope_Bessel J1(Jt1.l(), k1);
    const Envelope_Bessel J2(Jt2.l(), k2);
    res += Jt2.EJ()*augmented_integral(Jt1, Jt2);
    res -= -k2*k2*atomic_integral(J1, J2, at_meshes[cryst.atom_index(site)].r_back());
    return res;
}

GSL::Complex Simulation::H_element(const size_t i1, const size_t i2, const GSL::Vector& kp)
{
    const Augmented_spherical_wave w1 = basis_valence[i1];
    const Augmented_spherical_wave w2 = basis_valence[i2];
    Site_t<3> at_i = w1.center(), at_j = w2.center();

    GSL::Complex res = 0.;
    GSL::Vector tau_ij(at_i.pos() - at_j.pos());
    size_t l_i = Hs_m[cryst.atom_index(at_i)].get_index(w1.l(), w1.kappa(), w1.s());
    size_t l_j = Hs_m[cryst.atom_index(at_j)].get_index(w2.l(), w2.kappa(), w2.s());

    const double k1 = w1.kappa(), k2 = w2.kappa();

    if(at_i == at_j && w1.l().l == w2.l().l && w1.l().m == w2.l().m){
        res += XH1[at_i.index()][l_i][l_j];
    }

    if(std::abs(k1 - k2) < 1E-15){
        res += -GSL::pow_int(k2, 2)*
        B_m.dot(cryst, {w1.l().n - 1, 0}, w1.l(), w2.l(), k1, tau_ij, kp);
    }

    res += XH2[at_i.index()][l_i][Bs_m[cryst.atom_index(at_i)].get_index(w1.l(), k2, w1.s())]*
    B_m(cryst, {w2.l().n - 1, 0}, w1.l(), w2.l(), k2, tau_ij, kp);

    res += B_m(cryst, {std::min(4, w1.l().n - 1), 0}, w1.l(), w2.l(), k1, tau_ij, kp)*
        XH2[at_j.index()][l_j][Bs_m[cryst.atom_index(at_j)].get_index(w2.l(), k1, w1.s())];
    if(std::abs(k1 - k2) < 1E-10){
        res += B_m(cryst, {std::min(4, w1.l().n - 1), 0}, w1.l(), w2.l(), k1, tau_ij, kp);
    }

    for(auto& site : cryst.sites()){
	int npp = Bs_m[cryst.atom_index(site)].max_lm().n;
        for(int lpp = 0; lpp <= npp - 1; lpp++){
            for(int mpp = - lpp; mpp <= lpp; mpp++){
                size_t lpp1 = Bs_m[cryst.atom_index(site)].get_index({w1.l().n, lpp, -lpp}, k1, w1.s());
                size_t lpp2 = Bs_m[cryst.atom_index(site)].get_index({w2.l().n, lpp, -lpp}, k2, w2.s());
                GSL::Vector ti = site.pos() - at_i.pos();
                GSL::Vector tj = site.pos() - at_j.pos();
                res += B_m(cryst, {std::min(4, npp - 1), 0}, {w1.l().n, lpp, mpp}, w1.l(), k1, ti, kp).conjugate()*
                    XH3[site.index()][lpp1][lpp2]*
                    B_m(cryst, {std::min(4, npp - 1), 0}, {w2.l().n, lpp, mpp}, w2.l(), k2, tj, kp);
            }
        }
    }
    return res;
}

double Simulation::X_S1(const Augmented_Hankel& Ht1, const Augmented_Hankel& Ht2, const Site_t<3>& site)
{
    const double k1 = Ht1.kappa(), k2 = Ht2.kappa();
    const Envelope_Hankel H1(Ht1.l(), k1);
    const Envelope_Hankel H2(Ht2.l(), k2);
    return augmented_integral(Ht1, Ht2) + off_atomic_integral(H1, H2, at_meshes[cryst.atom_index(site)].r_back());
}

double Simulation::X_S2(const Augmented_Hankel& Ht1, const Augmented_Bessel& Jt2, const Site_t<3>& site)
{
    double res = 0;
    const double k1 = Ht1.kappa(), k2 = Jt2.kappa();
    const Envelope_Hankel H1(Ht1.l(), k1);
    const Envelope_Bessel J2(Jt2.l(), k2);
    if(std::abs(Ht1.EH() - Jt2.EJ()) < 1e-10){
        res += augmented_integral(Ht1, Ht1);
    }else{
        res += 1./(Ht1.EH() - Jt2.EJ());
    }
    res -= atomic_integral(H1, J2, at_meshes[cryst.atom_index(site)].r_back());
    if(std::abs(k1 - k2) > 1e-10){
        res += 1./(-k2*k2 + k1*k1);
    }
    return res;
}

double Simulation::X_S3(const Augmented_Bessel& Jt1, const Augmented_Bessel& Jt2, const Site_t<3>& site)
{
    const double k1 = Jt1.kappa(), k2 = Jt2.kappa();
    const Envelope_Bessel J1(Jt1.l(), k1);
    const Envelope_Bessel J2(Jt2.l(), k2);
    return augmented_integral(Jt1, Jt2) - atomic_integral(J1, J2, at_meshes[cryst.atom_index(site)].r_back());
}

GSL::Complex Simulation::S_element(const size_t i1, const size_t i2, const GSL::Vector& kp)
{
    const Augmented_spherical_wave w1 = basis_valence[i1];
    const Augmented_spherical_wave w2 = basis_valence[i2];
    const Site_t<3> at_i = w1.center(), at_j = w2.center();
    const double k1 = w1.kappa(), k2 = w2.kappa();

    GSL::Complex res = 0.;
    const GSL::Vector tau_ij(at_i.pos() - at_j.pos());
    const size_t l_i = Hs_m[cryst.atom_index(at_i)].get_index(w1.l(), k1, w1.s());
    const size_t l_j = Hs_m[cryst.atom_index(at_j)].get_index(w2.l(), k2, w2.s());

    if(at_i == at_j && lm {w1.l().l, w1.l().m} == lm {w2.l().l, w2.l().m}){
        res += XS1[at_i.index()][l_i][l_j];
    }

    if(std::abs(k1 - k2) < 1E-10){
        res += B_m.dot(cryst, {std::min(4, w1.l().n - 1), 0}, w1.l(), w2.l(), k1, tau_ij, kp);
    }

    res += XS2[at_i.index()][l_i][Bs_m[cryst.atom_index(at_i)].get_index(w1.l(), k2, w2.s())]*
        B_m(cryst, {std::min(4, w2.l().n - 1), 0}, w1.l(), w2.l(), k2, tau_ij, kp);

    res += B_m(cryst, {std::min(4, w1.l().n - 1), 0}, w1.l(), w2.l(), k1, tau_ij, kp)*
        XS2[at_j.index()][l_j][Bs_m[cryst.atom_index(at_j)].get_index(w2.l(), k1, w1.s())];

    for(const auto& site : cryst.sites()){
	int npp = Bs_m[cryst.atom_index(site)].max_lm().n;
        for(int lpp = 0; lpp <= npp - 1; lpp++){
            for(int mpp = - lpp; mpp <= lpp; mpp++){
                size_t lpp1 = Bs_m[cryst.atom_index(site)].get_index({w1.l().n, lpp, -lpp}, k1, w1.s());
                size_t lpp2 = Bs_m[cryst.atom_index(site)].get_index({w2.l().n, lpp, -lpp}, k2, w2.s());
                    GSL::Vector ti = site.pos() - at_i.pos();
                    GSL::Vector tj = site.pos() - at_j.pos();
                    res += B_m(cryst, {std::min(4, npp - 1), 0}, {w1.l().n, lpp, mpp}, w1.l(), k1, ti, kp).conjugate()*
                        XS3[site.index()][lpp1][lpp2]*
                        B_m(cryst, {std::min(4, npp - 1), 0}, {w2.l().n, lpp, mpp}, w2.l(), k2, tj, kp)*0;
            }
        }
    }
    return res;
}


void Simulation::set_up_X_matrices()
{
    for(size_t i = 0; i < basis_valence.size(); i++){
        Augmented_spherical_wave w1 = basis_valence[i];
	const double k1 = w1.kappa();
        Site_t<3> at_i = w1.center();
        size_t i1 = Hs_m[cryst.atom_index(at_i)].get_index(w1.l(), k1, w1.s());
        for(size_t j = 0; j < basis_valence.size(); j++){
            Augmented_spherical_wave w2 = basis_valence[j];
	    const double k2 = w2.kappa();
            Site_t<3> at_j = w2.center();
            size_t i2 = Hs_m[cryst.atom_index(at_j)].get_index(w2.l(), k2, w2.s());

            XH1[at_i.index()][i1][i2] = X_H1(Hs_m[cryst.atom_index(at_i)].get_function(w1.l(), k1, w1.s()), Hs_m[cryst.atom_index(at_j)].get_function(w2.l(), k2, w2.s()), at_i);
            XS1[at_i.index()][i1][i2] = X_S1(Hs_m[cryst.atom_index(at_i)].get_function(w1.l(), k1, w1.s()), Hs_m[cryst.atom_index(at_j)].get_function(w2.l(), k2, w2.s()), at_i);
        }

        for(Augmented_Bessel J2s : Bs_m[cryst.atom_index(at_i)]){
            if(w1.l().l != J2s.l().l){
                continue;
            }
                size_t ipp2 = Bs_m[cryst.atom_index(at_i)].get_index(J2s.l(), J2s.kappa(), J2s.s());
                XH2[at_i.index()][i1][ipp2] = X_H2(Hs_m[cryst.atom_index(at_i)].get_function(w1.l(), k1, w1.s()), J2s, at_i);
                XS2[at_i.index()][i1][ipp2] = X_S2(Hs_m[cryst.atom_index(at_i)].get_function(w1.l(), k1, w1.s()), J2s, at_i);
        }

    }
    for(auto& site : cryst.sites()){
        size_t n_s = cryst.atom_index(site);
        for(const Augmented_Bessel& J1s : Bs_m[n_s]){
            for(const Augmented_Bessel& J2s : Bs_m[n_s]){
                    if(J1s.l().l != J2s.l().l){
                        continue;
                    }
                    size_t ipp1 = Bs_m[n_s].get_index(J1s.l(), J1s.kappa(), J1s.s());
                    size_t ipp2 = Bs_m[n_s].get_index(J2s.l(), J2s.kappa(), J2s.s());
                    XH3[site.index()][ipp1][ipp2] = X_H3(J1s, J2s, site);
                    XS3[site.index()][ipp1][ipp2] = X_S3(J1s, J2s, site);
            }
        }
    }

    std::cout << "XS1 =\n";
    for(auto XS1_i : XS1){
        for(auto row : XS1_i){
            std::cout << row << "\n";
        }
        std::cout << "\n";
    }
    std::cout << "XS2 =\n";
    for(auto XS2_i : XS2){
        for(auto row : XS2_i){
            std::cout << row << "\n";
        }
        std::cout << "\n";
    }
    std::cout << "XS3 =\n";
    for(auto XS3_i : XS3){
        for(auto row :  XS3_i){
            std::cout << row << "\n";
        }
        std::cout << "\n";
    }
}

const GSL::Matrix_cx Simulation::set_H(const GSL::Vector& kp)
{
    size_t N = basis_valence.size();
    GSL::Matrix_cx H(N, N);
    for(size_t i = 0; i < N; i++){
        for(size_t j = 0; j <= i; j++){
        // for(size_t j = 0; j < N; j++){
            H[i][j] = H_element(i, j, kp);
            H[j][i] = H_element(i, j, kp).conjugate();
        }
    }
    return H;
}

const GSL::Matrix_cx Simulation::set_S(const GSL::Vector& kp)
{
    size_t N = basis_valence.size();
    GSL::Matrix_cx S(N, N);
    for(size_t i = 0; i < N; i++){
        for(size_t j = 0; j <= i; j++){
        // for(size_t j = 0; j < N; j++){

            S[i][j] = S_element(i, j, kp);
            S[j][i] = S_element(i, j, kp).conjugate();
        }
    }
    return S;
}

void Simulation::calc_eigen(const GSL::Vector& kp)
{
    const GSL::Matrix_cx H {set_H(kp)};
    const GSL::Matrix_cx S {set_S(kp)};
    size_t N = basis_valence.size()/2;
    GSL::Matrix_cx eigvecs_up(N, N);
    GSL::Vector eigvals_up(N);

    GSL::Matrix_cx::View H_up = H.view(0,0, N, N);
    GSL::Matrix_cx::View S_up = S.view(0, 0, N, N);;

    for(size_t i = 0; i < N; i++){
        for(size_t j = 0; j < N; j++){
            if(!(basis_valence[i].s() == UP && basis_valence[j].s() == UP)){
                std::string msg("Valence function ");
                if(basis_valence[i].s() != UP){
                    msg += std::to_string(i);
                }else{
                    msg += std::to_string(j);
                }
                msg += " is not spin UP.\n";
                throw std::runtime_error(msg + "FATAL ERROR!\n");
            }
        }
    }

    std::pair<GSL::Matrix_cx, GSL::Vector> tmp, overlap;
    overlap = GSL::hermitian_eigen(S_up);
    std::cout << "k-point " << kp << "\n";
    std::cout << " Overlap matrix\n";
    for(size_t i = 0 ; i < N; i++){
        std::cout << "  " << S_up[i] << "\n";
    }
    std::cout << "Eigenvalues of Overlap matrix\n" << overlap.second << "\n";
    try{
        tmp = GSL::general_hermitian_eigen(H_up, S_up);
        eigvecs_up = tmp.first;
        eigvals_up = tmp.second;
        std::cout << "\tEigenvalues (up):\n\t" << eigvals_up << " (Ry)\n";
        std::cout << "\tEigenvectors (up) :\n";
        for(auto row : eigvecs_up.T()){
            std::cout << "\t\t" << row << "\n";
        }
        std::cout << "\n";
        // k_eigenvals[kp] = GSL::sort_general_hermitian_eigen(eigvals_up, eigvecs_up);
        k_eigenvals[kp] = {eigvecs_up, eigvals_up};

    }catch (const std::runtime_error &e){
		std::cerr << e.what();
		// std::cerr << " Hamiltonian matrix\n";
		// for(size_t i = 0 ; i < N; i++){
		// 	std::cerr << "  " << H_up[i] << "\n";
		// }
        std::cerr << "Setting eigenvalues to "<< nan("") <<"\n";
        for(auto& eigval : eigvals_up){
            eigval = nan("");
        }
    }
    k_eigenvals[kp] = {eigvecs_up, eigvals_up};
}

void Simulation::print_eigvals(const K_mesh& k_mesh)
{
    std::ofstream of;
    of.open("bands.dat", std::ios::trunc);
    for(const auto& kp :  k_mesh.k_points){
        auto it = k_eigenvals.find(kp);
		if(it != k_eigenvals.end()){
				auto eigvals = it->second.second;
				of << kp[0] << " " << kp[1] << " " << kp[2] << " ";
				for(size_t i = 0; i < eigvals.dim(); i++){
						of << eigvals[i] << " ";
				}
				of << "\n";
		}
    }
    of.close();
}

double Simulation::canonical_band(const lm l, const double kappa, const spin s, const GSL::Vector& kp)
{
    auto h = Hs_m[0].get_function(l, kappa, s);
    auto j = Bs_m[0].get_function(l, kappa, s);
    const double eh = h.EH(), ej = j.EJ();

    const double shh = augmented_integral(h, h), sjh = std::abs(eh - ej) < 5e-4 ? augmented_integral(h, h) : 1./(eh - ej);

    lm lint = {4, 0};
    GSL::Complex B = B_m(this->cryst, lint, l, l, kappa, GSL::Vector(3), kp);
    if(std::abs(B.im()) > 1e-12 ){
        std::cerr << "Imaginary part of Bloch summed structure constant is too large!\n";
    }
    double tmp = sjh/shh*B.re();
    return (eh - ej*tmp)/(1 - tmp);
}
