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
    for(auto site : cryst.sites()){
        cryst.atom(site).set_MT(nn[site.index()][0].pos().norm<double>()/2);
        at_vol += GSL::pow_int(cryst.atom(site).get_MT(), 3);
        std::cout << "RMT = " << cryst.atom(site).get_MT() << "\n";
    }
    at_vol *= 4*M_PI/3;
    // Set AS radii, atomic sphere volumes should add up to the crystal volume
    for(auto site : cryst.sites()){
        cryst.atom(site).set_AS( std::cbrt(cryst.volume()/at_vol) *
        cryst.atom(site).get_MT());
    }
}

void Simulation::add_states(const Atom& center, const double kappa)
{
    std::set<lm> p;
    lm l{1, 0, 0};
    size_t nel = 0;
    size_t Z = center.get_Z();
    while(nel <= Z || l.l < l.n - 1){
        p.insert(l);
        l++;
        nel += 2;
    }
    for(; l.m < l.l; l++){
        p.insert(l);
    }
    p.insert(l);
    nel = 0;
    auto center_index = std::distance(cryst.atoms().begin(), std::find(cryst.atoms().begin(), cryst.atoms().end(), center));
    for(lm tmp : p){
        if(nel + static_cast<size_t>(2*tmp.n*tmp.n) >= Z){
            std::cout << "valence : " << tmp << std::endl;
            basis_valence.push_back(Augmented_spherical_wave(Hs_m, Bs_m, static_cast<size_t>(center_index), kappa, tmp, UP));
            basis_valence.push_back(Augmented_spherical_wave(Hs_m, Bs_m, static_cast<size_t>(center_index), kappa, tmp, DOWN));
        }else{
            std::cout << "core : " << tmp << std::endl;
            basis_core.push_back(Augmented_spherical_wave(Hs_m, Bs_m, static_cast<size_t>(center_index), kappa, tmp, UP));
            basis_core.push_back(Augmented_spherical_wave(Hs_m, Bs_m, static_cast<size_t>(center_index), kappa, tmp, DOWN));
        }
        nel += 2;
    }
}

void Simulation::set_up_basis()
{
    std::cout << "Setting up logarithmic meshes" << std::endl;
    // Set logarithmic meshes for the atoms and count number of electrons
    for(auto at : cryst.atoms()){
        // auto index = std::distance(cryst.atoms().begin(), std::find(cryst.atoms().begin(), cryst.atoms().end(), at));
        // at_meshes[static_cast<size_t>(index)] = Logarithmic_mesh(at.get_AS(),
        at_meshes.push_back(Logarithmic_mesh(at.get_AS(),
            500 + static_cast<size_t>(1000*std::ceil(
                static_cast<double>(at.get_Z())/118)) ));
    }
    set_up_augmented_functions();
    // Divide electrons into core and valence states
    for(const Atom& at: cryst.atoms()){
        for(auto kappa : kappas){
            add_states(at, kappa);
        }
    }
    std::cout << "# valence states " << basis_valence.size() << std::endl;

    std::cout << "Setting up core states." << std::endl;
    for(Augmented_spherical_wave& wave : basis_core){
        wave.core_state = true;
        // wave.set_up(pot);
    }
    std::cout << "Setting up valence states." << std::endl;
    for(Augmented_spherical_wave& wave: basis_valence){
        wave.core_state = false;
        // wave.set_up(pot);
    }
}

void Simulation::set_up_augmented_functions()
{
    std::cout << "Setting up augmented functions" << std::endl;
    lm lH{1, 0, 0};
    size_t nel = 0;
    for(auto& at : cryst.atoms()){
        auto at_index = std::distance(cryst.atoms().begin(), std::find(cryst.atoms().begin(), cryst.atoms().end(), at));
        Hs_m.push_back(Hankel_container());
        Bs_m.push_back(Bessel_container());
        while(nel <= at.get_Z() || lH.l > 0){
            for(auto kappa : kappas){
                for(auto s : {spin{UP}, spin{DOWN}}){
                    Hs_m.back().add_function(Augmented_Hankel(lH, kappa, s, at.pos,
                        at_meshes[static_cast<size_t>(at_index)]));
                }
            }
            nel += 2*static_cast<size_t>(2*lH.l + 1);
            lH += 2*lH.l + 1;
        }
        for(lm lB = {1, 0, 0}; lB != lH + lH.n*lH.n; lB += 2*lB.l + 1){
        // for(lm lB = {1, 0, 0}; lB < lH; lB += 2*lB.l + 1){
            for(auto kappa : kappas){
                for(auto s : {spin{UP}, spin{DOWN}}){
                    Bs_m.back().add_function(Augmented_Bessel(lB, kappa, s, at.pos,
                        at_meshes[static_cast<size_t>(at_index)]));
                }
            }
        }
        nel = 0;
        lH = {1, 0, 0};
        std::cout << "Hs_m contains " << Hs_m.back().size() << " elements\n";
        std::cout << "Bs_m contains " << Bs_m.back().size() << " elements\n";
    }
}

void Simulation::set_up_potential(std::function<double(const size_t, const double)> at_pot, const XC_FUN func)
{
    std::cout << "Setting up potential" << std::endl;
    size_t nel = 0;
    for(auto at : cryst.atoms()){
        nel += at.get_Z();
    }
    // Set up initial potential
    pot = Potential(cryst, at_meshes, at_pot);
    pot.set_xc_fun(func);
    pot.initial_pot(nel, cryst.volume());
}

void Simulation::init_augmented_functions()
{
    std::cout << "Initialising augmented Hankel functions" << std::endl;
    for(auto wf : basis_valence){
        if(wf.l.m != -wf.l.l){
            continue;
        }
        auto H = Hs_m[wf.center].get_function(wf.l, wf.kappa, wf.s);
        H.update(pot.sphere(wf.center), -GSL::pow_int(static_cast<double>(cryst.atoms()[wf.center].get_Z())/H.l.n, 2) + pot.MT0(), wf.core_state);
        H.EH() -= pot.MT0();
    }
    std::cout << "Initialising core states" << std::endl;
    for(auto wf : basis_core){
        if(wf.l.m != -wf.l.l){
            continue;
        }
        auto H = Hs_m[wf.center].get_function(wf.l, wf.kappa, wf.s);
        H.update(pot.sphere(wf.center), -GSL::pow_int(static_cast<double>(cryst.atoms()[wf.center].get_Z())/H.l.n, 2) + pot.MT0(), wf.core_state);
        H.EH() -= pot.MT0();
    }
    std::cout << "Initialising augmented Bessel functions" << std::endl;
    for(size_t at_i = 0; at_i < cryst.atoms().size(); at_i++){
        for(auto J : Bs_m[at_i]){
            if(J.l.m != -J.l.l){
                continue;
            }
            J.update(pot.sphere(at_i), -GSL::pow_int(static_cast<double>(cryst.atoms()[at_i].get_Z())/J.l.n, 2) + pot.MT0(), false);
            J.EJ() -= pot.MT0();
        }
    }
}

Simulation::Simulation(const Crystal_t<3, Atom>& crystal, const XC_FUN func, const std::vector<double> kappas_n,
 std::function<double(const size_t, const double)> at_pot)
 : kappas(kappas_n), Hs_m(), Bs_m(), cryst(crystal), pot(), n(), at_meshes(), basis_valence(), basis_core(), k_eigenvals(), XH1(),
   XS1(), XH2(), XS2(), XH3(), XS3()
{

    set_up_crystal();
    set_up_basis();
    set_up_potential(at_pot, func);
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
    XH1 = std::vector<GSL::Matrix>(cryst.atoms().size(), GSL::Matrix(nH, nH));
    XS1 = std::vector<GSL::Matrix>(cryst.atoms().size(), GSL::Matrix(nH, nH));

    XH2 = std::vector<GSL::Matrix>(cryst.atoms().size(), GSL::Matrix(nH, nH));
    XS2 = std::vector<GSL::Matrix>(cryst.atoms().size(), GSL::Matrix(nH, nH));

    XH3 = std::vector<GSL::Matrix>(cryst.atoms().size(), GSL::Matrix(nB, nB));
    XS3 = std::vector<GSL::Matrix>(cryst.atoms().size(), GSL::Matrix(nB, nB));
}

double Simulation::X_H1(const Augmented_Hankel& Ht1, const Augmented_Hankel& Ht2, const size_t& at)
{
    double res = 0;
    const Envelope_Hankel H1(at, Ht1.l, Ht1.kappa);
    const Envelope_Hankel H2(at, Ht2.l, Ht2.kappa);
    res += Ht2.EH()*augmented_integral(Ht1, Ht2);
    res += -Ht2.kappa*Ht2.kappa*off_atomic_integral(H1, H2, at_meshes[at].r_back());
    return res;
}

double Simulation::X_H2(const Augmented_Hankel& Ht1, const Augmented_Bessel& Jt2, const size_t& at)
{
    double res = 0;
    const Envelope_Hankel H1(at, Ht1.l, Ht1.kappa);
    const Envelope_Bessel J2(at, Jt2.l, Jt2.kappa);
    //res += Jt2.EJ*augmented_integral(Ht1, Jt2);
    if(std::abs(Ht1.EH() - Jt2.EJ()) < 1e-7){
        res += Jt2.EJ()*augmented_integral(Ht1, Ht1);
    }else{
        res += Jt2.EJ()/(Ht1.EH() - Jt2.EJ());
    }
    res -= -Jt2.kappa*Jt2.kappa*atomic_integral(H1, J2, at_meshes[at].r_back());
    if(Ht1.kappa*Ht1.kappa != Jt2.kappa*Jt2.kappa){
        res += -Jt2.kappa*Jt2.kappa/(-Jt2.kappa*Jt2.kappa + Ht1.kappa*Ht1.kappa);
    }

    return res;
}

double Simulation::X_H3(const Augmented_Bessel& Jt1, const Augmented_Bessel& Jt2, const size_t& at)
{
    double res = 0;
    const Envelope_Bessel J1(at, Jt1.l, Jt1.kappa);
    const Envelope_Bessel J2(at, Jt2.l, Jt2.kappa);
    res += Jt2.EJ()*augmented_integral(Jt1, Jt2);
    res -= -Jt2.kappa*Jt2.kappa*atomic_integral(J1, J2, at_meshes[at].r_back());
    return res;
}

GSL::Complex Simulation::H_element(const size_t i1, const size_t i2, const GSL::Vector& kp) const
{
    const Augmented_spherical_wave w1 = basis_valence[i1];
    const Augmented_spherical_wave w2 = basis_valence[i2];
    size_t at_i = w1.center, at_j = w2.center;
    lm l_low1 = Bs_m[at_i].back().l, l_low2 = Bs_m[at_j].back().l;
    l_low1.m *= -1;
    l_low2.m *= -1;
    GSL::Complex res(0,0);
    const GSL::Vector tau_ij(cryst.atoms()[at_i].pos - cryst.atoms()[at_j].pos);
    size_t l_i = Hs_m[at_i].get_index(w1.l, w1.kappa, w1.s);
    size_t l_j = Hs_m[at_j].get_index(w2.l, w2.kappa, w2.s);


    if(at_i == at_j && w1.l == w2.l){
        res += XH1[at_i][l_i][l_j];
    }

    Bloch_summed_structure_constant B1(l_low2.l + 1, w2.kappa, cryst, w1.l, w2.l);
    res += XH2[at_i][l_i][l_j]*B1(tau_ij, kp);

    B1 = Bloch_summed_structure_constant(l_low1.l + 1, w1.kappa, cryst, w1.l, w2.l);
    if(std::abs(w1.kappa*w1.kappa - w2.kappa*w2.kappa) < 1E-15){
        B1 = Bloch_summed_structure_constant(l_low1.l + 1, w1.kappa, cryst, w1.l, w2.l);
        // res += -w2.kappa*w2.kappa*B1.dot(tau_ij, kp);
        res += B1(tau_ij, kp);
    }

    res += B1(tau_ij, kp)*XH2[at_j][l_j][l_i];

    for(size_t at_m = 0; at_m < cryst.atoms().size(); at_m++){
        lm l_low3 = Bs_m[at_m].back().l;
        l_low3.m *= -1;
        for(Augmented_Bessel J1 : Bs_m[at_m]){
            for(Augmented_Bessel J2 : Bs_m[at_m]){
                if(J1.l == J2.l){
                    for(int dl = 0; dl < 2*J1.l.l + 1; dl++){
                        size_t lpp1 = Bs_m[at_m].get_index(J1.l + dl, J1.kappa, J1.s);
                        size_t lpp2 = Bs_m[at_m].get_index(J2.l + dl, J2.kappa, J2.s);
                        B1 = Bloch_summed_structure_constant(l_low3.l + 1, w1.kappa, cryst, J1.l + dl, w1.l);
                        Bloch_summed_structure_constant B2 = Bloch_summed_structure_constant(l_low3.l + 1, w2.kappa, cryst, J2.l + dl, w2.l);
                        res += B1(cryst.atoms()[at_m].pos - cryst.atoms()[at_i].pos, kp).conjugate()*
                        XH3[at_m][lpp1][lpp2]*
                        B2(J2.center - cryst.atoms()[at_j].pos, kp);
                    }
                }
            }
        }
    }
    return res;
}

double Simulation::X_S1(const Augmented_Hankel& Ht1, const Augmented_Hankel& Ht2, const size_t& at)
{
    double res = 0;
    const Envelope_Hankel H1(at, Ht1.l, Ht1.kappa);
    const Envelope_Hankel H2(at, Ht2.l, Ht2.kappa);
    res += augmented_integral(Ht1, Ht2);
    res += off_atomic_integral(H1, H2, at_meshes[at].r_back());
    return res;
}

double Simulation::X_S2(const Augmented_Hankel& Ht1, const Augmented_Bessel& Jt2, const size_t& at)
{
    double res = 0;
    const Envelope_Hankel H1(at, Ht1.l, Ht1.kappa);
    const Envelope_Bessel J2(at, Jt2.l, Jt2.kappa);
//    res += augmented_integral(Ht1, Jt2);
    if(std::abs(Ht1.EH() - Jt2.EJ()) < 1e-7){
        res += augmented_integral(Ht1, Ht1);
    }else{
        res += 1./(Ht1.EH() - Jt2.EJ());
    }
    res -= atomic_integral(H1, J2, at_meshes[at].r_back());
    if(Ht1.kappa*Ht1.kappa != Jt2.kappa*Jt2.kappa){
        res += 1./(-Jt2.kappa*Jt2.kappa + Ht1.kappa*Ht1.kappa);
    }
    return res;
}

double Simulation::X_S3(const Augmented_Bessel& Jt1, const Augmented_Bessel& Jt2, const size_t& at)
{
    double res = 0;
    const Envelope_Bessel J1(at, Jt1.l, Jt1.kappa);
    const Envelope_Bessel J2(at, Jt2.l, Jt2.kappa);
    res += augmented_integral(Jt1, Jt2);
    res -= atomic_integral(J1, J2, at_meshes[at].r_back());
    return res;
}

GSL::Complex Simulation::S_element(const size_t i1, const size_t i2, const GSL::Vector& kp) const
{
    const Augmented_spherical_wave w1 = basis_valence[i1];
    const Augmented_spherical_wave w2 = basis_valence[i2];
    size_t at_i = w1.center, at_j = w2.center;
    lm l_low1 = Bs_m[at_i].back().l, l_low2 = Bs_m[at_j].back().l;
    l_low1.m *= -1;
    l_low2.m *= -1;

    GSL::Complex res = 0.;
    const GSL::Vector tau_ij(cryst.atoms()[at_i].pos - cryst.atoms()[at_j].pos);
    size_t l_i = Hs_m[at_i].get_index(w1.l, w1.kappa, w1.s);
    size_t l_j = Hs_m[at_j].get_index(w2.l, w2.kappa, w2.s);

    if(at_i == at_j && l_i == l_j){
        res += XS1[at_i][l_i][l_j];
    }

    Bloch_summed_structure_constant B1 = Bloch_summed_structure_constant(l_low1.l + 1, w2.kappa, cryst, w1.l, w2.l);
    res += XS2[at_i][l_i][l_j]*B1(tau_ij, kp);

    if(std::abs(w1.kappa*w1.kappa - w2.kappa*w2.kappa) < 1E-15){
        B1 = Bloch_summed_structure_constant(l_low1.l + 1, w1.kappa, cryst, w1.l, w2.l);
        res += B1.dot(tau_ij, kp);
    }


    B1 = Bloch_summed_structure_constant(l_low2.l + 1, w1.kappa, cryst, w1.l, w2.l);
    res += B1(tau_ij, kp)*XS2[at_j][l_j][l_i];

    for(size_t at_m = 0; at_m < cryst.atoms().size(); at_m++){
        lm l_low3 = Bs_m[at_m].back().l;
        l_low3.m *= -1;
        for(Augmented_Bessel J1 : Bs_m[at_m]){
            for(Augmented_Bessel J2 : Bs_m[at_m]){
                if(J1.l == J2.l){
                    for(int dl = 0; dl < 2*J1.l.l + 1; dl++){
                        size_t lpp1 = Bs_m[at_m].get_index(J1.l + dl, J1.kappa, J1.s);
                        size_t lpp2 = Bs_m[at_m].get_index(J2.l + dl, J2.kappa, J2.s);
                        B1 = Bloch_summed_structure_constant(l_low3.l + 1, w1.kappa, cryst, J1.l + dl, w1.l);
                        Bloch_summed_structure_constant B2 = Bloch_summed_structure_constant(l_low3.l + 1, w2.kappa, cryst, J2.l + dl, w2.l);
                        res += B1(cryst.atoms()[at_m].pos - cryst.atoms()[at_i].pos, kp).conjugate()*
                        XS3[at_m][lpp1][lpp2]*
                        B2(J2.center - cryst.atoms()[at_j].pos, kp);
                    }
                }
            }
        }
    }
    return res;
}


void Simulation::set_up_X_matrices()
{
    for(size_t i = 0; i < basis_valence.size(); i++){
        for(size_t j = 0; j < basis_valence.size(); j++){

            Augmented_spherical_wave w1 = basis_valence[i];
            Augmented_spherical_wave w2 = basis_valence[j];
            if(w1.s != w2.s){
                continue;
            }
            size_t at_i = w1.center, at_j = w2.center;
            size_t i1 = Hs_m[at_i].get_index(w1.l, w1.kappa, w1.s);
            size_t i2 = Hs_m[at_j].get_index(w2.l, w2.kappa, w2.s);
            size_t ipp1 = 0;
            size_t ipp2 = 0;

            if(w1.center == w2.center && w1.l == w2.l && w1.s == w2.s){
                XH1[at_i][i1][i2] = X_H1(Hs_m[at_i].get_function(w1.l, w1.kappa, w1.s), Hs_m[at_j].get_function(w2.l, w2.kappa, w2.s), at_i);
                XS1[at_i][i1][i2] = X_S1(Hs_m[at_i].get_function(w1.l, w1.kappa, w1.s), Hs_m[at_j].get_function(w2.l, w2.kappa, w2.s), at_i);
            }

            // Find off center expansion on site w1.center of wave w2
            for(Augmented_Bessel J2s : Bs_m[at_i]){
                if(J2s.l == w1.l && J2s.s == w1.s){
                    ipp2 = Bs_m[at_i].get_index(J2s.l, J2s.kappa, J2s.s);
                    XH2[at_i][i1][ipp2] = X_H2(Hs_m[at_i].get_function(w1.l, w1.kappa, w1.s), J2s, at_i);
                    XS2[at_i][i1][ipp2] = X_S2(Hs_m[at_i].get_function(w1.l, w1.kappa, w1.s), J2s, at_i);
                }
            }

            // Find off center expansion on site s of waves w1 and w2
            for(size_t n_s = 0; n_s < cryst.atoms().size(); n_s++){
                for(Augmented_Bessel J1s : Bs_m[n_s]){
                    for(Augmented_Bessel J2s : Bs_m[n_s]){
                        if(J1s.l == J2s.l && J1s.s == J2s.s){
                            ipp1 = Bs_m[n_s].get_index(J1s.l, J1s.kappa, J1s.s);
                            ipp2 = Bs_m[n_s].get_index(J2s.l, J2s.kappa, J2s.s);
                            XH3[n_s][ipp1][ipp2] =  X_H3(J1s, J2s, n_s);
                            XS3[n_s][ipp1][ipp2] =  X_S3(J1s, J2s, n_s);
                        }
                    }
                }
            }
        }
    }
    std::cout << "XS1 = ";
    for(auto XS1_i : XS1){
        std::cout << XS1_i << std::endl;
    }
    std::cout << "XS2 = ";
    for(auto XS2_i : XS2){
        std::cout << XS2_i << std::endl;
    }
    std::cout << "XS3 = ";
    for(auto XS3_i : XS3){
        std::cout << XS3_i << std::endl;
    }
}

const GSL::Matrix_cx Simulation::set_H(const GSL::Vector& kp) const
{
    size_t N = basis_valence.size();
    GSL::Matrix_cx H(N, N);
    for(size_t i = 0; i < N; i++){
        for(size_t j = 0; j < N; j++){
            // auto tmp = H_element(i, j, kp);
            // std::cout << "H[" << i << "][" << j << "] = " << tmp << "\n";
            H[i][j] = H_element(i, j, kp);
        }
    }
    return H;
}

const GSL::Matrix_cx Simulation::set_S(const GSL::Vector& kp) const
{
    size_t N = basis_valence.size();
    GSL::Matrix_cx S(N, N);
    for(size_t i = 0; i < basis_valence.size(); i++){
        for(size_t j = 0; j < basis_valence.size(); j++){
//            auto tmp = S_element(i, j, kp);
//            std::cout << "S[" << i << "][" << j << "] = " << tmp << "\n";
            S[i][j] = S_element(i, j, kp);
        }
    }
    return S;
}

void Simulation::add_eigvals(const GSL::Vector& kp, const GSL::Vector& eigvals)
{
    std::pair<GSL::Vector, GSL::Vector> tmp(kp, eigvals);
    k_eigenvals.insert(tmp);
}

// const std::pair<GSL::Matrix_cx, GSL::Vector> Simulation::calc_eigen(const GSL::Vector& kp) const
void Simulation::calc_eigen(const GSL::Vector& kp) const
{
    const GSL::Matrix_cx H {set_H(kp)};
    const GSL::Matrix_cx S {set_S(kp)};
//    size_t N = basis_valence.size();
    size_t N = basis_valence.size()/2;
    GSL::Matrix_cx eigvecs_up(N, N);
    // GSL::Matrix_cx eigvecs_down(N, N);
    GSL::Vector eigvals_up(N);
    // GSL::Vector eigvals_down(N);

    GSL::Matrix_cx H_up(N, N);
    // GSL::Matrix_cx H_down(N, N);
    GSL::Matrix_cx S_up(N, N);
    // GSL::Matrix_cx S_down(N, N);

    for(size_t i = 0; i < N; i++){
        for(size_t j = 0; j < N; j++){
            if(!(basis_valence[2*i].s == UP && basis_valence[2*j].s == UP)){
                std::cout << "FATAL ERROR!\n";
            }
        //    H_up[i][j] = H[i][j];
        //    S_up[i][j] = S[i][j];
            H_up[i][j] = H[2*i][2*j];
            S_up[i][j] = S[2*i][2*j];
//            H_down[i][j] = H[2*i+1][2*j+1];
//            S_down[i][j] = S[2*i+1][2*j+1];
        }
    }

    std::pair<GSL::Matrix_cx, GSL::Vector> tmp, overlap;
    overlap = GSL::hermitian_eigen(S_up);
    try{
        tmp = GSL::general_hermitian_eigen(H_up, S_up);
        eigvecs_up = tmp.first;
        eigvals_up = tmp.second;
    }catch (const std::runtime_error &e){
	    std::cerr << e.what();
	    std::cout << " Hamiltonian matrix\n";
	    for(size_t i = 0 ; i < N; i++){
		    std::cout << "  " << H_up[i] << "\n";
	    }
	    std::cout << " Overlap matrix\n";
	    for(size_t i = 0 ; i < N; i++){
		    std::cout << "  " << S_up[i] << "\n";
	    }
        std::cout << "Eigenvalues of Overlap matrix\n" << overlap.second << "\n";
    }
    // for(size_t i = 0; i < N; i++){
    std::cout << "  Eigenvalues (up): " << eigvals_up << " (Ry)\n";
        // std::cout << "  Eigenvectors (up) :\n";
        // std::cout << "\t" << eigvecs_up.get_col(i) << "\n";

    // }
    std::cout << "\n\n";
}

void Simulation::print_eigvals(K_mesh& k_mesh)
{
    std::ofstream of;
    of.open("bands.dat", std::ios::trunc);
    for(const auto& kp :  k_mesh.k_points){
        const GSL::Vector eigvals {k_eigenvals[kp]};
        of << kp[0] << " " << kp[1] << " " << kp[2] << " ";
        for(size_t i = 0; i < eigvals.dim(); i++){
            of << eigvals[i] << " ";
        }
        of << "\n";
    }
    of.close();
}
