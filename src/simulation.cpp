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

Simulation::Simulation(const Crystal_t<3, Atom>& crystal, const XC_FUN func, const double kappa,
 std::function<double(const size_t, const double)> at_pot)
 : cryst(crystal), pot(), n(), basis_valence(), basis_core(), k_eigenvals(), XH1(),
   XS1(), XH2(), XS2(), XH3(), XS3()
{
    // Calculate MT radii
    double at_vol = 0;
    size_t l_low = 2;

    std::vector<std::vector<Site_t<3>>> nn = cryst.calc_nearest_neighbours();
    // for(size_t i = 0; i < cryst.sites().size(); i++){
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
    // Set logarithmic meshes for the atoms and count number of electrons
    size_t nel = 0;
    for(auto site : cryst.sites()){
        std::cout << cryst.atom(site).get_Z() << std::flush;
        cryst.atom(site).mesh = Logarithmic_mesh(cryst.atom(site).get_AS(),
            500 + static_cast<size_t>(1000*std::ceil(
                static_cast<double>(cryst.atom(site).get_Z())/118)) );
        nel += cryst.atom(site).get_Z();
    }
    // Set up initial potential
    pot = Potential(cryst, at_pot);
    pot.set_xc_fun(func);
    pot.initial_pot(nel, cryst.volume());
    // Divide electrons into core and valence states
    for(const Atom& at: cryst.atoms()){
        add_states(at, kappa);
        if(at.get_Z() > 20){
            l_low = 3;
        }
    }
    std::cout << "# valence states " << basis_valence.size() << std::endl;

    if(basis_core.size() > nel){
        throw std::runtime_error("Number of core states exceeds number of electrons in the system!");
    }
    std::cout << "Setting up core states." << std::endl;
    for(Augmented_spherical_wave& wave : basis_core){
        wave.core_state = true;
        wave.set_up(pot);
    }
    std::cout << "Setting up valence states." << std::endl;
    for(Augmented_spherical_wave& wave: basis_valence){
        wave.core_state = false;
        wave.set_up(pot);
    }

    // Reserve space for each atom ad each 0 <= l <= l_low, -l <= m <= l
    // XH1 = GSL::Matrix(cryst.atoms().size(),(l_low + 1)*(l_low + 1));
    XH1 = GSL::Matrix(cryst.atoms().size(),(l_low + 1)*(2*l_low + 1));
    // XS1 = GSL::Matrix(cryst.atoms().size(),(l_low + 1)*(l_low + 1));
    XS1 = GSL::Matrix(cryst.atoms().size(),(l_low + 1)*(2*l_low + 1));
    XH2 = GSL::Matrix(cryst.atoms().size(),(l_low + 1)*(l_low + 1));
    XS2 = GSL::Matrix(cryst.atoms().size(),(l_low + 1)*(l_low + 1));

    // Reserve space for each atom and each 0<= l <= l_int = l_low + 1, -l <= m <= l
    XH3 = GSL::Matrix(cryst.atoms().size(),(l_low + 2)*(l_low + 2));
    XS3 = GSL::Matrix(cryst.atoms().size(),(l_low + 2)*(l_low + 2));
}

struct nl_comp {
    bool operator() (const std::pair<int, int>& lhs, const std::pair<int, int>& rhs)
    {
        if(lhs.first + lhs.second == rhs.first + rhs.second){
            return lhs.first < rhs.first;
        }else{
            return (lhs.first + lhs.second) < (rhs.first + rhs.second);
        }
    }
};

void Simulation::add_states(const Atom& at, const double kappa)
{
    std::set<std::pair<int, int>, nl_comp> p({{1, 0}});
    int n_s = 1, ln = 0;
    size_t nel = 0;
    size_t Z = at.get_Z();
    while(nel < Z){
        for(int l = 0; l < n_s; l++){
            p.insert({n_s+1, l});
            nel += static_cast<size_t>(2*(2*l + 1));
        }
        p.insert({n_s+1, n_s});
        nel += static_cast<size_t>(2*(2*n_s + 1));
        n_s++;
    }
    nel = 0;
    for(std::pair<int, int> tmp : p){
        n_s = tmp.first;
        ln = tmp.second;
        if(nel + static_cast<size_t>(2*n_s*n_s) >= Z){
            std::cout << "valence : ( " << n_s << " " << ln << " )" << std::endl;
            for(int m = -ln; m <= ln; m++){
                basis_valence.push_back(Augmented_spherical_wave(kappa, n_s, lm {ln, m}, UP, at, cryst.atoms()));
//                basis_valence.push_back(Augmented_spherical_wave(kappa, n_s, lm {ln, m}, DOWN, at, cryst.atoms));
            }
        }else{
            std::cout << "core : ( " << n_s << " " << ln << " )" << std::endl;
            for(int m = -ln; m <= ln; m++){
                basis_core.push_back(Augmented_spherical_wave(kappa, n_s, lm {ln, m}, UP, at, cryst.atoms()));
//                basis_core.push_back(Augmented_spherical_wave(kappa, n_s, lm {ln, m}, DOWN, at, cryst.atoms));
            }
        }
        nel += static_cast<size_t>(2*(2*ln + 1));
    }
}

double X_H1(const Augmented_Hankel& Ht1, const Augmented_Hankel& Ht2, const Atom& at)
{
    double res = 0;
    const Envelope_Hankel H1(at, Ht1.l, Ht1.kappa);
    const Envelope_Hankel H2(at, Ht2.l, Ht2.kappa);
    res += Ht2.EH()*augmented_integral(Ht1, Ht2);
    res += -Ht2.kappa*Ht2.kappa*off_atomic_integral(H1, H2);
    return res;
}

double X_H2(const Augmented_Hankel& Ht1, const Augmented_Bessel& Jt2, const Atom& at)
{
    double res = 0;
    const Envelope_Hankel H1(at, Ht1.l, Ht1.kappa);
    const Envelope_Bessel J2(at, Jt2.l, Jt2.kappa);
    //res += Jt2.EJ*augmented_integral(Ht1, Jt2);
    res += Jt2.EJ()/(Ht1.EH() - Jt2.EJ());
    res -= -Jt2.kappa*Jt2.kappa*atomic_integral(H1, J2);
    if(Ht1.kappa*Ht1.kappa != Jt2.kappa*Jt2.kappa){
        res += -Jt2.kappa*Jt2.kappa/(-Jt2.kappa*Jt2.kappa + Ht1.kappa*Ht1.kappa);
    }

    return res;
}

double X_H3(const Augmented_Bessel& Jt1, const Augmented_Bessel& Jt2, const Atom& at)
{
    double res = 0;
    const Envelope_Bessel J1(at, Jt1.l, Jt1.kappa);
    const Envelope_Bessel J2(at, Jt2.l, Jt2.kappa);
    res += Jt2.EJ()*augmented_integral(Jt1, Jt2);
    res -= -Jt2.kappa*Jt2.kappa*atomic_integral(J1, J2);
    return res;
    // return 0;
}

GSL::Complex Simulation::H_element(const size_t i1, const size_t i2, const GSL::Vector& kp) const
{
    const Augmented_spherical_wave w1 = basis_valence[i1];
    const Augmented_spherical_wave w2 = basis_valence[i2];
    int l_low_1 = 2, l_low_2 = 2, l_low  = 2;
    if(w1.center.Z > 20){
        l_low_1 = 3;
    }
    if(w2.center.Z > 20){
        l_low_2 = 3;
    }
    GSL::Complex res(0,0);
    const GSL::Vector tau_ij(w1.center.pos - w2.center.pos);
    size_t at_i = 0, at_j = 0;
    for(size_t i = 0; i < cryst.atoms().size(); i++){
        if(cryst.atoms()[i] == w1.center){
            at_i = i;
        }
        if(cryst.atoms()[i] == w2.center){
            at_j = i;
        }
    }
    size_t l_i = static_cast<size_t>(w1.l.l*(l_low_1 + 1) + w1.l.m);
    size_t l_j = static_cast<size_t>(w2.l.l*(l_low_2 + 1) + w2.l.m);

    if(at_i == at_j && l_i == l_j){
        std::cout << "Hola!\t";
        res += XH1[at_i][l_i];
    }

    Bloch_summed_structure_constant B1 {l_low_2 + 1, w2.kappa, cryst, w1.l, w2.l};
    res += XH2[at_i][l_i]*B1(tau_ij, kp);

    if(std::abs(w1.kappa*w1.kappa - w2.kappa*w2.kappa) < 1E-15){
        B1 = Bloch_summed_structure_constant(l_low_1 + 1, w1.kappa, cryst, w1.l, w2.l);
        res += -w2.kappa*w2.kappa*B1.dot(tau_ij, kp);
    }

    B1 = Bloch_summed_structure_constant(l_low_1 + 1, w1.kappa, cryst, w1.l, w2.l);
    res += B1(tau_ij, kp)*XH2[at_j][l_j];
    if(std::abs(w1.kappa*w1.kappa - w2.kappa*w2.kappa) < 1E15){
        res += B1(tau_ij, kp);
    }

    for(size_t at_m = 0; at_m < cryst.atoms().size(); at_m++){
        if(cryst.atoms()[at_m].get_Z() > 20){
            l_low = 3;
        }
        for(Augmented_Bessel J1 : w1.J[at_m]){
            // for(Augmented_Bessel J2 : w2.J[at_m]){
                // if(J1.l == J2.l){
                    size_t ipp = static_cast<size_t>(J1.l.l*(J1.l.l + 1) + J1.l.m);
                    B1 = Bloch_summed_structure_constant(l_low + 1, w1.kappa, cryst, J1.l, w1.l);
                    Bloch_summed_structure_constant B2 = Bloch_summed_structure_constant(l_low + 1, w2.kappa, cryst, J1.l, w2.l);
                    res += B1(J1.center - w1.center.pos, kp).conjugate()*
                    XH3[at_m][ipp]*
                    B2(J1.center - w2.center.pos, kp);
                // }
            // }
        }
    }
    return res;
}

double X_S1(const Augmented_Hankel& Ht1, const Augmented_Hankel& Ht2, const Atom& at)
{
    double res = 0;
    const Envelope_Hankel H1(at, Ht1.l, Ht1.kappa);
    const Envelope_Hankel H2(at, Ht2.l, Ht2.kappa);
    res += augmented_integral(Ht1, Ht2);
    res += off_atomic_integral(H1, H2);
    return res;
}

double X_S2(const Augmented_Hankel& Ht1, const Augmented_Bessel& Jt2, const Atom& at)
{
    double res = 0;
    const Envelope_Hankel H1(at, Ht1.l, Ht1.kappa);
    const Envelope_Bessel J2(at, Jt2.l, Jt2.kappa);
//    res += augmented_integral(Ht1, Jt2);
    res += 1./(Ht1.EH() - Jt2.EJ());
    res -= atomic_integral(H1, J2);
    if(Ht1.kappa*Ht1.kappa != Jt2.kappa*Jt2.kappa){
        res += 1./(-Jt2.kappa*Jt2.kappa + Ht1.kappa*Ht1.kappa);
    }
    return res;
}

double X_S3(const Augmented_Bessel& Jt1, const Augmented_Bessel& Jt2, const Atom& at)
{
    double res = 0;
    const Envelope_Bessel J1(at, Jt1.l, Jt1.kappa);
    const Envelope_Bessel J2(at, Jt2.l, Jt2.kappa);
    res += augmented_integral(Jt1, Jt2);
    res -= atomic_integral(J1, J2);
    return res;
}

GSL::Complex Simulation::S_element(const size_t i1, const size_t i2, const GSL::Vector& kp) const
{
    const Augmented_spherical_wave w1 = basis_valence[i1];
    const Augmented_spherical_wave w2 = basis_valence[i2];
    int l_low_1 = 2, l_low_2 = 2, l_low  = 2;
    if(w1.center.Z > 20){
        l_low_1 = 3;
    }
    if(w2.center.Z > 20){
        l_low_2 = 3;
    }
    GSL::Complex res = 0.;
    const GSL::Vector tau_ij(w1.center.pos - w2.center.pos);
    size_t at_i = 0, at_j = 0;
    for(size_t i = 0; i < cryst.atoms().size(); i++){
        if(cryst.atoms()[i] == w1.center){
            at_i = i;
        }
        if(cryst.atoms()[i] == w2.center){
            at_j = i;
        }
    }
    size_t l_i = static_cast<size_t>(w1.l.l*(l_low_1 + 1) + w1.l.m);
    size_t l_j = static_cast<size_t>(w2.l.l*(l_low_2 + 1) + w2.l.m);

    if(at_i == at_j && l_i == l_j){
        res += XS1[at_i][l_i];
    }

    Bloch_summed_structure_constant B1 = Bloch_summed_structure_constant(l_low_1 + 1, w2.kappa, cryst, w1.l, w2.l);
    res += XS2[at_i][l_i]*B1(tau_ij, kp);

    if(std::abs(w1.kappa*w1.kappa - w2.kappa*w2.kappa) < 1E-15){
        B1 = Bloch_summed_structure_constant(l_low_1 + 1, w1.kappa, cryst, w1.l, w2.l);
        res += B1.dot(tau_ij, kp);
    }


    B1 = Bloch_summed_structure_constant(l_low_2 + 1, w1.kappa, cryst, w1.l, w2.l);
    res += B1(tau_ij, kp)*XS2[at_j][l_j];

    for(size_t at_m = 0; at_m < cryst.atoms().size(); at_m++){
        if(cryst.atoms()[at_m].get_Z() > 20){
            l_low = 3;
        }
        for(Augmented_Bessel J1 : w1.J[at_m]){
            for(Augmented_Bessel J2 : w2.J[at_m]){
                if(J1.l == J2.l){
                    size_t ipp = static_cast<size_t>(J1.l.l*(l_low + 1) + J1.l.m);
                    B1 = Bloch_summed_structure_constant(l_low + 1, w1.kappa, cryst, J1.l, w1.l);
                    Bloch_summed_structure_constant B2 = Bloch_summed_structure_constant(l_low + 1, w2.kappa, cryst, J2.l, w2.l);
                    res += B1(J1.center - w1.center.pos, kp).conjugate()*
	                   XS3[at_m][ipp]*
        	           B2(J2.center - w2.center.pos, kp);
                }
            }
        }
    }
    return res;
}


void Simulation::set_up_X_matrices()
{
    size_t at_i = 0;

    std::cout << "Setting up X matrices (for H and S)" << std::endl;
    for(size_t i = 0; i < basis_valence.size(); i++){
        for(size_t j = 0; j < basis_valence.size(); j++){

            Augmented_spherical_wave w1 = basis_valence[i];
            Augmented_spherical_wave w2 = basis_valence[j];
            int l_low_1 = 2, l_low_2 = 2, l_low  = 2;
            if(w1.center.Z > 20){
                l_low_1 = 3;
            }
            if(w2.center.Z > 20){
                l_low_2 = 3;
            }
            size_t l_i = static_cast<size_t>(w1.l.l*(l_low_1 + 1) + w1.l.m);

            for(size_t n_s = 0; n_s < cryst.atoms().size(); n_s++){
                if(cryst.atoms()[n_s] == w1.center){
                    at_i = n_s;
                }
            }

            if(w1.center == w2.center && w1.H.l == w2.H.l){
                XH1[at_i][l_i] = X_H1(w1.H, w2.H, w1.center);
                XS1[at_i][l_i] = X_S1(w1.H, w2.H, w1.center);
            }

            // Find off center expansion on site w1.center of wave w2
            for(Augmented_Bessel J2s : w2.J[at_i]){
                if(J2s.l == w1.H.l){
                    XH2[at_i][l_i] = X_H2(w1.H, J2s, w1.center);
                    XS2[at_i][l_i] = X_S2(w1.H, J2s, w1.center);
                }
            }

            // Find off center expansion on site s of waves w1 ans w2
            for(size_t n_s = 0; n_s < cryst.atoms().size(); n_s++){
                if(cryst.atoms()[n_s].get_Z() > 20){
                    l_low = 3;
                }
                for(Augmented_Bessel J1s : w1.J[n_s]){
                    for(Augmented_Bessel J2s : w2.J[n_s]){
                        if(J1s.l == J2s.l){
                            size_t ipp = static_cast<size_t>(
                                J1s.l.l*(l_low + 1) + J1s.l.m);
                            XH3[n_s][ipp]
                              =  X_H3(J1s, J2s, cryst.atoms()[n_s]);
                            XS3[n_s][ipp]
                              =  X_S3(J1s, J2s, cryst.atoms()[n_s]);
                        }
                    }
                }
            }

        }
    }
    std::cout << "XS1 = " << XS1 << std::endl;
    std::cout << "XS2 = " << XS2 << std::endl;
    std::cout << "XS3 = " << XS3 << std::endl;
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
            auto tmp = S_element(i, j, kp);
            std::cout << "S[" << i << "][" << j << "] = " << tmp << "\n";
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
    size_t N = basis_valence.size();
    // size_t N = basis_valence.size()/2;
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
            H_up[i][j] = H[i][j];
            S_up[i][j] = S[i][j];
            // H_up[i][j] = H[2*i][2*j];
            // S_up[i][j] = S[2*i][2*j];
//            H_down[i][j] = H[2*i+1][2*j+1];
//            S_down[i][j] = S[2*i+1][2*j+1];
        }
    }

    std::pair<GSL::Matrix_cx, GSL::Vector> tmp;
    tmp = GSL::hermitian_eigen(S_up);
    try{
        tmp = GSL::general_hermitian_eigen(H_up, S_up);
        eigvecs_up = tmp.first;
        eigvals_up = tmp.second;
    }catch (const std::runtime_error &e){
	    std::cerr << e.what();
	    std::cout << " Overlap matrix\n";
	    for(size_t i = 0 ; i < N; i++){
		    std::cout << "  " << S_up[i] << "\n";
	    }
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
