#include <gsl-asw/brillouin_zone_integration.h>
#include <gsl-asw/numerical-mesh.h>
#include <fstream>
#include <GSLpp/error.h>

const GSL::Matrix eigvals{{-1, 0, 1},{ -1.5, 0, 1.5},{-0.5, 0, 0.5},{-1, 0, 1}};
// const GSL::Matrix eigvals{{0},{0},{0},{0}};

const double kB = 6.33362083E-4; // Ry/K


int main()
{
    GSL::Error_handler e_handler;
    e_handler.off();
    double T = 300;
    Simple_sampler sbz(eigvals, 5e-5, 1e-2);
    Fermi_sampler fbz(eigvals, kB*T);
    Gaussian_sampler gbz(eigvals, kB*T);
    Lorentzian_sampler lbz(eigvals, kB*T);
    size_t N = 2;
    MethfesselPaxton_sampler mpbz(eigvals, kB*T, N);
    ColdSmearing_sampler csbz(eigvals, kB*T);

    Linear_mesh<1, double> mesh(-5, 5, 1001);

    std::ofstream out;
    out.open("brillouin_zone_integration_test.dat");
    out << "# Brillouin zone integration test plots\n";
    out << "# Peak at 0 should have double the intensity of the peaks at -1 and +1, who should have double intensity compared to peaks at -1.5, -0.5, 0.5 and, 1.5\n";
    out << "Energy\t";
    out << "\"Simple\"\t";
    out << "\"FermiSampler\"\t";
    out << "\"GaussianSsampler\"\t";
    out << "\"LorentzianSampler\"\t";
    out << "\"MethfesselPaxtonSampler_" << N << "\"\t";
    out << "\"ColdSmearingSampler\"\n";

    for(auto point : mesh){
        auto e = point.r();
        out << e << "\t" << sbz.dos(e) << "\t" << fbz.dos(e) << "\t" << gbz.dos(e)
        << "\t" << lbz.dos(e) << "\t" << mpbz.dos(e) << "\t" << csbz.dos(e) << "\n";
    }
    out.close();
    return 0 ;
}
