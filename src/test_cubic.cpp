#include <spherical_fun.h>
#include <GSLpp/vector.h>
#include <cmath>
#include <iostream>
#include <fstream>

template<class T>
int signum(T t)
{
    return (T(0) < t) - (t < T(0));
}

int main()
{
    std::fstream out_file("cubic_harmonic_test.dat", std::ios::out);
    double x_max = 1., y_max = 1., z_max = 1.;
    double dx = 0.08, dy = 0.08, dz = 0.08;
    for(int l = 0; l < 5; l++){
	    for(int m = -l; m <= l; m++){
		    for(double x = -x_max; x <= x_max; x += dx ){
		    for(double y = -y_max; y <= y_max; y += dy ){
		    for(double z = -z_max; z <= z_max; z += dz ){
			    GSL::Vector r{x, y, z},tmp;
			    tmp.copy(r);
			    r.normalize<double>();
			    GSL::Vector res = r*std::abs(cubic_harmonic(lm{l, m}, tmp).val);
			    out_file << res[0] << " " << res[1] << " " << res[2] << " " <<
                signum(cubic_harmonic(lm{l, m}, tmp).val) << "\n";
		    }}}
		    out_file << "\n\n";
	    }
    }
    out_file.close();
    return 0;
}
