#ifndef EWALD_INT_H
#define EWALD_INT_H

#include <vector>
#include "log_mesh.h"
#include "spherical_fun.h"

class Ewald_integral{
        double ewald_param;
        double kappa;
        double bar_ew_int(lm l, double r);
        double bar_comp_ew_int(lm l, double r);
    public:

        void set_ewald_param(double eta);
        void set_kappa(double kappa);

        Ewald_integral();
        std::vector<double> evaluate(lm l, Logarithmic_mesh &mesh);
        std::vector<double> evaluate_comp(lm l, Logarithmic_mesh &mesh);
        double ewald_int(lm l, double r);
        double comp_ewald_int(lm l, double r);
};

#endif //EWALD_INT_H
