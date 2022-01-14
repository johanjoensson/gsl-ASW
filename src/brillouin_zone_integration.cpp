/*
#include <gsl-asw/brillouin_zone_integration.h>

double simple_delta(const double x)
{
    return std::abs(x) < 1e-12 ? 1 : 0;
}

std::function<double(const double, const double)>
    Brillouin_sampling::delta_approx(const Sampling_method method) const
{
    switch (method)
    {
        case Simple: return simple_delta;
    }
}
*/
