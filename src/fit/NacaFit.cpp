#include "fit/NacaFit.hpp"

#include <cmath>
#include <stdexcept>

namespace BladeGen::Fit {

Profile::ProfileParams naca65_start(const Naca::Naca65Params& p)
{
    p.validate();
    auto q = Profile::ProfileParams::naca65_1210();
    const double ts = p.t_max_over_c / 0.10;
    const auto ref = Naca::equivalent_arc_angles(1.2), cur = Naca::equivalent_arc_angles(p.cl0);
    const double ws = ref[1] > 0.0 ? cur[1] / ref[1] : 1.0;
    // camber: the SS control point height is mean-line camber + half thickness
    const double camber_ref = Naca::a10(0.5, 1.2)[0], camber = Naca::a10(0.5, p.cl0)[0];
    const double half_ref = q.ss.ctrl2_theta - camber_ref;
    q.ss.ctrl2_theta = camber + half_ref * ts;
    q.le.beta_deg *= ws;
    q.te.beta_deg *= ws;
    q.le.radius = std::max(0.002, q.le.radius * ts * ts);
    q.thickness.t_max *= ts;
    q.thickness.t_le *= ts;
    q.thickness.t_te *= ts;
    if (p.variant == Naca::Variant::KovachSandercock) {
        q.te.radius = 0.008;  // the Kovach TE circle
        q.te.stretch_ab = 1.0;
        q.te.d_beta_ss_deg = q.te.d_beta_ps_deg = 6.0;
        q.thickness.t_te = 0.016;
    }
    return q;
}

FitOptions naca65_fit_options()
{
    FitOptions o;
    o.fairness_weight = 1.0;
    o.stages = 2;
    o.max_iterations = 40;
    o.n_edge = 50;
    o.n_body = 140;
    return o;
}

FitResult fit_naca65(const Naca::Naca65Section& s, const FitOptions& opt_in, int n_target)
{
    if (n_target < 20) throw std::invalid_argument("fit_naca65: n_target >= 20");
    const int ne = n_target / 4, nb = n_target - ne;
    const auto ss = s.sample_side(Profile::Side::Suction, ne, nb).xy;
    const auto ps = s.sample_side(Profile::Side::Pressure, ne, nb).xy;
    FitOptions opt = opt_in;
    const auto start = naca65_start(s.params());
    if (opt.params.empty()) opt.params = default_fit_params(start);
    return fit_profile(ss, ps, start, opt);
}

}  // namespace BladeGen::Fit
