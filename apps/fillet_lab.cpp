// fillet_lab: development probe of the rolling-ball fillet resolution.
#include "occt/Fillet.hpp"
#include "occt/OcctBlade.hpp"
#include <chrono>
#include <cstdio>
using namespace BladeGen;
int main()
{
    const Blade::Blade blade(Blade::BladeDefinition::demo_rotor());
    const auto sh = Occt::build_blade(blade, 0.0, 24, 80);
    for (auto [ns, rows] : {std::pair{720, 9}, std::pair{720, 15}, std::pair{1440, 9}, std::pair{1440, 15}}) {
        Occt::FilletOptions fo;
        fo.radius = 0.003;
        fo.stations = ns;
        fo.arc_rows = rows;
        const auto t0 = std::chrono::steady_clock::now();
        const auto r = Occt::rolling_ball_fillet(sh.blade, sh.blade_face, sh.hub, blade.definition().channel.hub(), fo);
        const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::printf("stations %4d rows %2d: ok=%d g1=%.3f deg ball=%.2e arc %.1f..%.1f deg  %.1f s\n", ns, rows, r.ok, r.max_g1_deg, r.max_ball_error, r.min_arc_deg, r.max_arc_deg, dt);
    }
}
