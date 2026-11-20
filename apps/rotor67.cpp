// rotor67: reconstruct NASA Rotor 67 with BladeGen.
//
//   rotor67 <data file> <output dir> [--fillet R_mm] [--no-solid] [--sections a,b,...]
//
// 1. reads the 14 published sections and the flowpath (io/Rotor67),
// 2. fits the BladeGen profile parameters to every section (fit/SectionFit),
// 3. builds the C2 NURBS blade through the fitted sections, trimmed on the hub and
//    on the casing contour offset by the tip clearance (TP-2879), optionally with
//    the 1.78 mm hub fillet, and writes STEP files and a parameter JSON.
#include "fit/SectionFit.hpp"
#include "io/BladeJson.hpp"
#include "io/Rotor67.hpp"
#include "occt/OcctBlade.hpp"
#include "profile/Quality.hpp"

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <future>
#include <string>
#include <thread>
#include <vector>

using namespace BladeGen;

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::fprintf(stderr, "usage: rotor67 <Rotor_14_Sections_and_Flowpath_EN.txt> <output dir> [--fillet R_mm] [--no-solid]\n");
        return 2;
    }
    const std::string data_file = argv[1];
    const std::filesystem::path out = argv[2];
    double fillet_mm = 1.78;
    bool solid = true;
    for (int i = 3; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--fillet" && i + 1 < argc) fillet_mm = std::stod(argv[++i]);
        else if (a == "--no-solid") solid = false;
    }
    std::filesystem::create_directories(out);
    try {
        const auto data = IO::read_rotor67(data_file);
        std::printf("NASA Rotor 67: %zu sections, %d blades, tip clearance %.4f mm\n", data.sections.size(), data.n_blades,
                    1e3 * data.tip_clearance);
        // ── fit every section (in parallel) ─────────────────────────────────────
        const auto n = data.sections.size();
        std::vector<std::future<Fit::SectionFitResult>> jobs;
        for (std::size_t k = 0; k < n; ++k) {
            jobs.push_back(std::async(std::launch::async, [&data, k, n] {
                const auto& s = data.sections[k];
                const auto loop = s.profile_conformal(161, 33);
                Fit::FitOptions opt;
                opt.fairness_weight = 0.3;
                return Fit::fit_section(loop, s.streamline, static_cast<double>(k) / static_cast<double>(n - 1),
                                        Profile::ProfileParams::cda_start(), opt);
            }));
        }
        std::vector<Fit::SectionFitResult> fits;
        std::printf("\n sec   chord'   stagger   RMS 2D     max 2D   | RMS 3D   max 3D  | quality        time\n");
        std::printf("       [-]      [deg]    [%%chord]   [%%chord] |  [um]     [um]   |\n");
        for (std::size_t k = 0; k < n; ++k) {
            fits.push_back(jobs[k].get());
            const auto& f = fits.back();
            const auto q = Profile::check_quality(Profile::Profile2D(f.fit.params));
            std::printf(" %2d   %7.4f  %7.2f   %8.4f   %8.4f  | %7.2f  %7.2f  | %d err %d warn  %5.1f s\n",
                        data.sections[k].index, f.frame.chord, Math::rad2deg(f.frame.stagger), 100 * f.fit.rms,
                        100 * f.fit.max, 1e6 * f.rms_3d, 1e6 * f.max_3d, q.count(Profile::Severity::Error),
                        q.count(Profile::Severity::Warning), f.fit.seconds);
        }
        // ── the blade definition, written as JSON ───────────────────────────────
        Blade::BladeDefinition def;
        def.name = "NASA_Rotor_67";
        def.n_blades = data.n_blades;
        def.tip_clearance = data.tip_clearance;
        def.interpolation = Blade::SpanInterpolation::Cubic;
        for (const auto& f : fits) def.stations.push_back(f.station);
        def.channel = Blade::MeridionalChannel(data.hub_curve(), data.shroud_curve());
        IO::write_blade_json(def, (out / "rotor67_bladegen.json").string());
        {
            std::ofstream csv(out / "rotor67_fit_errors.csv");
            csv << "section,rms_2d_pct_chord,max_2d_pct_chord,rms_3d_um,max_3d_um\n";
            for (std::size_t k = 0; k < n; ++k)
                csv << data.sections[k].index << "," << 100 * fits[k].fit.rms << "," << 100 * fits[k].fit.max << ","
                    << 1e6 * fits[k].rms_3d << "," << 1e6 * fits[k].max_3d << "\n";
        }
        if (!solid) return 0;
        // ── 3-D blade: C2 skin through the 14 fitted sections ────────────────────
        const Blade::Blade blade(def);
        std::vector<std::vector<Math::Vec3>> loops;
        for (const auto& s : blade.sections(60, 200)) loops.push_back(s.xyz);
        loops = Occt::extend_loops(loops, 1.0, 1.0);
        Occt::SolidOptions so;
        so.hub = data.hub_curve();
        so.casing = data.sections.back().streamline;  // TP-2879: tip from the casing contour of surface 14
        so.tip_clearance = data.tip_clearance;
        const auto inside = Occt::inside_point(blade.section(0.5, 60, 200));
        auto shapes = Occt::build_shapes(loops, so, inside);
        std::printf("\nblade solid: valid %d, volume %.3f cm^3, area %.2f cm^2, skin error %.2e m\n", shapes.blade_valid,
                    1e6 * shapes.blade_volume, 1e4 * shapes.blade_area, shapes.skin.max_error);
        // published sections against the solid's skin
        double worst = 0.0;
        for (const auto& f : fits)
            for (double d : Occt::distances(std::vector<Math::Vec3>(f.published_xyz.begin(), f.published_xyz.begin() + 40), shapes.blade_face))
                worst = std::max(worst, d);
        std::printf("published points to the NURBS skin: max %.2f um (first 40 points of every section)\n", 1e6 * worst);
        Occt::write_shape(shapes.blade, (out / "rotor67_blade.step").string());
        Occt::write_shape(shapes.blade_on_hub, (out / "rotor67_hub_blade.step").string());
        if (fillet_mm > 0.0) {
            so.fillet_radius = 1e-3 * fillet_mm;
            shapes = Occt::build_shapes(loops, so, inside);
            std::printf("hub fillet r = %.2f mm: %s (G1 %.3f deg, ball %.1e m)\n", fillet_mm, shapes.fillet_message.c_str(),
                        shapes.fillet_g1_deg, shapes.fillet_ball_error);
            if (shapes.fillet_ok) {
                Occt::write_shape(shapes.blade_on_hub, (out / "rotor67_hub_blade_fillet.step").string());
                Occt::write_shape(shapes.blade_with_fillet, (out / "rotor67_blade_fillet.step").string());
            }
        }
        std::printf("written to %s\n", out.string().c_str());
    } catch (const std::exception& e) {
        std::fprintf(stderr, "rotor67: %s\n", e.what());
        return 1;
    }
    return 0;
}
