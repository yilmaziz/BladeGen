// bladegen: command-line front end of the BladeGen engine.
//
//   bladegen profile <profile.json> [--csv curvature.csv]
//       quality report of a profile (BladeGen or blade2d JSON); optional curvature CSV
//   bladegen build <blade.json> <out.step|.stp|.igs|.stl|.brep> [--fillet R_mm] [--unit mm|m]
//       C2 NURBS blade, trimmed solid, optional hub fillet; writes the blade and hub+blade
//   bladegen fit <sides.csv> <out.json> [--ss x,y] [--ps x,y] [--negate-y] [--fairness w]
//       inverse design from published coordinates (unstaggered frame)
//   bladegen naca65 <cl0> <t_max/c> [--variant kovach_sandercock|emery_base] [--csv sides.csv] [--fit out.json]
//       parametric NACA 65-(cl0 A10) t/c section: quality report, side coordinates and curvature,
//       optionally the custom profile parameters fitted to it (NACA 65 + Custom)
//   bladegen export-sections <blade.json> --span s[,s...] [--out f.json] [--n-edge 40] [--n-body 160] [--channel-points 201]
//       versioned JSON (format "bladegen.sections") of the sections at the span fractions (0 = hub):
//       SS / PS points LE -> TE in (m', theta) and (z, r, theta), LE / TE, streamline, channel, pitch
//   bladegen export-sections naca65 <cl0> <t_max/c> --variant kovach_sandercock|emery_base --chord c
//            --stagger deg --solidity s [--stations 101] [--n-edge 40] [--out f.json]
//       the planar NACA 65 cascade section (axial x, tangential y), staggered, pitch = chord / solidity
//   bladegen demo <out dir>
//       writes the demo rotor definition and its STEP files
//   bladegen version
#include "fit/InverseFit.hpp"
#include "fit/NacaFit.hpp"
#include "naca/Naca65.hpp"
#include "io/BladeJson.hpp"
#include "io/SectionExport.hpp"
#include "occt/OcctBlade.hpp"
#include "profile/Quality.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace BladeGen;
namespace fs = std::filesystem;

namespace {

int usage()
{
    std::fprintf(stderr,
                 "usage:\n"
                 "  bladegen profile <profile.json> [--csv curvature.csv]\n"
                 "  bladegen build <blade.json> <out.step> [--fillet R_mm] [--unit mm|m]\n"
                 "  bladegen fit <sides.csv> <out.json> [--ss x,y] [--ps x,y] [--negate-y] [--fairness w]\n"
                 "  bladegen naca65 <cl0> <t_max/c> [--variant kovach_sandercock|emery_base] [--csv sides.csv] [--fit out.json]\n"
                 "  bladegen export-sections <blade.json> --span s[,s...] [--out f.json] [--n-edge 40] [--n-body 160]\n"
                 "                           [--channel-points 201]\n"
                 "  bladegen export-sections naca65 <cl0> <t_max/c> --variant kovach_sandercock|emery_base --chord c\n"
                 "                           --stagger deg --solidity s [--stations 101] [--n-edge 40] [--out f.json]\n"
                 "  bladegen demo <out dir>\n"
                 "  bladegen version\n");
    return 2;
}

std::map<std::string, std::string> options(int argc, char** argv, int from)
{
    std::map<std::string, std::string> o;
    for (int i = from; i < argc; ++i) {
        std::string a = argv[i];
        if (a.rfind("--", 0) != 0) throw std::invalid_argument("unexpected argument '" + a + "'");
        if (a == "--negate-y") o[a] = "1";
        else if (i + 1 < argc) o[a] = argv[++i];
        else throw std::invalid_argument("missing value for " + a);
    }
    return o;
}

int cmd_profile(int argc, char** argv)
{
    if (argc < 3) return usage();
    const auto opt = options(argc, argv, 3);
    const auto p = IO::read_profile_json(argv[2]);
    const Profile::Profile2D prof(p);
    const auto q = Profile::check_quality(prof);
    std::printf("connector %s\n%s\n", Profile::to_string(p.connector).c_str(), q.summary().c_str());
    for (const auto& j : q.junctions)
        std::printf("  %-20s gap %.1e  angle %.1e deg  dkappa %.1e\n", j.name.c_str(), j.gap, Math::rad2deg(j.tangent_angle),
                    j.curvature_jump);
    if (opt.count("--csv")) {
        std::ofstream f(opt.at("--csv"));
        f << "side,s_norm,mprime,theta,kappa,piece\n";
        for (int sd = 0; sd < 2; ++sd) {
            const auto& c = sd == 0 ? q.ss : q.ps;
            for (std::size_t i = 0; i < c.xy.size(); ++i)
                f << (sd == 0 ? "SS" : "PS") << "," << c.s[i] / c.s.back() << "," << c.xy[i].x() << "," << c.xy[i].y() << ","
                  << c.kappa[i] << "," << Profile::to_string(static_cast<Profile::Piece>(c.piece[i])) << "\n";
        }
    }
    return q.ok() ? 0 : 3;  // 3: built, but with geometry errors
}

void write_all(const Occt::BladeShapes& s, const fs::path& out, Occt::Unit unit)
{
    Occt::write_shape(s.blade, out.string(), unit);
    const fs::path hb = out.parent_path() / (out.stem().string() + "_hub_blade" + (s.fillet_ok ? "_fillet" : "") + out.extension().string());
    Occt::write_shape(s.blade_on_hub, hb.string(), unit);
    std::printf("wrote %s and %s\n", out.string().c_str(), hb.string().c_str());
}

int cmd_build(int argc, char** argv)
{
    if (argc < 4) return usage();
    const auto opt = options(argc, argv, 4);
    const auto def = IO::read_blade_json(argv[2]);
    const double fillet = opt.count("--fillet") ? 1e-3 * std::stod(opt.at("--fillet")) : 0.0;
    const auto unit = opt.count("--unit") && opt.at("--unit") == "m" ? Occt::Unit::Metre : Occt::Unit::Millimetre;
    const Blade::Blade blade(def);
    for (const auto& st : def.stations) {
        const auto q = Profile::check_quality(Profile::Profile2D(st.profile));
        if (!q.ok()) std::printf("station span %.3f: %s\n", st.span, q.summary().c_str());
    }
    const auto s = Occt::build_blade(blade, fillet);
    std::printf("%s: blade solid valid %d, volume %.6g m^3, skin interpolation error %.1e m\n", def.name.c_str(),
                s.blade_valid ? 1 : 0, s.blade_volume, s.skin.max_error);
    if (fillet > 0.0) std::printf("hub fillet: %s\n", s.fillet_message.c_str());
    fs::create_directories(fs::absolute(argv[3]).parent_path());
    write_all(s, argv[3], unit);
    return s.blade_valid ? 0 : 1;
}

std::vector<Math::Vec2> read_side(const std::string& path, const std::string& cols, bool negate)
{
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    const auto comma = cols.find(',');
    if (comma == std::string::npos) throw std::invalid_argument("columns must be given as x,y");
    const std::string cx = cols.substr(0, comma), cy = cols.substr(comma + 1);
    std::string line;
    std::vector<std::string> header;
    std::vector<Math::Vec2> out;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        std::vector<std::string> cells;
        std::stringstream s(line);
        std::string c;
        while (std::getline(s, c, ',')) cells.push_back(c);
        if (header.empty()) { header = cells; continue; }
        std::size_t ix = header.size(), iy = header.size();
        for (std::size_t k = 0; k < header.size(); ++k) {
            if (header[k] == cx) ix = k;
            if (header[k] == cy) iy = k;
        }
        if (ix == header.size() || iy == header.size()) throw std::invalid_argument("columns " + cols + " not in " + path);
        if (ix >= cells.size() || iy >= cells.size() || cells[ix].empty() || cells[iy].empty()) continue;
        const double x = std::stod(cells[ix]), y = std::stod(cells[iy]);
        if (std::isnan(x) || std::isnan(y)) continue;
        const Math::Vec2 q(x, negate ? -y : y);
        if (out.empty() || (q - out.back()).norm() > 0.0) out.push_back(q);
    }
    return out;
}

int cmd_fit(int argc, char** argv)
{
    if (argc < 4) return usage();
    const auto opt = options(argc, argv, 4);
    const bool neg = opt.count("--negate-y") > 0;
    const auto ss = read_side(argv[2], opt.count("--ss") ? opt.at("--ss") : "upper_x,upper_y", neg);
    const auto ps = read_side(argv[2], opt.count("--ps") ? opt.at("--ps") : "lower_x,lower_y", neg);
    Fit::FitOptions fo;
    fo.fairness_weight = opt.count("--fairness") ? std::stod(opt.at("--fairness")) : 1.0;
    const auto r = Fit::fit_profile(ss, ps, Profile::ProfileParams::cda_start(), fo);
    const auto q = Profile::check_quality(Profile::Profile2D(r.params));
    std::printf("%zu + %zu points: RMS %.4f %% chord, max %.4f %%, %d evaluations, %.1f s\n%s\n", ss.size(), ps.size(),
                100 * r.rms, 100 * r.max, r.evaluations, r.seconds, q.summary().c_str());
    for (const auto& b : r.at_bound) std::printf("  at its bound: %s\n", b.c_str());
    IO::write_profile_json(r.params, argv[3]);
    std::printf("wrote %s\n", argv[3]);
    return 0;
}

int cmd_naca65(int argc, char** argv)
{
    if (argc < 4) return usage();
    const auto opt = options(argc, argv, 4);
    Naca::Naca65Params p{std::stod(argv[2]), std::stod(argv[3]),
                         Naca::variant_from_string(opt.count("--variant") ? opt.at("--variant") : "kovach_sandercock")};
    const Naca::Naca65Section s(p);
    const auto& m = s.metadata();
    const auto q = Naca::check_quality(s);
    std::printf("%s (%s): t_max %.6f c at station %.4f, R_LE %.6f c, R_TE %.4f c\n", m.designation.c_str(),
                Naca::to_string(p.variant).c_str(), m.max_thickness, m.max_thickness_station, m.le_radius, m.te_radius);
    std::printf("equivalent arc: camber %.3f deg, chi_LE %.3f deg, chi_TE %.3f deg\n%s\n", Math::rad2deg(m.camber_angle),
                Math::rad2deg(m.chi_le), Math::rad2deg(m.chi_te), q.summary().c_str());
    if (opt.count("--csv")) {
        std::ofstream f(opt.at("--csv"));
        f << "side,x,y,kappa,piece\n";
        f.precision(12);
        for (int sd = 0; sd < 2; ++sd) {
            const auto& c = sd == 0 ? q.ss : q.ps;
            for (std::size_t i = 0; i < c.xy.size(); ++i)
                f << (sd == 0 ? "SS" : "PS") << "," << c.xy[i].x() << "," << c.xy[i].y() << "," << c.kappa[i] << ","
                  << Profile::to_string(static_cast<Profile::Piece>(c.piece[i])) << "\n";
        }
    }
    if (opt.count("--fit")) {
        const auto r = Fit::fit_naca65(s);
        std::printf("NACA 65 + Custom: RMS %.4f %% c, max %.4f %% c (%.1f s)\n", 100.0 * r.rms, 100.0 * r.max, r.seconds);
        IO::write_profile_json(r.params, opt.at("--fit"));
        std::printf("wrote %s\n", opt.at("--fit").c_str());
    }
    return q.ok() ? 0 : 3;
}

/// A whole-string number (std::stod accepts a trailing junk suffix).
double number(const std::string& s, const std::string& what)
{
    std::size_t used = 0;
    double v = 0.0;
    try {
        v = std::stod(s, &used);
    } catch (const std::exception&) {
        used = 0;
    }
    if (used == 0 || used != s.size() || !std::isfinite(v)) throw std::invalid_argument(what + ": not a number: '" + s + "'");
    return v;
}

int integer(const std::string& s, const std::string& what)
{
    const double v = number(s, what);
    if (v != std::floor(v) || std::abs(v) > 1e6) throw std::invalid_argument(what + ": not an integer: '" + s + "'");
    return static_cast<int>(v);
}

const std::string& required(const std::map<std::string, std::string>& opt, const std::string& key)
{
    const auto it = opt.find(key);
    if (it == opt.end()) throw std::invalid_argument("missing option " + key);
    return it->second;
}

void check_known(const std::map<std::string, std::string>& opt, const std::vector<std::string>& known)
{
    for (const auto& [k, v] : opt)
        if (std::find(known.begin(), known.end(), k) == known.end()) throw std::invalid_argument("unknown option " + k);
}

int emit_sections(const nlohmann::ordered_json& j, const std::map<std::string, std::string>& opt, const std::string& what)
{
    if (!opt.count("--out")) {
        std::cout << IO::dump_sections(j);
        return 0;
    }
    const fs::path out = opt.at("--out");
    if (out.has_parent_path()) fs::create_directories(fs::absolute(out).parent_path());
    IO::write_sections(j, out.string());
    std::printf("wrote %s (%s)\n", out.string().c_str(), what.c_str());
    return 0;
}

int cmd_export_sections(int argc, char** argv)
{
    if (argc < 3) return usage();
    if (std::string(argv[2]) == "naca65") {
        if (argc < 5) return usage();
        const auto opt = options(argc, argv, 5);
        check_known(opt, {"--variant", "--chord", "--stagger", "--solidity", "--stations", "--n-edge", "--out"});
        IO::CascadeExportOptions o;
        o.naca = {number(argv[3], "cl0"), number(argv[4], "t_max/c"), Naca::variant_from_string(required(opt, "--variant"))};
        o.chord = number(required(opt, "--chord"), "--chord");
        o.stagger_deg = number(required(opt, "--stagger"), "--stagger");
        o.solidity = number(required(opt, "--solidity"), "--solidity");
        if (opt.count("--stations")) o.n_stations = integer(opt.at("--stations"), "--stations");
        if (opt.count("--n-edge")) o.n_edge = integer(opt.at("--n-edge"), "--n-edge");
        const auto j = IO::export_naca65_cascade(o, BLADEGEN_VERSION);
        char buf[160];
        std::snprintf(buf, sizeof buf, "%s %s, %zu SS + %zu PS points, pitch %.6g m",
                      j["source"]["naca"]["designation"].get<std::string>().c_str(), Naca::to_string(o.naca.variant).c_str(),
                      j["suction_side"]["xy"].size(), j["pressure_side"]["xy"].size(), j["cascade"]["pitch"].get<double>());
        return emit_sections(j, opt, buf);
    }
    const auto opt = options(argc, argv, 3);
    check_known(opt, {"--span", "--n-edge", "--n-body", "--channel-points", "--out"});
    IO::DesignExportOptions o;
    o.spans.clear();
    std::stringstream ss(required(opt, "--span"));
    for (std::string item; std::getline(ss, item, ',');) o.spans.push_back(number(item, "--span"));
    if (opt.count("--n-edge")) o.n_edge = integer(opt.at("--n-edge"), "--n-edge");
    if (opt.count("--n-body")) o.n_body = integer(opt.at("--n-body"), "--n-body");
    if (opt.count("--channel-points")) o.n_channel = integer(opt.at("--channel-points"), "--channel-points");
    const fs::path design = argv[2];
    const auto def = IO::read_blade_json(design.string());
    // the file name only: no directories of the user's machine in the output
    const auto j = IO::export_design_sections(def, o, design.filename().string(), BLADEGEN_VERSION);
    char buf[160];
    std::snprintf(buf, sizeof buf, "%s: %zu section(s), %zu points per side, %d blades", def.name.c_str(),
                  j["sections"].size(), j["sections"][0]["suction_side"]["mprime_theta"].size(), def.n_blades);
    return emit_sections(j, opt, buf);
}

int cmd_demo(int argc, char** argv)
{
    if (argc < 3) return usage();
    const fs::path dir = argv[2];
    fs::create_directories(dir);
    const auto def = Blade::BladeDefinition::demo_rotor();
    IO::write_blade_json(def, (dir / "demo_rotor.json").string());
    const auto s = Occt::build_blade(Blade::Blade(def), 0.003);
    std::printf("demo rotor: blade valid %d, volume %.4f cm^3, fillet: %s\n", s.blade_valid ? 1 : 0, 1e6 * s.blade_volume,
                s.fillet_message.c_str());
    write_all(s, dir / "demo_rotor.step", Occt::Unit::Millimetre);
    return 0;
}

}  // namespace

int main(int argc, char** argv)
{
    if (argc < 2) return usage();
    const std::string cmd = argv[1];
    try {
        if (cmd == "profile") return cmd_profile(argc, argv);
        if (cmd == "build") return cmd_build(argc, argv);
        if (cmd == "fit") return cmd_fit(argc, argv);
        if (cmd == "naca65") return cmd_naca65(argc, argv);
        if (cmd == "export-sections") return cmd_export_sections(argc, argv);
        if (cmd == "demo") return cmd_demo(argc, argv);
        if (cmd == "version") {
            std::printf("BladeGen %s\n", BLADEGEN_VERSION);
            return 0;
        }
    } catch (const std::exception& e) {
        std::fprintf(stderr, "bladegen %s: %s\n", cmd.c_str(), e.what());
        return 1;
    }
    return usage();
}
