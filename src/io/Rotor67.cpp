#include "io/Rotor67.hpp"

#include "blade/Blade.hpp"
#include "math/Interpolation.hpp"
#include "math/Numerics.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>

namespace BladeGen::IO {

namespace {

constexpr double kCm = 0.01;

std::vector<std::string> split_tabs(std::string line)
{
    if (!line.empty() && line.back() == '\r') line.pop_back();
    std::vector<std::string> out;
    std::stringstream s(line);
    std::string cell;
    while (std::getline(s, cell, '\t')) out.push_back(cell);
    return out;
}

bool is_number(const std::string& s)
{
    if (s.empty()) return false;
    char* end = nullptr;
    std::strtod(s.c_str(), &end);
    return end && *end == '\0';
}

std::map<std::string, std::vector<std::vector<std::string>>> read_blocks(const std::string& path)
{
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + path);
    std::map<std::string, std::vector<std::vector<std::string>>> blocks;
    std::string line, name;
    bool first = true;
    while (std::getline(f, line)) {
        if (first && line.size() >= 3 && static_cast<unsigned char>(line[0]) == 0xEF) line = line.substr(3);  // BOM
        first = false;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.rfind("===", 0) == 0) {
            name = line;
            name.erase(std::remove(name.begin(), name.end(), '='), name.end());
            name.erase(0, name.find_first_not_of(' '));
            name.erase(name.find_last_not_of(' ') + 1);
            blocks[name];
        } else if (!name.empty() && line.find_first_not_of(" \t") != std::string::npos) {
            blocks[name].push_back(split_tabs(line));
        }
    }
    return blocks;
}

double num(const std::vector<std::string>& r, std::size_t i) { return std::stod(r.at(i)); }

}  // namespace

double Rotor67Section::theta(int surface, double mm) const { return (surface == 1 ? surf1 : surf2).eval(mm)[0]; }

std::vector<Vec2> Rotor67Section::profile_m_theta(int n_surface, int n_arc) const
{
    std::vector<Vec2> out;
    auto cosine = [](double a, double b, int n) {
        std::vector<double> v;
        for (int i = 0; i < n; ++i) v.push_back(a + (b - a) * 0.5 * (1.0 - std::cos(Math::kPi * i / (n - 1))));
        return v;
    };
    auto local_angle = [&](int j, int edge) {
        const EdgeCircle& c = edge == 0 ? le : te;
        const double mt = m_tan[j - 1][edge];
        return std::atan2(c.rc * (theta(j, mt) - c.thc), mt - c.mc);
    };
    auto arc = [&](int edge) {
        const EdgeCircle& c = edge == 0 ? le : te;
        double a, b;
        if (edge == 1) {  // TE, through phi = 0, from surface 1 to surface 2
            a = local_angle(1, 1);
            b = local_angle(2, 1);
        } else {          // LE, through phi = pi, from surface 2 to surface 1
            a = std::fmod(local_angle(2, 0) + 2.0 * Math::kPi, 2.0 * Math::kPi);
            b = std::fmod(local_angle(1, 0) + 2.0 * Math::kPi, 2.0 * Math::kPi);
        }
        std::vector<Vec2> p;
        for (int i = 1; i < n_arc - 1; ++i) {
            const double f = a + (b - a) * i / (n_arc - 1);
            p.emplace_back(c.mc + c.rho * std::cos(f), c.thc + c.rho * std::sin(f) / c.rc);
        }
        return p;
    };
    for (double mm : cosine(m_tan[0][0], m_tan[0][1], n_surface)) out.emplace_back(mm, theta(1, mm));
    for (const auto& q : arc(1)) out.push_back(q);
    for (double mm : cosine(m_tan[1][1], m_tan[1][0], n_surface)) out.emplace_back(mm, theta(2, mm));
    for (const auto& q : arc(0)) out.push_back(q);
    out.push_back(out.front());
    return out;
}

std::vector<Vec2> Rotor67Section::profile_conformal(int n_surface, int n_arc) const
{
    auto p = profile_m_theta(n_surface, n_arc);
    for (auto& q : p) q.x() = streamline.mprime_of_m(q.x());
    return p;
}

std::vector<Vec3> Rotor67Section::profile_xyz(int n_surface, int n_arc) const
{
    std::vector<Vec3> out;
    for (const auto& q : profile_m_theta(n_surface, n_arc)) {
        const Vec2 zr = streamline.at_m(q.x());
        out.push_back(Blade::cartesian(zr.x(), zr.y(), q.y()));
    }
    return out;
}

Blade::MeridionalCurve Rotor67Data::hub_curve() const { return Blade::MeridionalCurve::through_points(hub); }
Blade::MeridionalCurve Rotor67Data::shroud_curve() const { return Blade::MeridionalCurve::through_points(shroud); }

Rotor67Data read_rotor67(const std::string& path)
{
    auto blocks = read_blocks(path);
    for (const char* b : {"Notes", "Sections", "Edge circles", "Flowpath"})
        if (!blocks.count(b)) throw std::runtime_error(std::string("Rotor 67 file: block '") + b + "' missing");
    Rotor67Data d;
    for (const auto& r : blocks["Notes"]) {
        if (r.size() < 2 || !is_number(r[1])) continue;
        if (r[0] == "Number of rotor blades") d.n_blades = static_cast<int>(num(r, 1));
        else if (r[0] == "Tip clearance") d.tip_clearance = num(r, 1) * kCm;
        else if (r[0] == "Rotating hub start Z") d.rotating_hub[0] = num(r, 1) * kCm;
        else if (r[0] == "Rotating hub end Z") d.rotating_hub[1] = num(r, 1) * kCm;
    }
    std::map<int, Rotor67Section> secs;
    for (const auto& r : blocks["Sections"]) {
        if (r.size() < 10 || !is_number(r[0])) continue;
        auto& s = secs[static_cast<int>(num(r, 0))];
        s.index = static_cast<int>(num(r, 0));
        s.z.push_back(num(r, 2) * kCm);
        s.r.push_back(num(r, 3) * kCm);
        s.th1.push_back(num(r, 4));
        s.th2.push_back(num(r, 5));
        s.beta1.push_back(Math::deg2rad(num(r, 8)));
        s.beta2.push_back(Math::deg2rad(num(r, 9)));
    }
    for (const auto& r : blocks["Edge circles"]) {
        if (r.size() < 7 || !is_number(r[0])) continue;
        auto it = secs.find(static_cast<int>(num(r, 0)));
        if (it == secs.end()) continue;
        EdgeCircle c{num(r, 2) * kCm, num(r, 3) * kCm, num(r, 4) * kCm, num(r, 5) * kCm, num(r, 6)};
        (r[1].rfind("Leading", 0) == 0 ? it->second.le : it->second.te) = c;
    }
    for (const auto& r : blocks["Flowpath"]) {
        if (r.size() < 4 || !is_number(r[0])) continue;
        d.hub.emplace_back(num(r, 1) * kCm, num(r, 2) * kCm);
        d.shroud.emplace_back(num(r, 1) * kCm, num(r, 3) * kCm);
    }
    for (auto& [idx, s] : secs) {
        if (s.z.size() < 4) throw std::runtime_error("Rotor 67 section " + std::to_string(idx) + ": too few points");
        if (!(s.le.rho > 0.0) || !(s.te.rho > 0.0))
            throw std::runtime_error("Rotor 67 section " + std::to_string(idx) + ": edge circle missing");
        std::vector<Vec2> zr;
        for (std::size_t k = 0; k < s.z.size(); ++k) zr.emplace_back(s.z[k], s.r[k]);
        s.streamline = Blade::MeridionalCurve::through_points(zr);
        // meridional distance of the K points (they are the interpolation nodes)
        const auto tk = Math::chord_length_params<2>(zr);
        s.m.assign(tk.size(), 0.0);
        for (std::size_t k = 1; k < tk.size(); ++k) s.m[k] = s.m[k - 1] + Math::arc_length(s.streamline.curve(), tk[k - 1], tk[k]);
        auto surface = [&](const std::vector<double>& th, const std::vector<double>& beta) {
            std::vector<Math::Vec<1>> y;
            for (double v : th) y.push_back(Math::Vec<1>(v));
            const double s0 = std::tan(beta.front()) / s.r.front(), s1 = std::tan(beta.back()) / s.r.back();
            return Math::interpolate_cubic<1>(s.m, y, Math::EndCondition<1>::derivative(Math::Vec<1>(s0)),
                                              Math::EndCondition<1>::derivative(Math::Vec<1>(s1)));
        };
        s.surf1 = surface(s.th1, s.beta1);
        s.surf2 = surface(s.th2, s.beta2);
        // tangency: foot of the perpendicular from the circle centre to each surface
        for (int j = 1; j <= 2; ++j)
            for (int edge = 0; edge < 2; ++edge) {
                const EdgeCircle& c = edge == 0 ? s.le : s.te;
                const auto& f = j == 1 ? s.surf1 : s.surf2;
                auto g = [&](double mm) {
                    const auto dd = f.derivs(mm, 1);
                    return (mm - c.mc) + c.rc * (dd[0][0] - c.thc) * c.rc * dd[1][0];
                };
                const double lo = edge == 0 ? 0.0 : c.mc - c.rho, hi = edge == 0 ? c.mc + c.rho : s.m.back();
                double best = std::numeric_limits<double>::quiet_NaN(), err = 1e300;
                const int n = 400;
                double x0 = lo, g0 = g(lo);
                for (int i = 1; i <= n; ++i) {
                    const double x1 = lo + (hi - lo) * i / n, g1 = g(x1);
                    if ((g0 > 0.0) != (g1 > 0.0)) {
                        const double root = Math::brent_root(g, x0, x1, 1e-15);
                        const double dist = std::hypot(root - c.mc, c.rc * (f.eval(root)[0] - c.thc));
                        if (std::abs(dist - c.rho) < err) { err = std::abs(dist - c.rho); best = root; }
                    }
                    x0 = x1;
                    g0 = g1;
                }
                if (!std::isfinite(best))
                    throw std::runtime_error("Rotor 67 section " + std::to_string(idx) + ": no tangency point");
                s.m_tan[j - 1][edge] = best;
            }
        d.sections.push_back(std::move(s));
    }
    if (d.sections.empty()) throw std::runtime_error("Rotor 67 file: no sections");
    return d;
}

}  // namespace BladeGen::IO
