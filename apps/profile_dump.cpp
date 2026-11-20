// profile_dump: writes the side samples (m', theta, s, kappa, piece) of a built-in profile as CSV.
#include "profile/Profile2D.hpp"
#include "profile/Quality.hpp"
#include <cstdio>
#include <string>
using namespace BladeGen::Profile;
int main(int argc, char** argv)
{
    std::string which = argc > 1 ? argv[1] : "naca";
    auto prm = which == "cda" ? ProfileParams::cda_start() : ProfileParams::naca65_1210();
    if (argc > 2) prm.connector = connector_mode_from_string(argv[2]);
    for (int i = 3; i + 1 < argc; i += 2) set_param(prm, argv[i], std::stod(argv[i + 1]));
    const Profile2D p(prm);
    const auto q = check_quality(p);
    std::fprintf(stderr, "%s\n", q.summary().c_str());
    std::printf("side,x,y,s,kappa,piece\n");
    for (int sd = 0; sd < 2; ++sd) {
        const auto& c = sd == 0 ? q.ss : q.ps;
        for (std::size_t i = 0; i < c.xy.size(); ++i)
            std::printf("%d,%.12g,%.12g,%.12g,%.12g,%d\n", sd, c.xy[i].x(), c.xy[i].y(), c.s[i], c.kappa[i], c.piece[i]);
    }
}
