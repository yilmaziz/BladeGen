// bladegen_designer: the BladeGen desktop designer (Qt Widgets, MDI).
//
//   bladegen_designer [design.json] [--screenshot out.png] [--self-test out_dir]
//
// --screenshot  evaluates the design, renders the default layout and writes the window image
// --self-test   additionally cycles every view content over every layer, writes one image per
//               content to out_dir and exits non-zero when a view fails (CI, offscreen)
#include "DesignerWindow.hpp"

#include <QApplication>
#include <QCommandLineParser>
#include <QDir>
#include <QImage>
#include <QTimer>

#include <cmath>
#include <cstdio>

int main(int argc, char** argv)
{
    QApplication app(argc, argv);
    QApplication::setApplicationName("BladeGen Designer");
    QApplication::setApplicationVersion(BLADEGEN_VERSION);
    QApplication::setStyle("Fusion");
    QCommandLineParser cli;
    cli.setApplicationDescription("BladeGen designer: meridional channel, five layers, conformal 3-D blade");
    cli.addHelpOption();
    cli.addVersionOption();
    cli.addPositionalArgument("design", "design file (bladegen-design-2, or a v1.0.1 bladegen-design-1)");
    QCommandLineOption shot("screenshot", "write the window image and exit", "png");
    QCommandLineOption self("self-test", "render every view content and exit", "dir");
    cli.addOption(shot);
    cli.addOption(self);
    cli.process(app);

    Designer::DesignerWindow w;
    if (!cli.positionalArguments().isEmpty() && !w.openDesign(cli.positionalArguments().front(), cli.isSet(shot) || cli.isSet(self))) return 2;
    w.show();
    if (cli.isSet(self)) {
        const QString dir = cli.value(self);
        QDir().mkpath(dir);
        w.settle();
        const auto R = w.session().result;
        if (!R || !R->ok()) {
            std::fprintf(stderr, "self-test: the design did not evaluate\n");
            return 3;
        }
        auto views = w.views();
        if (views.isEmpty()) return 4;
        auto* v = views.front();
        int n = 0;
        for (const auto& grp : Designer::contentGroups())
            for (auto c : grp)
                for (int k = 0; k < static_cast<int>(w.session().design.layers.size()); ++k) {
                    v->layer = k;
                    v->setContent(c);
                    QApplication::processEvents();
                    const QImage img = v->image();
                    if (img.isNull() || img.width() < 10) {
                        std::fprintf(stderr, "self-test: empty image for %s\n", qPrintable(Designer::contentName(c)));
                        return 5;
                    }
                    if (k == 2) img.save(QString("%1/view_%2.png").arg(dir).arg(static_cast<int>(c), 2, 10, QChar('0')));
                    ++n;
                }
        // a geometry change must reach the views
        const double before = w.session().result->layers[2].beta_le;
        w.session().design.layers[2].stagger_deg += 3.0;
        w.submit();
        w.settle();
        const double after = w.session().result->layers[2].beta_le;
        if (std::abs(after - before - 3.0) > 1e-6) {
            std::fprintf(stderr, "self-test: a stagger change of 3 deg moved beta_LE by %g\n", after - before);
            return 7;
        }
        // dragging the free SS control point CP2 of the mid layer moves ss.ctrl2
        v->layer = 2;
        v->setContent(Designer::Content::ProfileSection);
        v->simulateDrag(Designer::ContentView::kSsHandle + 2, QPointF(0.48, 0.11));
        const auto& c = w.session().design.layers[2].custom.ss;
        if (std::abs(c.ctrl2_mprime - 0.48) > 1e-12 || std::abs(c.ctrl2_theta - 0.11) > 1e-12) {
            std::fprintf(stderr, "self-test: the control-point drag did not reach the design\n");
            return 8;
        }
        w.submit();
        w.settle();
        if (!w.session().result->ok()) {
            std::fprintf(stderr, "self-test: the dragged profile did not evaluate\n");
            return 9;
        }
        w.grab().save(dir + "/window.png");
        std::printf("self-test: %d view renders OK, geometry %.1f ms\n", n, 1000.0 * w.session().result->seconds);
        return 0;
    }
    if (cli.isSet(shot)) {
        w.settle();
        QApplication::processEvents();
        const bool ok = w.grab().save(cli.value(shot));
        return ok ? 0 : 6;
    }
    return QApplication::exec();
}
