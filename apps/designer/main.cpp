// bladegen_designer: the BladeGen desktop designer (Qt Widgets, MDI).
//
//   bladegen_designer [design.json] [--screenshot out.png] [--self-test out_dir]
//
// --screenshot  evaluates the design, renders the default layout and writes the window image
// --self-test   additionally cycles every view content over every layer, writes one image per
//               content to out_dir and exits non-zero when a view fails (CI, offscreen)
#include "DesignerWindow.hpp"
#include "OcctView.hpp"

#include <QApplication>
#include <QMdiSubWindow>
#include <QSet>
#include <QToolBox>
#include <QCommandLineParser>
#include <QDir>
#include <QImage>
#include <QTimer>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

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
        // the Blade Sections view in its three coordinate systems (maximised: readable images)
        auto* vsub = qobject_cast<QMdiSubWindow*>(v->parentWidget());
        if (vsub) vsub->showMaximized();
        for (const auto sc : {Designer::SectionCoords::MRTheta, Designer::SectionCoords::ZRTheta, Designer::SectionCoords::XYAxial}) {
            v->sections = sc;
            v->setContent(Designer::Content::BladeSections);
            QApplication::processEvents();
            const QImage img = v->image();
            if (img.isNull() || img.width() < 10) {
                std::fprintf(stderr, "self-test: empty Blade Sections image (coordinates %d)\n", static_cast<int>(sc));
                return 18;
            }
            img.save(QString("%1/sections_%2.png").arg(dir).arg(static_cast<int>(sc)));
            ++n;
        }
        v->sections = Designer::SectionCoords::MRTheta;
        if (vsub) {
            vsub->showNormal();
            w.gridLayout(6);
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
        // re-stagger of the whole blade: every layer's beta_LE moves by the angle
        {
            std::vector<double> b0;
            for (const auto& L : w.session().result->layers) b0.push_back(L.beta_le);
            w.session().design.restagger.angle_deg = 15.0;  // large: the camber / TE placement must follow
            w.submit();
            w.settle();
            const auto& Ls = w.session().result->layers;
            for (std::size_t k = 0; k < Ls.size(); ++k)
                if (!w.session().result->ok() || std::abs(Ls[k].beta_le - b0[k] - 15.0) > 1e-9) {
                    std::fprintf(stderr, "self-test: the re-stagger did not turn layer %zu by 15 deg\n", k);
                    return 12;
                }
            v->layer = 4;  // the camber line and the TE marker follow the re-staggered section
            v->setContent(Designer::Content::BladeToBlade);
            if (auto* s = qobject_cast<QMdiSubWindow*>(v->parentWidget())) s->showMaximized();
            QApplication::processEvents();
            v->image().save(dir + "/b2b_restagger_tip.png");
            if (auto* s = qobject_cast<QMdiSubWindow*>(v->parentWidget())) s->showNormal();
            w.gridLayout(6);
            w.session().design.restagger.angle_deg = 0.0;
            w.submit();
            w.settle();
        }
        // 3D Solid: the OCCT shapes are built; the viewer renders when OpenGL is available
        {
            v->setContent(Designer::Content::Solid3D);
            QString err = v->waitForSolid();
            if (!err.isEmpty()) {
                std::fprintf(stderr, "self-test: the 3-D solid was not built: %s\n", qPrintable(err));
                return 10;
            }
            auto* ov = v->occtView();
            if (ov->available()) {
                using P = Designer::OcctView::Part;
                using L = Designer::OcctView::Look;
                auto save = [&](const QString& name) {
                    const QImage img = v->image();
                    QSet<QRgb> colors;
                    for (int y = 0; y < img.height(); y += 3)
                        for (int x = 0; x < img.width(); x += 3) colors.insert(img.pixel(x, y));
                    img.save(dir + "/" + name);
                    return static_cast<int>(colors.size());
                };
                const int c1 = save("solid_single.png");
                v->setAnnulus(true);
                err = v->waitForSolid();
                if (!err.isEmpty()) {
                    std::fprintf(stderr, "self-test: the annulus was not built: %s\n", qPrintable(err));
                    return 11;
                }
                const int c2 = save("solid_annulus.png");
                ov->setLook(P::Hub, L::Wireframe);
                ov->setLook(P::Casing, L::Hidden);
                ov->setLook(P::Blade, L::Shaded);
                const int c3 = save("solid_annulus_wire.png");
                ov->setLook(P::Blade, L::Solid);
                ov->setLook(P::Hub, L::Solid);
                ov->setLook(P::Casing, L::Shaded);
                v->setAnnulus(false);
                (void)v->waitForSolid();
                // the exact OCCT solid
                v->setExact(true);
                err = v->waitForSolid();
                if (!err.isEmpty()) {
                    std::fprintf(stderr, "self-test: the exact OCCT solid was not built: %s\n", qPrintable(err));
                    return 15;
                }
                const int c4 = save("solid_exact.png");
                v->setExact(false);
                if (std::min({c1, c2, c3, c4}) < 50) {
                    std::fprintf(stderr, "self-test: the OCCT view rendered blank images (%d, %d, %d, %d colours)\n", c1, c2, c3, c4);
                    return 13;
                }
                std::printf("self-test: OCCT viewer renders (design mesh: single blade, annulus, wireframe; exact solid)\n");
            } else {
                v->setExact(true);  // the exact solid still builds without OpenGL
                err = v->waitForSolid();
                v->setExact(false);
                if (!err.isEmpty()) {
                    std::fprintf(stderr, "self-test: the exact OCCT solid was not built: %s\n", qPrintable(err));
                    return 15;
                }
                std::printf("self-test: OCCT viewer not available here (%s); mesh and exact solid were built\n", qPrintable(ov->unavailableReason()));
            }
            v->setContent(Designer::Content::BladeToBlade);
        }
        // the equal-size layouts: 4, 2 and 6 views
        for (const int k : {4, 2, 6}) {
            w.gridLayout(k);
            QApplication::processEvents();
            QCoreApplication::sendPostedEvents(nullptr, QEvent::DeferredDelete);
            QList<QMdiSubWindow*> subs;
            for (auto* cv : w.views())
                if (auto* s = qobject_cast<QMdiSubWindow*>(cv->parentWidget()); s && !s->isHidden()) subs << s;
            bool equal = subs.size() == k;
            for (auto* s : subs)
                equal = equal && std::abs(s->width() - subs.front()->width()) <= 1 && std::abs(s->height() - subs.front()->height()) <= 1;
            if (!equal) {
                std::fprintf(stderr, "self-test: the %d-view layout has %d views or unequal sizes\n", k, static_cast<int>(subs.size()));
                return 14;
            }
        }
        // a 3D Solid view opened after the layouts closed views still gets its viewer
        {
            auto* sv = w.addView(Designer::Content::Solid3D);
            QApplication::processEvents();
            const QString err = sv->waitForSolid();
            const bool first = v->occtView()->available();
            if (!err.isEmpty() || sv->occtView()->available() != first) {
                std::fprintf(stderr, "self-test: a new 3D Solid view failed (%s; viewer %s)\n", qPrintable(err),
                             qPrintable(sv->occtView()->unavailableReason()));
                return 16;
            }
            sv->image().save(dir + "/solid_new_view.png");
            // closing one 3D Solid view leaves the others working (they share one GL driver)
            auto* sv2 = w.addView(Designer::Content::Solid3D);
            (void)sv2->waitForSolid();
            if (auto* sub = qobject_cast<QMdiSubWindow*>(sv->parentWidget())) sub->close();
            QApplication::processEvents();
            QCoreApplication::sendPostedEvents(nullptr, QEvent::DeferredDelete);
            if (first) {
                const QImage img = sv2->image();
                QSet<QRgb> colors;
                for (int y = 0; y < img.height(); y += 3)
                    for (int x = 0; x < img.width(); x += 3) colors.insert(img.pixel(x, y));
                if (colors.size() < 50) {
                    std::fprintf(stderr, "self-test: a 3D Solid view stopped rendering after another one was closed\n");
                    return 17;
                }
            }
            w.gridLayout(6);
        }
        if (auto* box = w.findChild<QToolBox*>()) {  // the meridional page (re-stagger inputs)
            const int page = box->currentIndex();
            box->setCurrentIndex(0);
            QApplication::processEvents();
            box->grab().save(dir + "/panel_meridional.png");
            box->setCurrentIndex(page);
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
