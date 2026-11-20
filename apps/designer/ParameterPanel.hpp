// ParameterPanel: every input of the design, grouped.  Meridional: the eight cylindrical
// coordinates of the initial meridional configuration and the hub / casing law (Bell &
// Mehta or linear).  Layer: hub, 25 %, mid, 75 %, tip; the construction method and its
// knobs (NACA 65, Custom families, NACA 65 + Custom fit), the placement.  Blade: number of
// blades and the 3-D mesh resolution.
#pragma once

#include "Engine.hpp"

#include <QWidget>

#include <functional>
#include <vector>

class QButtonGroup;
class QCheckBox;
class QComboBox;
class QDoubleSpinBox;
class QLabel;
class QPushButton;
class QSlider;
class QSpinBox;
class QTabWidget;
class QGroupBox;

namespace Designer {

/// Slider + spin box for one scalar; emits changed() on user edits only.
class Knob : public QWidget
{
    Q_OBJECT
public:
    Knob(const QString& label, double lo, double hi, int decimals, QWidget* parent = nullptr, const QString& tip = {});
    void setValue(double v);  ///< silent (no changed())
    [[nodiscard]] double value() const;

signals:
    void changed(double v);

private:
    QSlider* slider_;
    QDoubleSpinBox* spin_;
    double lo_, hi_;
};

class ParameterPanel : public QWidget
{
    Q_OBJECT
public:
    explicit ParameterPanel(Session& s, QWidget* parent = nullptr);
    /// Set every widget from the session (after loading or switching the layer).
    void sync();
    void setFitStatus(const QString& t, bool running);

signals:
    void designChanged();
    void activeLayerChanged();
    void fitRequested();

private:
    QWidget* meridionalPage();
    QWidget* layerPage();
    QWidget* bladePage();
    [[nodiscard]] BladeGen::Design::Layer& layer();
    void changed();
    void updateEnabled();

    Session& s_;
    bool syncing_ = false;
    // meridional
    struct Wall
    {
        QDoubleSpinBox *z_le, *r_le, *z_te, *r_te;
        QPushButton *bell, *linear;
    } hub_{}, casing_{};
    QDoubleSpinBox *ext_in_ = nullptr, *ext_out_ = nullptr;
    // layer
    QButtonGroup* layers_ = nullptr;
    QComboBox* method_ = nullptr;
    Knob *stagger_ = nullptr, *theta_le_ = nullptr, *solidity_ = nullptr, *cl0_ = nullptr, *tc_ = nullptr;
    QComboBox* variant_ = nullptr;
    QPushButton* fit_ = nullptr;
    QLabel* fit_status_ = nullptr;
    QGroupBox *naca_box_ = nullptr, *custom_box_ = nullptr;
    QLabel* custom_hint_ = nullptr;
    std::vector<std::pair<std::string, Knob*>> custom_;
    std::vector<Knob*> ps_ctrl_;
    QComboBox* connector_ = nullptr;
    QCheckBox* enforce_c1_ = nullptr;
    // blade
    QSpinBox *n_blades_ = nullptr, *rows_ = nullptr, *n_body_ = nullptr, *n_edge_ = nullptr;
    QDoubleSpinBox* clustering_ = nullptr;
    QCheckBox* mirror_ = nullptr;
};

}  // namespace Designer
