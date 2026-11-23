// TurbolenzaWindow: the TURBOlenza window of BladeGen Designer, the blade-to-blade compressible
// potential-flow solver (Aungier 2003, Ch. 5; the C++ port of TURBOlenza v1).  A tree on the left,
// one page per node on the right:
//
//   Pre-processing   Blade     the layer's parametrization (BladeGen's Knob class): a COPY of one
//                              layer of the open design ("Take from layer"), edited here without
//                              changing the 3-D design; the profile preview
//                    Cascade   chord, solidity, stagger, inlet flow, gas; derived pitch, angle of
//                              attack, kappa_le / kappa_te, throat o/s; three blades with the throat
//                    Mesh      grid level (or explicit n_passage x n_eta), extensions; "Generate mesh"
//                    Solver    every CaseSettings solver input, cores, linear solver; Run, Stop
//   Run                        the log of the passes (also printed to the terminal) and the live
//                              residual plot
//   Post-processing  Residuals, Contours (pressure, Mach), Blade loading (Mach, p, Dunavant S):
//                              separate windows after a run; "Save results..."
//
// The window holds no numerics: it fills CaseSettings, calls section_from_layer (BladeGen ->
// BladeSection), make_cascade, build_case and run_case, and draws.  The solve runs on a worker
// thread; Stop sets run_case's cancel flag.  The residual plot is ONE widget fed by the one record
// list of the run: it lives on the Run page during the run and moves into the Residuals window
// after it (v1 drew the history twice: a live window and a replayed figure).
#pragma once

#include "Engine.hpp"
#include "TurbolenzaPlots.hpp"

#include "turbolenza/BladeGenSource.hpp"
#include "turbolenza/Case.hpp"

#include <QMainWindow>

#include <atomic>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <thread>
#include <vector>

class QCheckBox;
class QComboBox;
class QDoubleSpinBox;
class QGroupBox;
class QLabel;
class QLineEdit;
class QPlainTextEdit;
class QPushButton;
class QSpinBox;
class QStackedWidget;
class QTimer;
class QTreeWidget;
class QTreeWidgetItem;
class QVBoxLayout;

namespace Designer {

class Knob;
class PlotWidget;

class TurbolenzaWindow : public QMainWindow
{
public:
    explicit TurbolenzaWindow(Session& s, QWidget* parent = nullptr);
    ~TurbolenzaWindow() override;

    /// Copies layer k of the open design (and its stagger and 2-D solidity) into the window.
    void takeLayer(int k);
    /// Builds the case of the current settings and shows its grid; false (and the message on the
    /// Mesh page) when the section, cascade or grid cannot be built.
    bool generateMesh();
    void startRun();
    void stopRun();
    [[nodiscard]] bool running() const { return worker_.joinable(); }
    /// export_run + the PNGs of the post-processing plots into dir; the files written.
    QStringList saveResults(const QString& dir, bool overwrite);
    /// --turbolenza-self-test: Mid layer, very_coarse mesh, a `passes`-pass solve, post-processing
    /// views and the saved results in dir.  0 when everything worked.
    int selfTest(const QString& dir, int passes);

protected:
    void closeEvent(QCloseEvent* e) override;

private:
    struct Channel;       // worker thread -> GUI
    class Observer;
    class PostWindow;

    // pages (TurbolenzaPages.cpp)
    QWidget* bladePage();
    QWidget* cascadePage();
    QWidget* meshPage();
    QWidget* solverPage();
    QWidget* runPage();
    QWidget* postPage();
    void syncBlade();          ///< widgets <- layer_
    void updateBladeEnabled();
    void bladeChanged();       ///< layer_ edited: the section, cascade and mesh are stale
    void refreshProfile();
    void refreshCascade();
    void updateSolverEnabled();
    [[nodiscard]] BladeGen::Turbolenza::CaseSettings settings() const;
    void setCascadeDefaults(double stagger_deg, double solidity);

    // section / cascade
    const BladeGen::Turbolenza::BladeSection* section();  ///< built on demand; null + section_error_
    void showPage(int index);
    void buildTree();
    void buildMenus();

    // run
    void pollRun();
    void runFinished();
    void updateLivePlot(bool final);
    void openPostWindows();
    void showResiduals();
    void showContours();
    void showLoading();
    void saveDialog();
    void redrawContours();
    void redrawLoading();

    Session& s_;
    BladeGen::Design::Evaluator ev_;
    BladeGen::Design::Layer layer_;
    int source_layer_ = 2;
    bool syncing_ = false;
    std::optional<BladeGen::Turbolenza::BladeSection> section_;
    BladeGen::Turbolenza::SectionFitReport fit_;
    QString section_error_;
    std::optional<BladeGen::Turbolenza::Cascade> cascade_;
    std::optional<BladeGen::Turbolenza::Case> mesh_case_;

    // run state
    std::unique_ptr<Channel> ch_;
    std::thread worker_;
    std::vector<BladeGen::Turbolenza::PassRecord> records_;  ///< the one history of the run
    std::optional<BladeGen::Turbolenza::Case> run_case_;
    std::optional<BladeGen::Turbolenza::RunResult> result_;
    BladeGen::Turbolenza::CaseSettings run_settings_;
    QTimer* poll_ = nullptr;
    qint64 last_plot_ms_ = 0;

    // layout
    QTreeWidget* tree_ = nullptr;
    QStackedWidget* stack_ = nullptr;
    std::vector<QTreeWidgetItem*> items_;
    QLabel* status_ = nullptr;

    // blade page
    QComboBox *layer_combo_ = nullptr, *method_ = nullptr, *variant_ = nullptr, *connector_ = nullptr;
    Knob *cl0_ = nullptr, *tc_ = nullptr;
    QGroupBox *naca_box_ = nullptr, *custom_box_ = nullptr;
    std::vector<std::pair<std::string, Knob*>> custom_;
    std::vector<Knob*> ps_ctrl_;
    QCheckBox* enforce_c1_ = nullptr;
    QPushButton* fit_button_ = nullptr;
    QLabel* blade_info_ = nullptr;
    PlotWidget* profile_plot_ = nullptr;
    QTimer* profile_timer_ = nullptr;

    // cascade page
    QDoubleSpinBox *chord_ = nullptr, *solidity_ = nullptr, *stagger_ = nullptr, *beta1_ = nullptr, *mach1_ = nullptr;
    QDoubleSpinBox *pt_ = nullptr, *tt_ = nullptr, *gamma_ = nullptr, *gas_r_ = nullptr, *omega_ = nullptr, *mdot_ = nullptr;
    QCheckBox *mach_given_ = nullptr, *mdot_given_ = nullptr;
    QLabel* derived_ = nullptr;
    PlotWidget* cascade_plot_ = nullptr;
    QTimer* cascade_timer_ = nullptr;

    // mesh page
    QComboBox* level_ = nullptr;
    QSpinBox *n_passage_ = nullptr, *n_eta_ = nullptr;
    QDoubleSpinBox *ext_up_ = nullptr, *ext_down_ = nullptr;
    QLabel *mesh_info_ = nullptr, *mesh_error_ = nullptr;
    PlotWidget* grid_plot_ = nullptr;

    // solver page
    QLineEdit* case_name_ = nullptr;
    QComboBox *mode_ = nullptr, *exit_start_ = nullptr, *linear_ = nullptr;
    QCheckBox *fixed_is_kappa_ = nullptr, *vel_is_rho_ = nullptr, *stop_sonic_ = nullptr, *measure_ = nullptr;
    QDoubleSpinBox *fixed_angle_ = nullptr, *lambda_rho_ = nullptr, *lambda_k_ = nullptr, *exit_start_deg_ = nullptr;
    QLineEdit *tol_rho_ = nullptr, *tol_w_ = nullptr, *tol_k_ = nullptr;
    QSpinBox *max_rho_ = nullptr, *max_k_ = nullptr, *threads_ = nullptr;
    QPushButton *run_ = nullptr, *stop_ = nullptr, *run2_ = nullptr, *stop2_ = nullptr;

    // run page
    QPlainTextEdit* log_ = nullptr;
    QLabel* run_status_ = nullptr;
    PlotWidget* residuals_ = nullptr;  ///< the one residual plot
    QVBoxLayout* run_plot_box_ = nullptr;
    QLabel* run_plot_note_ = nullptr;

    // post-processing
    QLabel* report_ = nullptr;
    PostWindow *res_win_ = nullptr, *contour_win_ = nullptr, *loading_win_ = nullptr;
    PlotWidget *contour_plot_ = nullptr, *loading_plot_ = nullptr;
    QComboBox *contour_field_ = nullptr, *loading_q_ = nullptr;
    QCheckBox* isolines_ = nullptr;
    QSpinBox* bands_ = nullptr;
};

}  // namespace Designer
