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
//                      Mesh independence  the same case on the checked standard grid levels
//                              (run_mesh_study), a results table; Run mesh study, Stop, Save
//   Run                        the log of the passes (also printed to the terminal) and the live
//                              residual plot
//   Post-processing  Residuals, Contours (pressure, Mach), Blade loading (Mach, p, Dunavant S; a
//                              viscous run also its inviscid k = 0 loading), Boundary layer (H,
//                              theta/c, delta*/c, c_f, Re_theta with Aungier's Fig. 5-15; the
//                              coupling history):
//                              separate windows after a run; "Save results..."; Mesh study:
//                              beta_out, Mesh study: blade loading (after a mesh study)
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
#include "turbolenza/MeshStudy.hpp"

#include <QFutureWatcher>
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
class QTableWidget;
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
    /// Builds the case of the current settings on a worker thread and shows its grid (or the
    /// message on the Mesh page).  wait: block (with the event loop running) until it is built and
    /// return whether it was; otherwise return whether the build started.
    bool generateMesh(bool wait = false);
    void startRun();
    void stopRun();
    [[nodiscard]] bool running() const { return worker_.joinable(); }
    /// run_mesh_study of the current settings on the checked levels, on a worker thread.
    void startStudy();
    void stopStudy();
    [[nodiscard]] bool studyRunning() const { return study_worker_.joinable(); }
    /// write_mesh_study + the PNGs of the mesh-study plots into dir; the files written.
    QStringList saveStudy(const QString& dir, bool overwrite);
    /// export_run + the PNGs of the post-processing plots into dir; the files written.
    QStringList saveResults(const QString& dir, bool overwrite);
    /// --turbolenza-self-test: Mid layer, very_coarse mesh, a `passes`-pass solve, post-processing
    /// views and the saved results in dir; then the same case with the one-way boundary layer to
    /// convergence, its window and files in dir/viscous.  0 when everything worked.
    int selfTest(const QString& dir, int passes);

protected:
    void closeEvent(QCloseEvent* e) override;

private:
    struct Channel;       // worker thread -> GUI
    struct StudyChannel;  // mesh-study worker thread -> GUI
    struct MeshJob
    {
        std::optional<BladeGen::Turbolenza::Case> c;
        QString err;
    };
    class Observer;
    class PostWindow;

    // pages (TurbolenzaPages.cpp)
    QWidget* bladePage();
    QWidget* cascadePage();
    QWidget* meshPage();
    QWidget* solverPage();
    QWidget* runPage();
    QWidget* postPage();
    QWidget* studyPage();
    void syncBlade();          ///< widgets <- layer_
    void updateBladeEnabled();
    void bladeChanged();       ///< layer_ edited: the section, cascade and mesh are stale
    void refreshProfile();
    void refreshCascade();
    void setDerivedOmega(double omega);  ///< the derived omega field, without invalidating the cascade
    void meshFinished();
    void updateSolverEnabled();
    [[nodiscard]] BladeGen::Turbolenza::CaseSettings settings() const;
    void setCascadeDefaults(double stagger_deg, double solidity);

    // section / cascade
    const BladeGen::Turbolenza::BladeSection* section();  ///< built on demand; null + section_error_
    /// Not in v1: the layer copy on its stream surface (stream_surface_blade, the Cascade page's stagger as
    /// BladeGen's), built on demand; null + stream_error_.
    const BladeGen::Turbolenza::StreamSurfaceBlade* streamBlade();
    [[nodiscard]] bool layerMode() const;     ///< the Cascade page's stream surface is the BladeGen layer's
    [[nodiscard]] bool rotatingRow() const;   ///< the Cascade page's row rotates
    void updateRowEnabled();
    /// The case of settings st: build_case(st, section) or, on the layer's stream surface, build_case(st, blade).
    struct CaseInput
    {
        BladeGen::Turbolenza::BladeSection section;
        std::optional<BladeGen::Turbolenza::StreamSurfaceBlade> stream;
    };
    [[nodiscard]] std::optional<CaseInput> caseInput(QString& error);
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
    void showBoundaryLayer();
    void redrawBoundaryLayer();
    void setP0FromReynolds();          ///< Cascade page: Pt from a reference chord Reynolds number
    [[nodiscard]] bool referenceCase() const;  ///< the case of Aungier's Fig. 5-15 (beta1 45, M1 0.65, sigma 1, stagger 33)
    void saveDialog();
    void redrawContours();
    void redrawLoading();
    void updateRunButtons();   ///< a normal run and a mesh study exclude each other

    // mesh study
    void pollStudy();
    void studyFinished();
    void addStudyRow(const BladeGen::Turbolenza::MeshStudyEntry& e);
    void redrawStudyResiduals(const QString& status);
    void showStudyBetaOut();
    void showStudyLoading();
    void redrawStudyBetaOut();
    void redrawStudyLoading();
    void saveStudyDialog();

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
    QDoubleSpinBox *re_ref_ = nullptr, *m_ref_ = nullptr;
    QLabel* re_info_ = nullptr;
    QLabel* derived_ = nullptr;
    PlotWidget* cascade_plot_ = nullptr;
    QTimer* cascade_timer_ = nullptr;
    // cascade page: the stream surface and the row (not in v1)
    QComboBox *surface_ = nullptr, *row_ = nullptr;
    QDoubleSpinBox *rpm_ = nullptr, *alpha1_ = nullptr, *mach_abs_ = nullptr;
    QLabel* triangle_info_ = nullptr;
    PlotWidget *sheet_plot_ = nullptr, *triangle_plot_ = nullptr;
    std::optional<BladeGen::Turbolenza::StreamSurfaceBlade> stream_;
    QString stream_error_;

    // mesh page
    QComboBox* level_ = nullptr;
    QComboBox* edge_nodes_ = nullptr;  ///< first / last blade nodes: inside (default) or touching (experiment)
    QSpinBox *n_passage_ = nullptr, *n_eta_ = nullptr;
    QDoubleSpinBox *ext_up_ = nullptr, *ext_down_ = nullptr;
    QLabel *mesh_info_ = nullptr, *mesh_error_ = nullptr;
    QPushButton* gen_button_ = nullptr;
    QFutureWatcher<MeshJob>* mesh_watch_ = nullptr;
    PlotWidget* grid_plot_ = nullptr;

    // solver page
    QLineEdit* case_name_ = nullptr;
    QComboBox *mode_ = nullptr, *exit_start_ = nullptr, *linear_ = nullptr;
    QCheckBox *fixed_is_kappa_ = nullptr, *vel_is_rho_ = nullptr, *stop_sonic_ = nullptr, *measure_ = nullptr;
    QDoubleSpinBox *fixed_angle_ = nullptr, *lambda_rho_ = nullptr, *lambda_k_ = nullptr, *exit_start_deg_ = nullptr;
    QLineEdit *tol_rho_ = nullptr, *tol_w_ = nullptr, *tol_k_ = nullptr;
    QSpinBox *max_rho_ = nullptr, *max_k_ = nullptr, *threads_ = nullptr;
    QPushButton *run_ = nullptr, *stop_ = nullptr, *run2_ = nullptr, *stop2_ = nullptr;
    QGroupBox* bl_box_ = nullptr;
    QCheckBox* viscous_ = nullptr;
    QComboBox *bl_coupling_ = nullptr, *bl_cf_ = nullptr, *bl_wake_ = nullptr;
    QDoubleSpinBox *bl_omega_ = nullptr, *bl_re_tr_ = nullptr, *bl_hk_ = nullptr;
    QLineEdit* bl_tol_ = nullptr;
    QSpinBox* bl_max_ = nullptr;

    // run page
    QPlainTextEdit* log_ = nullptr;
    QLabel* run_status_ = nullptr;
    PlotWidget* residuals_ = nullptr;  ///< the one residual plot
    QVBoxLayout* run_plot_box_ = nullptr;
    QLabel* run_plot_note_ = nullptr;

    // post-processing
    QLabel* report_ = nullptr;
    PostWindow *res_win_ = nullptr, *contour_win_ = nullptr, *loading_win_ = nullptr, *bl_win_ = nullptr;
    PlotWidget *contour_plot_ = nullptr, *loading_plot_ = nullptr;
    PlotWidget *bl_plot_ = nullptr, *bl_change_plot_ = nullptr, *bl_beta_plot_ = nullptr;
    QComboBox *contour_field_ = nullptr, *loading_q_ = nullptr, *bl_q_ = nullptr, *bl_x_ = nullptr;
    QCheckBox *loading_inviscid_ = nullptr, *loading_dunavant_ = nullptr, *bl_overlay_ = nullptr;
    QCheckBox* isolines_ = nullptr;
    QSpinBox* bands_ = nullptr;

    // mesh study
    std::unique_ptr<StudyChannel> sch_;
    std::thread study_worker_;
    BladeGen::Turbolenza::CaseSettings study_settings_;
    std::optional<BladeGen::Turbolenza::MeshStudyResult> study_;  ///< the levels finished so far / the last study
    std::string study_level_;                                     ///< the level of study_records_
    std::vector<BladeGen::Turbolenza::PassRecord> study_records_;  ///< the current level's history
    qint64 study_plot_ms_ = 0;
    std::vector<QCheckBox*> study_check_;  ///< one per kGridLevels entry
    QPushButton *study_run_ = nullptr, *study_stop_ = nullptr, *study_save_ = nullptr;
    QLabel* study_status_ = nullptr;
    QTableWidget* study_table_ = nullptr;
    PlotWidget* study_residuals_ = nullptr;
    PostWindow *study_beta_win_ = nullptr, *study_load_win_ = nullptr;
    PlotWidget *study_beta_plot_ = nullptr, *study_load_plot_ = nullptr;
    QComboBox *study_x_ = nullptr, *study_q_ = nullptr;
};

}  // namespace Designer
