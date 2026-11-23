// DesignerWindow: the BladeGen designer main window (MDI).  Views are sub-windows whose
// content is picked from their right-click menu; the parameter panel and the issue list
// are docks.  Geometry runs on the Engine thread; every view redraws when it finishes.  The
// TURBOlenza menu opens the blade-to-blade solver in its own window (turbolenza/TurbolenzaWindow).
#pragma once

#include "ContentView.hpp"
#include "Engine.hpp"

#include <QMainWindow>

#include <QFutureWatcher>

class QMdiArea;
class QTreeWidget;
class QLabel;

namespace Designer {

class ParameterPanel;
class TurbolenzaWindow;

class DesignerWindow : public QMainWindow
{
    Q_OBJECT
public:
    DesignerWindow();
    ~DesignerWindow() override;

    /// quiet: report a failure on stderr instead of a message box (command line).
    bool openDesign(const QString& path, bool quiet = false);
    ContentView* addView(Content c);
    void defaultLayout();
    /// n views (6, 4 or 2) in equal cells: 3 x 2, 2 x 2, 2 x 1.  Views beyond n close (the
    /// active one is kept); missing ones are opened.  The grid follows the window size until
    /// Tile or Cascade is chosen.
    void gridLayout(int n);
    /// Waits for the engine and redraws (screenshots, tests).
    void settle();
    [[nodiscard]] Session& session() { return s_; }
    [[nodiscard]] QList<ContentView*> views() const;
    /// Sends the session design to the engine.
    void submit();
    /// Opens (or raises) the TURBOlenza window (menu TURBOlenza > TURBOlenza...).
    TurbolenzaWindow* openTurbolenza();

protected:
    void closeEvent(QCloseEvent* e) override;
    bool eventFilter(QObject* o, QEvent* e) override;

private:
    void buildMenus();
    void onResult();
    void refreshViews(bool refit = false);
    void updateIssues();
    void newDesign(bool naca);
    void openDialog();
    bool save(bool as);
    void exportStep();
    void exportObj();
    void exportCsv();
    void exportImage();
    void meshPointsDialog();
    void startFit();
    void updateTitle();
    void arrangeGrid();

    Session s_;
    Engine engine_;
    QMdiArea* mdi_;
    ParameterPanel* panel_;
    QTreeWidget* issues_;
    QLabel* status_;
    QFutureWatcher<QString>* fit_watch_;
    bool dirty_ = false;
    int grid_n_ = 0;  ///< the equal-size grid in use (0: none)
    TurbolenzaWindow* turbolenza_ = nullptr;
};

}  // namespace Designer
