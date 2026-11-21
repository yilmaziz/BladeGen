// DesignerWindow: the BladeGen designer main window (MDI).  Views are sub-windows whose
// content is picked from their right-click menu; the parameter panel and the issue list
// are docks.  Geometry runs on the Engine thread; every view redraws when it finishes.
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
    /// Waits for the engine and redraws (screenshots, tests).
    void settle();
    [[nodiscard]] Session& session() { return s_; }
    [[nodiscard]] QList<ContentView*> views() const;
    /// Sends the session design to the engine.
    void submit();

protected:
    void closeEvent(QCloseEvent* e) override;

private:
    void buildMenus();
    void onResult();
    void refreshViews(bool refit = false);
    void updateIssues();
    void newDesign();
    void openDialog();
    bool save(bool as);
    void exportStep();
    void exportObj();
    void exportCsv();
    void exportImage();
    void meshPointsDialog();
    void startFit();
    void updateTitle();

    Session s_;
    Engine engine_;
    QMdiArea* mdi_;
    ParameterPanel* panel_;
    QTreeWidget* issues_;
    QLabel* status_;
    QFutureWatcher<QString>* fit_watch_;
    bool dirty_ = false;
};

}  // namespace Designer
