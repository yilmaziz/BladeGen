// Engine: evaluates designs on a worker thread; a newer request replaces a pending one,
// so knobs stay responsive.  The window reads latest() when finished() arrives.
#pragma once

#include "design/BladeDesign.hpp"

#include <QObject>

#include <condition_variable>
#include <memory>
#include <mutex>
#include <optional>
#include <thread>

namespace Designer {

struct Session
{
    BladeGen::Design::BladeDesign design = BladeGen::Design::BladeDesign::defaults();
    std::shared_ptr<const BladeGen::Design::DesignResult> result;
    int active_layer = 2;
    QString path;  ///< current design file
};

class Engine : public QObject
{
    Q_OBJECT
public:
    Engine();
    ~Engine() override;

    void submit(const BladeGen::Design::BladeDesign& d);
    [[nodiscard]] std::shared_ptr<const BladeGen::Design::DesignResult> latest() const;
    [[nodiscard]] BladeGen::Design::Evaluator& evaluator() { return ev_; }
    [[nodiscard]] bool busy() const;
    /// Blocks until every submitted design is evaluated (tests, screenshots).
    void wait();

signals:
    void finished();

private:
    void run();

    BladeGen::Design::Evaluator ev_;
    mutable std::mutex mu_;
    std::condition_variable cv_;
    std::optional<BladeGen::Design::BladeDesign> job_;
    std::shared_ptr<const BladeGen::Design::DesignResult> latest_;
    bool working_ = false, stop_ = false;
    std::thread thread_;
};

}  // namespace Designer
