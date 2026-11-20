#include "Engine.hpp"

namespace Designer {

Engine::Engine() : thread_([this] { run(); }) {}

Engine::~Engine()
{
    {
        std::lock_guard<std::mutex> g(mu_);
        stop_ = true;
    }
    cv_.notify_all();
    thread_.join();
}

void Engine::submit(const BladeGen::Design::BladeDesign& d)
{
    {
        std::lock_guard<std::mutex> g(mu_);
        job_ = d;
    }
    cv_.notify_all();
}

std::shared_ptr<const BladeGen::Design::DesignResult> Engine::latest() const
{
    std::lock_guard<std::mutex> g(mu_);
    return latest_;
}

bool Engine::busy() const
{
    std::lock_guard<std::mutex> g(mu_);
    return working_ || job_.has_value();
}

void Engine::wait()
{
    std::unique_lock<std::mutex> l(mu_);
    cv_.wait(l, [this] { return !working_ && !job_.has_value(); });
}

void Engine::run()
{
    for (;;) {
        BladeGen::Design::BladeDesign d;
        {
            std::unique_lock<std::mutex> l(mu_);
            cv_.wait(l, [this] { return stop_ || job_.has_value(); });
            if (stop_) return;
            d = std::move(*job_);
            job_.reset();
            working_ = true;
        }
        auto r = std::make_shared<const BladeGen::Design::DesignResult>(ev_.evaluate(d));
        {
            std::lock_guard<std::mutex> g(mu_);
            latest_ = std::move(r);
            working_ = false;
        }
        cv_.notify_all();
        emit finished();
    }
}

}  // namespace Designer
