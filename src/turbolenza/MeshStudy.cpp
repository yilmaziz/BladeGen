// turbolenza/MeshStudy.cpp: the mesh-independence study on the run layer (see MeshStudy.hpp).
// No numerics: build_case, run_case and export_run of Case.hpp, once per grid level.
#include "turbolenza/MeshStudy.hpp"

#include <chrono>
#include <cmath>
#include <exception>
#include <filesystem>
#include <format>
#include <fstream>
#include <numbers>
#include <stdexcept>
#include <utility>

namespace BladeGen::Turbolenza {

namespace {

/// Forwards every pass of one level's solve to the study's callback.
class PassForwarder : public RunObserver
{
public:
    PassForwarder(const MeshStudyCallbacks& cb, const std::string& level) : cb_(cb), level_(level) {}
    void update(const PassRecord& r) override
    {
        if (cb_.pass) cb_.pass(level_, r);
    }

private:
    const MeshStudyCallbacks& cb_;
    const std::string& level_;
};

std::string num(double v)
{
    return std::isfinite(v) ? std::format("{:.10g}", v) : std::string();
}

std::string quoted(const std::string& s)
{
    if (s.find_first_of(",\"\r\n") == std::string::npos) return s;
    std::string out = "\"";
    for (const char ch : s) {
        if (ch == '"') out += '"';
        out += ch;
    }
    return out + "\"";
}

}  // namespace

MeshStudyResult run_mesh_study(const CaseSettings& base, const BladeSection& section, const std::vector<std::string>& levels,
                               const MeshStudyCallbacks& callbacks, const std::atomic<bool>* cancel,
                               const StreamSurfaceBlade* stream)
{
    MeshStudyResult out;
    for (std::size_t k = 0; k < levels.size(); ++k) {
        if (cancel && cancel->load()) {
            out.cancelled = true;
            break;
        }
        const auto t0 = std::chrono::steady_clock::now();
        MeshStudyEntry e;
        e.level = levels[k];
        CaseSettings st = base;
        st.grid_level = levels[k];
        st.n_passage.reset();
        st.n_eta.reset();
        try {
            Case c = stream ? build_case(st, *stream) : build_case(st, section);
            e.n_m = c.grid.n_m();
            e.n_eta = c.grid.n_eta();
            e.total_nodes = e.n_m * e.n_eta;
            if (c.standard_info) {
                e.dm_over_chord = c.standard_info->dm_over_chord;
                e.total_nodes = c.standard_info->total_nodes;
            }
            if (callbacks.level_started) callbacks.level_started(e.level, static_cast<int>(k), e.n_m, e.n_eta);
            PassForwarder obs(callbacks, e.level);
            RunResult r = run_case(c, &obs, cancel);
            e.converged = r.converged;
            e.passes = r.iterations;
            e.error = r.error;
            if (r.solution) e.beta_out_deg = r.solution->exit_angle * 180.0 / std::numbers::pi;
            if (r.summary) e.mach_max = r.summary->mach_max;
            else if (!r.records.empty()) e.mach_max = r.records.back().mach_max;
            e.c = std::move(c);
            e.result = std::move(r);
        } catch (const std::exception& ex) {
            e.converged = false;
            e.error = std::string("the case cannot be run: ") + ex.what();
        }
        if (!e.converged && e.error.empty()) e.error = "not converged";
        e.wall_time_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (callbacks.level_finished) callbacks.level_finished(e);
        out.levels.push_back(std::move(e));
        if (cancel && cancel->load()) {
            out.cancelled = k + 1 < levels.size();
            break;
        }
    }
    return out;
}

std::vector<std::string> write_mesh_study(const MeshStudyResult& result, const std::string& directory, bool overwrite)
{
    namespace fs = std::filesystem;
    const fs::path dir(directory);
    const std::string table = "mesh_study.csv";
    if (!overwrite && fs::exists(dir / table))
        throw std::runtime_error(std::format("{} exists: write_mesh_study overwrites only when asked", (dir / table).string()));
    fs::create_directories(dir);
    std::vector<std::string> files;
    {
        std::ofstream csv(dir / table, std::ios::binary);
        if (!csv) throw std::runtime_error("cannot write " + (dir / table).string());
        csv << "level,n_m,n_eta,nodes,dm_over_chord,converged,passes,beta_out_deg,mach_max,wall_s,error\r\n";
        for (const MeshStudyEntry& e : result.levels)
            csv << std::format("{},{},{},{},{},{},{},{},{},{},{}\r\n", quoted(e.level), e.n_m, e.n_eta, e.total_nodes, num(e.dm_over_chord),
                               e.converged ? "true" : "false", e.passes, num(e.beta_out_deg), num(e.mach_max), num(e.wall_time_s),
                               quoted(e.error));
        if (!csv) throw std::runtime_error("cannot write " + (dir / table).string());
    }
    files.push_back(table);
    for (const MeshStudyEntry& e : result.levels) {
        if (!e.c || !e.result) continue;
        for (const std::string& f : export_run(*e.result, *e.c, (dir / e.level).string(), overwrite))
            files.push_back((fs::path(e.level) / f).generic_string());
    }
    return files;
}

}  // namespace BladeGen::Turbolenza
