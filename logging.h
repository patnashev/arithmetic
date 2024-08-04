#pragma once

#include <memory>
#include <vector>
#include <map>
#include "inputnum.h"
#include "file.h"

class Progress
{
public:
    Progress() { }

    void reset() { _cur_stage = 0; _cur_progress = 0; _time_total = 0; _time_stage = 0; _time_op = 0; _timer = 0; _op_count = 0; }
    void add_stage(double cost) { _total_cost += cost; _costs.push_back(cost); }
    void next_stage() { _cur_stage++; update(0, 0); }
    void skip_stage() { _cur_stage++; }
    void update(double progress, int op_count);
    void time_init(double elapsed);
    void set_parent(Progress* parent) { _parent = parent; }

    std::vector<double>& costs() { return _costs; }
    double progress_stage() { return _cur_progress; }
    double progress_total() { if (_cur_stage >= _costs.size()) return 1; double cost = 0; for (int i = 0; i < _cur_stage; i++) cost += _costs[i]; return (cost + _costs[_cur_stage]*_cur_progress)/_total_cost; }
    double cost_total() { return _total_cost; }
    double time_total() { return _time_total + _time_stage; }
    double time_op() { return _time_op; }
    int op_count() { return _op_count; }
    int num_stages() { return (int)_costs.size(); }
    int cur_stage() { return _cur_stage; }
    std::string& param(const std::string& name) { return _params[name]; }
    int param_int(const std::string& name) { std::string& val = _params[name]; return std::strtol(val.data(), nullptr, 10); }
    double param_double(const std::string& name) { std::string& val = _params[name]; return std::strtod(val.data(), nullptr); }
    std::map<std::string, std::string>& params() { return _params; }

private:
    std::vector<double> _costs;
    double _total_cost = 0;
    int _cur_stage = 0;
    double _cur_progress = 0;
    double _time_total = 0;
    double _time_stage = 0;
    double _time_op = 0;
    double _timer = 0;
    int _op_count = 0;
    Progress* _parent = nullptr;
    std::map<std::string, std::string> _params;
};

class Logging
{
public:
    static const int LEVEL_DEBUG = 0;
    static const int LEVEL_INFO = 1;
    static const int LEVEL_PROGRESS = 2;
    static const int LEVEL_WARNING = 3;
    static const int LEVEL_ERROR = 4;
    static const int LEVEL_RESULT = 5;
    int level_result_success = LEVEL_WARNING;
    int level_result_not_success = LEVEL_INFO;

public:
    Logging(int level = LEVEL_INFO) : _level(level) { }
    virtual ~Logging() { }

    void debug(const char* message...);
    void info(const char* message...);
    void warning(const char* message...);
    void error(const char* message...);
    void result(bool success, const char* message...);

    virtual void report(const std::string& message, int level);
    virtual void report_progress();
    virtual void report_param(const std::string& name, const std::string& value) { progress().param(name) = value; }
    virtual void report_param(const std::string& name, int value) { progress().param(name) = std::to_string(value); }
    virtual void report_param(const std::string& name, double value);
    virtual void report_factor(InputNum& input, const arithmetic::Giant& f);

    virtual void result_save(const std::string& message);

    virtual bool state_save_flag() { return false; }
    virtual void state_save() { }
    virtual void progress_save();
    virtual void heartbeat() { }

    virtual void file_progress(File* file_progress);
    void file_result(const std::string& filename) { _file_result = filename; }
    void file_factor(const std::string& filename) { _file_factor = filename; }
    void file_log(const std::string& filename) { _file_log = filename; }

    int level() { return _level; }
    Progress& progress() { return _progress; }
    const std::string& prefix() { return _prefix; }
    void set_prefix(const std::string& prefix) { _prefix = prefix; }

protected:
    int _level;
    Progress _progress;
    File* _file_progress = nullptr;
    std::string _file_result = "result.txt";
    std::string _file_factor = "factors.txt";
    std::string _file_log;
    std::string _prefix;
    bool _print_prefix = true;
    bool _overwrite_line = false;
};

class SubLogging : public Logging
{
public:
    SubLogging(Logging& parent, int level = LEVEL_WARNING) : Logging(level), _parent(parent) { progress().set_parent(&parent.progress()); }

    virtual void report(const std::string& message, int level) override { _parent.report(_prefix.empty() ? message : _prefix + message, level); }
    virtual void report_param(const std::string& name, const std::string& value) override { Logging::report_param(name, value); _parent.report_param(name, value); }
    virtual void report_param(const std::string& name, int value) override { Logging::report_param(name, value); _parent.report_param(name, value); }
    virtual void report_param(const std::string& name, double value) override { Logging::report_param(name, value); _parent.report_param(name, value); }
    virtual void report_factor(InputNum& input, const arithmetic::Giant& f) override { _parent.report_factor(input, f); }
    virtual void progress_save() override { Logging::progress_save(); _parent.progress_save(); }
    virtual void heartbeat() override { _parent.heartbeat(); }

protected:
    Logging& _parent;
};
