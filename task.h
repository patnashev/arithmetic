#pragma once

#include <chrono>
#include <stdexcept>
#include <memory>
#include "arithmetic.h"
#include "inputnum.h"
#include "file.h"
#include "logging.h"

class TaskRestartException : public std::exception
{
public:
    TaskRestartException() { }
};

class TaskAbortException : public std::exception
{
public:
    TaskAbortException() { }
};

class TaskState
{
public:
    TaskState(char type) : _type(type), _written(false), _iteration(0) { }
    virtual ~TaskState() { }
    
    void set(int iteration) { _written = false; _iteration = iteration; }
    virtual bool read(Reader& reader);
    virtual void write(Writer& writer);
    void set_written() { _written = true; }

    char type() { return _type; }
    char version() { return 0; }
    bool is_written() { return _written; }
    int iteration() { return _iteration; }

protected:
    char _type;
    bool _written;
    int _iteration;
};

template<class State>
State* read_state(File* file)
{
    State* state = nullptr;
    if (file != nullptr)
    {
        state = new State();
        if (!file->read(*state))
        {
            delete state;
            state = nullptr;
        }
    }
    return state;
}

class Task
{
public:
    static int MULS_PER_STATE_UPDATE;
    static int DISK_WRITE_TIME;
    static int PROGRESS_TIME;

public:
    Task() { }
    virtual ~Task() { }

    static void abort() { _abort_flag = true; }
    static void abort_reset() { _abort_flag = false; }
    static bool abort_flag() { return _abort_flag; }

    virtual void run();

    arithmetic::GWArithmetic& gw() { return *_gw; }
    TaskState* state() { return _state.get(); }
    int iterations() { return _iterations; }
    int state_update_period() { return _state_update_period; }
    bool is_last(int iteration) { return iteration + 1 - (_state ? _state->iteration() : 0) >= _state_update_period || iteration + 1 == _iterations || abort_flag() || _logging->state_save_flag(); }
    virtual double progress() { return _state ? _state->iteration()/(double)iterations() : 0.0; }
    int ops() { return (int)(_op_count + (_gw != nullptr ? _gwstate->ops() - _op_base : 0)); }

protected:
    virtual void init(arithmetic::GWState* gwstate, File* file, TaskState* state, Logging* logging, int iterations);
    virtual void setup() = 0;
    virtual void execute() = 0;
    virtual void reinit_gwstate() = 0;
    virtual void release() = 0;

    void check();
    void commit_setup();
    template<class TState, class... Args>
    void commit_execute(int iteration, Args&&... args)
    {
        if (iteration - (_state ? _state->iteration() : 0) >= state_update_period() || iteration == iterations() || abort_flag() || _logging->state_save_flag())
        {
            check();
            set_state<TState>(iteration, std::forward<Args>(args)...);
        }
    }
    template<class TState, class... Args>
    void set_state(int iteration, Args&&... args)
    {
        if (!_tmp_state || dynamic_cast<TState*>(_tmp_state.get()) == nullptr)
            _tmp_state.reset(new TState());
        static_cast<TState*>(_tmp_state.get())->set(iteration, std::forward<Args>(args)...);
        _tmp_state.swap(_state);
        on_state();
    }
    template<class TState>
    void reset_state()
    {
        _state.reset(new TState());
        on_state();
    }
    void on_state();
    virtual void write_state();

protected:
    bool _error_check = false;
    arithmetic::GWState* _gwstate = nullptr;
    arithmetic::GWArithmetic* _gw = nullptr;
    std::unique_ptr<TaskState> _state;
    std::unique_ptr<TaskState> _tmp_state;
    int _iterations = 0;
    int _state_update_period = MULS_PER_STATE_UPDATE;
    File* _file;
    Logging* _logging;
    std::chrono::system_clock::time_point _last_write;
    std::chrono::system_clock::time_point _last_progress;
    static bool _abort_flag;
    int _restart_count = 0;
    int _restart_op = 0;
    double _op_count = 0;
    double _op_base = 0;
};

class InputTask : public Task
{
public:
    InputTask() : Task()
    {
    }

    void set_error_check(bool near, bool check);
    double timer() { return _timer; }

private:
    using Task::init;
protected:
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, TaskState* state, Logging* logging, int iterations);
    void reinit_gwstate() override;
    virtual void done();

protected:
    InputNum* _input = nullptr;
    double _timer = 0;
    bool _error_check_near = true;
    bool _error_check_forced = false;
};
