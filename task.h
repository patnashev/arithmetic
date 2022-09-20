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
    TaskState(char type) : _type(type), _iteration(0) { }
    virtual ~TaskState() { }
    
    void set(int iteration) { _iteration = iteration; }
    virtual bool read(Reader& reader);
    virtual void write(Writer& writer);

    int iteration() { return _iteration; }
    char type() { return _type; }
    char version() { return 0; }

protected:
    char _type;
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
    Task(bool error_check) : _error_check(error_check){ }
    virtual ~Task() { }

    virtual void run();
    static void abort() { _abort_flag = true; }
    static void abort_reset() { _abort_flag = false; }

    arithmetic::GWArithmetic& gw() { return *_gw; }
    TaskState* state() { return _state.get(); }
    int iterations() { return _iterations; }
    int state_update_period() { return _state_update_period; }
    bool is_last(int iteration) { return iteration + 1 == _iterations || (iteration + 1)%_state_update_period == 0; }
    static bool abort_flag() { return _abort_flag; }

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
        if (iteration%state_update_period() == 0 || iteration == iterations() || abort_flag())
        {
            check();
            set_state<TState>(iteration, std::forward<Args>(args)...);
        }
    }
    template<class TState, class... Args>
    void set_state(int iteration, Args&&... args)
    {
        if (!_tmp_state)
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

protected:
    bool _error_check;
    arithmetic::GWState* _gwstate = nullptr;
    arithmetic::GWArithmetic* _gw = nullptr;
    std::unique_ptr<TaskState> _state;
    int _iterations = 0;
    int _state_update_period = MULS_PER_STATE_UPDATE;
    File* _file;
    Logging* _logging;
    std::chrono::system_clock::time_point _last_write;
    std::chrono::system_clock::time_point _last_progress;
    static bool _abort_flag;
private:
    int _restart_op = 0;
    std::unique_ptr<TaskState> _tmp_state;
};
