
#include <stdlib.h>
#include "gwnum.h"
#include "cpuid.h"
#include "task.h"
#include "exception.h"
#include "file.h"

using namespace arithmetic;

bool TaskState::read(Reader& reader)
{
    if (!reader.read(_iteration))
    {
        _iteration = 0;
        return false;
    }
    return true;
}

void TaskState::write(Writer& writer)
{
    writer.write(_iteration);
}

bool Task::_abort_flag = false;
int Task::MULS_PER_STATE_UPDATE = 20000;
int Task::DISK_WRITE_TIME = 300;
int Task::PROGRESS_TIME = 10;

void Task::init(arithmetic::GWState* gwstate, File* file, TaskState* state, Logging* logging, int iterations)
{
    _gwstate = gwstate;
    _iterations = iterations;
    _gw = nullptr;
    _file = file;
    _state.reset(state);
    _logging = logging;
    _op_count = 0;
    logging->progress().update(0, ops());
    _last_write = std::chrono::system_clock::now();
    _last_progress = std::chrono::system_clock::now();
}

void Task::run()
{
    ReliableGWArithmetic* reliable;
    std::unique_ptr<arithmetic::GWArithmetic> arithmetics[2];

    int setup_restart_count = 0;
    _restart_count = 0;
    _restart_op = 0;
    int* restart_count[2] { &setup_restart_count, &_restart_count };
    while (true)
    {
        _op_base = _gwstate->ops();
        int i;
        for (i = 0; i < 2; )
        {
            if (!arithmetics[i])
                arithmetics[i].reset(_error_check ? new ReliableGWArithmetic(*_gwstate) : new GWArithmetic(*_gwstate));
            _gw = arithmetics[i].get();
            reliable = dynamic_cast<ReliableGWArithmetic*>(_gw);
            try
            {
                if (i == 0)
                    setup();
                if (i == 1)
                    execute();
                *restart_count[i] = 0;
                i++;
                continue;
            }
            catch (const TaskRestartException&)
            {
                if (!reliable)
                {
                    release();
                    _gw = nullptr;
                    _op_count = _gwstate->ops() - _op_base;
                    _op_base = _gwstate->ops();
                    _error_check = true;
                    arithmetics[0].reset();
                    arithmetics[1].reset();
                }
            }
            catch (const ArithmeticException& e)
            {
                _logging->warning("ArithmeticException: %s\n", e.what());
            }
            catch (...)
            {
                release();
                throw;
            }
            //GWASSERT(0);
            _logging->warning("Arithmetic error, restarting at %.1f%%.\n", 100.0*progress());
            if (reliable && reliable->restart_flag() && !reliable->failure_flag())
            {
                std::string opstr = "pass " + std::to_string(i) + " suspicious ops:";
                for (auto it = reliable->suspect_ops().begin(); it != reliable->suspect_ops().end(); it++)
                    if (*it >= (i == 1 ? _restart_op : 0))
                        opstr += " " + std::to_string(*it);
                opstr += ".\n";
                _logging->debug(opstr.data());
                reliable->restart(i == 1 ? _restart_op : 0);
                i = 0;
                continue;
            }
            if ((!reliable || !reliable->failure_flag()) && *restart_count[i] < 2)
            {
                if (reliable)
                    reliable->restart(i == 1 ? _restart_op : 0);
                (*restart_count[i])++;
                i = 0;
                continue;
            }
            break;
        }
        if (_gw != nullptr)
            _op_count = _gwstate->ops() - _op_base;
        if (i == 2)
            break;

        if (_gw != nullptr)
            release();
        _gw = nullptr;
        if (_gwstate->next_fft_count >= 5)
        {
            _logging->error("Maximum FFT increment reached.\n");
            throw TaskAbortException();
        }
        if (*restart_count[i] >= 5)
        {
            _logging->error("Too many restarts, aborting.\n");
            throw TaskAbortException();
        }
        (*restart_count[i])++;

        _gwstate->next_fft_count++;
        reinit_gwstate();
        _restart_op = 0;
        arithmetics[0].reset();
        arithmetics[1].reset();
    }
    _tmp_state.reset();
    if (_gw != nullptr)
        release();
    _gw = nullptr;
}

void Task::check()
{
    if (_error_check)
    {
        ReliableGWArithmetic* reliable = dynamic_cast<ReliableGWArithmetic*>(_gw);
        if (reliable->restart_flag() || reliable->failure_flag())
            throw TaskRestartException();
    }
}

void Task::on_state()
{
    bool state_save_flag = _logging->state_save_flag();
    if (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - _last_write).count() >= DISK_WRITE_TIME || abort_flag() || state_save_flag)
    {
        _logging->debug("saving state to disk.\n");
        _logging->progress().update(progress(), ops());
        _logging->state_save();
        write_state();
        _last_write = std::chrono::system_clock::now();
        _logging->progress().update(progress(), ops());
        _logging->progress_save();
    }
    if (abort_flag())
        throw TaskAbortException();
    if (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - _last_progress).count() >= PROGRESS_TIME)
    {
        _logging->progress().update(progress(), ops());
        _logging->report_progress();
        _last_progress = std::chrono::system_clock::now();
    }
    if (_error_check && _gw != nullptr)
    {
        ReliableGWArithmetic* reliable = dynamic_cast<ReliableGWArithmetic*>(_gw);
        _restart_op = reliable->op();
    }
    _logging->heartbeat();
}

void Task::write_state()
{
    if (_file != nullptr)
    {
        if (_state && !_state->is_written())
            _file->write(*_state);
        if (!_state)
            _file->clear();
    }
}

void Task::commit_setup()
{
    check();
    if (_error_check)
    {
        ReliableGWArithmetic* reliable = dynamic_cast<ReliableGWArithmetic*>(_gw);
        reliable->reset();
    }
}

void InputTask::init(InputNum* input, arithmetic::GWState* gwstate, File* file, TaskState* state, Logging* logging, int iterations)
{
    Task::init(gwstate, file, state, logging, iterations);
    _input = input;
    _timer = getHighResTimer();
    _error_check = _error_check_near ? gwnear_fft_limit(gwstate->gwdata(), 1) == TRUE : _error_check_forced;
}

void InputTask::reinit_gwstate()
{
    _gwstate->done();
    _input->setup(*_gwstate);
    std::string prefix = _logging->prefix();
    _logging->set_prefix("");
    _logging->warning("Restarting using %s\n", _gwstate->fft_description.data());
    _logging->set_prefix(prefix);
    _logging->report_param("fft_desc", _gwstate->fft_description);
    _logging->report_param("fft_len", _gwstate->fft_length);
    _error_check = _error_check_near ? gwnear_fft_limit(_gwstate->gwdata(), 1) == TRUE : _error_check_forced;
}

void InputTask::done()
{
    _timer = (getHighResTimer() - _timer)/getHighResTimerFrequency();
    _logging->progress().update(1, ops());
}

void InputTask::set_error_check(bool near, bool check)
{
    _error_check_near = near;
    _error_check_forced = check;
    if (_gwstate != nullptr)
        _error_check = _error_check_near ? gwnear_fft_limit(_gwstate->gwdata(), 1) == TRUE : _error_check_forced;
}
