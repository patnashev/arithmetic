
#include <stdlib.h>
#include "gwnum.h"
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
    logging->progress().update(0, (int)gwstate->handle.fft_count/2);
    _last_write = std::chrono::system_clock::now();
    _last_progress = std::chrono::system_clock::now();
}

void Task::run()
{
    ReliableGWArithmetic* reliable;
    std::unique_ptr<arithmetic::GWArithmetic> gw_setup(_error_check ? new ReliableGWArithmetic(*_gwstate) : new GWArithmetic(*_gwstate));
    std::unique_ptr<arithmetic::GWArithmetic> gw_execute(_error_check ? new ReliableGWArithmetic(*_gwstate) : new GWArithmetic(*_gwstate));

    _restart_op = 0;
    int restart_count = 0;
    while (true)
    {
        _gw = gw_setup.get();
        reliable = dynamic_cast<ReliableGWArithmetic*>(_gw);
        if (reliable)
            reliable->reset();
        int setup_restart_count = 0;
        bool setup_failed = false;
        while (true)
        {
            try
            {
                setup();
                break;
            }
            catch (const TaskRestartException&)
            {
                if (!reliable)
                {
                    release();
                    _error_check = true;
                    gw_setup.reset(new ReliableGWArithmetic(*_gwstate));
                    gw_execute.reset(new ReliableGWArithmetic(*_gwstate));
                    _gw = gw_setup.get();
                    reliable = dynamic_cast<ReliableGWArithmetic*>(_gw);
                }
            }
            catch (const ArithmeticException& e)
            {
                _logging->error("ArithmeticException: %s\n", e.what());
            }
            catch (...)
            {
                release();
                throw;
            }
            //GWASSERT(0);
            _logging->error("Arithmetic error, restarting at %.1f%%.\n", state() ? 100.0*state()->iteration()/iterations() : 0.0);
            if (reliable && reliable->restart_flag() && !reliable->failure_flag())
            {
                reliable->restart();
                continue;
            }
            if ((!reliable || !reliable->failure_flag()) && setup_restart_count < 5)
            {
                setup_restart_count++;
                continue;
            }
            setup_failed = true;
            break;
        }

        _gw = gw_execute.get();
        reliable = dynamic_cast<ReliableGWArithmetic*>(_gw);

        if (!setup_failed)
        {
            try
            {
                execute();
                break;
            }
            catch (const TaskRestartException&)
            {
                if (!reliable)
                {
                    release();
                    _error_check = true;
                    gw_setup.reset(new ReliableGWArithmetic(*_gwstate));
                    gw_execute.reset(new ReliableGWArithmetic(*_gwstate));
                }
            }
            catch (const ArithmeticException& e)
            {
                _logging->error("ArithmeticException: %s\n", e.what());
            }
            catch (...)
            {
                release();
                throw;
            }
            //GWASSERT(0);
            _logging->error("Arithmetic error, restarting at %.1f%%.\n", state() ? 100.0*state()->iteration()/iterations() : 0.0);
            if (reliable && reliable->restart_flag() && !reliable->failure_flag())
            {
                reliable->restart(_restart_op);
                continue;
            }
            if ((!reliable || !reliable->failure_flag()) && restart_count < 5)
            {
                restart_count++;
                continue;
            }
        }

        release();
        if (_gwstate->next_fft_count >= 5)
        {
            _logging->error("Maximum FFT increment reached.\n");
            throw TaskAbortException();
        }
        _gwstate->next_fft_count++;
        reinit_gwstate();
        restart_count = 0;
        _restart_op = 0;
        if (reliable)
            reliable->reset();
    }
    release();
    _tmp_state.reset();
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
    if (_file != nullptr && _state && (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - _last_write).count() >= DISK_WRITE_TIME || abort_flag()))
    {
        _logging->progress().update(_state->iteration()/(double)iterations(), (int)_gwstate->handle.fft_count/2);
        _file->write(*_state);
        _last_write = std::chrono::system_clock::now();
    }
    if (abort_flag())
        throw TaskAbortException();
    if (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - _last_progress).count() >= PROGRESS_TIME)
    {
        _logging->progress().update(_state ? _state->iteration()/(double)iterations() : 0.0, (int)_gwstate->handle.fft_count/2);
        _logging->report_progress();
        _last_progress = std::chrono::system_clock::now();
    }
    if (_error_check && _gw != nullptr)
    {
        ReliableGWArithmetic* reliable = dynamic_cast<ReliableGWArithmetic*>(_gw);
        _restart_op = reliable->op();
    }
}

void Task::commit_setup()
{
    check();
    if (_error_check)
    {
        ReliableGWArithmetic* reliable = dynamic_cast<ReliableGWArithmetic*>(_gw);
        _restart_op = reliable->op();
    }
}
