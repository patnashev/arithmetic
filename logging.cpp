
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <iostream>
#include "gwnum.h"
#include "cpuid.h"
#include "logging.h"

void Progress::update(double progress, int op_count)
{
    _cur_progress = progress;
    if (_timer == 0)
        _timer = getHighResTimer();
    double elapsed = (getHighResTimer() - _timer)/getHighResTimerFrequency();
    _timer = getHighResTimer();
    _time_total += elapsed;
    _time_stage = 0;
    if (_op_count < op_count)
        _time_op = elapsed/(op_count - _op_count);
    _op_count = op_count;
    Progress* cur = this;
    Progress* parent = _parent;
    while (parent != nullptr)
    {
        parent->_cur_progress = cur->progress_total();
        parent->_time_stage += elapsed;
        cur = parent;
        parent = cur->_parent;
    }
}

void Progress::time_init(double elapsed)
{
    _time_total = elapsed;
    _timer = getHighResTimer();
}

void Logging::debug(const char* message...)
{
    if (_level > LEVEL_DEBUG)
        return;
    char buf[1000];
    va_list args;
    va_start(args, message);
    vsnprintf(buf, 1000, message, args);
    va_end(args);
    report(buf, LEVEL_DEBUG);
}

void Logging::info(const char* message...)
{
    if (_level > LEVEL_INFO)
        return;
    char buf[1000];
    va_list args;
    va_start(args, message);
    vsnprintf(buf, 1000, message, args);
    va_end(args);
    report(buf, LEVEL_INFO);
}

void Logging::warning(const char* message...)
{
    if (_level > LEVEL_WARNING)
        return;
    char buf[1000];
    va_list args;
    va_start(args, message);
    vsnprintf(buf, 1000, message, args);
    va_end(args);
    report(buf, LEVEL_WARNING);
    result_save(buf);
}

void Logging::error(const char* message...)
{
    if (_level > LEVEL_ERROR)
        return;
    char buf[1000];
    va_list args;
    va_start(args, message);
    vsnprintf(buf, 1000, message, args);
    va_end(args);
    report(buf, LEVEL_ERROR);
    result_save(buf);
}

void Logging::result(bool success, const char* message...)
{
    if (_level > (success ? LEVEL_WARNING : LEVEL_INFO))
        return;
    char buf[1000];
    va_list args;
    va_start(args, message);
    vsnprintf(buf, 1000, message, args);
    va_end(args);
    report(buf, LEVEL_RESULT);
}

void Logging::report(const std::string& message, int level)
{
    if (_overwrite_line)
        std::cout << "\r                                                                      \r";
    if (_print_prefix)
        std::cout << _prefix;
    std::cout << message;
    if (!_file_log.empty() && level != LEVEL_PROGRESS)
    {
        FILE *fp = fopen(_file_log.data(), "a");
        if (fp)
        {
            if (_print_prefix)
                fwrite(_prefix.data(), 1, _prefix.length(), fp);
            fwrite(message.data(), 1, message.length(), fp);
            fclose(fp);
        }
    }
    _overwrite_line = (level == LEVEL_PROGRESS);
    _print_prefix = !message.empty() && (level == LEVEL_PROGRESS || message.back() == '\n');
}

void Logging::report_progress()
{
    if (_level > LEVEL_PROGRESS)
        return;
    char buf[100];
    if (progress().num_stages() > 1 && progress().progress_stage() > 0 && progress().progress_stage() < 1)
    {
        if (progress().time_op() > 0)
            snprintf(buf, 100, "%.1f%% stage / %.1f%% total, time per op: %.3f ms. ", progress().progress_stage()*100, progress().progress_total()*100, progress().time_op()*1000);
        else
            snprintf(buf, 100, "%.1f%% stage / %.1f%% total. ", progress().progress_stage()*100, progress().progress_total()*100);
        report(buf, LEVEL_PROGRESS);
    }
    else if (progress().num_stages() > 0)
    {
        if (progress().time_op() > 0)
            snprintf(buf, 100, "%.1f%% done, time per op: %.3f ms. ", progress().progress_total()*100, progress().time_op()*1000);
        else
            snprintf(buf, 100, "%.1f%% done. ", progress().progress_total()*100);
        report(buf, LEVEL_PROGRESS);
    }
}

void Logging::report_factor(InputNum& input, const arithmetic::Giant& f)
{
    if (_file_factor.empty())
        return;
    std::string str = f.to_string();
    result(true, "found factor %s\n", str.data());
    result_save(input.input_text() + " found factor " + str + ", time: " + std::to_string((int)progress().time_total()) + " s.\n");
    FILE *fp = fopen(_file_factor.data(), "a");
    if (fp)
    {
        fprintf(fp, "%s | %s\n", str.data(), input.input_text().data());
        fclose(fp);
    }
}

void Logging::result_save(const std::string& message)
{
    if (_file_result.empty())
        return;
    FILE *fp = fopen(_file_result.data(), "a");
    if (fp)
    {
        fwrite(_prefix.data(), 1, _prefix.length(), fp);
        fwrite(message.data(), 1, message.length(), fp);
        fclose(fp);
    }
}

void Logging::report_param(const std::string& name, double value)
{
    char buf[32];
    snprintf(buf, 32, "%.16g", value);
    progress().param(name) = buf;
}

void Logging::file_progress(File* file_progress)
{
    _file_progress = file_progress;
    if (file_progress == nullptr)
        return;
    progress().params().clear();
    std::string st;
    std::unique_ptr<TextReader> reader(file_progress->get_textreader());
    while (reader->read_textline(st))
    {
        auto i = st.find('=');
        if (i == std::string::npos)
            continue;
        progress().param(st.substr(0, i)) = st.substr(i + 1);
    }
    progress().time_init(progress().param_double("time_total"));
}

void Logging::progress_save()
{
    if (_file_progress == nullptr)
        return;
    report_param("time_total", progress().time_total());
    std::unique_ptr<Writer> writer(_file_progress->get_writer());
    for (auto& pair : progress().params())
        if (!pair.second.empty())
        {
            writer->write_text(pair.first);
            writer->write_text("=");
            writer->write_textline(pair.second);
        }
    _file_progress->commit_writer(*writer);
}
