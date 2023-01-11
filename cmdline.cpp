
#include <cctype>
#include <cmath>
#include <charconv>
#include <stdlib.h>
#include "cmdline.h"

using namespace arithmetic;

void CmdLine::parse(int argc, char *argv[])
{
    std::vector<CmdObject*> objects;
    for (auto it = _objects.begin(); it != _objects.end(); it++)
        objects.push_back(it->get());
    int cur = 1;
    while (cur < argc)
    {
        int last = cur;
        bool prev_flag = false;
        for (auto it = objects.begin(); it != objects.end(); )
        {
            CmdPrev* cmd_prev = dynamic_cast<CmdPrev*>(*it);
            if (cmd_prev != nullptr)
            {
                if (prev_flag)
                    cmd_prev->on_prev();
                it++;
            }
            else
            {
                int prev = cur;
                (*it)->parse(argc, argv, cur);
                prev_flag = prev != cur;
                if (prev_flag && dynamic_cast<CmdExclusive*>(*it))
                    it = objects.erase(it);
                else
                    it++;
            }
        }
        if (cur == last)
        {
            _default_code(argv[cur]);
            cur++;
        }
    }
}

template<class T>
bool CmdObject::parse_number(const char* str, T& value)
{
    auto [ptr, ec] { std::from_chars(str, str + strlen(str), value) };
    if (ec != std::errc())
        return false;
    if (*ptr && value != 0)
    {
        T max = std::numeric_limits<T>::max()/(value < 0 ? -value : value);
        switch (*ptr)
        {
        case 'k':
        case 'K':
            if (max < 1000)
                return false;
            value = (T)(value*1000);
            ptr++;
            break;
        case 'M':
            if (max < 1000000)
                return false;
            value = (T)(value*1000000);
            ptr++;
            break;
        case 'G':
            if (max < 1000000000)
                return false;
            value = (T)(value*1000000000);
            ptr++;
            break;
        case 'T':
            if (max < 1000000000000ULL)
                return false;
            value = (T)(value*1000000000000ULL);
            ptr++;
            break;
        case 'P':
            if (max < 1000000000000000ULL)
                return false;
            value = (T)(value*1000000000000000ULL);
            ptr++;
            break;
        }
    }
    while (*ptr && std::isspace(*ptr))
        ptr++;
    return *ptr == 0;
}

template bool CmdObject::parse_number(const char* str, char& value);
template bool CmdObject::parse_number(const char* str, int8_t& value);
template bool CmdObject::parse_number(const char* str, uint8_t& value);
template bool CmdObject::parse_number(const char* str, int16_t& value);
template bool CmdObject::parse_number(const char* str, uint16_t& value);
template bool CmdObject::parse_number(const char* str, int32_t& value);
template bool CmdObject::parse_number(const char* str, uint32_t& value);
template bool CmdObject::parse_number(const char* str, int64_t& value);
template bool CmdObject::parse_number(const char* str, uint64_t& value);
template bool CmdObject::parse_number(const char* str, float& value);
template bool CmdObject::parse_number(const char* str, double& value);
template bool CmdObject::parse_number(const char* str, long double& value);


void CmdKeyValue::parse(int argc, char *argv[], int& cur)
{
    if (cur >= argc)
        return;
    if (strncmp(argv[cur], _name.data(), _name.size()) != 0)
        return;
    char* val_s = argv[cur] + _name.size();
    int inc = 1;
    if (_delim == ' ')
    {
        if (*val_s != 0)
            return;
        if (cur + 1 >= argc)
            return;
        val_s = argv[cur + 1];
        inc = 2;
    }
    else if (_delim != 0)
    {
        if (*val_s != _delim)
            return;
        val_s++;
    }
    if (!set_value(val_s))
        return;
    cur += inc;
}

void CmdKey::parse(int argc, char *argv[], int& cur)
{
    if (cur >= argc)
        return;
    if (strncmp(argv[cur], _name.data(), _name.size()) != 0)
        return;
    char* val_s = argv[cur] + _name.size();
    if (*val_s != 0)
        return;
    on_key();
    cur++;
}

void CmdGroup::parse(int argc, char *argv[], int& cur)
{
    if (!_name.empty())
    {
        if (cur >= argc)
            return;
        if (strncmp(argv[cur], _name.data(), _name.size()) != 0)
            return;
        char* val_s = argv[cur] + _name.size();
        if (*val_s != 0)
            return;
        cur++;
    }
    std::vector<CmdObject*> objects;
    for (auto it = _objects.begin(); it != _objects.end(); it++)
        objects.push_back(it->get());
    int last = cur - 1;
    while (last != cur)
    {
        last = cur;
        bool prev_flag = false;
        for (auto it = objects.begin(); it != objects.end(); )
        {
            CmdPrev* cmd_prev = dynamic_cast<CmdPrev*>(*it);
            if (cmd_prev != nullptr)
            {
                if (prev_flag)
                    cmd_prev->on_prev();
                it++;
            }
            else
            {
                int prev = cur;
                (*it)->parse(argc, argv, cur);
                prev_flag = prev != cur;
                if (prev_flag && dynamic_cast<CmdExclusive*>(*it))
                    it = objects.erase(it);
                else
                    it++;
            }
        }
    }
}

void CmdKeyList::parse(int argc, char *argv[], int& cur)
{
    if (cur >= argc)
        return;
    if (strncmp(argv[cur], _name.data(), _name.size()) != 0)
        return;
    char* val_s = argv[cur] + _name.size();
    int local_cur = cur + 1;
    if (_delim == ' ')
    {
        if (*val_s != 0)
            return;
        if (local_cur >= argc)
            return;
        val_s = argv[local_cur];
        local_cur++;
    }
    else if (_delim != 0)
    {
        if (*val_s != _delim)
            return;
        val_s++;
    }
    if (_list_delim != ' ')
    {
        int i, j, k;
        for (i = 0, k = 0; val_s[i]; i++)
            if (val_s[i] == _list_delim)
                k++;
        if (k != (int)_objects.size() - 1)
            return;
        for (i = 0, j = 0, k = 0; ; i++)
            if (!val_s[i] || val_s[i] == _list_delim)
            {
                std::string val_str(val_s + j, i - j);
                if (!_objects[k]->set_value(val_str.data()))
                    return;
                k++;
                if (!val_s[i])
                    break;
                j = i + 1;
            }
    }
    else
    {
        int k = 0;
        while (true)
        {
            if (!_objects[k]->set_value(val_s))
                return;
            k++;
            if (k >= _objects.size())
                break;
            if (local_cur >= argc)
                return;
            val_s = argv[local_cur];
            local_cur++;
        }
    }
    cur = local_cur;
}

void CmdExclusive::parse(int argc, char *argv[], int& cur)
{
    int last = cur;
    for (auto it = _objects.begin(); it != _objects.end(); it++)
    {
        (*it)->parse(argc, argv, cur);
        if (cur != last)
            break;
    }
}

void CmdExclusive::CmdCase::parse(int argc, char *argv[], int& cur)
{
    std::vector<CmdObject*> objects;
    for (auto it = _objects.begin(); it != _objects.end(); it++)
        objects.push_back(it->get());
    int local_cur = cur;
    int last = cur - 1;
    while (last != local_cur)
    {
        last = local_cur;
        bool prev_flag = false;
        for (auto it = objects.begin(); it != objects.end(); )
        {
            CmdPrev* cmd_prev = dynamic_cast<CmdPrev*>(*it);
            if (cmd_prev != nullptr)
            {
                if (prev_flag)
                {
                    cmd_prev->on_prev();
                    it = objects.erase(it);
                }
                else
                    it++;
            }
            else
            {
                int prev = local_cur;
                (*it)->parse(argc, argv, local_cur);
                prev_flag = prev != local_cur;
                if (prev_flag)
                    it = objects.erase(it);
                else
                    it++;
            }
        }
    }
    if (!objects.empty())
        return;
    if (_optional != nullptr)
        ((CmdObject*)_optional.get())->parse(argc, argv, local_cur);
    cur = local_cur;
}
