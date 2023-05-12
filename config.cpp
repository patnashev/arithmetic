
#include <cctype>
#include <cmath>
#include <stdlib.h>
#include <string.h>
#include "config.h"

using namespace arithmetic;

void ConfigObject::on_parsed()
{
    if (_parsed)
        _parsed->on_key();
}

void Config::parse_args(int argc, char *argv[])
{
    int cur = 1;
    while (cur < argc)
    {
        int last = cur;
        for (auto it = _objects.begin(); it != _objects.end(); it++)
            (*it)->parse_args(argc, argv, cur);
        if (cur == last)
        {
            _default_code(argv[cur]);
            cur++;
        }
    }
}

void Config::parse_ini(File& file)
{
    bool not_eof = true;
    std::unique_ptr<TextReader> reader(file.get_textreader());

    std::string st;
    std::string key;
    std::string value;
    std::string group = "";
    std::map<std::string, std::map<std::string, std::string>> values;
    values[""];

    while (not_eof)
    {
        while ((not_eof = reader->read_textline(st)) && (st.empty() || st[0] != '['))
        {
            if (st.empty() || st[0] == ';' || st[0] == '#')
                continue;
            size_t i = st.find('=');
            if (i == 0)
            {
                _default_code(st.data());
                continue;
            }
            if (i == std::string::npos)
                values[group][st] = "";
            else
                values[group][st.substr(0, i)] = st.substr(i + 1);
        }
        if (not_eof)
        {
            if (st[st.length() - 1] != ']' || st == "[]")
            {
                _default_code(st.data());
                continue;
            }
            group = "-" + st.substr(1, st.length() - 2);
            for (auto it = group.begin(); it != group.end(); it++)
                if (*it == '\\')
                    values[group.substr(0, it - group.begin())];
            values[group];
        }
    }

    ConfigGroup::parse_ini(values, "");
    for (auto& v : values[""])
        _default_code(v.first.data());
    values.erase("");

    for (auto& s : values)
        if (s.second.empty())
            printf("Unknown section %s.\n", s.first.data());
        else
            for (auto& v : s.second)
                printf("Unknown option %s\\%s.\n", s.first.data(), v.first.data());
}

#ifdef CHARCONV
#include <charconv>
template<class T>
const char* from_chars(const char* str, T& value)
{
    auto [ptr, ec] { std::from_chars(str, str + strlen(str), value) };
    if (ec != std::errc())
        return nullptr;
    return ptr;
}
#else
#include <cerrno>
template<typename  T, typename std::enable_if_t<std::is_integral<T>::value && std::is_unsigned<T>::value, bool> = true>
const char* from_chars(const char* str, T& value)
{
    char* str_end;
    unsigned long long res = std::strtoull(str, &str_end, 10);
    if (errno == ERANGE)
        return nullptr;
    if (res > std::numeric_limits<T>::max())
        return nullptr;
    value = (T)res;
    return str_end;
}
template<class T, typename std::enable_if_t<!std::is_floating_point<T>::value && std::is_signed<T>::value, bool> = true>
const char* from_chars(const char* str, T& value)
{
    char* str_end;
    long long res = std::strtoll(str, &str_end, 10);
    if (errno == ERANGE)
        return nullptr;
    if (res < std::numeric_limits<T>::min() || res > std::numeric_limits<T>::max())
        return nullptr;
    value = (T)res;
    return str_end;
}
template<typename  T, typename std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
const char* from_chars(const char* str, T& value)
{
    char* str_end;
    double res = std::strtod(str, &str_end);
    if (errno == ERANGE)
        return nullptr;
    if (res < std::numeric_limits<T>::lowest() || res > std::numeric_limits<T>::max())
        return nullptr;
    if (!std::isnormal(res))
        return nullptr;
    value = (T)res;
    return str_end;
}
#endif

template<class T>
bool ConfigObject::parse_number(const char* str, T& value_ref)
{
    T value;
    const char* ptr = from_chars(str, value);
    if (ptr == nullptr)
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
    if (*ptr != 0)
        return false;
    value_ref = value;
    return true;
}

template bool ConfigObject::parse_number(const char* str, char& value);
template bool ConfigObject::parse_number(const char* str, int8_t& value);
template bool ConfigObject::parse_number(const char* str, uint8_t& value);
template bool ConfigObject::parse_number(const char* str, int16_t& value);
template bool ConfigObject::parse_number(const char* str, uint16_t& value);
template bool ConfigObject::parse_number(const char* str, int32_t& value);
template bool ConfigObject::parse_number(const char* str, uint32_t& value);
template bool ConfigObject::parse_number(const char* str, int64_t& value);
template bool ConfigObject::parse_number(const char* str, uint64_t& value);
template bool ConfigObject::parse_number(const char* str, float& value);
template bool ConfigObject::parse_number(const char* str, double& value);
template bool ConfigObject::parse_number(const char* str, long double& value);


void ConfigKeyValue::parse_args(int argc, char *argv[], int& cur)
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
    else
    {
        if (*val_s == 0)
            return;
    }
    if (!set_value(val_s))
        return;
    cur += inc;
    on_parsed();
}

void ConfigKeyValue::parse_ini(std::map<std::string, std::map<std::string, std::string>>& values, const std::string& path)
{
    if (values.at(path).find(_name) == values.at(path).end())
        return;
    if (!set_value(values.at(path).at(_name).data()))
        return;
    values.at(path).erase(_name);
    on_parsed();
}

void ConfigKey::parse_args(int argc, char *argv[], int& cur)
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
    on_parsed();
}

void ConfigKey::parse_ini(std::map<std::string, std::map<std::string, std::string>>& values, const std::string& path)
{
    if (values.at(path).find(_name) == values.at(path).end())
        return;
    on_key();
    values.at(path).erase(_name);
    on_parsed();
}

void ConfigGroup::parse_args(int argc, char *argv[], int& cur)
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
        on_parsed();
        for (auto it = _objects.begin(); it != _objects.end(); it++)
        {
            ConfigExclusive* ex = dynamic_cast<ConfigExclusive*>(it->get());
            if (ex != nullptr)
                ex->reset();
        }
    }
    int last = cur - 1;
    while (last != cur)
    {
        last = cur;
        for (auto it = _objects.begin(); it != _objects.end(); it++)
            (*it)->parse_args(argc, argv, cur);
    }
}


void ConfigGroup::parse_ini(std::map<std::string, std::map<std::string, std::string>>& values, const std::string& path)
{
    std::string local_path = path;
    if (!_name.empty())
    {
        if (!local_path.empty())
            local_path += "\\";
        local_path += _name;
        if (values.find(local_path) == values.end())
            return;
        on_parsed();
        for (auto it = _objects.begin(); it != _objects.end(); it++)
        {
            ConfigExclusive* ex = dynamic_cast<ConfigExclusive*>(it->get());
            if (ex != nullptr)
                ex->reset();
        }
    }
    for (auto it = _objects.begin(); it != _objects.end(); it++)
        (*it)->parse_ini(values, local_path);
    if (!_name.empty() && values.at(local_path).empty())
        values.erase(local_path);
}

void ConfigKeyList::parse_args(int argc, char *argv[], int& cur)
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
        if (_fixed_size && k != (int)_objects.size() - 1)
            return;
        for (i = 0, j = 0, k = 0; ; i++)
            if (!val_s[i] || val_s[i] == _list_delim)
            {
                std::string val_str(val_s + j, i - j);
                if (!_objects[k%_objects.size()]->set_value(val_str.data()))
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
    if (cur != local_cur)
    {
        cur = local_cur;
        on_parsed();
    }
}

void ConfigKeyList::parse_ini(std::map<std::string, std::map<std::string, std::string>>& values, const std::string& path)
{
    if (values.at(path).find(_name) == values.at(path).end())
        return;
    const char* val_s = values.at(path).at(_name).data();
    int i, j, k;
    for (i = 0, k = 0; val_s[i]; i++)
        if (val_s[i] == _list_delim)
            k++;
    if (_fixed_size && k != (int)_objects.size() - 1)
        return;
    for (i = 0, j = 0, k = 0; ; i++)
        if (!val_s[i] || val_s[i] == _list_delim)
        {
            std::string val_str(val_s + j, i - j);
            if (!_objects[k%_objects.size()]->set_value(val_str.data()))
                return;
            k++;
            if (!val_s[i])
                break;
            j = i + 1;
        }
    values.at(path).erase(_name);
    on_parsed();
}

void ConfigExclusive::parse_args(int argc, char *argv[], int& cur)
{
    if (_case != nullptr)
    {
        _case->optional()->parse_args(argc, argv, cur);
        return;
    }

    int last = cur;
    for (auto it = _objects.begin(); it != _objects.end(); it++)
    {
        (*it)->parse_args(argc, argv, cur);
        if (_case != nullptr)
        {
            on_parsed();
            _case->optional()->parse_args(argc, argv, cur);
            return;
        }
    }
}

void ConfigExclusive::parse_ini(std::map<std::string, std::map<std::string, std::string>>& values, const std::string& path)
{
    for (auto it = _objects.begin(); it != _objects.end(); it++)
    {
        (*it)->parse_ini(values, path);
        if (_case != nullptr)
        {
            on_parsed();
            _case->optional()->parse_ini(values, path);
            return;
        }
    }
}

void ConfigExclusive::ConfigCase::parse_args(int argc, char *argv[], int& cur)
{
    std::vector<ConfigObject*> objects;
    for (auto it = _objects.begin(); it != _objects.end(); it++)
        objects.push_back(it->get());
    int local_cur = cur;
    int last = cur - 1;
    while (last != local_cur)
    {
        last = local_cur;
        for (auto it = objects.begin(); it != objects.end(); )
        {
            int prev = local_cur;
            (*it)->parse_args(argc, argv, local_cur);
            if (prev != local_cur)
                it = objects.erase(it);
            else
                it++;
        }
    }
    if (!objects.empty())
        return;
    if (cur != local_cur)
    {
        cur = local_cur;
        _exclusive->_case = this;
        on_parsed();
    }
}

void ConfigExclusive::ConfigCase::parse_ini(std::map<std::string, std::map<std::string, std::string>>& values, const std::string& path)
{
    bool parsed = true;
    std::map<std::string, std::string> copies;
    for (auto it = _objects.begin(); it != _objects.end(); it++)
    {
        ConfigGroup* group = dynamic_cast<ConfigGroup*>(it->get());
        if (group != nullptr)
        {
            if (values.find(path + "\\" + group->name()) == values.end())
            {
                parsed = false;
                break;
            }
        }
        else
        {
            if (values.at(path).find((*it)->name()) == values.at(path).end())
            {
                parsed = false;
                break;
            }
            copies[(*it)->name()] = values.at(path).at((*it)->name());
            (*it)->parse_ini(values, path);
            if (values.at(path).find((*it)->name()) != values.at(path).end())
            {
                parsed = false;
                break;
            }
        }
    }

    if (!parsed)
    {
        for (auto& c : copies)
            values.at(path)[c.first] = c.second;
        return;
    }
    for (auto it = _objects.begin(); it != _objects.end(); it++)
    {
        ConfigGroup* group = dynamic_cast<ConfigGroup*>(it->get());
        if (group != nullptr)
            group->parse_ini(values, path);
    }
    _exclusive->_case = this;
    on_parsed();
}
