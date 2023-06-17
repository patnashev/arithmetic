#pragma once

#include <memory>
#include <vector>
#include <string>
#include <map>
#include <optional>
#include <functional>
#include <climits>
#include "giant.h"
#include "file.h"

template<class V>
class Enum
{
public:
    Enum& add(const std::string& key, const V& value) { values[key] = value; return *this; }

public:
    std::map<std::string, V> values;
};

class ConfigKey;

class ConfigObject
{
public:
    ConfigObject(const std::string& name) : _name(name) { }
    virtual ~ConfigObject() { }

    const std::string& name() const { return _name; }

    template<class T>
    static bool parse_number(const char* str, T& value);
    void set_parsed(ConfigKey* parsed) { _parsed.reset(parsed); }

    virtual void parse_args(int argc, char *argv[], int& cur) = 0;
    virtual void parse_ini(std::map<std::string, std::map<std::string, std::string>>& values, const std::string& path) = 0;

protected:
    virtual void on_parsed();

protected:
    std::string _name;
    std::unique_ptr<ConfigKey> _parsed;
};

class ConfigKeyValue : public ConfigObject
{
public:
    ConfigKeyValue(const std::string& name, char delim) : ConfigObject(name), _delim(delim) { }

    void parse_args(int argc, char *argv[], int& cur) override;
    void parse_ini(std::map<std::string, std::map<std::string, std::string>>& values, const std::string& path) override;

protected:
    virtual bool set_value(const char* value_s) = 0;

protected:
    char _delim;
};

class ConfigKeyValueCode : public ConfigKeyValue
{
public:
    ConfigKeyValueCode(const std::string& name, char delim, std::function<bool(const char*)> func) : ConfigKeyValue(name, delim), _func(func) { }

protected:
    bool set_value(const char* value_s) { return _func(value_s); }

protected:
    std::function<bool(const char*)> _func;
};

template<class V>
class ConfigKeyValueString : public ConfigKeyValue
{
public:
    ConfigKeyValueString(const std::string& name, char delim, V& value) : ConfigKeyValue(name, delim), _value(value) { }

protected:
    bool set_value(const char* value_s) { _value = value_s; return true; }

protected:
    V& _value;
};

template<class T, class V>
class ConfigKeyValueNumber : public ConfigKeyValue
{
public:
    ConfigKeyValueNumber(const std::string& name, char delim, T& value, V min, V max) : ConfigKeyValue(name, delim), _value(value), _min(min), _max(max) { }

protected:
    bool set_value(const char* value_s) override { V val; if (!parse_number(value_s, val) || val < _min || val > _max) return false; _value = val; return true; }

protected:
    T& _value;
    V _min;
    V _max;
};

template<class T, class V>
class ConfigKeyValueEnum : public ConfigKeyValue
{
public:
    ConfigKeyValueEnum(const std::string& name, char delim, T& value, const Enum<V>& enm) : ConfigKeyValue(name, delim), _value(value), _enum(enm) { }

protected:
    bool set_value(const char* value_s) { auto it = _enum.values.find(value_s); if (it == _enum.values.end()) return false; _value = it->second; return true; }

protected:
    T& _value;
    const Enum<V>& _enum;
};


class ConfigKey : public ConfigObject
{
    friend class ConfigObject;
public:
    ConfigKey(const std::string& name) : ConfigObject(name) { }

    void parse_args(int argc, char *argv[], int& cur) override;
    void parse_ini(std::map<std::string, std::map<std::string, std::string>>& values, const std::string& path) override;

protected:
    virtual void on_key() { }
};


template<class T, class V>
class ConfigKeyCheck : public ConfigKey
{
public:
    ConfigKeyCheck(const std::string& name, T& dst, V value) : ConfigKey(name), _dst(dst), _value(value) { }

protected:
    void on_key() override { _dst = _value; }

protected:
    T& _dst;
    V _value;
};

class ConfigKeyCode : public ConfigKey
{
public:
    ConfigKeyCode(const std::string& name, std::function<void()> func) : ConfigKey(name), _func(func) { }

protected:
    void on_key() override { _func(); }

protected:
    std::function<void()> _func;
};

class ConfigKeyList : public ConfigObject
{
public:
    class ConfigListValue
    {
    public:
        virtual ~ConfigListValue() { }
        virtual bool set_value(const char* value_s) = 0;
    };
    template<class V>
    class ConfigListValueString : public ConfigListValue
    {
    public:
        ConfigListValueString(V& value) : _value(value) { }

        bool set_value(const char* value_s) override { _value = value_s; return true; }

    protected:
        V& _value;
    };
    template<class T, class V>
    class ConfigListValueNumber : public ConfigListValue
    {
    public:
        ConfigListValueNumber(T& value, V min, V max) : _value(value), _min(min), _max(max) { }

        bool set_value(const char* value_s) override { V val; if (!parse_number(value_s, val) || val < _min || val >= _max) return false; _value = val; return true; }

    protected:
        T& _value;
        V _min;
        V _max;
    };
    class ConfigListValueCode : public ConfigListValue
    {
    public:
        ConfigListValueCode(std::function<bool(const char*)> func) : _func(func) { }

        bool set_value(const char* value_s) { return _func(value_s); }

    protected:
        std::function<bool(const char*)> _func;
    };

public:
    ConfigKeyList(const std::string& name, char delim, char list_delim, bool fixed_size) : ConfigObject(name), _delim(delim), _list_delim(list_delim), _fixed_size(fixed_size){ }

    template<class T>
    T* add(T* object) { _objects.emplace_back(object);  return object; }

    void parse_args(int argc, char *argv[], int& cur) override;
    void parse_ini(std::map<std::string, std::map<std::string, std::string>>& values, const std::string& path) override;

protected:
    char _delim;
    char _list_delim;
    bool _fixed_size;
    std::vector<std::unique_ptr<ConfigListValue>> _objects;
};

class ConfigGroup : public ConfigObject
{
public:
    ConfigGroup() : ConfigObject("") { }
    ConfigGroup(const std::string& name) : ConfigObject(name) { }

    template<class T>
    T* add(T* object) { _objects.emplace_back(object);  return object; }
    ConfigObject* last() { return _objects[_objects.size() - 1].get(); }

    void parse_args(int argc, char *argv[], int& cur) override;
    void parse_ini(std::map<std::string, std::map<std::string, std::string>>& values, const std::string& path) override;

protected:
    std::vector<std::unique_ptr<ConfigObject>> _objects;
};

class ConfigExclusive : public ConfigObject
{
public:
    class ConfigCase : public ConfigGroup
    {
        friend class ConfigExclusive;
    public:
        ConfigCase(ConfigExclusive* exclusive) : ConfigGroup(), _exclusive(exclusive) { }

        void parse_args(int argc, char *argv[], int& cur) override;
        void parse_ini(std::map<std::string, std::map<std::string, std::string>>& values, const std::string& path) override;

        ConfigGroup* optional() { return &_optional; }

    protected:
        ConfigExclusive* _exclusive;
        ConfigGroup _optional;
    };

public:
    ConfigExclusive() : ConfigObject("") { }

    ConfigCase* add(ConfigCase* object) { _objects.emplace_back(object);  return object; }
    void reset() { _case = nullptr; }

    void parse_args(int argc, char *argv[], int& cur) override;
    void parse_ini(std::map<std::string, std::map<std::string, std::string>>& values, const std::string& path) override;

protected:
    std::vector<std::unique_ptr<ConfigCase>> _objects;
    ConfigCase* _case = nullptr;
};



template<class R>
class ConfigKeyListSetup
{
public:
    ConfigKeyListSetup(ConfigKeyList* list, R* parent) : _list(list), _parent(parent) { }

    ConfigKeyListSetup& value_code(std::function<bool(const char*)> func) { _list->add(new ConfigKeyList::ConfigListValueCode(func)); return *this; }
    ConfigKeyListSetup& value_string(std::string& value) { _list->add(new ConfigKeyList::ConfigListValueString(value)); return *this; }
    ConfigKeyListSetup& value_string(std::optional<std::string>& value) { _list->add(new ConfigKeyList::ConfigListValueString(value)); return *this; }
    template<class T>
    ConfigKeyListSetup& value_number(T& value, T min = std::numeric_limits<T>::lowest(), T max = std::numeric_limits<T>::max()) { _list->add(new ConfigKeyList::ConfigListValueNumber<T, T>(value, min, max)); return *this; }
    template<class T>
    ConfigKeyListSetup& value_number(std::optional<T>& value, T min = std::numeric_limits<T>::lowest(), T max = std::numeric_limits<T>::max()) { _list->add(new ConfigKeyList::ConfigListValueNumber<std::optional<T>, T>(value, min, max)); return *this; }

    R& end() { return (R&)*_parent; }

protected:
    ConfigKeyList* _list;
    R* _parent;
};

template<class R>
class ConfigKeyGroupSetup;
template<class R>
class ConfigExclusiveSetup;

template<class P, class R>
class ConfigGroupSetup
{
public:
    ConfigGroupSetup(ConfigGroup* group, R* parent) : _group(group), _parent(parent) { }

    P& ignore(const std::string& key) { _group->add(new ConfigKey(key)); return (P&)*this; }
    P& value_code(const std::string& key, char delim, std::function<bool(const char*)> func) { _group->add(new ConfigKeyValueCode(key, delim, func)); return (P&)*this; }
    P& value_string(const std::string& key, char delim, std::string& value) { _group->add(new ConfigKeyValueString<std::string>(key, delim, value)); return (P&)*this; }
    P& value_string(const std::string& key, char delim, std::optional<std::string>& value) { _group->add(new ConfigKeyValueString<std::optional<std::string>>(key, delim, value)); return (P&)*this; }
    template<class T>
    P& value_number(const std::string& key, char delim, T& value, T min = std::numeric_limits<T>::lowest(), T max = std::numeric_limits<T>::max()) { _group->add(new ConfigKeyValueNumber<T, T>(key, delim, value, min, max)); return (P&)*this; }
    template<class T>
    P& value_number(const std::string& key, char delim, std::optional<T>& value, T min = std::numeric_limits<T>::lowest(), T max = std::numeric_limits<T>::max()) { _group->add(new ConfigKeyValueNumber<std::optional<T>, T>(key, delim, value, min, max)); return (P&)*this; }
    template<class T, class V>
    P& value_enum(const std::string& key, char delim, T& value, const Enum<V>& enm) { _group->add(new ConfigKeyValueEnum<T, V>(key, delim, value, enm)); return (P&)*this; }
    template<class T, class V>
    P& check(const std::string& key, T& dst, V value) { _group->add(new ConfigKeyCheck<T, V>(key, dst, value)); return (P&)*this; }
    P& check_code(const std::string& key, std::function<void()> func) { _group->add(new ConfigKeyCode(key, func)); return (P&)*this; }
    template<class T, class V>
    P& on_check(T& dst, V value) { _group->last()->set_parsed(new ConfigKeyCheck<T, V>("", dst, value)); return (P&)*this; }
    P& on_code(std::function<void()> func) { _group->last()->set_parsed(new ConfigKeyCode("", func)); return (P&)*this; }
    ConfigKeyGroupSetup<P> group(const std::string& key) { return ConfigKeyGroupSetup<P>(_group->add(new ConfigGroup(key)), (P*)this); }
    ConfigKeyListSetup<P> list(const std::string& key, char delim, char list_delim, bool fixed_size = true) { return ConfigKeyListSetup<P>(_group->add(new ConfigKeyList(key, delim, list_delim, fixed_size)), (P*)this); }
    ConfigExclusiveSetup<P> exclusive() { return ConfigExclusiveSetup<P>(_group->add(new ConfigExclusive()), (P*)this); }

    R& end() { return (R&)*_parent; }

protected:
    ConfigGroup* _group;
    R* _parent;
};

template<class R>
class ConfigKeyGroupSetup : public ConfigGroupSetup<ConfigKeyGroupSetup<R>, R>
{
public:
    ConfigKeyGroupSetup(ConfigGroup* group, R* parent) : ConfigGroupSetup<ConfigKeyGroupSetup<R>, R>(group, parent) { }
};

template<class R>
class ConfigExclusiveSetup : public ConfigExclusive
{
public:
    class ConfigOptionalSetup : public ConfigGroupSetup<ConfigOptionalSetup, ConfigExclusiveSetup>
    {
    public:
        ConfigOptionalSetup(ConfigGroup* optional, ConfigExclusiveSetup* parent) : ConfigGroupSetup<ConfigOptionalSetup, ConfigExclusiveSetup>(optional, parent) { }
    };
    class ConfigCaseSetup : public ConfigGroupSetup<ConfigCaseSetup, ConfigExclusiveSetup>
    {
    public:
        ConfigCaseSetup(ConfigCase* ex_case, ConfigExclusiveSetup* parent) : ConfigGroupSetup<ConfigCaseSetup, ConfigExclusiveSetup>(ex_case, parent) { }

        ConfigOptionalSetup optional() { return ConfigOptionalSetup(((ConfigCase*)this->_group)->optional(), this->_parent); }
    };

public:
    ConfigExclusiveSetup(ConfigExclusive* ex, R* parent) : _ex(ex), _parent(parent) { }

    ConfigCaseSetup ex_case() { return ConfigCaseSetup(_ex->add(new ConfigCase(_ex)), this); }

    R& end() { return (R&)*_parent; }

protected:
    ConfigExclusive* _ex;
    R* _parent;
};

class Config : public ConfigGroup, public ConfigGroupSetup<Config, std::nullptr_t>
{
public:
    Config() : ConfigGroupSetup(this, nullptr) { }

    Config& default_code(std::function<void(const char*)> func) { _default_code = func; return *this; }
    void parse_args(int argc, char *argv[]);
    void parse_ini(File& file);

private:
    using ConfigObject::parse_args;
    using ConfigObject::parse_ini;
    std::function<void(const char*)> _default_code;
};
