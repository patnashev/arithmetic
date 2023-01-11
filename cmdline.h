#pragma once

#include <memory>
#include <vector>
#include <string>
#include <map>
#include <optional>
#include <functional>
#include "giant.h"

template<class V>
class Enum
{
public:
    Enum& add(const std::string& key, const V& value) { values[key] = value; return *this; }

public:
    std::map<std::string, V> values;
};

class CmdGroup;
class CmdLine;
class CmdExclusive;

class CmdObject
{
    friend class CmdGroup;
    friend class CmdLine;
    friend class CmdExclusive;
public:
    CmdObject(const std::string& name) : _name(name) { }
    virtual ~CmdObject() { }

    template<class T>
    static bool parse_number(const char* str, T& value);

protected:
    virtual void parse(int argc, char *argv[], int& cur) = 0;

protected:
    std::string _name;
};

class CmdKeyValue : public CmdObject
{
public:
    CmdKeyValue(const std::string& name, char delim) : CmdObject(name), _delim(delim) { }

protected:
    void parse(int argc, char *argv[], int& cur) override;
    virtual bool set_value(const char* value_s) = 0;

protected:
    char _delim;
};

class CmdKeyValueCode : public CmdKeyValue
{
public:
    CmdKeyValueCode(const std::string& name, char delim, std::function<bool(const char*)> func) : CmdKeyValue(name, delim), _func(func) { }

protected:
    bool set_value(const char* value_s) { return _func(value_s); }

protected:
    std::function<bool(const char*)> _func;
};

template<class V>
class CmdKeyValueString : public CmdKeyValue
{
public:
    CmdKeyValueString(const std::string& name, char delim, V& value) : CmdKeyValue(name, delim), _value(value) { }

protected:
    bool set_value(const char* value_s) { _value = value_s; return true; }

protected:
    V& _value;
};

template<class T, class V>
class CmdKeyValueNumber : public CmdKeyValue
{
public:
    CmdKeyValueNumber(const std::string& name, char delim, T& value, V min, V max) : CmdKeyValue(name, delim), _value(value), _min(min), _max(max) { }

protected:
    bool set_value(const char* value_s) override { V val; if (!parse_number(value_s, val) || val < _min || val > _max) return false; _value = val; return true; }

protected:
    T& _value;
    V _min;
    V _max;
};

template<class T, class V>
class CmdKeyValueEnum : public CmdKeyValue
{
public:
    CmdKeyValueEnum(const std::string& name, char delim, T& value, const Enum<V>& enm) : CmdKeyValue(name, delim), _value(value), _enum(enm) { }

protected:
    bool set_value(const char* value_s) { auto it = _enum.values.find(value_s); if (it == _enum.values.end()) return false; _value = it->second; return true; }

protected:
    V& _value;
    const Enum<V>& _enum;
};


class CmdKey : public CmdObject
{
public:
    CmdKey(const std::string& name) : CmdObject(name) { }

protected:
    void parse(int argc, char *argv[], int& cur) override;
    virtual void on_key() { }
};


template<class T, class V>
class CmdKeyCheck : public CmdKey
{
public:
    CmdKeyCheck(const std::string& name, T& dst, V value) : CmdKey(name), _dst(dst), _value(value) { }

protected:
    void on_key() override { _dst = _value; }

protected:
    T& _dst;
    V _value;
};

class CmdKeyCode : public CmdKey
{
public:
    CmdKeyCode(const std::string& name, std::function<void()> func) : CmdKey(name), _func(func) { }

protected:
    void on_key() override { _func(); }

protected:
    std::function<void()> _func;
};

class CmdPrev : public CmdObject
{
public:
    CmdPrev() : CmdObject("") { }

    void parse(int argc, char *argv[], int& cur) override { }
    virtual void on_prev() = 0;
};

template<class T, class V>
class CmdPrevCheck : public CmdPrev
{
public:
    CmdPrevCheck(T& dst, V value) : _dst(dst), _value(value) { }

    void on_prev() override { _dst = _value; }

protected:
    T& _dst;
    V _value;
};

class CmdPrevCode : public CmdPrev
{
public:
    CmdPrevCode(std::function<void()> func) : _func(func) { }

    void on_prev() override { _func(); }

protected:
    std::function<void()> _func;
};

class CmdKeyList : public CmdObject
{
public:
    class CmdListValue
    {
    public:
        virtual bool set_value(const char* value_s) = 0;
    };
    template<class V>
    class CmdListValueString : public CmdListValue
    {
    public:
        CmdListValueString(V& value) : _value(value) { }

        bool set_value(const char* value_s) override { _value = value_s; return true; }

    protected:
        V& _value;
    };
    template<class T, class V>
    class CmdListValueNumber : public CmdListValue
    {
    public:
        CmdListValueNumber(T& value, V min, V max) : _value(value), _min(min), _max(max) { }

        bool set_value(const char* value_s) override { V val; if (!parse_number(value_s, val) || val < _min || val >= _max) return false; _value = val; return true; }

    protected:
        T& _value;
        V _min;
        V _max;
    };
    class CmdListValueCode : public CmdListValue
    {
    public:
        CmdListValueCode(std::function<bool(const char*)> func) : _func(func) { }

        bool set_value(const char* value_s) { return _func(value_s); }

    protected:
        std::function<bool(const char*)> _func;
    };

public:
    CmdKeyList(const std::string& name, char delim, char list_delim) : CmdObject(name), _delim(delim), _list_delim(list_delim) { }

    template<class T>
    T* add(T* object) { _objects.emplace_back(object);  return object; }

protected:
    void parse(int argc, char *argv[], int& cur) override;

protected:
    char _delim;
    char _list_delim;
    std::vector<std::unique_ptr<CmdListValue>> _objects;
};

class CmdGroup : public CmdObject
{
public:
    CmdGroup() : CmdObject("") { }
    CmdGroup(const std::string& name) : CmdObject(name) { }

    template<class T>
    T* add(T* object) { _objects.emplace_back(object);  return object; }

protected:
    void parse(int argc, char *argv[], int& cur) override;

protected:
    std::vector<std::unique_ptr<CmdObject>> _objects;
};

class CmdExclusive : public CmdObject
{
public:
    class CmdCase : public CmdGroup
    {
        friend class CmdExclusive;
    protected:
        void parse(int argc, char *argv[], int& cur) override;

    public:
        CmdGroup* optional() { if (!_optional) _optional.reset(new CmdGroup()); return _optional.get(); }

    protected:
        std::unique_ptr<CmdGroup> _optional;
    };

public:
    CmdExclusive() : CmdObject("") { }

    CmdCase* add(CmdCase* object) { _objects.emplace_back(object);  return object; }

protected:
    void parse(int argc, char *argv[], int& cur) override;

protected:
    std::vector<std::unique_ptr<CmdCase>> _objects;
};



template<class R>
class CmdKeyListConf
{
public:
    CmdKeyListConf(CmdKeyList* list, R* parent) : _list(list), _parent(parent) { }

    CmdKeyListConf& value_code(std::function<bool(const char*)> func) { _list->add(new CmdKeyList::CmdListValueCode(func)); return *this; }
    CmdKeyListConf& value_string(std::string& value) { _list->add(new CmdKeyList::CmdListValueString(value)); return *this; }
    CmdKeyListConf& value_string(std::optional<std::string>& value) { _list->add(new CmdKeyList::CmdListValueString(value)); return *this; }
    template<class T>
    CmdKeyListConf& value_number(T& value, T min = std::numeric_limits<T>::lowest(), T max = std::numeric_limits<T>::max()) { _list->add(new CmdKeyList::CmdListValueNumber(value, min, max)); return *this; }
    template<class T>
    CmdKeyListConf& value_number(std::optional<T>& value, T min = std::numeric_limits<T>::lowest(), T max = std::numeric_limits<T>::max()) { _list->add(new CmdKeyList::CmdListValueNumber(value, min, max)); return *this; }

    R& end() { return (R&)*_parent; }

protected:
    CmdKeyList* _list;
    R* _parent;
};

template<class R>
class CmdKeyGroupConf;
template<class R>
class CmdExclusiveConf;

template<class P, class R>
class CmdGroupConf
{
public:
    CmdGroupConf(CmdGroup* group, R* parent) : _group(group), _parent(parent) { }

    P& ignore(const std::string& key) { _group->add(new CmdKey(key)); return (P&)*this; }
    P& value_code(const std::string& key, char delim, std::function<bool(const char*)> func) { _group->add(new CmdKeyValueCode(key, delim, func)); return (P&)*this; }
    P& value_string(const std::string& key, char delim, std::string& value) { _group->add(new CmdKeyValueString(key, delim, value)); return (P&)*this; }
    P& value_string(const std::string& key, char delim, std::optional<std::string>& value) { _group->add(new CmdKeyValueString(key, delim, value)); return (P&)*this; }
    template<class T>
    P& value_number(const std::string& key, char delim, T& value, T min = std::numeric_limits<T>::lowest(), T max = std::numeric_limits<T>::max()) { _group->add(new CmdKeyValueNumber(key, delim, value, min, max)); return (P&)*this; }
    template<class T>
    P& value_number(const std::string& key, char delim, std::optional<T>& value, T min = std::numeric_limits<T>::lowest(), T max = std::numeric_limits<T>::max()) { _group->add(new CmdKeyValueNumber(key, delim, value, min, max)); return (P&)*this; }
    template<class T, class V>
    P& value_enum(const std::string& key, char delim, T& value, const Enum<V>& enm) { _group->add(new CmdKeyValueEnum(key, delim, value, enm)); return (P&)*this; }
    template<class T, class V>
    P& check(const std::string& key, T& dst, V value) { _group->add(new CmdKeyCheck(key, dst, value)); return (P&)*this; }
    P& check_code(const std::string& key, std::function<void()> func) { _group->add(new CmdKeyCode(key, func)); return (P&)*this; }
    template<class T, class V>
    P& prev_check(T& dst, V value) { _group->add(new CmdPrevCheck(dst, value)); return (P&)*this; }
    P& prev_code(std::function<void()> func) { _group->add(new CmdPrevCode(func)); return (P&)*this; }
    CmdKeyGroupConf<P> group(const std::string& key) { return CmdKeyGroupConf<P>(_group->add(new CmdGroup(key)), (P*)this); }
    CmdKeyListConf<P> list(const std::string& key, char delim, char list_delim) { return CmdKeyListConf<P>(_group->add(new CmdKeyList(key, delim, list_delim)), (P*)this); }
    CmdExclusiveConf<P> exclusive() { return CmdExclusiveConf<P>(_group->add(new CmdExclusive()), (P*)this); }

    R& end() { return (R&)*_parent; }

protected:
    CmdGroup* _group;
    R* _parent;
};

template<class R>
class CmdKeyGroupConf : public CmdGroupConf<CmdKeyGroupConf<R>, R>
{
public:
    CmdKeyGroupConf(CmdGroup* group, R* parent) : CmdGroupConf<CmdKeyGroupConf<R>, R>(group, parent) { }
};

template<class R>
class CmdExclusiveConf : public CmdExclusive
{
public:
    class CmdOptionalConf : public CmdGroupConf<CmdOptionalConf, CmdExclusiveConf>
    {
    public:
        CmdOptionalConf(CmdGroup* optional, CmdExclusiveConf* parent) : CmdGroupConf<CmdOptionalConf, CmdExclusiveConf>(optional, parent) { }
    };
    class CmdCaseConf : public CmdGroupConf<CmdCaseConf, CmdExclusiveConf>
    {
    public:
        CmdCaseConf(CmdCase* ex_case, CmdExclusiveConf* parent) : CmdGroupConf<CmdCaseConf, CmdExclusiveConf>(ex_case, parent) { }

        CmdOptionalConf optional() { return CmdOptionalConf(((CmdCase*)this->_group)->optional(), this->_parent); }
    };

public:
    CmdExclusiveConf(CmdExclusive* ex, R* parent) : _ex(ex), _parent(parent) { }

    CmdCaseConf ex_case() { return CmdCaseConf(_ex->add(new CmdCase()), this); }

    R& end() { return (R&)*_parent; }

protected:
    CmdExclusive* _ex;
    R* _parent;
};

class CmdLine : public CmdGroup, public CmdGroupConf<CmdLine, nullptr_t>
{
public:
    CmdLine() : CmdGroupConf(this, nullptr) { }

    CmdLine& default_code(std::function<void(const char*)> func) { _default_code = func; return *this; }
    void parse(int argc, char *argv[]);

private:
    std::function<void(const char*)> _default_code;
};
