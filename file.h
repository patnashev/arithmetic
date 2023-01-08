#pragma once

#include <memory>
#include <vector>
#include "arithmetic.h"

class Writer
{
public:
    Writer() { }
    Writer(std::vector<char>& buffer) : _buffer(buffer) { }
    Writer(std::vector<char>&& buffer) : _buffer(std::move(buffer)) { }

    void write(int32_t value);
    void write(uint32_t value);
    void write(uint64_t value);
    void write(double value);
    void write(const std::string& value);
    void write(const arithmetic::Giant& value);
    void write(const char* ptr, int count);
    void write_text(const char* ptr);
    void write_text(const std::string& value);
    void write_textline(const char* ptr);
    void write_textline(const std::string& value);

    std::vector<char>& buffer() { return _buffer; }
    std::vector<char> hash();
    std::string hash_str();

private:
    std::vector<char> _buffer;
};

class Reader
{
public:
    Reader(char format_version, char type, char version, char *data, int size, int pos) : _format_version(format_version), _type(type), _version(version), _data(data), _size(size), _pos(pos) { }

    bool read(int32_t& value);
    bool read(uint32_t& value);
    bool read(uint64_t& value);
    bool read(double& value);
    bool read(std::string& value);
    bool read(arithmetic::Giant& value);

    char type() { return _type; }
    char version() { return _version; }

private:
    char _format_version;
    char _type;
    char _version;
    char* _data;
    int _size;
    int _pos = 0;
};

class TextReader
{
public:
    TextReader(char *data, int size, int pos) : _data(data), _size(size), _pos(pos) { }

    bool read_textline(std::string& value);

private:
    char* _data;
    int _size;
    int _pos = 0;
};

class TaskState;

class File
{
public:
    static const uint32_t MAGIC_NUM = 0x9f2b3cd4;
    static int FILE_APPID;

public:
    File(const std::string& filename, uint32_t fingerprint) : _filename(filename), _fingerprint(fingerprint), _hash_filename(filename + ".md5") { }
    virtual ~File() { }
    File& operator = (File&&) = default;

    virtual File* add_child(const std::string& name, uint32_t fingerprint);
    const std::vector<std::unique_ptr<File>>& children() { return _children; }

    bool read(TaskState& state);
    void write(TaskState& state);
    void write_text(const std::string& value);
    void write_textline(const std::string& value);

    virtual Reader* get_reader();
    virtual TextReader* get_textreader();
    virtual Writer* get_writer();
    virtual Writer* get_writer(char type, char version);

    virtual void read_buffer();
    virtual void commit_writer(Writer& writer);
    virtual void free_buffer();
    virtual void clear(bool recursive = false);

    static uint32_t unique_fingerprint(uint32_t fingerprint, const std::string& unique_id);

    void set_filename(std::string& filename, std::string& hash_filename) { _filename = filename; _hash_filename = hash_filename; }
    std::string& filename() { return _filename; }
    uint32_t fingerprint() { return _fingerprint; }
    std::vector<char>& buffer() { return _buffer; }
    bool hash = true;
    int appid = FILE_APPID;

protected:
    std::string _filename;
    std::string _hash_filename;
    uint32_t _fingerprint;
    std::vector<char> _buffer;
    std::vector<std::unique_ptr<File>> _children;
};

class FileEmpty : public File
{
public:
    FileEmpty() : File("", 0) { }

    File* add_child(const std::string& /*name*/, uint32_t /*fingerprint*/) override { return _children.emplace_back(new FileEmpty()).get(); }
    void read_buffer() override { }
    void commit_writer(Writer& /*writer*/) override { }
};
