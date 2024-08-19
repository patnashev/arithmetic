
#define _FILE_OFFSET_BITS 64
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <cctype>
#include <cerrno>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "container.h"
#include "md5.h"
#ifdef _WIN32
#include <io.h>
#include "windows.h"
#else
#include <unistd.h>
#endif

namespace container
{
    void skip_space(const char*& buffer, size_t& count)
    {
        while (count > 0 && std::isspace(*buffer))
            buffer++, count--;
    }

    size_t len_string(const char* buffer, size_t count)
    {
        int utf8state = 0;
        int utf16state = 0;
        size_t res = 1;
        buffer++, count--;
        for (; count > 0; buffer++, count--, res++)
            if ((*buffer & 0x80) != 0 || utf8state != 0 || utf16state != 0)
            {
                if ((*buffer & 0xC0) == 0x80)
                {
                    if (utf8state == 0) // missing lead
                        return 0;
                    if (utf8state < 0)
                    {
                        if (utf8state == -3 && (*buffer & 0xF0) != 0x80) // > 0x10FFFF
                            return 0;
                        if (utf8state == -2 && (*buffer & 0xE0) != 0x80) // surrogate
                        {
                            //return 0;
                            if (utf16state == 0 && (*buffer & 0x10) != 0x00)
                                return 0;
                            if (utf16state == 1 && (*buffer & 0x10) != 0x10)
                                return 0;
                            utf16state = ((utf16state + 1) & 1);
                        }
                        utf8state = -utf8state;
                    }
                    utf8state--;
                    continue;
                }
                if (utf8state != 0) // missing continuation
                    return 0;
                if ((*buffer & 0xF8) == 0xF8 || ((*buffer & 0xF4) == 0xF4 && (*buffer & 3) != 0)) // > 0x10FFFF
                    return 0;
                else if ((*buffer & 0xF4) == 0xF4)
                    utf8state = -3;
                else if ((*buffer & 0xF8) == 0xF0)
                    utf8state = 3;
                else if (*buffer == (char)0xED)
                    utf8state = -2;
                else if ((*buffer & 0xF0) == 0xE0)
                    utf8state = 2;
                else if ((*buffer & 0xE0) == 0xC0)
                    utf8state = 1;
                else // missing surrogate
                    return 0;
                if (utf16state == 1 && utf8state != -2) // missing surrogate
                    return 0;
            }
            else if (*buffer == '\"')
                return res + 1;
            else if (*buffer == '\\')
            {
                buffer++, count--, res++;
                if (count == 0)
                    return 0;
                switch (*buffer)
                {
                case 'u':
                    if (count < 5)
                        break;
                    if (buffer[1] == 'd' || buffer[1] == 'D') // surrogate
                    {
                        if (buffer[2] != '8' && buffer[2] != '9' && buffer[2] != 'a' && buffer[2] != 'A' && buffer[2] != 'b' && buffer[2] != 'B')
                            break;
                        for (size_t i = 0; i < 4; i++, buffer++, count--, res++)
                            if (!std::isxdigit(buffer[1]))
                                break;
                        buffer++, count--, res++;
                        if (count < 6)
                            break;
                        if (*buffer != '\\')
                            break;
                        buffer++, count--, res++;
                        if (*buffer != 'u' || (buffer[1] != 'd' && buffer[1] != 'D'))
                            break;
                        if (buffer[2] != 'c' && buffer[2] != 'C' && buffer[2] != 'd' && buffer[2] != 'D' && buffer[2] != 'e' && buffer[2] != 'E' && buffer[2] != 'f' && buffer[2] != 'F')
                            break;
                    }
                    for (size_t i = 0; i < 4; i++, buffer++, count--, res++)
                        if (!std::isxdigit(buffer[1]))
                            break;
                    continue;
                case 'b':
                case 'f':
                case 'n':
                case 'r':
                case 't':
                case '\"':
                case '\\':
                case '/':
                    continue;
                }
                return 0;
            }
        return 0;
    }

    std::string parse_string(const char* buffer, size_t count)
    {
        count--;
        buffer++, count--;
        std::string res;
        res.reserve(count);
        for (; count > 0; buffer++, count--)
            if (*buffer == '\\')
            {
                buffer++, count--;
                switch (*buffer)
                {
                case 'b': res += '\b'; break;
                case 'f': res += '\f'; break;
                case 'n': res += '\n'; break;
                case 'r': res += '\r'; break;
                case 't': res += '\t'; break;
                case '\"': res += '\"'; break;
                case '\\': res += '\\'; break;
                case '/': res += '/'; break;
                case 'u':
                    std::string su(buffer + 1, 4);
                    int u16 = std::stoi(su, nullptr, 16);
                    if (u16 <= 0x7F)
                        res += (char)u16;
                    else if (u16 <= 0x7FF)
                    {
                        res += (char)(0xC0 | (u16 >> 6));
                        res += (char)(0x80 | (u16 & 0x3F));
                    }
                    else if (u16 >= 0xD800 && u16 <= 0xDBFF)
                    {
                        int u32 = 0x10000 + ((u16 & 0x3FF) << 10);
                        buffer += 6, count -= 6;
                        std::string su2(buffer + 1, 4);
                        u32 += (std::stoi(su2, nullptr, 16) & 0x3FF);
                        res += (char)(0xF0 | (u32 >> 18));
                        res += (char)(0x80 | ((u32 >> 12) & 0x3F));
                        res += (char)(0x80 | ((u32 >> 6) & 0x3F));
                        res += (char)(0x80 | (u32 & 0x3F));
                    }
                    else
                    {
                        res += (char)(0xE0 | (u16 >> 12));
                        res += (char)(0x80 | ((u16 >> 6) & 0x3F));
                        res += (char)(0x80 | (u16 & 0x3F));
                    }
                    buffer += 4, count -= 4;
                }
            }
            else if (*buffer == (char)0xED && (buffer[1] & 0xE0) != 0x80) // surrogate
            {
                buffer++, count--;
                int u32 = 0x10000 + ((*buffer & 0xF) << 16);
                buffer++, count--;
                u32 += ((*buffer & 0x3F) << 10);
                buffer++, count--;
                buffer++, count--;
                u32 += ((*buffer & 0xF) << 6);
                buffer++, count--;
                u32 += (*buffer & 0x3F);
                res += (char)(0xF0 | (u32 >> 18));
                res += (char)(0x80 | ((u32 >> 12) & 0x3F));
                res += (char)(0x80 | ((u32 >> 6) & 0x3F));
                res += (char)(0x80 | (u32 & 0x3F));
            }
            else
                res += *buffer;
        return res;
    }

    bool JSON::Node::parse(const char*& buffer, size_t& count)
    {
        skip_space(buffer, count);
        if (count == 0)
            return false;

        if (*buffer == '[')
        {
            auto arr = std::vector<Node>();
            while (true)
            {
                if (!arr.empty() && *buffer == ']')
                    break;
                buffer++, count--;
                skip_space(buffer, count);
                if (count == 0)
                    return false;
                if (arr.empty() && *buffer == ']')
                    break;
                if (!arr.emplace_back().parse(buffer, count))
                    return false;
                skip_space(buffer, count);
                if (count == 0 || (*buffer != ',' && *buffer != ']'))
                    return false;
            }
            buffer++, count--;
            _value = std::move(arr);
        }
        else if (*buffer == '{')
        {
            auto keys = std::vector<std::string>();
            auto vals = std::vector<Node>();
            while (true)
            {
                if (!keys.empty() && *buffer == '}')
                    break;
                buffer++, count--;
                skip_space(buffer, count);
                if (count == 0)
                    return false;
                if (keys.empty() && *buffer == '}')
                    break;

                if (*buffer != '\"')
                    return false;
                size_t len = len_string(buffer, count);
                if (len == 0)
                    return false;
                std::string key = parse_string(buffer, len);
                auto it = std::lower_bound(keys.begin(), keys.end(), key);
                size_t index = it - keys.begin();
                keys.insert(it, std::move(key));
                buffer += len, count -= len;

                skip_space(buffer, count);
                if (count == 0)
                    return false;
                if (*buffer != ':')
                    return false;
                buffer++, count--;
                skip_space(buffer, count);
                if (count == 0)
                    return false;

                if (!vals.emplace(vals.begin() + index)->parse(buffer, count))
                    return false;
                skip_space(buffer, count);
                if (count == 0 || (*buffer != ',' && *buffer != '}'))
                    return false;
            }
            buffer++, count--;
            _value = std::pair(keys, vals);
        }
        else if (*buffer == '\"')
        {
            size_t len = len_string(buffer, count);
            if (len == 0)
                return false;
            _value = std::pair<const char*, size_t>(buffer, len);
            buffer += len, count -= len;
        }
        else
        {
            size_t len;
            for (len = 0; len < count && !std::isspace(buffer[len]) && buffer[len] != ',' && buffer[len] != ']' && buffer[len] != '}'; len++);
            if (len == 0)
                return false;
            if (len == 4 && strncmp(buffer, "null", 4) == 0)
                _value = nullptr;
            else
            {
                _value = std::pair<const char*, size_t>(buffer, len);
                if (!is_bool() && !is_int() && !is_double())
                    return false;
            }
            buffer += len, count -= len;
        }

        return true;
    }

    bool JSON::Node::is_bool() const
    {
        if (type() != Type::Value)
            return false;
        auto [buffer, count] { value_buffer() };
        if (buffer[0] == '\"' && buffer[count - 1] == '\"')
            buffer++, count -= 2;
        if (count == 4 && strncmp(buffer, "true", 4) == 0)
            return true;
        if (count == 5 && strncmp(buffer, "false", 4) == 0)
            return true;
        return false;
    }

    bool JSON::Node::is_int() const
    {
        if (type() != Type::Value)
            return false;
        auto [buffer, count] { value_buffer() };
        if (buffer[0] == '\"' && buffer[count - 1] == '\"')
            buffer++, count -= 2;
        if (*buffer == '-')
            buffer++, count--;
        if (count == 0 || count > 18)
            return false;
        for (; count > 0 && std::isdigit(*buffer); buffer++, count--);
        return count == 0;
    }

    bool JSON::Node::is_double() const
    {
        if (type() != Type::Value)
            return false;
        auto [buffer, count] { value_buffer() };
        if (buffer[0] == '\"' && buffer[count - 1] == '\"')
            buffer++, count -= 2;
        if (*buffer == '-')
            buffer++, count--;
        if (count == 0)
            return false;
        for (; count > 0 && std::isdigit(*buffer); buffer++, count--);
        if (count > 0 && *buffer == '.')
            buffer++, count--;
        for (; count > 0 && std::isdigit(*buffer); buffer++, count--);
        if (count > 0 && (*buffer == 'e' || *buffer == 'E'))
        {
            buffer++, count--;
            if (count == 0)
                return false;
            if (*buffer == '-' || *buffer == '+')
                buffer++, count--;
            if (count == 0)
                return false;
            for (; count > 0 && std::isdigit(*buffer); buffer++, count--);
        }
        return count == 0;
    }

    bool JSON::Node::is_string() const
    {
        if (type() != Type::Value)
            return false;
        auto [buffer, count] { value_buffer() };
        return buffer[0] == '\"' && buffer[count - 1] == '\"';
    }

    bool JSON::Node::value_bool() const
    {
        auto [buffer, count] { value_buffer() };
        if (buffer[0] == '\"' && buffer[count - 1] == '\"')
            buffer++, count -= 2;
        return count == 4;
    }

#ifdef CHARCONV
#include <charconv>
#endif
    int64_t JSON::Node::value_int() const
    {
        auto [buffer, count] { value_buffer() };
        if (buffer[0] == '\"' && buffer[count - 1] == '\"')
            buffer++, count -= 2;
#ifdef CHARCONV
        int64_t res;
        auto [ptr, ec] { std::from_chars(buffer, buffer + count, res) };
        if (ec != std::errc())
            throw std::out_of_range(std::string(buffer, count));
#else
        char* str_end;
        long long res = std::strtoll(buffer, &str_end, 10);
        if (errno == ERANGE || (str_end - buffer) != count)
            throw std::out_of_range(std::string(buffer, count));
#endif
        return res;
    }

    double JSON::Node::value_double() const
    {
        auto [buffer, count] { value_buffer() };
        if (buffer[0] == '\"' && buffer[count - 1] == '\"')
            buffer++, count -= 2;
        char* str_end;
        double res = std::strtod(buffer, &str_end);
        if (errno == ERANGE || !std::isnormal(res) || (str_end - buffer) != count)
            throw std::out_of_range(std::string(buffer, count));
        return res;
    }

    std::string JSON::Node::value_string() const
    {
        auto [buffer, count] { value_buffer() };
        if (buffer[0] != '\"' || buffer[count - 1] != '\"')
            throw std::bad_variant_access();
        return parse_string(buffer, count);
    }

    void encode_string(const std::string& src, std::string& dst)
    {
        dst.reserve(dst.size() + src.size() + 10);
        dst += '\"';
        for (auto it = src.begin(); it != src.end(); it++)
            switch (*it)
            {
            case '\b': dst += "\\b"; break;
            case '\f': dst += "\\f"; break;
            case '\n': dst += "\\n"; break;
            case '\r': dst += "\\r"; break;
            case '\t': dst += "\\t"; break;
            case '\"': dst += "\\\""; break;
            case '\\': dst += "\\\\"; break;
            default:
                if ((*it > 0 && *it <= 0x1f) || *it == 0x7f)
                {
                    dst += "\\u0000";
                    snprintf(&*(dst.end() - 2), 3, "%02x", (int)*it);
                }
                else
                    dst += *it;
            }
        dst += '\"';
    }

    void JSON::Node::to_string(std::string& res) const
    {
        switch (type())
        {
        case Type::Null:
            res += "null";
            break;
        case Type::Value:
        {
            auto [buffer, count] { value_buffer() };
            res.insert(res.end(), buffer, buffer + count);
            break;
        }
        case Type::Array:
        {
            res += '[';
            auto& vals = values();
            for (int i = 0; i < vals.size(); i++)
            {
                if (i > 0)
                    res += ',';
                vals[i].to_string(res);
            }
            res += ']';
            break;
        }
        case Type::Object:
        {
            res += '{';
            auto& mkeys = keys();
            auto& mvals = values();
            for (int i = 0; i < mkeys.size(); i++)
            {
                if (i > 0)
                    res += ',';
                encode_string(mkeys[i], res);
                res += ':';
                mvals[i].to_string(res);
            }
            res += '}';
            break;
        }
        }
    }

    bool JSON::parse_node(Node& node)
    {
        const char* buffer = _buffers.back().data();
        size_t count = _buffers.back().size();
        bool success = node.parse(buffer, count);
        if (success)
        {
            skip_space(buffer, count);
            success = (count == 0);
        }
        if (!success)
        {
            node._value = nullptr;
            _buffers.pop_back();
        }
        return success;
    }

    void JSON::set_value(Node& node, std::nullptr_t value)
    {
        set_json(node, "null");
    }

    void JSON::set_value(Node& node, bool value)
    {
        set_json(node, value ? "true" : "false");
    }

    void JSON::set_value(Node& node, int64_t value)
    {
        set_json(node, std::to_string(value));
    }

    void JSON::set_value(Node& node, double value)
    {
        char buf[32];
        snprintf(buf, 32, "%.16g", value);
        set_json(node, buf);
    }

    void JSON::set_value(Node& node, const std::string& value)
    {
        std::string json;
        encode_string(value, json);
        set_json(node, json);
    }



    void MemoryStream::set_position(int64_t value)
    {
        if (value < 0)
            throw std::out_of_range("Invalid position.");
        if (value > (int64_t)_data.size())
            _pos = (int64_t)_data.size();
        else
            _pos = value;
    }

    void MemoryStream::set_length(int64_t value)
    {
        if (value == -1)
            _data.resize(_pos);
        else
        {
            _data.resize(value);
            if (_pos > value)
                _pos = value;
        }
        _length = value;
    }

    size_t MemoryStream::read(char* buffer, size_t count)
    {
        if ((size_t)_pos >= _data.size())
            return 0;
        if ((size_t)_pos + count > _data.size())
            count = _data.size() - (size_t)_pos;
        memcpy(buffer, _data.data() + _pos, count);
        _pos += count;
        return count;
    }

    void MemoryStream::write(const char* buffer, size_t count)
    {
        if (_length != -1 && (size_t)_pos >= _data.size())
            return;
        if (_length == -1 && (size_t)_pos + count > _data.size())
            _data.resize((size_t)_pos + count);
        if ((size_t)_pos + count > _data.size())
            count = _data.size() - (size_t)_pos;
        memcpy(_data.data() + _pos, buffer, count);
        _pos += count;
    }

    void LimitedReadStream::set_position(int64_t value)
    {
        if (value < _pos)
            throw std::out_of_range("Invalid position.");
        if (value > _length)
            value = _length;
        int64_t pos = _stream.position();
        _stream.set_position(pos + value - _pos);
        _pos += _stream.position() - pos;
    }

    size_t LimitedReadStream::read(char* buffer, size_t count)
    {
        if (_pos >= _length)
            return 0;
        if (_pos + (int64_t)count > _length)
            count = (size_t)(_length - _pos);
        count = _stream.read(buffer, count);
        _pos += count;
        return count;
    }

    void FileStream::open(const std::string& filename)
    {
        if (_stream)
            close();
        _destroy = true;
        _stream = fopen(filename.data(), _read && _write ? "a+b" : _read ? "rb" : "wb");
        if (!_stream)
        {
            _length = 0;
            return;
        }
        if (_read)
        {
            if (fseek((FILE*)_stream, 0, SEEK_END) == 0)
            {
#ifdef _WIN32
                _length = _ftelli64((FILE*)_stream);
                _fseeki64((FILE*)_stream, 0, SEEK_SET);
#else
                _length = ftello((FILE*)_stream);
                fseeko((FILE*)_stream, 0, SEEK_SET);
#endif
            }
            else
                throw std::runtime_error("File random access failed.");
        }
    }

    void FileStream::open(void* filestream, bool destroy)
    {
        if (_stream)
            close();
        _destroy = destroy;
        _stream = filestream;
        if (!_stream)
        {
            _length = 0;
            return;
        }
        if (_read)
        {
#ifdef _WIN32
            int64_t start = _ftelli64((FILE*)_stream);
#else
            int64_t start = ftello((FILE*)_stream);
#endif
            if (fseek((FILE*)_stream, 0, SEEK_END) == 0)
            {
#ifdef _WIN32
                _length = _ftelli64((FILE*)_stream) - start;
                _fseeki64((FILE*)_stream, start, SEEK_SET);
#else
                _length = ftello((FILE*)_stream) - start;
                fseeko((FILE*)_stream, start, SEEK_SET);
#endif
            }
            else
                throw std::runtime_error("File random access failed.");
        }
    }

    void FileStream::set_position(int64_t value)
    {
        if (!_stream)
            return;
        if (value < 0 || value > _length)
            throw std::out_of_range("Invalid file position.");
        if (
#ifdef _WIN32
            _fseeki64((FILE*)_stream, value - _pos, SEEK_CUR)
#else
            fseeko((FILE*)_stream, value - _pos, SEEK_CUR)
#endif
            != 0)
            throw std::runtime_error("File random access failed.");
        _pos = value;
    }

    void FileStream::set_length(int64_t value)
    {
        if (!_stream)
            return;
        if (value == -1 && _pos < _length)
        {
            if (fflush((FILE*)_stream) != 0)
                throw std::runtime_error("File flush failed.");
#ifdef _WIN32
            int64_t cur = _ftelli64((FILE*)_stream);
            _chsize_s(_fileno((FILE*)_stream), cur);
            _fseeki64((FILE*)_stream, cur, SEEK_SET);
#else
            int64_t cur = ftello((FILE*)_stream);
            ftruncate(fileno((FILE*)_stream), cur);
            fseeko((FILE*)_stream, cur, SEEK_SET);
#endif
        }
        if (value != -1 && ((_length != -1 && _length != value) || (_length == -1 && _pos > value)))
        {
            if (fflush((FILE*)_stream) != 0)
                throw std::runtime_error("File flush failed.");
#ifdef _WIN32
            int64_t start = _ftelli64((FILE*)_stream) - _pos;
            _chsize_s(_fileno((FILE*)_stream), start + value);
#else
            int64_t start = ftello((FILE*)_stream) - _pos;
            ftruncate(fileno((FILE*)_stream), start + value);
#endif
            if (_pos > value)
                _pos = value;
#ifdef _WIN32
            _fseeki64((FILE*)_stream, start + _pos, SEEK_SET);
#else
            fseeko((FILE*)_stream, start + _pos, SEEK_SET);
#endif
        }
        _length = value;
    }

    size_t FileStream::read(char* buffer, size_t count)
    {
        if (!_stream)
            return 0;
        if (_length != -1 && _pos >= _length)
            return 0;
        if (_length != -1 && _pos + (int64_t)count > _length)
            count = (size_t)(_length - _pos);
        count = fread(buffer, 1, count, (FILE*)_stream);
        _pos += count;
        return count;
    }

    bool FileStream::readline(std::string& res, int max_chars)
    {
        if (!_stream)
            return false;
        int64_t start = _pos;
        while (true)
        {
            int c;
            if ((c = fgetc((FILE*)_stream)) == EOF)
                return start < _pos;
            _pos++;
            if (c == '\r')
            {
                if ((c = fgetc((FILE*)_stream)) != EOF)
                    _pos++;
                if (c == '\n')
                    return true;
                else if (max_chars != -1 && _pos - start > max_chars)
                    break;
                else
                    res += '\r';
                if (c == EOF)
                    return true;
            }
            else if (max_chars != -1 && _pos - start > max_chars)
                break;
            if (c == 0)
                break;
            res += (char)c;
        }
        set_position(start);
        return false;
    }

    void FileStream::write(const char* buffer, size_t count)
    {
        if (!_stream)
            return;
        if (_length != -1 && _pos >= _length)
            return;
        if (_length != -1 && _pos + (int64_t)count > _length)
            count = (size_t)(_length - _pos);
        size_t written = fwrite(buffer, 1, count, (FILE*)_stream);
        _pos += written;
        if (written != count)
            throw std::runtime_error("File write failed.");
    }

    void FileStream::flush()
    {
        if (!_stream)
            return;
        if (fflush((FILE*)_stream) != 0)
            throw std::runtime_error("File flush failed.");
    }

    void FileStream::close()
    {
        if (!_stream)
            return;
        int error;
        if (_destroy)
            error = fclose((FILE*)_stream);
        else
            error = fflush((FILE*)_stream);
        _stream = nullptr;
        if (error != 0)
            throw std::runtime_error("File close failed.");
    }

    int ChunkedWriteStream::MAX_BUFFER_SIZE = 16777200;

    void ChunkedWriteStream::set_length(int64_t value)
    {
        if (_buffer_size > (size_t)value)
            _buffer_size = (size_t)value;
    }

    void ChunkedWriteStream::write(const char* buffer, size_t count)
    {
        if (_buffer.size() >= _buffer_size)
            write_chunk();
        while (_buffer.size() + count >= _buffer_size)
        {
            size_t left = _buffer_size - _buffer.size();
            _buffer.insert(_buffer.end(), buffer, buffer + left);
            write_chunk();
            buffer += left, count -= left;
        }
        _buffer.reserve(_buffer_size);
        if (count > 0)
            _buffer.insert(_buffer.end(), buffer, buffer + count);
    }

    void ChunkedWriteStream::write_chunk()
    {
        if (_buffer.empty())
            return;
        std::string json = "[" + std::to_string(_buffer.size()) + "," + std::to_string(_id) + (_codec_json.empty() ? "" : ",") + (_codec_json.empty() ? "" : _codec_json) + "]\r\n";
        _stream.write(json.data(), json.size());
        if (_container && _id > 0)
        {
            auto& chunk = _container->_streams[_id].chunks.emplace_back();
            chunk.offset = _container->_stream->position();
            chunk.size = (int64_t)_buffer.size();
            if (!_codec_json.empty())
                chunk.codec.reset(new JSON(std::move(_codec_json)));
        }
        _stream.write(_buffer.data(), _buffer.size());

        _codec_json.clear();
        _buffer.clear();
    }

    void ChunkedWriteStream::flush()
    {
        write_chunk();
        _stream.flush();
    }

    void ChunkedWriteStream::close()
    {
        write_chunk();
        std::vector<char>().swap(_buffer);
        //std::string json = "[0," + std::to_string(_id) + "]\r\n";
        //_stream.write(json.data(), json.size());
    }

    char ctobase64(char val)
    {
        if (val < 26)
            return 'A' + val;
        if (val < 52)
            return 'a' - 26 + val;
        if (val < 62)
            return '0' - 52 + val;
        if (val == 62)
            return '+';
        if (val == 63)
            return '/';
        return 0;
    }

    char base64toc(char c)
    {
        if (c >= 'A' && c <= 'Z')
            return c - 'A';
        if (c >= 'a' && c <= 'z')
            return c - 'a' + 26;
        if (c >= '0' && c <= '9')
            return c - '0' + 52;
        if (c == '+')
            return 62;
        if (c == '/')
            return 63;
        if (c == '=')
            return 64;
        return -1;
    }

    void Base64CoderStream::set_length(int64_t value)
    {
        _length = value;
    }

    void Base64CoderStream::write(const char* buffer, size_t count)
    {
        if (_length != -1 && _pos >= _length)
            return;
        if (_length != -1 && _pos + (int64_t)count > _length)
            count = (size_t)(_length - _pos);
        _buffer.reserve((count + 2)/3*4 + (count + 53)/54*2);
        for (const char* end = buffer + count; buffer != end; buffer++)
            if ((_state & 3) == 0)
            {
                _buffer.push_back(ctobase64(((unsigned char)*buffer >> 2)));
                _next = ((unsigned char)*buffer & 3) << 4;
                _state++;
            }
            else if ((_state & 3) == 1)
            {
                _buffer.push_back(ctobase64(_next + ((unsigned char)*buffer >> 4)));
                _next = ((unsigned char)*buffer & 15) << 2;
                _state++;
            }
            else if ((_state & 3) == 2)
            {
                _buffer.push_back(ctobase64(_next + ((unsigned char)*buffer >> 6)));
                _buffer.push_back(ctobase64(((unsigned char)*buffer & 63)));
                _state += 2;
                if (_state >= 72)
                {
                    _buffer.push_back('\r');
                    _buffer.push_back('\n');
                    _state = 0;
                }
            }
        _stream.write(_buffer.data(), _buffer.size());
        _buffer.clear();
        _pos += count;
    }

    void Base64CoderStream::write_last()
    {
        if (_state == 0)
            return;
        if ((_state & 3) == 1)
        {
            _buffer.push_back(ctobase64(_next));
            _buffer.push_back('=');
            _buffer.push_back('=');
        }
        if ((_state & 3) == 2)
        {
            _buffer.push_back(ctobase64(_next));
            _buffer.push_back('=');
        }
        _state = 0;
        _buffer.push_back('\r');
        _buffer.push_back('\n');
        _stream.write(_buffer.data(), _buffer.size());
        _buffer.clear();
    }

    void Base64CoderStream::flush()
    {
        write_last();
        _stream.flush();
    }

    void Base64CoderStream::close()
    {
        write_last();
    }

    void Base64DecoderReadStream::set_position(int64_t value)
    {
        if (value < _pos)
            throw std::out_of_range("Invalid position.");
        read(nullptr, value - _pos);
    }

    size_t Base64DecoderReadStream::read(char* buffer, size_t count)
    {
        size_t total = 0;
        while (count > 0)
        {
            if (_buffer_pos == _buffer.size())
            {
                _buffer_pos = 0;
                _buffer.resize((count + 2)/3*4 + (count + 53)/54*2);
                _buffer.resize(_stream.read(_buffer.data(), _buffer.size()));
                if (_buffer.size() == 0)
                {
                    if ((_state & 3) != 0)
                        throw container_error("Base64 codec error.");
                    break;
                }
            }
            char c = _buffer[_buffer_pos++];
            if (std::isspace(c))
                continue;
            if ((c = base64toc(c)) == -1)
                throw container_error("Base64 invalid character.");
            if (c == 64 && (_state == 0 || _state == 1))
                throw container_error("Base64 invalid character.");
            if (_state == 0)
            {
                if (buffer)
                    _next = c << 2;
                _state++;
            }
            else if (_state == 1)
            {
                if (buffer)
                {
                    *buffer = _next + (c >> 4);
                    buffer++, total++, count--;
                    _next = (c & 15) << 4;
                }
                _state++;
            }
            else if (_state == 2 && c == 64)
                _state = 4;
            else if (_state == 2)
            {
                if (buffer)
                {
                    *buffer = _next + (c >> 2);
                    buffer++, total++, count--;
                    _next = (c & 3) << 6;
                }
                _state++;
            }
            else if (_state == 3 && c == 64)
                _state = 0;
            else if (_state == 3)
            {
                if (buffer)
                {
                    *buffer = _next + c;
                    buffer++, total++, count--;
                }
                _state = 0;
            }
            else if (_state == 4 && c == 64)
                _state = 0;
            else if (_state == 4)
                throw container_error("Base64 invalid character.");
        }
        _pos += total;
        return total;
    }

    Packer::Packer(const std::string& filename) : _file(new FileStream(filename, false, true)), _stream(_file.get())
    {
        _stream->write("[\"PK\",1,0]\r\n", 12);
    }

    Packer::Packer(WriteStream* stream) : _stream(stream)
    {
        _stream->write("[\"PK\",1,0]\r\n", 12);
    }

    Packer::Packer(FileContainer& container) : _stream(container._stream), _next_id(container._next_stream_id), _index(std::move(container._index)), _container(&container)
    {
        if (!container._stream || container._stream->is_closed() || !container._stream->can_write())
        {
            _container->_index = std::move(_index);
            throw std::runtime_error("Can't write to the stream.");
        }

        int64_t pos = container._last_pos;

        for (auto stream = container._streams.begin(); stream != container._streams.end(); stream++)
        {
            if (stream->second.files.empty())
                continue;
            int64_t length = 0;
            for (auto it = stream->second.files.begin(); it != stream->second.files.end(); it++)
                if (length < (*it)->offset + (*it)->size)
                    length = (*it)->offset + (*it)->size;
            int64_t cur = 0;
            for (auto it = stream->second.chunks.begin(); it != stream->second.chunks.end(); it++)
            {
                if (length >= 0)
                {
                    if (it->codec && !it->codec->root().is_null())
                        length = -1;
                    else if (cur >= length)
                        break;
                    cur += it->size;
                }
                if (pos < it->offset + it->size)
                    pos = it->offset + it->size;
            }
        }

        try
        {
            container._stream->set_position(pos);
            container._stream->set_length(-1);

            if (container._error == container_error::EMPTY)
                _stream->write("[\"\xD0\x9A\",1,0]\r\n", 12);
            container._error = container_error::OK;
        }
        catch (const std::exception& e)
        {
            on_error(e);
        }

        container._cur_file = nullptr;
        container._cur_reader.reset();
    }

    void Packer::on_error(const std::exception& e)
    {
        _stream = nullptr;
        _file.reset();
        if (_container && _container->_stream)
        {
            _container->_filesize = 0;
            _container->_stream = nullptr;
            _container->_file.reset();
        }
        throw e;
    }

    void Packer::close()
    {
        if (!_stream)
            return;
        try
        {
            while (!_writers.empty())
                _writers.back()->close();
            _stream->write("[0,0]\r\n", 7);
            _file.reset();
        }
        catch (const std::exception& e)
        {
            on_error(e);
        }
        _stream = nullptr;
        if (_container && _container->_stream)
        {
            _container->_stream->flush();
            _container->_stream->set_length(_container->_stream->position());
            _container->_filesize = _container->_stream->length();
            _container->_index = std::move(_index);
            _container->_next_stream_id = _next_id;
        }
    }

    class Packer::Writer::Index : public WriteStream
    {
    public:
        Index(Packer::Writer& writer) : _writer(writer) { _json.parse("[]"); _cur = 0; }

        FileDesc& file() { return *_file; }
        JSON& json() { return _json; }

        JSON::Node& node(FileDesc* file)
        {
            _file = file;
            auto it = std::lower_bound(_files.begin(), _files.end(), file->name);
            _cur = it - _files.begin();
            if (it == _files.end() || *it != file->name)
            {
                _files.insert(it, file->name);
                JSON::Node& res = _json.root().insert(_cur);
                _json.set_value(res.insert("file"), file->name);
                return res;
            }
            return _json.root()[_cur];
        }

        template<class T>
        void set_value(const std::string& key, T value)
        {
            if (_cur >= _files.size())
                node(_file);
            _json.set_value(_json.root()[_cur].insert(key), value);
        }
        template<class T>
        bool set_json(const std::string& key, T&& json)
        {
            if (_cur >= _files.size())
                node(_file);
            return _json.set_json(_json.root()[_cur].insert(key), std::forward<T>(json));
        }

        void write_index()
        {
            if (_files.empty())
                return;
            int i;
            std::string index = "[";
            for (i = 0; i < _json.root().size(); i++)
            {
                if (i > 0)
                    index += ",\r\n";
                _json.root()[i].to_string(index);
            }
            index += "]\r\n";
            std::string chunk = "[" + std::to_string(index.size()) + ",0"/* + (_codec_json.empty() ? "" : ",") + (_codec_json.empty() ? "" : _codec_json)*/ + "]\r\n";
            _writer._packer._stream->write(chunk.data(), chunk.size());
            _writer._packer._stream->write(index.data(), index.size());
            _json.parse("[]");
            _files.clear();
            _cur = 0;
            if (_writer._packer._container != nullptr)
                _writer._packer._container->_last_pos = _writer._packer._container->_stream->position();
        }

        void set_length(int64_t value) override { }

        void write(const char* buffer, size_t count) override
        {
            if (!_writer._packer._stream)
                throw std::runtime_error("Invalid packer.");
            try
            {
                write_index();
                _writer._packer._stream->write(buffer, count);
            }
            catch (const std::exception& e)
            {
                _writer._packer.on_error(e);
            }
        }

        void flush() override
        {
            if (!_writer._packer._stream)
                throw std::runtime_error("Invalid packer.");
            try
            {
                write_index();
                _writer._packer._stream->flush();
            }
            catch (const std::exception& e)
            {
                _writer._packer.on_error(e);
            }
        }

        void close() override
        {
            if (!_writer._packer._stream)
                throw std::runtime_error("Invalid packer.");
            try
            {
                write_index();
            }
            catch (const std::exception& e)
            {
                _writer._packer.on_error(e);
            }
        }

    private:
        Packer::Writer& _writer;
        JSON _json;
        std::vector<std::string> _files;
        size_t _cur;
        FileDesc* _file;
    };

    class Packer::Writer::Stream : public WriteStream
    {
    public:
        Stream(Packer::Writer& writer, WriteStream* stream, int64_t length) : _writer(writer), _stream(stream), _length(length)
        {
            if (_writer._packer.md5)
                MD5Init(&md5_context);
        }
        ~Stream() { if (!closed()) close(); }

        bool closed() { return _stream == nullptr; }

        void set_length(int64_t value) override
        {
            if (closed())
                return;
            if (_length == value)
                return;
            if (_length >= 0)
                throw std::runtime_error("File size already set.");
            _length = value;
            Packer::Writer::Index* index = dynamic_cast<Packer::Writer::Index*>(_writer._streams.front().get());
            index->set_value("size", _length);
            index->file().size = _length;
        }

        void write(const char* buffer, size_t count) override
        {
            if (closed())
                return;
            if (_length != -1 && _pos + (int64_t)count > _length)
                count = (size_t)(_length - _pos);
            if (_writer._packer.md5)
                MD5Update(&md5_context, (unsigned char *)buffer, (unsigned int)count);
            _stream->write(buffer, count);
            _pos += count;
            if (_length != -1 && _pos == _length)
                close();
        }

        void flush() override
        {
            if (closed())
                return;
            _stream->flush();
        }

        void close() override
        {
            if (closed())
                return;
            if (_length == -1)
                set_length(_pos);
            if (_pos != _length)
                throw container_error("File size mismatch.");
            if (_writer._packer.md5)
            {
                std::vector<unsigned char> digest(16);
                MD5Final(digest.data(), &md5_context);
                char md5hash[33];
                for (int i = 0; i < 16; i++)
                    snprintf(md5hash + i*2, 3, "%02x", digest[i]);

                Packer::Writer::Index* index = dynamic_cast<Packer::Writer::Index*>(_writer._streams.front().get());
                if (!index->file().md5.empty() && index->file().md5 != md5hash)
                    throw container_error("MD5 hash mismatch.");
                index->file().md5 = md5hash;
                index->set_value("md5", index->file().md5);
            }

            _stream = nullptr;
            _writer._length += _length;
            _writer._streams.pop_back();
        }

    private:
        Packer::Writer& _writer;
        WriteStream* _stream;
        int64_t _length = -1;
        int64_t _pos = 0;
        MD5_CTX md5_context;
    };

    Packer::Writer::Writer(Packer& packer, int64_t id) : _packer(packer), _id(id)
    {
        if (!packer._stream)
            throw std::runtime_error("Invalid packer.");
        _streams.emplace_back(new Packer::Writer::Index(*this));
        ChunkedWriteStream* cs;
        _streams.emplace_back(cs = new ChunkedWriteStream(id, *_streams.back()));
        if (packer._container)
            cs->set_container(packer._container);
    }

    void Packer::Writer::add_codec(const std::string& codec)
    {
        if (codec == "base64")
            _streams.emplace_back(new Base64CoderStream(*_streams.back()));
        else
            throw std::invalid_argument("Codec not found.");
        ChunkedWriteStream* cs = dynamic_cast<ChunkedWriteStream*>(_streams[1].get());
        cs->set_codec("\"" + codec + "\"");
    }

    void Packer::Writer::clear_codecs()
    {

    }

    void Packer::Writer::add_file(const std::string& filename, const char* buffer, size_t count)
    {
        if (!_packer._stream)
            throw std::runtime_error("Invalid packer.");
        WriteStream* stream = add_file(filename, count);
        if (stream == nullptr)
            return;
        stream->write(buffer, count);
        stream->close();
    }

    WriteStream* Packer::Writer::add_file(const std::string& filename, int64_t size)
    {
        if (!_packer._stream)
            throw std::runtime_error("Invalid packer.");
        FileDesc& file = _packer.index()[filename];
        if (_packer._container && file.stream > 0)
        {
            auto& stream = _packer._container->_streams[file.stream];
            for (auto it = stream.files.begin(); it != stream.files.end(); )
                if ((*it) == &file)
                    it = stream.files.erase(it);
                else
                    it++;
            for (auto it = _packer._container->_corrupted.begin(); it != _packer._container->_corrupted.end(); )
                if (it->name == filename)
                    it = _packer._container->_corrupted.erase(it);
                else
                    it++;
        }
        file.reset();
        file.size = size;
        return add_file(&file);
    }

    WriteStream* Packer::Writer::add_file(FileDesc* file)
    {
        if (_streams.empty())
            throw std::runtime_error("Writer closed.");
        if (!_streams.empty() && dynamic_cast<Packer::Writer::Stream*>(_streams.back().get()) != nullptr)
            throw std::runtime_error("Another write operation in progress.");

        Packer::Writer::Index* index = dynamic_cast<Packer::Writer::Index*>(_streams.front().get());
        index->node(file);
        if (file->size >= 0)
            index->set_value("size", file->size);
        if (file->size == 0)
            return nullptr;
        index->set_value("stream", _id);
        file->stream = _id;
        if (_length > 0)
            index->set_value("offset", _length);
        file->offset = _length;
        if (!file->md5.empty())
            index->set_value("md5", file->md5);

        if (_packer._container)
            _packer._container->_streams[_id].files.push_back(file);

        return _streams.emplace_back(new Packer::Writer::Stream(*this, _streams.back().get(), file->size)).get();
    }

    void Packer::Writer::close()
    {
        while (!_streams.empty())
        {
            _streams.back()->close();
            _streams.pop_back();
        }
        _id = 0;
        for (auto it = _packer._writers.begin(); it != _packer._writers.end(); it++)
            if (it->get() == this)
            {
                _packer._writers.erase(it);
                break;
            }
    }


    class FileContainer::Reader::RawReadStream : public ReadStream
    {
    public:
        RawReadStream(FileStream& stream, const std::vector<Chunk>::iterator& chunk, const std::vector<Chunk>::iterator& end) : _stream(stream), _chunk(chunk), _end(end)
        {
            for (auto it = _chunk; it != end; it++)
                _length += it->size;
            if (_length > 0)
                _stream.set_position(_chunk->offset);
        }

        const std::vector<Chunk>::iterator& chunk() { return _chunk; }
        const std::vector<Chunk>::iterator& end() { return _end; }
        int64_t length() { return _length; }
        int64_t position() { return _pos; }

        void set_position(int64_t value)
        {
            if (value < _pos)
                throw std::out_of_range("Invalid position.");
            if (value >= _length)
            {
                _chunk = _end;
                _pos = _length;
                _chunk_pos = 0;
            }
            for (_chunk_pos += value - _pos; _chunk != _end && _chunk_pos >= _chunk->size; _chunk++)
                _chunk_pos -= _chunk->size;
            if (_chunk != _end)
                _stream.set_position(_chunk->offset + _chunk_pos);
            _pos = value;
        }

        size_t read(char* buffer, size_t count)
        {
            if (_pos + (int64_t)count > _length)
                count = (size_t)(_length - _pos);
            size_t res = count;
            while (res > 0)
            {
                size_t to_read = res;
                if (to_read > (size_t)(_chunk->size - _chunk_pos))
                    to_read = (size_t)(_chunk->size - _chunk_pos);
                if (_stream.read(buffer, to_read) != to_read)
                    throw container_error("Unexpected end of stream.");
                buffer += to_read;
                _pos += to_read;
                _chunk_pos += to_read;
                res -= to_read;
                if (_chunk_pos == _chunk->size)
                {
                    _chunk++;
                    _chunk_pos = 0;
                    if (_chunk != _end)
                        _stream.set_position(_chunk->offset);
                }
            }
            return count;
        }

    private:
        FileStream& _stream;
        std::vector<Chunk>::iterator _chunk;
        std::vector<Chunk>::iterator _end;
        int64_t _length = 0;
        int64_t _pos = 0;
        int64_t _chunk_pos = 0;
    };

    void FileContainer::Reader::set_position(int64_t value)
    {
        if (value < _pos)
            throw std::out_of_range("Invalid position.");
        int64_t skip = value - _pos;
        while (skip > 0)
        {
            int64_t pos = _streams.back()->position();
            _streams.back()->set_position(pos + skip);
            skip -= _streams.back()->position() - pos;
            if (skip > 0)
            {
                auto it = dynamic_cast<RawReadStream*>(_streams.front().get())->chunk();
                if (it != dynamic_cast<RawReadStream*>(_streams.front().get())->end())
                    throw container_error("Codec error.");
                init_streams(it);
                if (_streams.empty())
                    break;
            }
        }
        _pos = value - skip;
    }

    size_t FileContainer::Reader::read(char* buffer, size_t count)
    {
        size_t to_read = count;
        while (to_read > 0 && !_streams.empty())
        {
            size_t res = _streams.back()->read(buffer, to_read);
            buffer += res;
            _pos += res;
            to_read -= res;
            if (to_read > 0)
            {
                auto it = dynamic_cast<RawReadStream*>(_streams.front().get())->chunk();
                if (it != dynamic_cast<RawReadStream*>(_streams.front().get())->end())
                    throw container_error("Codec error.");
                init_streams(it);
            }
        }
        return count - to_read;
    }

    void FileContainer::Reader::init_streams(const std::vector<Chunk>::iterator& begin)
    {
        while (!_streams.empty())
            _streams.pop_back();
        std::vector<Chunk>& chunks = _container._streams[_stream_id].chunks;
        if (begin == chunks.end())
            return;
        auto end = begin;
        for (end++; end != chunks.end() && (!end->codec || end->codec->root().is_null()); end++);
        _streams.emplace_back(new RawReadStream(*_container._stream, begin, end));
        if (begin->codec && (begin->codec->root().is_string() || begin->codec->root().is_object()))
            add_codec(begin->codec->root());
        if (begin->codec && begin->codec->root().is_array())
            for (auto& codec : begin->codec->root().values())
                add_codec(codec);
    }

    void FileContainer::Reader::add_codec(const JSON::Node& codec)
    {
        if (codec.is_string() && codec.value_string() == "base64")
            _streams.emplace_back(new Base64DecoderReadStream(*_streams.back()));
    }

    class FileContainer::FileReader : public ReadStream
    {
    public:
        FileReader(FileDesc* file, Reader& reader) : _file(file), _reader(reader)
        {
            try
            {
                if (_reader.position() != file->offset)
                    _reader.set_position(file->offset);
            }
            catch (const container_error&)
            {
                _reader.container().on_corrupted(_file, _reader.stream_id());
                throw;
            }
            if (!_file->md5.empty())
                MD5Init(&md5_context);
        }

        int64_t length() override { return _file->size; }
        int64_t position() override { return _pos; }

        void set_position(int64_t value) override
        {
            if (value < _pos)
                throw std::out_of_range("Invalid position.");
            if (value > _file->size)
                value = _file->size;
            int64_t pos = _reader.position();
            try
            {
                _reader.set_position(pos + value - _pos);
            }
            catch (const container_error&)
            {
                _reader.container().on_corrupted(_file, _reader.stream_id());
                throw;
            }
            _pos += _reader.position() - pos;
        }

        size_t read(char* buffer, size_t count) override
        {
            if (_pos >= _file->size)
                return 0;
            if (_pos + (int64_t)count > _file->size)
                count = (size_t)(_file->size - _pos);
            try
            {
                if (_reader.read(buffer, count) != count)
                    throw container_error("Unexpected end of stream.");
            }
            catch (const container_error&)
            {
                _reader.container().on_corrupted(_file, _reader.stream_id());
                throw;
            }
            _pos += count;

            if (!_file->md5.empty())
            {
                MD5Update(&md5_context, (unsigned char *)buffer, (unsigned int)count);
                if (_pos == _file->size)
                {
                    std::vector<unsigned char> digest(16);
                    MD5Final(digest.data(), &md5_context);
                    char md5hash[33];
                    for (int i = 0; i < 16; i++)
                        snprintf(md5hash + i*2, 3, "%02x", digest[i]);

                    if (_file->md5 != md5hash)
                    {
                        _reader.container().on_corrupted(_file, _reader.stream_id());
                        throw container_error("MD5 hash mismatch.");
                    }
                }
            }

            return count;
        }

    private:
        FileDesc* _file;
        Reader& _reader;
        int64_t _pos = 0;
        MD5_CTX md5_context;
    };

    int FileContainer::open()
    {
        std::string st;
        JSON json;
        int error = container_error::OK;

        if (!_stream || _stream->is_closed())
            return container_error::EMPTY;

        _last_pos = _stream->position();
        if (!_stream->can_read() ||
            !_stream->readline(st, 32) ||
            !json.parse(st) ||
            !json.root().is_array() || json.root().size() < 3 ||
            !json.root()[0].is_string() || (json.root()[0].value_string() != "PK" && json.root()[0].value_string() != "7f" && json.root()[0].value_string() != "\xD0\x9A") ||
            !json.root()[1].is_int() || json.root()[1].value_int() != 1)
            return container_error::EMPTY;

        while (true)
        {
            st.clear();
            if (!_stream->readline(st, 1023) ||
                !json.parse(st) ||
                !json.root().is_array() || json.root().size() < 2 ||
                !json.root()[0].is_int() ||
                !json.root()[1].is_int())
                return container_error::UNEXPECTED_END;

            int64_t chunk_size = json.root()[0].value_int();
            int64_t next_chunk_pos = _stream->position() + chunk_size;
            if (next_chunk_pos > _stream->length())
                return container_error::UNEXPECTED_END;

            int64_t stream_id = json.root()[1].value_int();
            if (stream_id == 0)
            {
                if (chunk_size == 0)
                    break;
                st.resize((size_t)chunk_size);
                _stream->read(st.data(), st.size());
                if (_stream->position() != next_chunk_pos ||
                    !json.parse(st))
                    return container_error::INDEX_CORRUPTED;
                if (json.root().is_array())
                    for (auto& file : json.root().values())
                        if (file.is_object() && file.contains("file") && file["file"].is_string())
                        {
                            auto& entry = _index[file["file"].value_string()];
                            if (file.contains("stream") && file["stream"].is_int())
                            {
                                stream_id = file["stream"].value_int();
                                if (_next_stream_id <= stream_id)
                                    _next_stream_id = stream_id + 1;
                                if (entry.stream > 0)
                                {
                                    auto& stream = _streams[entry.stream];
                                    for (auto it = stream.files.begin(); it != stream.files.end(); )
                                        if (*it == &entry)
                                            it = stream.files.erase(it);
                                        else
                                            it++;
                                }
                                entry.reset();
                                entry.stream = stream_id;
                                if (stream_id > 0)
                                    _streams[entry.stream].files.push_back(&entry);
                            }
                            if (file.contains("size") && file["size"].is_int())
                                entry.size = file["size"].value_int();
                            if (file.contains("offset") && file["offset"].is_int())
                                entry.offset = file["offset"].value_int();
                            if (file.contains("md5") && file["md5"].is_string())
                                entry.md5 = file["md5"].value_string();
                        }
                _last_pos = next_chunk_pos;
            }
            else
            {
                if (_next_stream_id <= stream_id)
                    _next_stream_id = stream_id + 1;
                auto& stream = _streams[stream_id];
                stream.chunks.emplace_back();
                stream.chunks.back().offset = _stream->position();
                stream.chunks.back().size = chunk_size;
                if (json.root().size() >= 3)
                    stream.chunks.back().codec.reset(new JSON(std::move(json), std::move(json.root()[2])));
                _stream->set_position(next_chunk_pos);
            }
        }

        for (auto&[id, stream] : _streams)
        {
            std::sort(stream.files.begin(), stream.files.end(), [](FileDesc*& a, FileDesc*& b) { return a->offset < b->offset; });
            int64_t offset = 0;
            for (auto it = stream.files.begin(); it != stream.files.end(); )
                if ((*it)->stream != id || (*it)->size < 0 || ((*it)->size > 0 && (*it)->offset < offset))
                {
                    on_corrupted(*it, 0);
                    it = stream.files.erase(it);
                    error = container_error::INDEX_CORRUPTED;
                }
                else
                {
                    if ((*it)->offset > 16777215)
                        error = container_error::INDEX_CORRUPTED;
                    if ((*it)->size > 0)
                        offset = (*it)->offset + (*it)->size;
                    it++;
                }
        }

        if (_stream->position() < _stream->length())
            return container_error::EXCESS_DATA;

        return error;
    }

    void FileContainer::on_corrupted(FileDesc* file, int64_t stream_id)
    {
        if (_cur_file == file)
            _cur_file = nullptr;
        if (stream_id > 0)
        {
            auto& stream = _streams[stream_id];
            for (auto it = stream.files.begin(); it != stream.files.end(); )
                if (*it == file)
                    it = stream.files.erase(it);
                else
                    it++;
        }
        auto ind_it = _index.lower_bound(*file);
        _corrupted.push_back(std::move(*file));
        _index.erase(ind_it);
    }

    void FileContainer::read_file(FileDesc* file)
    {
        if (file->data)
        {
            try
            {
                file->data->set_position(0);
                return;
            }
            catch (const container_error&)
            {
            }
        }
        if (file->size == 0)
        {
            file->data.reset(new MemoryStream("", 0));
            return;
        }
        if (_cur_file != nullptr)
        {
            if (dynamic_cast<FileReader*>(_cur_file->data.get()) != nullptr)
                _cur_file->data.reset();
            _cur_file = nullptr;
        }
        if (_cur_reader && (_cur_reader->stream_id() != file->stream || _cur_reader->position() > file->offset))
            _cur_reader.reset();
        if (!_cur_reader)
            _cur_reader.reset(new Reader(*this, file->stream));
        auto& stream = _streams[file->stream];
        auto pred = [](FileDesc*& a, const int64_t& b) { return a->offset < b; };
        for (auto it = std::lower_bound(stream.files.begin(), stream.files.end(), _cur_reader->position(), pred); it != stream.files.end(); )
            if ((*it) == file)
            {
                file->data.reset(new FileReader(file, *_cur_reader));
                _cur_file = file;
                return;
            }
            else
            {
                try
                {
                    std::unique_ptr<FileReader> tmp(new FileReader(*it, *_cur_reader));
                    (*it)->data.reset(new MemoryStream((ReadStream&)*tmp, (size_t)(*it)->size));
                    it++;
                }
                catch (const container_error&)
                {
                    it = std::lower_bound(stream.files.begin(), stream.files.end(), _cur_reader->position(), pred);
                }
            }
    }

    void FileContainer::close()
    {
        _stream = nullptr;
        _file.reset();
    }

    void FileContainer::reopen(bool read, bool write)
    {
        if (_filename.empty())
            return;
        if (_stream != nullptr)
        {
            if (_stream->can_read() == read && _stream->can_write() == write)
                return;
            close();
        }
        _file.reset(new FileStream(_filename, read, write));
        if (_file->is_closed())
        {
            _file.reset();
            return;
        }
        _stream = _file.get();
        if (_filesize != _stream->length())
        {
            _filesize = _stream->length();
            _error = open();
        }
    }
}
