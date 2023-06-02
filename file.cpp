
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "gwnum.h"
#include "file.h"
#include "md5.h"
#include "inputnum.h"
#include "task.h"
#ifdef _WIN32
#include "windows.h"
#endif

int File::FILE_APPID = 1;

void Writer::write(int32_t value)
{
    _buffer.insert(_buffer.end(), (char*)&value, 4 + (char*)&value);
}

void Writer::write(uint32_t value)
{
    _buffer.insert(_buffer.end(), (char*)&value, 4 + (char*)&value);
}

void Writer::write(uint64_t value)
{
    _buffer.insert(_buffer.end(), (char*)&value, 8 + (char*)&value);
}

void Writer::write(double value)
{
    _buffer.insert(_buffer.end(), (char*)&value, (char*)(&value + 1));
}

void Writer::write(const std::string& value)
{
    int32_t len = (int)value.size();
    _buffer.insert(_buffer.end(), (char*)&len, 4 + (char*)&len);
    _buffer.insert(_buffer.end(), (char*)value.data(), (char*)(value.data() + value.size()));
}

void Writer::write(const arithmetic::Giant& value)
{
    int32_t len = value.size();
    if (value < 0)
        len *= -1;
    _buffer.insert(_buffer.end(), (char*)&len, 4 + (char*)&len);
    _buffer.insert(_buffer.end(), (char*)value.data(), (char*)(value.data() + value.size()));
}

void Writer::write(const arithmetic::SerializedGWNum& value)
{
    int32_t len = value.size();
    _buffer.insert(_buffer.end(), (char*)&len, 4 + (char*)&len);
    _buffer.insert(_buffer.end(), (char*)value.data(), (char*)(value.data() + value.size()));
}

void Writer::write(const char* ptr, int count)
{
    _buffer.insert(_buffer.end(), ptr, ptr + count);
}

void Writer::write_text(const char* ptr)
{
    size_t count = strlen(ptr);
    _buffer.insert(_buffer.end(), ptr, ptr + count);
}

void Writer::write_text(const std::string& value)
{
    _buffer.insert(_buffer.end(), (char*)value.data(), (char*)(value.data() + value.length()));
}

void Writer::write_textline(const char* ptr)
{
    write_text(ptr);
    _buffer.insert(_buffer.end(), {'\r', '\n'});
}

void Writer::write_textline(const std::string& value)
{
    write_text(value);
    _buffer.insert(_buffer.end(), { '\r', '\n' });
}

std::vector<char> Writer::hash()
{
    std::vector<char> digest(16);
    MD5_CTX context;
    MD5Init(&context);
    MD5Update(&context, (unsigned char *)_buffer.data(), (unsigned int)_buffer.size());
    MD5Final((unsigned char *)digest.data(), &context);
    return digest;
}

std::string Writer::hash_str()
{
    char md5hash[33];
    md5_raw_input(md5hash, (unsigned char *)_buffer.data(), (unsigned int)_buffer.size());
    return md5hash;
}

bool Reader::read(int32_t& value)
{
    if (_size < _pos + 4)
        return false;
    value = *(int32_t*)(_data + _pos);
    _pos += 4;
    return true;
}

bool Reader::read(uint32_t& value)
{
    if (_size < _pos + 4)
        return false;
    value = *(uint32_t*)(_data + _pos);
    _pos += 4;
    return true;
}

bool Reader::read(uint64_t& value)
{
    if (_size < _pos + 8)
        return false;
    value = *(uint64_t*)(_data + _pos);
    _pos += 8;
    return true;
}

bool Reader::read(double& value)
{
    if (_size < _pos + 4)
        return false;
    value = *(double*)(_data + _pos);
    _pos += sizeof(double);
    return true;
}

bool Reader::read(std::string& value)
{
    if (_size < _pos + 4)
        return false;
    int len = *(int32_t*)(_data + _pos);
    _pos += 4;
    if (_size < _pos + len)
        return false;
    value.insert(value.end(), _data + _pos, _data + _pos + len);
    _pos += len;
    return true;
}

bool Reader::read(arithmetic::Giant& value)
{
    if (_size < _pos + 4)
        return false;
    int len = *(int32_t*)(_data + _pos);
    _pos += 4;
    if (_size < _pos + abs(len)*4)
        return false;
    value.arithmetic().init((uint32_t*)(_data + _pos), abs(len), value);
    if (len < 0)
        value.arithmetic().neg(value, value);
    _pos += abs(len)*4;
    return true;
}

bool Reader::read(arithmetic::SerializedGWNum& value)
{
    if (_size < _pos + 4)
        return false;
    int len = *(int32_t*)(_data + _pos);
    _pos += 4;
    if (_size < _pos + len*sizeof(uint32_t))
        return false;
    value.init((uint32_t*)(_data + _pos), len);
    _pos += len*sizeof(uint32_t);
    return true;
}

bool TextReader::read_textline(std::string& value)
{
    int i;
    if (_pos >= _size)
        return false;
    for (i = _pos + 1; i < _size && _data[i - 1] != '\n'; i++);
    if (i - 2 >= _pos && _data[i - 2] == '\r' && _data[i - 1] == '\n')
        value = std::string(_data + _pos, i - 2 - _pos);
    else if (_data[i - 1] == '\n')
        value = std::string(_data + _pos, i - 1 - _pos);
    else
        value = std::string(_data + _pos, i - _pos);
    _pos = i;
    return true;
}

File* File::add_child(const std::string& name, uint32_t fingerprint)
{
    _children.emplace_back(new File(_filename + "." + name, fingerprint));
    _children.back()->hash = hash;
    return _children.back().get();
}

uint32_t File::unique_fingerprint(uint32_t fingerprint, const std::string& unique_id)
{
    unsigned char digest[16];
    MD5_CTX context;
    MD5Init(&context);
    MD5Update(&context, (unsigned char *)&fingerprint, 4);
    MD5Update(&context, (unsigned char *)unique_id.data(), (unsigned int)unique_id.size());
    MD5Final(digest, &context);

    return *(uint32_t*)digest;
}

Writer* File::get_writer()
{
    Writer* writer = new Writer(std::move(_buffer));
    writer->buffer().clear();
    return writer;
}

Writer* File::get_writer(char type, char version)
{
    Writer* writer = get_writer();
    writer->write(MAGIC_NUM);
    writer->write(appid + (0 << 8) + ((unsigned char)type << 16) + ((unsigned char)version << 24));
    if (_fingerprint != 0)
        writer->write(_fingerprint);
    return writer;
}

void File::read_buffer()
{
    if (!_buffer.empty())
        return;

    FILE* fd = fopen(_filename.data(), "rb");
    if (!fd)
    {
        _buffer.clear();
        return;
    }
    fseek(fd, 0L, SEEK_END);
    int filelen = ftell(fd);
    _buffer.resize(filelen);
    fseek(fd, 0L, SEEK_SET);
    filelen = (int)fread(_buffer.data(), 1, filelen, fd);
    fclose(fd);
    if (filelen != _buffer.size())
    {
        _buffer.clear();
        return;
    }

    if (hash)
    {
        std::string md5_filename = _filename + ".md5";
        fd = fopen(md5_filename.data(), "rb");
        if (fd)
        {
            char md5hash[33];
            fread(md5hash, 1, 32, fd);
            fclose(fd);
            md5hash[32] = 0;
            std::string saved_hash(md5hash);
            if (!saved_hash.empty())
            {
                md5_raw_input(md5hash, (unsigned char *)_buffer.data(), (unsigned int)_buffer.size());
                if (saved_hash != md5hash)
                {
                    _buffer.clear();
                    return;
                }
            }
        }
    }
}

Reader* File::get_reader()
{
    read_buffer();

    if (_buffer.size() < 8)
        return nullptr;
    if (*(uint32_t*)_buffer.data() != MAGIC_NUM)
        return nullptr;
    if (_buffer[4] != appid)
        return nullptr;
    if (_fingerprint != 0 && _buffer.size() < 12)
        return nullptr;
    if (_fingerprint != 0 && *(uint32_t*)(_buffer.data() + 8) != _fingerprint)
        return nullptr;

    return new Reader(_buffer[5], _buffer[6], _buffer[7], _buffer.data(), (int)_buffer.size(), _fingerprint != 0 ? 12 : 8);
}

TextReader* File::get_textreader()
{
    read_buffer();
    return new TextReader(_buffer.data(), (int)_buffer.size(), 0);
}

bool writeThrough(char *filename, const void *buffer, size_t count)
{
    bool ret = true;
#ifdef _WIN32
    HANDLE hFile = CreateFileA(filename, GENERIC_WRITE, FILE_SHARE_READ, NULL, CREATE_ALWAYS, FILE_FLAG_WRITE_THROUGH, NULL);
    if (hFile == INVALID_HANDLE_VALUE)
        return false;
    DWORD written;
    if (!WriteFile(hFile, buffer, (DWORD)count, &written, NULL) || written != count)
        ret = false;
    CloseHandle(hFile);
#else
    FILE *fd = fopen(filename, "wb");
    if (!fd)
        return false;
    if (fwrite(buffer, 1, count, fd) != count)
        ret = false;
    fclose(fd);
#endif
    return ret;
}

void File::commit_writer(Writer& writer)
{
    std::string new_filename = _filename + ".new";
    if (!writeThrough(new_filename.data(), writer.buffer().data(), writer.buffer().size()))
    {
        remove(new_filename.data());
        return;
    }
    remove(_filename.data());
    rename(new_filename.data(), _filename.data());

    if (hash)
    {
        std::string hash = writer.hash_str();
        writeThrough(_hash_filename.data(), hash.data(), 32);
    }

    _buffer = std::move(writer.buffer());
}

void File::free_buffer()
{
    std::vector<char>().swap(_buffer);
}

void File::clear(bool recursive)
{
    remove(_filename.data());
    std::string md5_filename = _filename + ".md5";
    remove(md5_filename.data());
    std::vector<char>().swap(_buffer);

    if (recursive)
        for (auto it = _children.begin(); it != _children.end(); it++)
            (*it)->clear(true);
}

bool File::read(TaskState& state)
{
    std::unique_ptr<Reader> reader(get_reader());
    if (!reader)
        return false;
    if (reader->type() != state.type())
        return false;
    if (!state.read(*reader))
        return false;
    state.set_written();
    return true;
}

void File::write(TaskState& state)
{
    std::unique_ptr<Writer> writer(get_writer(state.type(), state.version()));
    state.write(*writer);
    commit_writer(*writer);
    state.set_written();
}

void File::write_text(const std::string& value)
{
    std::unique_ptr<Writer> writer(get_writer());
    writer->write_text(value);
    commit_writer(*writer);
}

void File::write_textline(const std::string& value)
{
    std::unique_ptr<Writer> writer(get_writer());
    writer->write_textline(value);
    commit_writer(*writer);
}
