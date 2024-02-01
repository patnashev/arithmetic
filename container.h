#pragma once

#include <memory>
#include <vector>
#include <string>
#include <variant>
#include <functional>
#include <stdexcept>
#include <deque>
#include <set>
#include <map>
#include <algorithm>

namespace container
{
    // JSON

    class JSON
    {
    public:
        class Node
        {
        public:
            enum Type : size_t { Null = 0, Value = 1, Array = 2, Object = 3 };

        public:
            Node() : _value(std::in_place_index<Type::Null>) { }

            Type type() const { return (Type)_value.index(); }
            bool is_null() const { return _value.index() == Type::Null; }
            bool is_bool() const;
            bool is_int() const;
            bool is_double() const;
            bool is_string() const;
            bool is_array() const { return _value.index() == Type::Array; }
            bool is_object() const { return _value.index() == Type::Object; }

            bool value_bool() const;
            int64_t value_int() const;
            double value_double() const;
            std::string value_string() const;
            const std::vector<Node>& values() const { if (is_array()) return array(); else if (is_object()) return flat_map().second; else throw std::bad_variant_access(); }
            size_t size() const { return values().size(); }
            bool empty() const { return values().empty(); }
            Node& operator [] (size_t index) const { return (Node&)values()[index]; }
            const std::vector<std::string>& keys() const { return flat_map().first; }
            bool contains(const std::string& key) const { auto& _keys = flat_map().first; auto it = std::lower_bound(_keys.begin(), _keys.end(), key); return it != _keys.end() && *it == key; }
            Node& operator [] (const std::string& key) const { auto& _keys = flat_map().first; auto it = std::lower_bound(_keys.begin(), _keys.end(), key); if (it == _keys.end() || *it != key) throw std::out_of_range(key); return (Node&)flat_map().second[it - _keys.begin()]; }

            void to_string(std::string& res) const;

            Node& push_back() { if (is_null()) _value = std::vector<Node>(); return array().emplace_back(); }
            Node& insert(size_t index) { if (is_null()) _value = std::vector<Node>(); return *array().emplace(array().begin() + index); }
            template<class T>
            Node& insert(T&& key) { if (is_null()) { _value = std::pair<std::vector<std::string>, std::vector<Node>>(); } auto& _keys = flat_map().first; auto& _values = flat_map().second; auto it = std::lower_bound(_keys.begin(), _keys.end(), key); size_t i = it - _keys.begin(); if (it == _keys.end() || *it != key) { _keys.insert(it, std::forward<T>(key)); _values.insert(_values.begin() + i, Node()); }; return (Node&)_values[i]; }

        private:
            friend class JSON;
            bool parse(const char*& buffer, size_t& count);
            std::vector<Node>& array() { return std::get<Type::Array>(_value); }
            std::pair<std::vector<std::string>, std::vector<Node>>& flat_map() { return std::get<Type::Object>(_value); }

            std::pair<const char*, size_t> value_buffer() const { return std::get<Type::Value>(_value); }
            const std::vector<Node>& array() const { return std::get<Type::Array>(_value); }
            const std::pair<std::vector<std::string>, std::vector<Node>>& flat_map() const { return std::get<Type::Object>(_value); }

        private:
            std::variant<nullptr_t, std::pair<const char*, size_t>, std::vector<Node>, std::pair<std::vector<std::string>, std::vector<Node>>> _value;
        };

    public:
        JSON() { }
        template<class T>
        JSON(T&& json) { parse(std::forward<T>(json)); }
        JSON(JSON&& a, Node&& root) : _buffers(std::move(a._buffers)), _root(std::move(root)) { }

        template<class T>
        bool parse(T&& json)
        {
            _buffers.clear();
            return set_json(_root, std::forward<T>(json));
        }

        std::string to_string()
        {
            std::string res;
            _root.to_string(res);
            return res;
        }

        Node& root() { return _root; }

        void set_value(Node& node, nullptr_t value);
        void set_value(Node& node, bool value);
        void set_value(Node& node, int64_t value);
        void set_value(Node& node, double value);
        void set_value(Node& node, const std::string& value);
        template<class T>
        bool set_json(Node& node, T&& json)
        {
            _buffers.push_back(std::forward<T>(json));
            return parse_node(node);
        }

    private:
        bool parse_node(Node& node);

    private:
        std::deque<std::string> _buffers;
        Node _root;
    };

    // Streams

    class ReadStream
    {
    public:
        virtual ~ReadStream() { }
        virtual int64_t length() = 0;
        virtual int64_t position() = 0;
        virtual void set_position(int64_t value) = 0;
        virtual size_t read(char* buffer, size_t count) = 0;
    };

    class WriteStream
    {
    public:
        virtual ~WriteStream() { }
        virtual void set_length(int64_t value) = 0;
        virtual void write(const char* buffer, size_t count) = 0;
        virtual void flush() = 0;
        virtual void close() = 0;
    };

    class MemoryStream : public ReadStream, public WriteStream
    {
    public:
        MemoryStream(bool read = true, bool write = false) : _read(read), _write(write) { }
        template<class T>
        MemoryStream(T&& data, bool read = true, bool write = false) : _data(std::forward<T>(data)), _read(read), _write(write), _length((int64_t)_data.size()) { }
        MemoryStream(const char* data, size_t count, bool read = true, bool write = false) : _data(data, data + count), _read(read), _write(write), _length((int64_t)_data.size()) { }
        MemoryStream(ReadStream& stream, size_t count) : _data(count), _read(true), _write(false) { _length = (int64_t)stream.read(_data.data(), count); _data.resize((size_t)_length); }

        char* data() { return _data.data(); }
        bool can_read() { return _read; }
        int64_t length() override { return _length; }
        int64_t position() override { return _pos; }
        void set_position(int64_t value) override;
        size_t read(char* buffer, size_t count) override;

        bool can_write() { return _write; }
        void set_length(int64_t value) override;
        void write(const char* buffer, size_t count) override;
        void flush() override { }
        void close() override { _data.clear(); }

    private:
        std::vector<char> _data;
        bool _read;
        bool _write;
        int64_t _length = -1;
        int64_t _pos = 0;
    };

    class LimitedReadStream : public ReadStream
    {
    public:
        LimitedReadStream(ReadStream& stream, size_t count) : _stream(stream), _length((int64_t)count) { }

        int64_t length() override { return _length; }
        int64_t position() override { return _pos; }
        void set_position(int64_t value) override;
        size_t read(char* buffer, size_t count) override;

    private:
        ReadStream& _stream;
        int64_t _length;
        int64_t _pos = 0;
    };

    class FileStream : public ReadStream, public WriteStream
    {
    public:
        FileStream(const std::string& filename, bool read = true, bool write = false);
        FileStream(void* filestream, bool destroy = false, bool read = true, bool write = false);
        FileStream(const FileStream&) = delete;

        bool can_read() { return _read; }
        int64_t length() override { return _length; }
        int64_t position() override { return _pos; }
        void set_position(int64_t value) override;
        size_t read(char* buffer, size_t count) override;
        bool readline(std::string& res, int max_chars = -1);

        bool can_write() { return _write; }
        void set_length(int64_t value) override;
        void write(const char* buffer, size_t count) override;
        void flush() override;
        void close() override;

    private:
        void* _stream;
        bool _destroy;
        bool _read;
        bool _write;
        int64_t _length = -1;
        int64_t _pos = 0;
    };

    // Container

    struct FileDesc
    {
        std::string name;
        int64_t size = -1;
        int64_t stream = -1;
        int64_t offset = 0;
        std::string md5;

        std::unique_ptr<ReadStream> data;

        template<class T>
        FileDesc(T&& name_) : name(std::forward<T>(name_)) { }

        bool operator < (const FileDesc& right) const { return name < right.name; }
        bool operator < (const std::string& right) const { return name < right; }

        void reset()
        {
            size = -1;
            stream = -1;
            offset = 0;
            md5.clear();
            data.reset();
        }
    };

    class Index : public std::set<FileDesc, std::less<>>
    {
    public:
        Index() { }

        bool contains(const std::string& filename) { auto it = lower_bound<std::string>(filename); return it != end() && it->name == filename; }
        template<class T>
        FileDesc& operator [](T&& filename) { auto it = lower_bound<std::string>(filename); if (it == end() || it->name != filename) { it = emplace_hint(it, std::forward<T>(filename)); } return (FileDesc&)*it; }
    };

    class container_error : public std::runtime_error
    {
    public:
        static const int OK = 0;
        static const int EMPTY = 1;
        static const int EXCESS_DATA = 2;
        static const int UNEXPECTED_END = 3;
        static const int INDEX_CORRUPTED = 4;

    public:
        container_error() : std::runtime_error("Container error.") { }
        container_error(const std::string& message) : std::runtime_error(message) { }
    };

    class FileContainer;

    class ChunkedWriteStream : public WriteStream
    {
    public:
        static int MAX_BUFFER_SIZE; // 16777200
    public:
        ChunkedWriteStream(int64_t id, WriteStream& stream, size_t buffer_size = MAX_BUFFER_SIZE) : _id(id), _stream(stream), _buffer_size(buffer_size) { }

        template<class T>
        void set_codec(T&& json) { _codec_json = std::forward<T>(json); }
        void set_ñontainer(FileContainer* ñontainer) { _ñontainer = ñontainer; }

        void set_length(int64_t value) override;
        void write(const char* buffer, size_t count) override;
        void flush() override;
        void close() override;
        void write_chunk();

    private:
        int64_t _id;
        WriteStream& _stream;
        size_t _buffer_size;
        std::vector<char> _buffer;
        std::string _codec_json;
        FileContainer* _ñontainer;
    };

    class Base64CoderStream : public WriteStream
    {
    public:
        Base64CoderStream(WriteStream& stream) : _stream(stream) { }

        void set_length(int64_t value) override;
        void write(const char* buffer, size_t count) override;
        void flush() override;
        void close() override;

    private:
        void write_last();

    private:
        WriteStream& _stream;
        std::vector<char> _buffer;
        int64_t _length = -1;
        int64_t _pos = 0;
        char _next = 0;
        char _state = 0;
    };

    class Base64DecoderReadStream : public ReadStream
    {
    public:
        Base64DecoderReadStream(ReadStream& stream) : _stream(stream) { }

        int64_t length() override { return -1; }
        int64_t position() override { return _pos; }
        void set_position(int64_t value) override;
        size_t read(char* buffer, size_t count) override;

    private:
        ReadStream& _stream;
        std::vector<char> _buffer;
        size_t _buffer_pos = 0;
        int64_t _pos = 0;
        char _next = 0;
        char _state = 0;
    };

    class Packer
    {
    public:
        class Writer
        {
        public:
            Writer(Packer& packer, int64_t id);
            Writer(const Writer&) = delete;
            ~Writer() { if (_id) close(); }

            void add_codec(const std::string& codec);
            void clear_codecs();

            void add_file(const std::string& filename, const char* buffer, size_t count);
            WriteStream* add_file(const std::string& filename, int64_t size = -1);

            void close();

        private:
            class Index; friend class Index;
            class Stream; friend class Stream;
            WriteStream* add_file(FileDesc* file);

        private:
            Packer& _packer;
            int64_t _id;
            std::vector<std::unique_ptr<WriteStream>> _streams;
            int64_t _length = 0;
        };
        friend class Writer;

    public:
        Packer(const std::string& filename);
        Packer(WriteStream* stream);
        Packer(FileContainer& ñontainer);
        Packer(const Packer&) = delete;
        ~Packer() { if (_stream) close(); }

        Writer& add_writer() { _writers.emplace_back(new Writer(*this, _next_id++)); return *_writers.back(); }
        void close();

        Index& index() { return _index; }
        bool md5 = false;

    protected:
        std::unique_ptr<FileStream> _file;
        WriteStream* _stream;
        int64_t _next_id = 1;
        Index _index;
        std::vector<std::unique_ptr<Writer>> _writers;
        FileContainer* _ñontainer = nullptr;
    };

    class Unpacker : public WriteStream
    {
    public:
        Unpacker(std::function<WriteStream*(FileDesc&)> on_file) : _on_file(on_file) { }

        Index& index() { return _index; }

    private:
        std::function<WriteStream*(FileDesc&)> _on_file;
        Index _index;
    };

    class FileContainer
    {
    private:
        struct Chunk
        {
            int64_t offset;
            int64_t size;
            std::unique_ptr<JSON> codec;
        };
        struct ChunkStream
        {
            std::vector<Chunk> chunks;
            std::vector<FileDesc*> files;
        };

        class Reader : public ReadStream
        {
        public:
            class RawReadStream;

        public:
            Reader(FileContainer& ñontainer, int64_t stream_id) : _ñontainer(ñontainer), _stream_id(stream_id) { init_streams(_ñontainer._streams[_stream_id].chunks.begin()); }

            FileContainer& ñontainer() { return _ñontainer; }
            int64_t stream_id() { return _stream_id; }
            int64_t length() override { return -1; }
            int64_t position() override { return _pos; }

            void set_position(int64_t value) override;
            size_t read(char* buffer, size_t count) override;

        private:
            void init_streams(const std::vector<Chunk>::iterator& begin);
            void add_codec(const JSON::Node& codec);

        private:
            FileContainer& _ñontainer;
            int64_t _stream_id;
            std::vector<std::unique_ptr<ReadStream>> _streams;
            int64_t _pos = 0;
        };

    public:
        FileContainer(const std::string& filename, bool read = true, bool write = false) : _file(new FileStream(filename, read, write)), _stream(_file.get()) { _error = open(); }
        FileContainer(FileStream* stream) : _stream(stream) { _error = open(); }
        ~FileContainer() { close(); }

        int error() { return _error; }
        Index& index() { return _index; }

        FileDesc* read_file(const std::string& filename) { if (!_index.contains(filename)) return nullptr; FileDesc* file = &_index[filename]; read_file(file); return file; }
        void read_file(FileDesc* file);
        void close();

    private:
        class FileReader; friend class FileReader;
        friend class Packer;
        friend class ChunkedWriteStream;

        int open();
        void on_corrupted(FileDesc* file, int64_t stream_id);

    private:
        std::unique_ptr<FileStream> _file;
        FileStream* _stream;

        int _error = container_error::OK;
        Index _index;
        std::vector<FileDesc> _corrupted;

        std::map<int64_t, ChunkStream> _streams;
        int64_t _next_stream_id = 1;
        int64_t _last_pos;

        std::unique_ptr<Reader> _cur_reader;
        FileDesc* _cur_file = nullptr;
    };
}
