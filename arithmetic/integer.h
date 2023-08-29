#pragma once

#include <vector>
#include <iterator>
#include <memory>

namespace arithmetic
{
    uint32_t gcd(uint32_t a, uint32_t b);
    inline int gcd(int a, int b) { return (int)gcd((uint32_t)a, (uint32_t)b); }

    bool is_prime(uint32_t a);
    inline bool is_prime(int a) { return is_prime((uint32_t)a); }

    uint32_t phi(uint32_t N);
    inline int phi(int N) { return (int)phi((uint32_t)N); }

    int jacobi(uint32_t a, uint32_t b);
    inline int jacobi(int a, int b) { return jacobi((uint32_t)a, (uint32_t)b); }

    int kronecker(uint32_t a, uint32_t b);
    inline int kronecker(int a, int b) { return kronecker((uint32_t)a, (uint32_t)b); }

    class PrimeIterator;

    class PrimeList
    {
        friend class PrimeIterator;

    public:
        PrimeList(int max);

        void sieve_range(int start, int end, std::vector<int>& list);

        size_t size() const { return primes.size(); }
        int operator[] (size_t pos) const { return primes[pos]; }

        PrimeIterator begin();

        static PrimeList& primes_16bit() { if (!_list65536) _list65536.reset(new PrimeList(65536)); return *_list65536; }

    private:
        std::vector<int> primes;

        static std::unique_ptr<PrimeList> _list65536;
    };

    class PrimeIterator : public std::iterator<std::input_iterator_tag, int, int, const int*, int>
    {
    public:
        PrimeIterator(PrimeList& list) : _list(list) { }
        PrimeIterator(const PrimeIterator& it) : _list(it._list), _cur(it._cur), _range(it._range), _range_pos(it._range_pos) { }
        PrimeIterator& operator=(const PrimeIterator& it)
        {
            _list = it._list; _cur = it._cur; _range = it._range; _range_pos = it._range_pos;
            return *this;
        }

        static PrimeIterator get() { return PrimeIterator(PrimeList::primes_16bit()); }
        static const PrimeIterator& max();

        void sieve_range(uint64_t start, uint64_t end, std::vector<uint64_t>& list);

        PrimeIterator& operator++();
        void operator++(int);
        PrimeIterator& operator+=(int offset);
        bool operator==(const PrimeIterator& other) const { return _cur == other._cur; }
        bool operator!=(const PrimeIterator& other) const { return !(*this == other); }
        int operator*() const { return _cur < _list.size() ? _list[_cur] : _range[_cur - _range_pos]; }
        size_t pos() const { return _cur; }

    private:
        PrimeList& _list;
        size_t _cur = 0;
        std::vector<int> _range;
        size_t _range_pos = 0;
    };
}
