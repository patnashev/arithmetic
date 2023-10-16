
#include <cmath>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include "gwnum.h"
#include "cpuid.h"
#include "giant.h"
#include "arithmetic.h"
#include "exception.h"
#include "integer.h"
#ifdef GMP
#include <gmp.h>
#ifdef _WIN32
#include "windows.h"
#include "tchar.h"
#endif
#endif

namespace arithmetic
{
    GiantsArithmetic* _defaultGiantsArithmetic = nullptr;
    void GiantsArithmetic::init_default_arithmetic()
    {
        if (_defaultGiantsArithmetic != nullptr)
            return;
#ifdef GMP
#ifdef _WIN32
        if (LoadLibrary(_T("libgmp-gw1.dll")) != NULL)
#else
        if (true)
#endif
            _defaultGiantsArithmetic = new GMPArithmetic();
        else
#endif
            _defaultGiantsArithmetic = new GiantsArithmetic();
    }
    GiantsArithmetic& GiantsArithmetic::default_arithmetic()
    {
        if (_defaultGiantsArithmetic == nullptr)
            init_default_arithmetic();
        return *_defaultGiantsArithmetic;
    }
    GiantsArithmetic* GiantsArithmetic::alloc_gwgiants(void* gwdata, int capacity)
    {
#ifdef GMP
        if (_defaultGiantsArithmetic == nullptr || dynamic_cast<GMPArithmetic*>(_defaultGiantsArithmetic) != nullptr)
            return new GWGMPArithmetic(gwdata, capacity);
#endif
        return new GWGiantsArithmetic(gwdata, capacity);
    }

#define giant(x) ((::giant)&(x)._capacity)

    void GiantsArithmetic::alloc(Giant& a)
    {
        a._capacity = 0;
        a._size = 0;
        a._data = nullptr;
        alloc(a, 0);
    }

    void GiantsArithmetic::free(Giant& a)
    {
        ::free(a._data);
        a._capacity = 0;
        a._size = 0;
        a._data = nullptr;
    }

    void GiantsArithmetic::alloc(Giant& a, int capacity)
    {
        if (capacity < this->capacity())
            capacity = this->capacity();
        if (a._capacity >= capacity)
            return;
        a._capacity = capacity;
        a._data = (uint32_t*)realloc(a._data, capacity*sizeof(uint32_t));
    }

    void GiantsArithmetic::copy(const Giant& a, Giant& res)
    {
        if (a.empty())
        {
            free(res);
            return;
        }
        if (a._size == 0)
        {
            init(0, res);
            return;
        }
        alloc(res, abs(a._size));
        memcpy(res._data, a._data, abs(a._size)*sizeof(uint32_t));
        res._size = a._size;
    }

    void GiantsArithmetic::move(Giant&& a, Giant& res)
    {
        if (!res.empty())
            free(res);
        res._capacity = a._capacity;
        res._size = a._size;
        res._data = a._data;
        a._capacity = 0;
        a._size = 0;
        a._data = nullptr;
    }

    void GiantsArithmetic::init(int32_t a, Giant& res)
    {
        alloc(res, 1);
        itog(a, giant(res));
    }

    void GiantsArithmetic::init(uint32_t a, Giant& res)
    {
        alloc(res, 1);
        ultog(a, giant(res));
    }

    void GiantsArithmetic::init(const std::string& a, Giant& res)
    {
        alloc(res, ((int)a.length() + 8)/9);
        ctog(a.data(), giant(res));
    }

    void GiantsArithmetic::init(uint32_t* data, int size, Giant& res)
    {
        if (size == 0)
        {
            init(0, res);
            return;
        }
        alloc(res, size);
        memcpy(res._data, data, size*sizeof(uint32_t));
        res._size = size;
        while (res._size && !res.data()[res._size - 1])
            res._size--;
    }

    void GiantsArithmetic::init(const GWNum& a, Giant& res)
    {
        int capacity = a.arithmetic().state().giants->capacity();
        res.arithmetic().alloc(res, capacity);
        if (gwtogiant(a.arithmetic().gwdata(), *a, giant(res)) < 0)
            throw InvalidFFTDataException();
    }

    std::string GiantsArithmetic::to_string(const Giant& a)
    {
        if (a.empty())
            return "";
        std::vector<char> buffer(abs(a._size)*10 + 10);
        if (a._size == 0)
            buffer[0] = '0';
        else if (a._size > 0)
            gtoc(giant(a), buffer.data(), (int)buffer.size());
        else
        {
            ((Giant&)a)._size = -a._size;
            buffer[0] = '-';
            gtoc(giant(a), buffer.data() + 1, (int)buffer.size() - 1);
            ((Giant&)a)._size = -a._size;
        }
        return std::string(buffer.data());
    }

    std::string Giant::to_res64() const
    {
        std::string res(16, '0');
        if (size() > 1)
            snprintf(res.data(), 17, "%08X%08X", data()[1], data()[0]);
        else if (size() > 0)
            snprintf(res.data() + 8, 9, "%08X", data()[0]);
        return res;
    }

    void GiantsArithmetic::to_GWNum(const Giant& a, GWNum& res)
    {
        if (a._size >= 0)
            gianttogw(res.arithmetic().gwdata(), giant(a), *res);
        else
        {
            Giant tmp(*this);
            add((Giant&)a, res.arithmetic().N(), tmp);
            gianttogw(res.arithmetic().gwdata(), giant(tmp), *res);
        }
    }

    int Giant::size() const
    {
#ifdef GMP
        if (dynamic_cast<GMPArithmetic*>(&arithmetic()) != nullptr)
        {
            int size = abs(_size)*GMP_NUMB_BITS/32;
            while (size > 0 && data()[size - 1] == 0)
                size--;
            return size;
        }
#endif
        return abs(_size);
    }

    int Giant::capacity() const
    {
#ifdef GMP
        if (dynamic_cast<GMPArithmetic*>(&arithmetic()) != nullptr)
            return _capacity*GMP_NUMB_BITS/32;
#endif
        return _capacity;
    }

    int GiantsArithmetic::cmp(const Giant& a, const Giant& b)
    {
        return gcompg(giant(a), giant(b));
    }

    int GiantsArithmetic::cmp(const Giant& a, int32_t b)
    {
        if (a._size > 1 || (a._size == 1 && ::bitlen(giant(a)) == 32))
            return 1;
        if (a._size < -1 || (a._size == -1 && ::bitlen(giant(a)) == 32))
            return -1;
        if (a._size == 0)
            return 0 > b ? 1 : 0 < b ? -1 : 0;
        int32_t val = (a._size == -1 ? -1 : 1)*(int32_t)a.data()[0];
        return val > b ? 1 : val < b ? -1 : 0;
    }

    void GiantsArithmetic::shiftleft(Giant& a, int b, Giant& res)
    {
        alloc(res, abs(a._size) + (b + 31)/32);
        if (res._data != a._data)
            copy(a, res);
        gshiftleft(b, giant(res));
    }

    void GiantsArithmetic::shiftright(Giant& a, int b, Giant& res)
    {
        int size = abs(a._size) - b/32;
        alloc(res, size);
        if (res._data != a._data)
        {
            memcpy(res.data(), a.data() + b/32, size*4);
            res._size = size*(a._size < 0 ? -1 : 1);
            gshiftright(b%32, giant(res));
        }
        else
            gshiftright(b, giant(res));
    }

    int GiantsArithmetic::bitlen(const Giant& a)
    {
        return ::bitlen(giant(a));
    }

    bool GiantsArithmetic::bit(const Giant& a, int b)
    {
        if (b/32 >= abs(a._size))
            return 0;
        return bitval(giant(a), b) != 0;
    }

    void GiantsArithmetic::substr(const Giant& a, int offset, int count, Giant& res)
    {
        int size = (count + 31 + offset%32)/32;
        if (size > abs(a._size) - offset/32)
            size = abs(a._size) - offset/32;
        if (size <= 0)
        {
            res = 0;
            return;
        }
        alloc(res, size);
        memcpy(res.data(), a.data() + offset/32, size*4);
        if ((offset%32 + count)/32 == size - 1)
            res.data()[size - 1] &= (1 << (offset + count)%32) - 1;
        for (; size > 0 && res.data()[size - 1] == 0; size--);
        res._size = size;
        gshiftright(offset%32, giant(res));
    }

    double GiantsArithmetic::log2(const Giant& a)
    {
        if (a._size == 0)
            return 0;
        if (a.size() == 1)
            return std::log2(a.data()[0]);
        if (a.size() == 2)
            return std::log2(((uint64_t*)a._data)[0]);
        return bitlen(a) - 1;
    }

    size_t GiantsArithmetic::digits(const Giant& a, int base)
    {
        if (a == 0)
            return 1;
        if (base == 2)
            return bitlen(a);
        if (a.size() == 1)
            return (size_t)std::floor(std::log2(a.data()[0])/std::log2(base)) + 1;
        int log = (int)std::floor((bitlen(a) - 1)/std::log2(base));
        Giant tmp(*this);
        tmp = base;
        power(tmp, log, tmp);
        for (; tmp <= a; tmp *= base, log++);
        return log;
    }

    void GiantsArithmetic::add(Giant& a, Giant& b, Giant& res)
    {
        int size = abs(a._size) + 1;
        if (size < abs(b._size) + 1)
            size = abs(b._size) + 1;
        alloc(res, size);
        if (res._data != a._data)
            copy(a, res);
        if (a._data != b._data)
            addg(giant(b), giant(res));
        else
            gshiftleft(1, giant(res));
    }

    void GiantsArithmetic::add(Giant& a, int32_t b, Giant& res)
    {
        alloc(res, abs(a._size) + 1);
        if (res._data != a._data)
            copy(a, res);
        sladdg(b, giant(res));
    }

    void GiantsArithmetic::sub(Giant& a, Giant& b, Giant& res)
    {
        int capacity = abs(a._size) + 1;
        if (capacity < abs(b._size) + 1)
            capacity = abs(b._size) + 1;
        alloc(res, capacity);
        if (res._data != a._data)
            copy(a, res);
        subg(giant(b), giant(res));
    }

    void GiantsArithmetic::sub(Giant& a, int32_t b, Giant& res)
    {
        alloc(res, abs(a._size) + 1);
        if (res._data != a._data)
            copy(a, res);
        sladdg(-b, giant(res));
    }

    void GiantsArithmetic::neg(Giant& a, Giant& res)
    {
        if (res._data != a._data)
            copy(a, res);
        res._size = -res._size;
    }

    void GiantsArithmetic::mul(Giant& a, Giant& b, Giant& res)
    {
        alloc(res, abs(a._size) + abs(b._size));
        if (res._data != a._data)
            copy(a, res);
        if (a._data != b._data)
        {
            if ((abs(res._size) > 10000 && abs(b._size) > 10 && abs(res._size) > 4*abs(b._size)) || (abs(b._size) > 10000 && abs(res._size) > 10 && abs(b._size) > 4*abs(res._size)))
            {
                setmulmode(FFT_MUL);
                mulg(giant(b), giant(res));
                setmulmode(AUTO_MUL);
            }
            else
                mulg(giant(b), giant(res));
        }
        else
            squareg(giant(res));
    }

    void GiantsArithmetic::mul(Giant& a, int32_t b, Giant& res)
    {
        alloc(res, abs(a._size) + 1);
        if (res._data != a._data)
            copy(a, res);
        imulg(b, giant(res));
    }

    void GiantsArithmetic::mul(Giant& a, uint32_t b, Giant& res)
    {
        alloc(res, abs(a._size) + 1);
        if (res._data != a._data)
            copy(a, res);
        ulmulg(b, giant(res));
    }

    void GiantsArithmetic::div(Giant& a, Giant& b, Giant& res)
    {
        alloc(res, abs(a._size));
        if (res._data != a._data)
            copy(a, res);
        divg(giant(b), giant(res));
    }

    void GiantsArithmetic::div(Giant& a, int32_t b, Giant& res)
    {
        alloc(res, abs(a._size));
        if (res._data != a._data)
            copy(a, res);
        dbldivg(b, giant(res));
    }

    void GiantsArithmetic::div(Giant& a, uint32_t b, Giant& res)
    {
        alloc(res, abs(a._size));
        if (res._data != a._data)
            copy(a, res);
        dbldivg(b, giant(res));
    }

    void GiantsArithmetic::mod(Giant& a, Giant& b, Giant& res)
    {
        alloc(res, abs(a._size));
        if (res._data != a._data)
            copy(a, res);
        modg(giant(b), giant(res));
    }

    void GiantsArithmetic::mod(Giant& a, uint32_t b, uint32_t& res)
    {
        uint64_t res64 = 0;
        for (int i = abs(a._size) - 1; i >= 0; i--)
        {
            res64 <<= 32;
            res64 += a.data()[i];
            res64 %= b;
        }
        res = (uint32_t)res64;
    }

    void GiantsArithmetic::gcd(Giant& a, Giant& b, Giant& res)
    {
        if (a == 0)
        {
            copy(b, res);
            return;
        }
        if (b == 0)
        {
            copy(a, res);
            return;
        }
        alloc(res, abs(a._size));
        if (res._data != a._data)
            copy(a, res);
        gcdg(giant(b), giant(res));
    }

    void GiantsArithmetic::inv(Giant& a, Giant& n, Giant& res)
    {
        alloc(res, abs(n._size));
        if (res._data != a._data)
            copy(a, res);
        if (res == 0)
            throw NoInverseException(res);
        invg(giant(n), giant(res));
        if (res._size < 0)
            throw NoInverseException(res);
    }

    void GiantsArithmetic::power(Giant& a, int32_t b, Giant& res)
    {
        if (b == 0)
        {
            init(1, res);
            return;
        }
        int capacity = abs(a._size);
        if (capacity == 1)
            capacity = (int)(std::log2(a.data()[0])*b/32) + 1;
        else
            capacity *= b;
        alloc(res, capacity);
        if (res._data != a._data)
            copy(a, res);
        ::power(giant(res), b);
    }

    void GiantsArithmetic::powermod(Giant& a, Giant& b, Giant& n, Giant& res)
    {
        alloc(res, abs(n._size));
        if (res._data != a._data)
            copy(a, res);
        ::powermodg(giant(res), giant(b), giant(n));
    }

    extern "C"
    {
        struct mt_state {
            unsigned long mt[624]; /* the array for the state vector */
            int mti;
        };

        void init_genrand(struct mt_state *x, unsigned long s);

        unsigned long genrand_int32(struct mt_state *x);
    }

    GiantsArithmetic::~GiantsArithmetic()
    {
        struct mt_state *state = (struct mt_state *)_rnd_state;
        if (state != nullptr)
            delete state;
    }

    void init_by_array(struct mt_state *x, uint32_t init_key[], int key_length)
    {
        const int N = 624;
        int i, j, k;
        init_genrand(x, 19650218UL);
        i = 1; j = 0;
        k = (N>key_length ? N : key_length);
        for (; k; k--) {
            x->mt[i] = (x->mt[i] ^ ((x->mt[i-1] ^ (x->mt[i-1] >> 30)) * 1664525UL))
                + init_key[j] + j; /* non linear */
            x->mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
            i++; j++;
            if (i>=N) { x->mt[0] = x->mt[N-1]; i = 1; }
            if (j>=key_length) j = 0;
        }
        for (k = N-1; k; k--) {
            x->mt[i] = (x->mt[i] ^ ((x->mt[i-1] ^ (x->mt[i-1] >> 30)) * 1566083941UL))
                - i; /* non linear */
            x->mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
            i++;
            if (i>=N) { x->mt[0] = x->mt[N-1]; i = 1; }
        }

        x->mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
    }

    void GiantsArithmetic::rnd_seed(Giant& a)
    {
        if (_rnd_state == nullptr)
            _rnd_state = new struct mt_state;
        struct mt_state *state = (struct mt_state *)_rnd_state;
        init_by_array(state, a.data(), a.size());
    }

    void GiantsArithmetic::rnd(Giant& res, int bits)
    {
        struct mt_state *state = (struct mt_state *)_rnd_state;
        if (_rnd_state == nullptr)
        {
            _rnd_state = state = new struct mt_state;
            double seed = getHighResTimer();
            init_by_array(state, (uint32_t*)&seed, 2);
        }
        res.arithmetic().alloc(res, (bits + 31)/32);
        int i;
        for (i = 0; i < (bits + 31)/32; i++)
            res.data()[i] = genrand_int32(state);
        if ((bits & 31) != 0)
            res.data()[i - 1] &= (1 << (bits & 31)) - 1;
        while (i > 0 && res.data()[i - 1] == 0)
            i--;
        res._size = i;
    }

    int GiantsArithmetic::kronecker(Giant& a, Giant& b)
    {
        GWASSERT(0);
        return 0;
    }

    int GiantsArithmetic::kronecker(Giant& a, uint32_t b)
    {
        if (!a.bit(0) && !(b & 1))
            return 0;
        int i;
        int res = 1;
        for (i = 0; b != 0 && !(b & 1); i++, b >>= 1);
        uint32_t amod8 = a.data()[0] & 7;
        if ((i & 1) && (amod8 == 3 || amod8 == 5))
            res *= -1;
        res *= arithmetic::kronecker(a%b, b);
        if (abs(res) > 1)
            res = 0;
        return res;
    }

    int GiantsArithmetic::kronecker(uint32_t a, Giant& b)
    {
        if (!(a & 1) && !b.bit(0))
            return 0;
        int i, j;
        int res = 1;
        for (i = 0; i < b.bitlen() && !b.bit(i); i++);
        if ((i & 1) && ((a & 7) == 3 || (a & 7) == 5))
            res *= -1;
        for (j = 0; a != 0 && !(a & 1); j++, a >>= 1);
        if (i > 0)
        {
            Giant tmp(b);
            tmp >>= i;
            res *= arithmetic::kronecker(tmp%a, a);
            if ((a & 2) && tmp.bit(1))
                res *= -1;
        }
        else
        {
            uint32_t bmod8 = b.data()[0] & 7;
            if ((j & 1) && (bmod8 == 3 || bmod8 == 5))
                res *= -1;
            res *= arithmetic::kronecker(b%a, a);
            if ((a & 2) && b.bit(1))
                res *= -1;
        }
        if (abs(res) > 1)
            res = 0;
        return res;
    }

    void GWGiantsArithmetic::alloc(Giant& a)
    {
        a._capacity = capacity();
        a._size = 0;
        a._data = popg(&((gwhandle*)_gwdata)->gdata, capacity())->n;
    }

    void GWGiantsArithmetic::free(Giant& a)
    {
        pushg(&((gwhandle*)_gwdata)->gdata, 1);
        a._capacity = 0;
        a._size = 0;
        a._data = nullptr;
    }

    void GWGiantsArithmetic::alloc(Giant& a, int capacity)
    {
        if (!a.empty())
            return;
        alloc(a);
    }

    void GWGiantsArithmetic::init(const GWNum& a, Giant& res)
    {
        res.arithmetic().alloc(res, capacity());
        if (gwtogiant((gwhandle*)_gwdata, *a, giant(res)) < 0)
            throw InvalidFFTDataException();
    }

    void GWGiantsArithmetic::to_GWNum(const Giant& a, GWNum& res)
    {
        if (a._size >= 0)
            gianttogw((gwhandle*)_gwdata, giant(a), *res);
        else
        {
            Giant tmp(*this);
            add((Giant&)a, res.arithmetic().N(), tmp);
            gianttogw((gwhandle*)_gwdata, giant(tmp), *res);
        }
    }

    void GWGiantsArithmetic::inv(Giant& a, Giant& n, Giant& res)
    {
        alloc(res, abs(n._size));
        if (res._data != a._data)
            copy(a, res);
        if (res == 0)
            throw NoInverseException(res);
        invgi(&((gwhandle*)_gwdata)->gdata, 0, giant(n), giant(res));
        if (res._size < 0)
            throw NoInverseException(res);
    }

#ifdef GMP
#define mpz(x) ((mpz_ptr)&(x)._capacity)

    void* realloc_func(void* block, size_t old_size, size_t size) { return realloc(block, size); }
    void free_func(void* block, size_t) { free(block); }
    GMPArithmetic::GMPArithmetic()
    {
        mp_set_memory_functions(malloc, realloc_func, free_func);
    }

    std::string GMPArithmetic::version()
    {
        return std::to_string(__GNU_MP_VERSION) + "." + std::to_string(__GNU_MP_VERSION_MINOR) + "." + std::to_string(__GNU_MP_VERSION_PATCHLEVEL);
    }

    void GMPArithmetic::alloc(Giant& a)
    {
        a._capacity = 0;
        a._size = 0;
        a._data = nullptr;
        alloc(a, 0);
    }

    void GMPArithmetic::free(Giant& a)
    {
        ::free(a._data);
        a._capacity = 0;
        a._size = 0;
        a._data = nullptr;
    }

    void GMPArithmetic::alloc(Giant& a, int capacity)
    {
        if (capacity < this->capacity())
            capacity = this->capacity();
        capacity = (capacity*32 + GMP_NUMB_BITS - 1)/GMP_NUMB_BITS;
        if (a._capacity >= capacity)
            return;
        a._capacity = capacity;
        a._data = (uint32_t*)realloc(a._data, capacity*GMP_NUMB_BITS/8);
    }

    void GMPArithmetic::copy(const Giant& a, Giant& res)
    {
        if (a.empty())
        {
            free(res);
            return;
        }
        if (a._size == 0)
        {
            init(0, res);
            return;
        }
        int size = a.size();
        res.arithmetic().alloc(res, size);
        memcpy(res._data, a._data, size*4);
        if (GMP_NUMB_BITS == 64 && (size & 1))
            res.data()[size] = 0;
        res._size = (size*32 + GMP_NUMB_BITS - 1)/GMP_NUMB_BITS;
    }

    void GMPArithmetic::move(Giant&& a, Giant& res)
    {
        if (!res.empty())
            free(res);
        res._capacity = a._capacity;
        res._size = a._size;
        res._data = a._data;
        a._capacity = 0;
        a._size = 0;
        a._data = nullptr;
    }

    void GMPArithmetic::init(int32_t a, Giant& res)
    {
        alloc(res, 1);
        mpz_set_si(mpz(res), a);
    }

    void GMPArithmetic::init(uint32_t a, Giant& res)
    {
        alloc(res, 1);
        mpz_set_ui(mpz(res), a);
    }

    void GMPArithmetic::init(const std::string& a, Giant& res)
    {
        alloc(res, ((int)a.length() + 8)/9 + 1);
        mpz_set_str(mpz(res), a.data(), 10);
    }

    void GMPArithmetic::init(uint32_t* data, int size, Giant& res)
    {
        if (size == 0)
        {
            init(0, res);
            return;
        }
        alloc(res, size);
        mpz_import(mpz(res), size, -1, 4, 0, 0, data);
    }

    void GMPArithmetic::init(const GWNum& a, Giant& res)
    {
        int capacity = a.arithmetic().state().giants->capacity();
        res.arithmetic().alloc(res, capacity);
        capacity = res._capacity;
        res._capacity = res.capacity();
        if (gwtogiant(a.arithmetic().gwdata(), *a, giant(res)) < 0)
            throw InvalidFFTDataException();
        if (GMP_NUMB_BITS == 64 && (res._size & 1))
            res.data()[res._size] = 0;
        res._size = (res._size*32 + GMP_NUMB_BITS - 1)/GMP_NUMB_BITS;
        res._capacity = capacity;
    }

    std::string GMPArithmetic::to_string(const Giant& a)
    {
        if (a.empty())
            return "";
        std::vector<char> buffer(a.size()*10 + 10);
        if (a._size == 0)
            buffer[0] = '0';
        else
            mpz_get_str(buffer.data(), 10, mpz(a));
        return std::string(buffer.data());
    }

    void GMPArithmetic::to_GWNum(const Giant& a, GWNum& res)
    {
        giantstruct g;
        if (a._size >= 0)
        {
            g.n = (uint32_t*)a._data;
            g.sign = a.size();
            gianttogw(res.arithmetic().gwdata(), &g, *res);
        }
        else
        {
            Giant tmp(*this);
            add((Giant&)a, res.arithmetic().N(), tmp);
            g.n = tmp.data();
            g.sign = tmp.size();
            gianttogw(res.arithmetic().gwdata(), &g, *res);
        }
    }

    int GMPArithmetic::cmp(const Giant& a, const Giant& b)
    {
        return mpz_cmp(mpz(a), mpz(b));
    }

    int GMPArithmetic::cmp(const Giant& a, int32_t b)
    {
        if (a._size == 0)
            return 0 > b ? 1 : 0 < b ? -1 : 0;
        if (b == 0)
            return a._size > 0 ? 1 : a._size < 0 ? -1 : 0;
        return mpz_cmp_si(mpz(a), b);
    }

    void GMPArithmetic::shiftleft(Giant& a, int b, Giant& res)
    {
        mpz_mul_2exp(mpz(res), mpz(a), b);
    }

    void GMPArithmetic::shiftright(Giant& a, int b, Giant& res)
    {
        mpz_fdiv_q_2exp(mpz(res), mpz(a), b);
    }

    int GMPArithmetic::bitlen(const Giant& a)
    {
        if (a._size == 0)
            return 0;
        return (int)mpz_sizeinbase(mpz(a), 2);
    }

    bool GMPArithmetic::bit(const Giant& a, int b)
    {
        return mpz_tstbit(mpz(a), b) != 0;
    }

    void GMPArithmetic::substr(const Giant& a, int offset, int count, Giant& res)
    {
        int size = (count + GMP_NUMB_BITS - 1 + offset%GMP_NUMB_BITS)/GMP_NUMB_BITS;
        if (size > abs(a._size) - offset/GMP_NUMB_BITS)
            size = abs(a._size) - offset/GMP_NUMB_BITS;
        if (size <= 0)
        {
            res = 0;
            return;
        }
        alloc(res, size*(GMP_NUMB_BITS/32));
        memcpy(res._data, (char*)a._data + offset/GMP_NUMB_BITS*(GMP_NUMB_BITS/8), size*(GMP_NUMB_BITS/8));
        if (GMP_NUMB_BITS == 32 && (offset%GMP_NUMB_BITS + count)/GMP_NUMB_BITS == size - 1)
            ((uint32_t*)res._data)[size - 1] &= (1 << (offset + count)%GMP_NUMB_BITS) - 1;
        if (GMP_NUMB_BITS == 64 && (offset%GMP_NUMB_BITS + count)/GMP_NUMB_BITS == size - 1)
            ((uint64_t*)res._data)[size - 1] &= (1ULL << (offset + count)%GMP_NUMB_BITS) - 1;
        res._size = size; // with leading zeroes
        res._size = (res.size()*32 + GMP_NUMB_BITS - 1)/GMP_NUMB_BITS;
        mpz_fdiv_q_2exp(mpz(res), mpz(res), offset%GMP_NUMB_BITS);
    }

    void GMPArithmetic::add(Giant& a, Giant& b, Giant& res)
    {
        mpz_add(mpz(res), mpz(a), mpz(b));
    }

    void GMPArithmetic::add(Giant& a, int32_t b, Giant& res)
    {
        if (b >= 0)
            mpz_add_ui(mpz(res), mpz(a), b);
        else
            mpz_sub_ui(mpz(res), mpz(a), -b);
    }

    void GMPArithmetic::sub(Giant& a, Giant& b, Giant& res)
    {
        mpz_sub(mpz(res), mpz(a), mpz(b));
    }

    void GMPArithmetic::sub(Giant& a, int32_t b, Giant& res)
    {
        if (b <= 0)
            mpz_add_ui(mpz(res), mpz(a), -b);
        else
            mpz_sub_ui(mpz(res), mpz(a), b);
    }

    void GMPArithmetic::neg(Giant& a, Giant& res)
    {
        if (res._data != a._data)
            copy(a, res);
        res._size = -res._size;
    }

    void GMPArithmetic::mul(Giant& a, Giant& b, Giant& res)
    {
        mpz_mul(mpz(res), mpz(a), mpz(b));
    }

    void GMPArithmetic::mul(Giant& a, int32_t b, Giant& res)
    {
        mpz_mul_si(mpz(res), mpz(a), b);
    }

    void GMPArithmetic::mul(Giant& a, uint32_t b, Giant& res)
    {
        mpz_mul_ui(mpz(res), mpz(a), b);
    }

    void GMPArithmetic::div(Giant& a, Giant& b, Giant& res)
    {
        mpz_fdiv_q(mpz(res), mpz(a), mpz(b));
    }

    void GMPArithmetic::div(Giant& a, int32_t b, Giant& res)
    {
        mpz_fdiv_q_ui(mpz(res), mpz(a), abs(b));
        if (b < 0)
            res._size = -res._size;
    }

    void GMPArithmetic::div(Giant& a, uint32_t b, Giant& res)
    {
        mpz_fdiv_q_ui(mpz(res), mpz(a), b);
    }

    void GMPArithmetic::mod(Giant& a, Giant& b, Giant& res)
    {
        mpz_fdiv_r(mpz(res), mpz(a), mpz(b));
    }

    void GMPArithmetic::mod(Giant& a, uint32_t b, uint32_t& res)
    {
        res = mpz_fdiv_ui(mpz(a), b);
    }

    void GMPArithmetic::gcd(Giant& a, Giant& b, Giant& res)
    {
        if (a == 0)
        {
            copy(b, res);
            return;
        }
        if (b == 0)
        {
            copy(a, res);
            return;
        }
        mpz_gcd(mpz(res), mpz(a), mpz(b));
    }

    void GMPArithmetic::inv(Giant& a, Giant& n, Giant& res)
    {
        Giant factor(*this);
        mpz_gcdext(mpz(factor), mpz(res), NULL, mpz(a), mpz(n));
        if (factor != 1)
            throw NoInverseException(factor);
        if (res < 0)
            res += n;
    }

    void GMPArithmetic::power(Giant& a, int32_t b, Giant& res)
    {
        if (b == 0)
        {
            init(1, res);
            return;
        }
        mpz_pow_ui(mpz(res), mpz(a), b);
    }

    void GMPArithmetic::powermod(Giant& a, Giant& b, Giant& n, Giant& res)
    {
        mpz_powm(mpz(res), mpz(a), mpz(b), mpz(n));
    }

    void GMPArithmetic::rnd(Giant& res, int bits)
    {
        GiantsArithmetic::rnd(res, bits);
        if (GMP_NUMB_BITS == 64 && (res._size & 1))
            res.data()[res._size] = 0;
        res._size = (res._size*32 + GMP_NUMB_BITS - 1)/GMP_NUMB_BITS;
    }

    int GMPArithmetic::kronecker(Giant& a, Giant& b)
    {
        return mpz_kronecker(mpz(a), mpz(b));
    }

    int GMPArithmetic::kronecker(Giant& a, uint32_t b)
    {
        return mpz_kronecker_ui(mpz(a), b);
    }

    int GMPArithmetic::kronecker(uint32_t a, Giant& b)
    {
        return mpz_ui_kronecker(a, mpz(b));
    }

    void GWGMPArithmetic::alloc(Giant& a)
    {
        a._capacity = (capacity()*32 + GMP_NUMB_BITS - 1)/GMP_NUMB_BITS;
        a._size = 0;
        a._data = popg(&((gwhandle*)_gwdata)->gdata, a._capacity*(GMP_NUMB_BITS/32))->n;
    }

    void GWGMPArithmetic::free(Giant& a)
    {
        pushg(&((gwhandle*)_gwdata)->gdata, 1);
        a._capacity = 0;
        a._size = 0;
        a._data = nullptr;
    }

    void GWGMPArithmetic::alloc(Giant& a, int capacity)
    {
        if (!a.empty())
            return;
        alloc(a);
    }

    void GWGMPArithmetic::init(const GWNum& a, Giant& res)
    {
        res.arithmetic().alloc(res, capacity());
        int capacity = res._capacity;
        res._capacity = res.capacity();
        if (gwtogiant((gwhandle*)_gwdata, *a, giant(res)) < 0)
            throw InvalidFFTDataException();
        if (GMP_NUMB_BITS == 64 && (res._size & 1))
            res.data()[res._size] = 0;
        res._size = (res._size*32 + GMP_NUMB_BITS - 1)/GMP_NUMB_BITS;
        res._capacity = capacity;
    }

    void GWGMPArithmetic::to_GWNum(const Giant& a, GWNum& res)
    {
        giantstruct g;
        if (a._size >= 0)
        {
            g.n = (uint32_t*)a._data;
            g.sign = a.size();
            gianttogw((gwhandle*)_gwdata, &g, *res);
        }
        else
        {
            Giant tmp(*this);
            add((Giant&)a, res.arithmetic().N(), tmp);
            g.n = tmp.data();
            g.sign = tmp.size();
            gianttogw((gwhandle*)_gwdata, &g, *res);
        }
    }
#endif
}
