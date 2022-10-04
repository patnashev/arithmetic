
#include <cmath>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include "gwnum.h"
#include "cpuid.h"
#include "giant.h"
#include "arithmetic.h"
#include "exception.h"

namespace arithmetic
{
    GiantsArithmetic _defaultGiantsArithmetic;
    GiantsArithmetic& GiantsArithmetic::default_arithmetic()
    {
        return _defaultGiantsArithmetic;
    }

    void GiantsArithmetic::alloc(Giant& a)
    {
        a._capacity = 0;
        a._giant = nullptr;
    }

    void GiantsArithmetic::free(Giant& a)
    {
        ::free(a._giant);
        a._capacity = 0;
        a._giant = nullptr;
    }

    void GiantsArithmetic::alloc(Giant& a, int capacity)
    {
        giant tmp = allocgiant(capacity);
        tmp->sign = 0;
        if (a._giant != nullptr)
        {
            gtog(a._giant, tmp);
            ::free(a._giant);
        }
        a._capacity = capacity;
        a._giant = tmp;
    }

    void GiantsArithmetic::copy(const Giant& a, Giant& res)
    {
        if (a.empty())
        {
            free(res);
            return;
        }
        if (a._giant->sign == 0)
        {
            init(0, res);
            return;
        }
        if (res._giant == nullptr || res._capacity < abs(a._giant->sign))
            alloc(res, abs(a._giant->sign));
        gtog(a._giant, res._giant);
    }

    void GiantsArithmetic::move(Giant&& a, Giant& res)
    {
        if (res._giant != nullptr)
            free(res);
        res._capacity = a._capacity;
        res._giant = a._giant;
        a._capacity = 0;
        a._giant = nullptr;
    }

    void GiantsArithmetic::init(int32_t a, Giant& res)
    {
        if (res._giant == nullptr || res._capacity < 1)
            alloc(res, 1);
        itog(a, res._giant);
    }

    void GiantsArithmetic::init(uint32_t a, Giant& res)
    {
        if (res._giant == nullptr || res._capacity < 1)
            alloc(res, 1);
        ultog(a, res._giant);
    }

    void GiantsArithmetic::init(const std::string& a, Giant& res)
    {
        if (res._giant == nullptr || res._capacity < ((int)a.length() + 8)/9)
            alloc(res, ((int)a.length() + 8)/9);
        ctog(a.data(), res._giant);
    }

    void GiantsArithmetic::init(uint32_t* data, int size, Giant& res)
    {
        if (size == 0)
        {
            init(0, res);
            return;
        }
        if (res._giant == nullptr || res._capacity < size)
            alloc(res, size);
        memcpy(res._giant->n, data, size*4);
        res._giant->sign = size;
        while (res._giant->sign && !res._giant->n[res._giant->sign - 1])
            res._giant->sign--;
    }

    void GiantsArithmetic::init(const GWNum& a, Giant& res)
    {
        int capacity = a.arithmetic().state().giants.capacity();
        if (res._giant == nullptr || res._capacity < capacity)
            res.arithmetic().alloc(res, capacity);
        if (gwtogiant(a.arithmetic().gwdata(), *a, res._giant) < 0)
            throw InvalidFFTDataException();
    }

    void GiantsArithmetic::to_GWNum(const Giant& a, GWNum& res)
    {
        if (a._giant->sign >= 0)
            gianttogw(res.arithmetic().gwdata(), a._giant, *res);
        else
        {
            Giant tmp(*this);
            add((Giant&)a, res.arithmetic().N(), tmp);
            gianttogw(res.arithmetic().gwdata(), tmp._giant, *res);
        }
    }

    std::string Giant::to_string() const
    {
        if (_giant == nullptr)
            return "";
        std::vector<char> buffer(abs(_giant->sign)*10 + 10);
        //std::iterator<char> x = buffer.begin();
        if (_giant->sign ==  0)
            buffer[0] = '0';
        else if (_giant->sign >  0)
            gtoc(_giant, buffer.data(), (int)buffer.size());
        else
        {
            _giant->sign = -_giant->sign;
            buffer[0] = '-';
            gtoc(_giant, buffer.data() + 1, (int)buffer.size() - 1);
            _giant->sign = -_giant->sign;
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

    void GiantsArithmetic::add(Giant& a, Giant& b, Giant& res)
    {
        int size = abs(a._giant->sign) + 1;
        if (size < abs(b._giant->sign) + 1)
            size = abs(b._giant->sign) + 1;
        if (res._giant == nullptr || res._capacity < size)
            alloc(res, size);
        if (res._giant != a._giant)
            copy(a, res);
        if (a._giant != b._giant)
            addg(b._giant, res._giant);
        else
            gshiftleft(1, res._giant);
    }

    void GiantsArithmetic::add(Giant& a, int32_t b, Giant& res)
    {
        if (res._giant == nullptr || res._capacity < abs(a._giant->sign) + 1)
            alloc(res, abs(a._giant->sign) + 1);
        if (res._giant != a._giant)
            copy(a, res);
        sladdg(b, res._giant);
    }

    void GiantsArithmetic::sub(Giant& a, Giant& b, Giant& res)
    {
        int capacity = abs(a._giant->sign) + 1;
        if (capacity < abs(b._giant->sign) + 1)
            capacity = abs(b._giant->sign) + 1;
        if (res._giant == nullptr || res._capacity < capacity)
            alloc(res, capacity);
        if (res._giant != a._giant)
            copy(a, res);
        subg(b._giant, res._giant);
    }

    void GiantsArithmetic::sub(Giant& a, int32_t b, Giant& res)
    {
        if (res._giant == nullptr || res._capacity < abs(a._giant->sign) + 1)
            alloc(res, abs(a._giant->sign) + 1);
        if (res._giant != a._giant)
            copy(a, res);
        sladdg(-b, res._giant);
    }

    void GiantsArithmetic::neg(Giant& a, Giant& res)
    {
        if (res._giant != a._giant)
            copy(a, res);
        res._giant->sign = -res._giant->sign;
    }

    void GiantsArithmetic::mul(Giant& a, Giant& b, Giant& res)
    {
        if (res._giant == nullptr || res._capacity < abs(a._giant->sign) + abs(b._giant->sign))
            alloc(res, abs(a._giant->sign) + abs(b._giant->sign));
        if (res._giant != a._giant)
            copy(a, res);
        if (a._giant != b._giant)
        {
            if ((abs(res._giant->sign) > 10000 && abs(b._giant->sign) > 10 && abs(res._giant->sign) > 4*abs(b._giant->sign)) || (abs(b._giant->sign) > 10000 && abs(res._giant->sign) > 10 && abs(b._giant->sign) > 4*abs(res._giant->sign)))
            {
                setmulmode(FFT_MUL);
                mulg(b._giant, res._giant);
                setmulmode(AUTO_MUL);
            }
            else
                mulg(b._giant, res._giant);
        }
        else
            squareg(res._giant);
    }

    void GiantsArithmetic::mul(Giant& a, int32_t b, Giant& res)
    {
        if (res._giant == nullptr || res._capacity < abs(a._giant->sign) + 1)
            alloc(res, abs(a._giant->sign) + 1);
        if (res._giant != a._giant)
            copy(a, res);
        imulg(b, res._giant);
    }

    void GiantsArithmetic::mul(Giant& a, uint32_t b, Giant& res)
    {
        if (res._giant == nullptr || res._capacity < abs(a._giant->sign) + 1)
            alloc(res, abs(a._giant->sign) + 1);
        if (res._giant != a._giant)
            copy(a, res);
        ulmulg(b, res._giant);
    }

    void GiantsArithmetic::shiftleft(Giant& a, int b, Giant& res)
    {
        if (res._giant == nullptr || res._capacity < abs(a._giant->sign) + (b + 31)/32)
            alloc(res, abs(a._giant->sign) + (b + 31)/32);
        if (res._giant != a._giant)
            copy(a, res);
        gshiftleft(b, res._giant);
    }

    void GiantsArithmetic::shiftright(Giant& a, int b, Giant& res)
    {
        if (res._giant == nullptr || res._capacity < abs(a._giant->sign) - b/32)
            alloc(res, abs(a._giant->sign) - b/32);
        if (res._giant != a._giant)
            copy(a, res);
        gshiftright(b, res._giant);
    }

    int GiantsArithmetic::bitlen(const Giant& a)
    {
        return ::bitlen(a._giant);
    }

    bool GiantsArithmetic::bit(const Giant& a, int b)
    {
        return bitval(a._giant, b) != 0;
    }

    void GiantsArithmetic::substr(const Giant& a, int offset, int count, Giant& res)
    {
        int size = (count + 31 + offset%32)/32;
        if (size > abs(a._giant->sign) - offset/32)
            size = abs(a._giant->sign) - offset/32;
        if (res._giant == nullptr || res._capacity < size)
            alloc(res, size);
        memcpy(res._giant->n, a._giant->n + offset/32, size*4);
        if ((offset%32 + count)/32 == size - 1)
            res._giant->n[size - 1] &= (1 << (offset + count)%32) - 1;
        for (; size > 0 && res._giant->n[size - 1] == 0; size--);
        res._giant->sign = size;
        gshiftright(offset%32, res._giant);
    }

    int GiantsArithmetic::cmp(const Giant& a, const Giant& b)
    {
        return gcompg(a._giant, b._giant);
    }

    double GiantsArithmetic::log2(const Giant& a)
    {
        if (a._giant == nullptr || a._giant->sign == 0)
            return 0;
        if (abs(a._giant->sign) == 1)
            return std::log2(a._giant->n[0]);
        if (abs(a._giant->sign) == 2)
            return std::log2(((uint64_t*)a._giant->n)[0]);
        return ::bitlen(a._giant);
    }

    int GiantsArithmetic::cmp(const Giant& a, int32_t b)
    {
        if (a._giant->sign > 1 || (a._giant->sign == 1 && ::bitlen(a._giant) == 32))
            return 1;
        if (a._giant->sign < -1 || (a._giant->sign == -1 && ::bitlen(a._giant) == 32))
            return -1;
        if (a._giant->sign == 0)
            return 0 > b ? 1 : 0 < b ? -1 : 0;
        return ((int)a._giant->n[0]) > b ? 1 : ((int)a._giant->n[0]) < b ? -1 : 0;
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
        if (res._giant == nullptr || res._capacity < abs(a._giant->sign))
            alloc(res, abs(a._giant->sign));
        if (res._giant != a._giant)
            copy(a, res);
        gcdg(b._giant, res._giant);
    }

    void GiantsArithmetic::inv(Giant& a, Giant& n, Giant& res)
    {
        if (res._giant == nullptr || res._capacity < abs(n._giant->sign))
            alloc(res, abs(n._giant->sign));
        if (res._giant != a._giant)
            copy(a, res);
        invg(n._giant, res._giant);
        if (res._giant->sign < 0)
            throw NoInverseException(res);
    }

    void GiantsArithmetic::div(Giant& a, Giant& b, Giant& res)
    {
        if (res._giant == nullptr || res._capacity < abs(a._giant->sign))
            alloc(res, abs(a._giant->sign));
        if (res._giant != a._giant)
            copy(a, res);
        divg(b._giant, res._giant);
    }

    void GiantsArithmetic::div(Giant& a, int32_t b, Giant& res)
    {
        if (res._giant == nullptr || res._capacity < abs(a._giant->sign))
            alloc(res, abs(a._giant->sign));
        if (res._giant != a._giant)
            copy(a, res);
        dbldivg(b, res._giant);
    }

    void GiantsArithmetic::div(Giant& a, uint32_t b, Giant& res)
    {
        if (res._giant == nullptr || res._capacity < abs(a._giant->sign))
            alloc(res, abs(a._giant->sign));
        if (res._giant != a._giant)
            copy(a, res);
        dbldivg(b, res._giant);
    }

    void GiantsArithmetic::mod(Giant& a, Giant& b, Giant& res)
    {
        if (res._giant == nullptr || res._capacity < abs(a._giant->sign))
            alloc(res, abs(a._giant->sign));
        if (res._giant != a._giant)
            copy(a, res);
        modg(b._giant, res._giant);
    }

    void GiantsArithmetic::mod(Giant& a, uint32_t b, uint32_t& res)
    {
        uint64_t res64 = 0;
        for (int i = abs(a._giant->sign) - 1; i >= 0; i--)
        {
            res64 <<= 32;
            res64 += a._giant->n[i];
            res64 %= b;
        }
        res = (uint32_t)res64;
    }

    void GiantsArithmetic::power(Giant& a, int32_t b, Giant& res)
    {
        if (b == 0)
        {
            init(1, res);
            return;
        }
        int capacity = abs(a._giant->sign);
        if (capacity == 1)
            capacity = (int)(std::log2(a._giant->n[0])*b/32) + 1;
        else
            capacity *= b;
        if (res._giant == nullptr || res._capacity < capacity)
            alloc(res, capacity);
        if (res._giant != a._giant)
            copy(a, res);
        ::power(res._giant, b);
    }

    void GiantsArithmetic::powermod(Giant& a, Giant& b, Giant& n, Giant& res)
    {
        if (res._giant == nullptr || res._capacity < abs(n._giant->sign))
            alloc(res, abs(n._giant->sign));
        if (res._giant != a._giant)
            copy(a, res);
        ::powermodg(res._giant, b._giant, n._giant);
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
        if (_rnd_state == nullptr)
        {
            Giant tmp(*this, 2);
            *(double*)tmp._giant->n = getHighResTimer();
            tmp._giant->sign = 2;
            rnd_seed(tmp);
        }
        struct mt_state *state = (struct mt_state *)_rnd_state;
        if (res._giant == nullptr || res._capacity < (bits + 31)/32)
            alloc(res, (bits + 31)/32);
        int i;
        for (i = 0; i < (bits + 31)/32; i++)
            res._giant->n[i] = genrand_int32(state);
        if ((bits & 31) != 0)
            res._giant->n[i - 1] &= (1 << (bits & 31)) - 1;
        while (i > 0 && res._giant->n[i - 1] == 0)
            i--;
        res._giant->sign = i;
    }

    void GWGiantsArithmetic::alloc(Giant& a)
    {
        a._capacity = _capacity;
        a._giant = popg(&((gwhandle*)_gwdata)->gdata, a._capacity);
    }

    void GWGiantsArithmetic::free(Giant& a)
    {
        pushg(&((gwhandle*)_gwdata)->gdata, 1);
        a._giant = nullptr;
    }

    void GWGiantsArithmetic::alloc(Giant& a, int capacity)
    {
        if (a._giant != nullptr)
            return;
        alloc(a);
    }

    void GWGiantsArithmetic::init(const GWNum& a, Giant& res)
    {
        if (res._giant == nullptr || res._capacity < capacity())
            res.arithmetic().alloc(res, capacity());
        if (gwtogiant((gwhandle*)_gwdata, *a, res._giant) < 0)
            throw InvalidFFTDataException();
    }

    void GWGiantsArithmetic::to_GWNum(const Giant& a, GWNum& res)
    {
        if (a._giant->sign >= 0)
            gianttogw((gwhandle*)_gwdata, a._giant, *res);
        else
        {
            Giant tmp(*this);
            add((Giant&)a, res.arithmetic().N(), tmp);
            gianttogw((gwhandle*)_gwdata, tmp._giant, *res);
        }
    }

    void GWGiantsArithmetic::inv(Giant& a, Giant& n, Giant& res)
    {
        if (res._giant == nullptr || res._capacity < abs(n._giant->sign))
            alloc(res, abs(n._giant->sign));
        if (res._giant != a._giant)
            copy(a, res);
        invgi(&((gwhandle*)_gwdata)->gdata, 0, n._giant, res._giant);
        if (res._giant->sign < 0)
            throw NoInverseException(res);
    }
}
