#pragma once

#include <memory>
#include <vector>

#include "arithmetic.h"
#include "polymult.h"

namespace arithmetic
{
    class Poly;

    class PolyMult
    {
    public:
        static int L3_CACHE_MB;

    public:
        PolyMult(GWArithmetic& gw, int max_threads = 1);
        ~PolyMult();

        void alloc(Poly& a, int size);
        void free(Poly& a);
        void copy(const Poly& a, Poly& res);
        void move(Poly&& a, Poly& res);
        void fft(const Poly& a, Poly& res);
        void init(bool monic, Poly& res);
        void init(const GWNum& a, bool monic, Poly& res);
        void init(GWNum&& a, bool monic, Poly& res);
        void init(gwnum* data, size_t size, bool freeable, bool monic, Poly& res);
        void mul(Poly& a, Poly& b, Poly& res, int options);
        void mul(Poly&& a, Poly&& b, Poly& res, int options);
        void mul_split(Poly& a, Poly& b, Poly& res_lo, Poly& res_hi, int size, int options);
        void mul_split(Poly&& a, Poly&& b, Poly& res_lo, Poly& res_hi, int size, int options);
        void mul_range(Poly& a, Poly& b, Poly& res, int offset, int count, int options);
        void mul_twohalf(Poly& a, Poly& b, Poly& c, Poly& res1, Poly& res2, int half, int options);
        void mul_twohalf(Poly&& a, Poly& b, Poly& c, Poly& res1, Poly& res2, int half, int options);
        void fma_range(Poly& a, Poly& b, Poly& fma, Poly& res, int offset, int count, int options);
        void preprocess(Poly& a, Poly& res, int size, int options);
        void preprocess_and_mul(Poly& a, Poly& b, Poly& res, int size, int options);
        void reciprocal(Poly& a, Poly& res, int options);
        void shiftleft(Poly& a, int b, Poly& res);
        void shiftright(Poly& a, int b, Poly& res);
        void convert(const Poly& a, PolyMult& pm_res, Poly& res);
        void insert(GWNum&& a, Poly& res, size_t pos);
        GWNum remove(Poly& a, size_t pos);

        void set_threads(int threads);

        GWArithmetic& gw() const { return _gw; }
        int max_output() const { return _max_output; }
        pmhandle* pmdata() { return &_pmdata; }
        static int max_polymult_output(GWState& state);

    private:
        void poly_seize(Poly& a, Poly& res, Poly& to_free, int size);

    private:
        GWArithmetic& _gw;
        int _max_output;
        pmhandle _pmdata;
    };

    class Poly
    {
        friend class PolyMult;

    public:
        Poly(PolyMult& pm) : _pm(pm), _cache(nullptr), _cache_size(0), _monic(false)
        {
        }
        Poly(PolyMult& pm, int size, bool monic) : _pm(pm), _cache(nullptr), _cache_size(0), _monic(monic)
        {
            pm.alloc(*this, size);
        }
        virtual ~Poly()
        {
            pm().free(*this);
        }
        Poly(const Poly& a) : _pm(a.pm()), _cache(nullptr)
        {
            pm().copy(a, *this);
        }
        Poly(Poly&& a) noexcept : _pm(a.pm()), _cache(nullptr)
        {
            pm().move(std::move(a), *this);
        }

        Poly& operator = (const Poly& a)
        {
            pm().copy(a, *this);
            return *this;
        }
        Poly& operator = (Poly&& a) noexcept
        {
            pm().move(std::move(a), *this);
            return *this;
        }

        Poly& operator <<= (int a)
        {
            pm().shiftleft(*this, a, *this);
            return *this;
        }
        friend Poly operator << (Poly& a, int b)
        {
            Poly res(a.pm());
            res.pm().shiftleft(a, b, res);
            return res;
        }
        friend Poly operator << (Poly&& a, int b)
        {
            Poly res(std::move(a));
            res <<= b;
            return res;
        }
        Poly& operator >>= (int a)
        {
            pm().shiftright(*this, a, *this);
            return *this;
        }
        friend Poly operator >> (Poly& a, int b)
        {
            Poly res(a.pm());
            res.pm().shiftright(a, b, res);
            return res;
        }
        friend Poly operator >> (Poly&& a, int b)
        {
            Poly res(std::move(a));
            res >>= b;
            return res;
        }

        GWNum eval(GWNum& x);
        Poly reciprocal(int precision, int options);

        PolyMult& pm() const { return _pm; }
        bool monic() const { return _monic; }
        bool preprocessed() const { return _cache != nullptr; }
        int degree() const { return (int)size() - (_monic ? 0 : 1); }
        size_t size() const { return _cache != nullptr ? _cache_size : _poly.size(); }
        gwnum* data() { return _cache != nullptr ? _cache : _poly.data(); }
        bool empty() const { return _cache == nullptr && _poly.empty(); }
        const GWNumWrapper at(size_t pos) { return GWNumWrapper(pm().gw(), data()[pos]); }
        void push_back(GWNum&& a) { pm().insert(std::move(a), *this, size()); }
        GWNum pop_back() { return pm().remove(*this, size() - 1); }

    private:
        PolyMult& _pm;
        std::vector<gwnum> _poly;
        gwnum* _cache;
        int _cache_size;
        bool _monic;
        bool _freeable = true;
    };
}
