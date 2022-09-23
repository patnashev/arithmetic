#pragma once

#include "field.h"

namespace arithmetic
{
    class Giant;
    class GWNum;

    class GiantsArithmetic : public FieldArithmetic<Giant>
    {
        friend class Giant;

    public:
        GiantsArithmetic() { }
        virtual ~GiantsArithmetic();

        virtual void alloc(Giant& a) override;
        virtual void alloc(Giant& a, int capacity);
        virtual void free(Giant& a) override;
        virtual void copy(const Giant& a, Giant& res) override;
        virtual void move(Giant&& a, Giant& res) override;
        virtual void init(int32_t a, Giant& res) override;
        virtual void init(uint32_t a, Giant& res);
        virtual void init(const std::string& a, Giant& res) override;
        virtual void init(uint32_t* data, int size, Giant& res);
        virtual void init(const GWNum& a, Giant& res);
        virtual void to_GWNum(const Giant& a, GWNum& res);
        virtual int cmp(const Giant& a, const Giant& b) override;
        virtual int cmp(const Giant& a, int32_t b) override;
        virtual void add(Giant& a, Giant& b, Giant& res) override;
        virtual void add(Giant& a, int32_t b, Giant& res) override;
        virtual void sub(Giant& a, Giant& b, Giant& res) override;
        virtual void sub(Giant& a, int32_t b, Giant& res) override;
        virtual void neg(Giant& a, Giant& res) override;
        virtual void mul(Giant& a, Giant& b, Giant& res) override;
        virtual void mul(Giant& a, int32_t b, Giant& res) override;
        virtual void mul(Giant& a, uint32_t b, Giant& res);
        virtual void div(Giant& a, Giant& b, Giant& res) override;
        virtual void div(Giant& a, int32_t b, Giant& res) override;
        virtual void div(Giant& a, uint32_t b, Giant& res);
        virtual void gcd(Giant& a, Giant& b, Giant& res) override;
        virtual void inv(Giant& a, Giant& n, Giant& res) override;
        virtual void mod(Giant& a, Giant& b, Giant& res) override;
        virtual void mod(Giant& a, uint32_t b, uint32_t& res);
        virtual void shiftleft(Giant& a, int b, Giant& res);
        virtual void shiftright(Giant& a, int b, Giant& res);
        virtual int bitlen(const Giant& a);
        virtual bool bit(const Giant& a, int b);
        virtual double log2(const Giant& a);
        virtual void power(Giant& a, int32_t b, Giant& res);
        virtual void powermod(Giant& a, Giant& b, Giant& n, Giant& res);
        virtual void rnd_seed(Giant& a);
        virtual void rnd(Giant& res, int bits);

        static GiantsArithmetic& default_arithmetic();

    private:
        void* _rnd_state = nullptr;
    };

    class GWGiantsArithmetic : public GiantsArithmetic
    {
        friend class GWState;

    public:
        GWGiantsArithmetic(void *gwdata) : _gwdata(gwdata) { }
        virtual ~GWGiantsArithmetic() { }

        virtual void alloc(Giant& a) override;
        virtual void alloc(Giant& a, int capacity) override;
        virtual void free(Giant& a) override;
        virtual void init(const GWNum& a, Giant& res);
        virtual void to_GWNum(const Giant& a, GWNum& res);
        virtual void inv(Giant& a, Giant& n, Giant& res) override;

        int capacity() const { return _capacity; }

    private:
        void* _gwdata;
        int _capacity = 0;
    };

    class Giant : public FieldElement<GiantsArithmetic, Giant>
    {
        friend class GiantsArithmetic;
        friend class GWGiantsArithmetic;

    public:
        Giant() : FieldElement<GiantsArithmetic, Giant>(GiantsArithmetic::default_arithmetic())
        {
            arithmetic().alloc(*this);
        }
        Giant(GiantsArithmetic& arithmetic) : FieldElement<GiantsArithmetic, Giant>(arithmetic)
        {
            arithmetic.alloc(*this);
        }
        Giant(GiantsArithmetic& arithmetic, int capacity) : FieldElement<GiantsArithmetic, Giant>(arithmetic)
        {
            arithmetic.alloc(*this, capacity);
        }
        virtual ~Giant()
        {
            if (_giant != nullptr)
                arithmetic().free(*this);
        }
        Giant(const Giant& a) : FieldElement<GiantsArithmetic, Giant>(a.arithmetic())
        {
            arithmetic().copy(a, *this);
        }
        Giant(Giant&& a) noexcept : FieldElement<GiantsArithmetic, Giant>(a.arithmetic())
        {
            arithmetic().move(std::move(a), *this);
        }

        Giant& operator = (const Giant& a)
        {
            arithmetic().copy(a, *this);
            return *this;
        }
        Giant& operator = (Giant&& a) noexcept
        {
            arithmetic().move(std::move(a), *this);
            return *this;
        }
        using FieldElement<GiantsArithmetic, Giant>::operator=;
        virtual Giant& operator = (const GWNum& a)
        {
            arithmetic().init(a, *this);
            return *this;
        };

        virtual std::string to_string() const override;
        virtual std::string to_res64() const;
        virtual void to_GWNum(GWNum& res) const { arithmetic().to_GWNum(*this, res); };

        uint32_t* data() { return _giant->n; }
        const uint32_t* data() const { return _giant->n; }
        int size() const { return abs(_giant->sign); }
        int capacity() const { return _capacity; }
        bool empty() const { return _giant == nullptr; }

        Giant& operator = (uint32_t a)
        {
            arithmetic().init(a, *this);
            return *this;
        }
        Giant& operator *= (uint32_t a)
        {
            arithmetic().mul(*this, a, *this);
            return *this;
        }
        using FieldElement<GiantsArithmetic, Giant>::operator*=;
        Giant& operator <<= (int a)
        {
            arithmetic().shiftleft(*this, a, *this);
            return *this;
        }
        friend Giant operator << (Giant& a, int b)
        {
            Giant res(a.arithmetic());
            res.arithmetic().shiftleft(a, b, res);
            return res;
        }
        friend Giant operator << (Giant&& a, int b)
        {
            Giant res(std::move(a));
            res <<= b;
            return res;
        }
        Giant& operator >>= (int a)
        {
            arithmetic().shiftright(*this, a, *this);
            return *this;
        }
        friend Giant operator >> (Giant& a, int b)
        {
            Giant res(a.arithmetic());
            res.arithmetic().shiftright(a, b, res);
            return res;
        }
        friend Giant operator >> (Giant&& a, int b)
        {
            Giant res(std::move(a));
            res >>= b;
            return res;
        }
        Giant& operator /= (Giant&& a)
        {
            *this /= a;
            return *this;
        }
        Giant& operator /= (uint32_t a)
        {
            arithmetic().div(*this, a, *this);
            return *this;
        }
        using FieldElement<GiantsArithmetic, Giant>::operator/=;
        friend Giant operator / (Giant& a, Giant& b)
        {
            Giant res(a.arithmetic());
            res.arithmetic().div(a, b, res);
            return res;
        }
        friend Giant operator / (Giant&& a, Giant& b)
        {
            Giant res(std::move(a));
            res /= ((void*)&a == (void*)&b) ? res : b;
            return res;
        }
        friend Giant operator / (Giant& a, Giant&& b)
        {
            return a / b;
        }
        friend Giant operator / (Giant&& a, Giant&& b)
        {
            return std::move(a) / b;
        }
        friend Giant operator / (Giant& a, uint32_t b)
        {
            Giant res(a.arithmetic());
            res.arithmetic().div(a, b, res);
            return res;
        }
        friend Giant operator / (Giant&& a, uint32_t b)
        {
            Giant res(std::move(a));
            res /= b;
            return res;
        }
        friend uint32_t operator % (Giant& a, uint32_t b)
        {
            uint32_t res;
            a.arithmetic().mod(a, b, res);
            return res;
        }
        friend uint32_t operator % (Giant&& a, uint32_t b)
        {
            uint32_t res;
            a.arithmetic().mod(a, b, res);
            return res;
        }
        int bitlen() const { return arithmetic().bitlen(*this); }
        bool bit(int b) const { return arithmetic().bit(*this, b); }
        friend double log2(Giant& a)
        {
            return a.arithmetic().log2(a); 
        }
        Giant& power(int32_t a)
        {
            arithmetic().power(*this, a, *this);
            return *this;
        }
        friend Giant power(Giant& a, int32_t b)
        {
            Giant res(a.arithmetic());
            res.arithmetic().power(a, b, res);
            return res;
        }
        friend Giant power(Giant&& a, int32_t b)
        {
            Giant res(std::move(a));
            res.power(b);
            return res;
        }
        static Giant rnd(int bits)
        {
            Giant res;
            res.arithmetic().rnd(res, bits);
            return res;
        }

    private:
        giant _giant = nullptr;
        int _capacity = 0;
    };
}
