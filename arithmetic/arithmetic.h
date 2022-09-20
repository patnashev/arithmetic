#pragma once

#include <unordered_set>
#include <mutex>

#include "field.h"
#include "giant.h"

namespace arithmetic
{
    class GWNum;

    class GWState
    {
    public:
        GWState() : giants(&handle)
        {
            gwinit(&handle);
        }
        GWState(GWState& state) : giants(&handle)
        {
            clone(state);
        }
        ~GWState()
        {
            done();
        }

        void setup(int k, int b, int n, int c);
        void setup(const Giant& g);
        void setup(int bitlen);
        void clone(GWState& state);
        void done();

        gwhandle* gwdata() { return &handle; }
        Giant popg() { return Giant(giants); }
        int max_polymult_output();

        int thread_count = 1;
        int next_fft_count = 0;
        double safety_margin = 0;
        int maxmulbyconst = 3;
        bool will_error_check = false;
        bool large_pages = false;
        bool force_general_mod = false;
        bool polymult = false;
        int spin_threads = 1;
        Giant known_factors;

        void copy(const GWState& a)
        {
            thread_count = a.thread_count;
            next_fft_count = a.next_fft_count;
            safety_margin = a.safety_margin;
            maxmulbyconst = a.maxmulbyconst;
            will_error_check = a.will_error_check;
            large_pages = a.large_pages;
            force_general_mod = a.force_general_mod;
            polymult = a.polymult;
            spin_threads = a.spin_threads;
            known_factors = a.known_factors;
        }

        gwhandle handle;
        GWGiantsArithmetic giants;
        Giant *N = nullptr;
        uint32_t fingerprint;
        std::string fft_description;
        int fft_length;
        int bit_length;
    };

    class CarefulGWArithmetic;
    class ReliableGWArithmetic;

    class GWArithmetic : public FieldArithmetic<GWNum>
    {
        friend class GWNum;
    public:
        using Element = GWNum;

    protected:
        GWArithmetic(GWState& state, CarefulGWArithmetic *careful) : _state(state), _careful(careful) { }
    public:
        GWArithmetic(GWState& state);
        virtual ~GWArithmetic();

        virtual void alloc(GWNum& a) override;
        virtual void free(GWNum& a) override;
        virtual void copy(const GWNum& a, GWNum& res) override;
        virtual void move(GWNum&& a, GWNum& res) override;
        virtual void init(int32_t a, GWNum& res) override;
        virtual void init(const std::string& a, GWNum& res) override;
        virtual int cmp(const GWNum& a, const GWNum& b) override;
        virtual int cmp(const GWNum& a, int32_t b) override;
        virtual int cmp(const GWNum& a, const Giant& b);
        virtual void add(GWNum& a, GWNum& b, GWNum& res) override;
        virtual void add(GWNum& a, int32_t b, GWNum& res) override;
        virtual void add(GWNum& a, GWNum& b, GWNum& res, int options);
        virtual void sub(GWNum& a, GWNum& b, GWNum& res) override;
        virtual void sub(GWNum& a, int32_t b, GWNum& res) override;
        virtual void sub(GWNum& a, GWNum& b, GWNum& res, int options);
        virtual void neg(GWNum& a, GWNum& res) override;
        virtual void mul(GWNum& a, GWNum& b, GWNum& res) override;
        virtual void mul(GWNum& a, int32_t b, GWNum& res) override;
        virtual void mul(GWNum& a, GWNum& b, GWNum& res, int options);
        virtual void div(GWNum& a, GWNum& b, GWNum& res) override;
        virtual void div(GWNum& a, int32_t b, GWNum& res) override;
        virtual void gcd(GWNum& a, GWNum& b, GWNum& res) override;
        virtual void inv(GWNum& a, GWNum& n, GWNum& res) override;
        virtual void inv(GWNum& a, GWNum& res);
        virtual void mod(GWNum& a, GWNum& b, GWNum& res) override;

        virtual void addsub(GWNum& a, GWNum& b, GWNum& res1, GWNum& res2, int options);
        virtual void addmul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options);
        virtual void submul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options);
        virtual void muladd(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options);
        virtual void mulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options);
        virtual void mulmuladd(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options);
        virtual void mulmulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options);
        virtual void fft(GWNum& a, GWNum& res);
        virtual void unfft(GWNum& a, GWNum& res);

        void square(GWNum& a, GWNum& res, int options) { mul(a, a, res, options); }
        void setmulbyconst(int32_t a) { gwsetmulbyconst(gwdata(), a); }
        void setaddin(int32_t a) { gwsetaddin(gwdata(), a); }

        GWState& state() { return _state; }
        gwhandle* gwdata() { return &_state.handle; }
        Giant popg() { return Giant(_state.giants); }
        Giant& N() { GWASSERT(_state.N != nullptr); return *_state.N; }
        CarefulGWArithmetic& carefully() { return *_careful; }

    protected:
        GWState& _state;
        CarefulGWArithmetic* _careful = nullptr;
    };

    class CarefulGWArithmetic : public GWArithmetic
    {
    public:
        CarefulGWArithmetic(GWState& state) : GWArithmetic(state, this) { }
        virtual ~CarefulGWArithmetic() { }

        using GWArithmetic::add;
        using GWArithmetic::sub;
        using GWArithmetic::mul;
        virtual void add(GWNum& a, int32_t b, GWNum& res) override;
        virtual void sub(GWNum& a, int32_t b, GWNum& res) override;
        virtual void add(GWNum& a, GWNum& b, GWNum& res, int options) override;
        virtual void sub(GWNum& a, GWNum& b, GWNum& res, int options) override;
        virtual void mul(GWNum& a, GWNum& b, GWNum& res, int options) override;
        virtual void addsub(GWNum& a, GWNum& b, GWNum& res1, GWNum& res2, int options) override;
        virtual void addmul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options) override;
        virtual void submul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options) override;
        virtual void muladd(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options) override;
        virtual void mulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options) override;
        virtual void mulmuladd(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options) override;
        virtual void mulmulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options) override;
    };

    class ReliableGWArithmetic : public GWArithmetic
    {
    public:
        ReliableGWArithmetic(GWState& state) : GWArithmetic(state) { }
        virtual ~ReliableGWArithmetic() { }

        virtual void mul(GWNum& a, GWNum& b, GWNum& res, int options) override;
        virtual void addmul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options) override;
        virtual void submul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options) override;
        virtual void muladd(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options) override;
        virtual void mulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options) override;
        virtual void mulmuladd(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options) override;
        virtual void mulmulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options) override;

        void reset();
        void restart(int op = 0);

        int op() { return _op; }
        bool restart_flag() { return _restart_flag; }
        bool failure_flag() { return _failure_flag; }

    private:
        int _op = 0;
        std::unordered_set<int> _suspect_ops;
        bool _restart_flag = false;
        bool _failure_flag = false;
        double _max_roundoff = 0.4;
    };

    class ThreadSafeGWArithmetic : public GWArithmetic
    {
    public:
        ThreadSafeGWArithmetic(GWState& state) : GWArithmetic(state) { }
        virtual ~ThreadSafeGWArithmetic() { }

        virtual void alloc(GWNum& a) override;
        virtual void free(GWNum& a) override;

    private:
        std::mutex _mutex;
    };


    class GWNum : public FieldElement<GWArithmetic, GWNum>
    {
        friend class GWArithmetic;
        friend class ThreadSafeGWArithmetic;
        friend class PolyMult;
        friend class GWNumWrapper;
    public:
        using Arithmetic = GWArithmetic;

    protected:
        GWNum(GWArithmetic& arithmetic, gwnum a) : FieldElement<GWArithmetic, GWNum>(arithmetic), _gwnum(a) { }
    public:
        GWNum(GWArithmetic& arithmetic) : FieldElement<GWArithmetic, GWNum>(arithmetic)
        {
            arithmetic.alloc(*this);
        }
        virtual ~GWNum()
        {
            if (_gwnum != nullptr)
                arithmetic().free(*this);
        }
        GWNum(const GWNum& a) : FieldElement<GWArithmetic, GWNum>(a.arithmetic())
        {
            arithmetic().copy(a, *this);
        }
        GWNum(GWNum&& a) noexcept : FieldElement<GWArithmetic, GWNum>(a.arithmetic())
        {
            arithmetic().move(std::move(a), *this);
        }

        virtual GWNum& operator = (const GWNum& a)
        {
            arithmetic().copy(a, *this);
            return *this;
        }
        virtual GWNum& operator = (GWNum&& a) noexcept
        {
            arithmetic().move(std::move(a), *this);
            return *this;
        }
        using FieldElement<GWArithmetic, GWNum>::operator=;

        virtual std::string to_string() const override;
        GWNum& operator = (const Giant& a)
        {
            arithmetic().state().giants.to_GWNum(a, *this);
            return *this;
        }

        friend gwnum operator *(const GWNum& a) { return a._gwnum; }

        GWNum& operator /= (GWNum&& a)
        {
            *this *= a.inv();
            return *this;
        }
        using FieldElement<GWArithmetic, GWNum>::operator/=;
        friend GWNum operator / (GWNum& a, GWNum& b)
        {
            GWNum t(b);
            t.inv();
            return a * std::move(t);
        }
        friend GWNum operator / (GWNum& a, GWNum&& b)
        {
            return a * std::move(b.inv());
        }
        friend GWNum operator / (GWNum&& a, GWNum& b)
        {
            GWNum res(std::move(a));
            res /= ((void*)&a == (void*)&b) ? res : b;
            return res;
        }
        friend GWNum operator / (GWNum&& a, GWNum&& b)
        {
            return std::move(a) * b.inv();
        }
        friend GWNum operator / (int32_t a, GWNum& b)
        {
            GWNum t(b);
            t.inv();
            return a * std::move(t);
        }
        friend GWNum operator / (int32_t a, GWNum&& b)
        {
            return a * std::move(b.inv());
        }
        friend Giant gcd(GWNum& a, Giant& b)
        {
            Giant res(b.arithmetic());
            res = a;
            res.gcd(b);
            return res;
        }
        GWNum& inv()
        {
            arithmetic().inv(*this, *this);
            return *this;
        }
        using FieldElement<GWArithmetic, GWNum>::inv;
        friend GWNum inv(GWNum& a)
        {
            GWNum t(a);
            t.inv();
            return t;
        }
        friend GWNum inv(GWNum&& a)
        {
            return a.inv();
        }

    protected:
        gwnum _gwnum = nullptr;
    };

    class GWNumWrapper : public GWNum
    {
    public:
        GWNumWrapper(GWArithmetic& gw, gwnum a) : GWNum(gw, a) { }
        GWNumWrapper(const GWNum& a) : GWNum(a.arithmetic(), a._gwnum) { }
        ~GWNumWrapper() { _gwnum = nullptr; }
        using GWNum::operator=;
        GWNum& operator = (GWNum&& a) noexcept override { *this = a; return *this; }
        operator GWNum&&() = delete; // No effect, just a remainder
    };
}

int gwconvert(
    gwhandle *gwdata_s,	/* Handle initialized by gwsetup */
    gwhandle *gwdata_d,	/* Handle initialized by gwsetup */
    gwnum	s,
    gwnum	d);
