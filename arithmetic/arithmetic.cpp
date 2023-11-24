
#include <vector>
#include <stdlib.h>
#include "gwnum.h"
#include "cpuid.h"
#include "arithmetic.h"
#include "exception.h"

gwnum convert_factor = NULL;

namespace arithmetic
{
    void GWState::init()
    {
        if (N)
            throw ArithmeticException(); // call done() first
        gwset_num_threads(gwdata(), thread_count);
        gwset_larger_fftlen_count(gwdata(), next_fft_count);
        gwset_safety_margin(gwdata(), safety_margin);
        gwset_maxmulbyconst(gwdata(), maxmulbyconst);
        if (will_error_check)
            gwset_will_error_check(gwdata());
        if (large_pages)
            gwset_use_large_pages(gwdata());
        if (force_mod_type == 1 || force_mod_type == 2)
            gwdata()->force_general_mod = force_mod_type;
        else
            gwdata()->force_general_mod = 0;
        if (polymult_safety_margin > 0)
        {
            gwset_using_polymult(gwdata());
            gwset_polymult_safety_margin(gwdata(), (float)polymult_safety_margin);
        }
        gwset_use_spin_wait(gwdata(), spin_threads);
        if (instructions == "SSE2")
        {
            gwdata()->cpu_flags &= ~(CPU_AVX | CPU_FMA3 | CPU_AVX512F);
            gwdata()->cpu_flags |= CPU_SSE2;
        }
        else if (instructions == "AVX")
        {
            gwdata()->cpu_flags &= ~(CPU_FMA3 | CPU_AVX512F);
            gwdata()->cpu_flags |= CPU_AVX;
        }
        else if (instructions == "FMA3")
        {
            gwdata()->cpu_flags &= ~(CPU_AVX512F);
            gwdata()->cpu_flags |= CPU_AVX | CPU_FMA3;
        }
        else if (instructions == "AVX512F")
        {
            gwdata()->cpu_flags |= CPU_AVX512F;
        }
        if (information_only)
            gwset_information_only(gwdata());
    }

    void GWState::setup(uint64_t k, uint64_t b, int n, int c)
    {
        init();
        if (gwsetup(gwdata(), (double)k, (uint32_t)b, n, c))
            throw ArithmeticException();
        bit_length = (int)gwdata()->bit_length;
        if (gwdata()->GENERAL_MOD)
            bit_length /= 2;
        giants.reset(GiantsArithmetic::alloc_gwgiants(gwdata(), (bit_length >> 5) + 10));
        N.reset(new Giant());
        *N = (uint32_t)k*power(std::move(*N = (uint32_t)b), n) + c;
        if (gwdata()->GENERAL_MOD)
            bit_length = N->bitlen();
        fingerprint = *N%3417905339UL;
        if (!known_factors.empty() && known_factors > 1 && *N%known_factors != 0)
            throw ArithmeticException();
        if (!known_factors.empty() && known_factors > 1)
            *N /= known_factors;
        char buf[200];
        gwfft_description(gwdata(), buf);
        fft_description = buf;
        fft_length = gwfftlen(gwdata());
    }

    void GWState::setup(const Giant& g)
    {
        init();
        if (gwsetup_general_mod(gwdata(), g.data(), g.size()))
            throw ArithmeticException();
        bit_length = (int)gwdata()->bit_length;
        if (gwdata()->GENERAL_MOD)
            bit_length /= 2;
        giants.reset(GiantsArithmetic::alloc_gwgiants(gwdata(), (bit_length >> 5) + 10));
        N.reset(new Giant());
        *N = g;
        if (gwdata()->GENERAL_MOD)
            bit_length = N->bitlen();
        fingerprint = *N%3417905339UL;
        if (!known_factors.empty() && known_factors > 1 && *N%known_factors != 0)
            throw ArithmeticException();
        if (!known_factors.empty() && known_factors > 1)
            *N /= known_factors;
        char buf[200];
        gwfft_description(gwdata(), buf);
        fft_description = buf;
        fft_length = gwfftlen(gwdata());
    }

    void GWState::setup(int bitlen)
    {
        init();
        if (gwsetup_without_mod(gwdata(), bitlen))
            throw ArithmeticException();
        bit_length = (int)gwdata()->bit_length;
        giants.reset(GiantsArithmetic::alloc_gwgiants(gwdata(), (bit_length >> 5) + 10));
        N.reset();
        fingerprint = 0;
        char buf[200];
        gwfft_description(gwdata(), buf);
        fft_description = buf;
        fft_length = gwfftlen(gwdata());
    }

    void GWState::clone(GWState& state)
    {
        if (N)
            throw ArithmeticException(); // call done() first
        copy(state);
        gwclone(&handle, &state.handle);
        bit_length = state.bit_length;
        giants.reset(GiantsArithmetic::alloc_gwgiants(gwdata(), state.giants->capacity()));
        N.reset(new Giant());
        *N = *state.N;
        fingerprint = state.fingerprint;
        fft_description = state.fft_description;
        fft_length = state.fft_length;
        mod_gwstate.reset();
    }

    void GWState::done()
    {
        mod_gwstate.reset();
        N.reset();
        giants.reset();
        fft_description.clear();
        fft_length = 0;
        gwdone(&handle);
        gwinit(&handle);
        convert_factor = NULL;
    }

    void GWState::mod(arithmetic::Giant& a, arithmetic::Giant& res)
    {
        if (!mod_gwstate)
        {
            mod_gwstate.reset(new GWState());
            mod_gwstate->force_mod_type = 2;
            if (gwdata()->GENERAL_MOD)
                mod_gwstate->setup(square(*N));
            else
                mod_gwstate->setup(*N);
        }

        GWArithmetic gw(*mod_gwstate);
        GWNum X(gw);
        if (known_factors > *N)
        {
            mod_gwstate->mod(a, res);
            X = res;
        }
        else
            X = a;
        gw.fft(X, X);
        res = X;
    }

    double GWState::ops()
    {
        return gw_get_fft_count(gwdata())*(gwdata()->GENERAL_MMGW_MOD ? 1.0/7.5 : gwdata()->GENERAL_MOD ? 1.0/6 : 1.0/2);
    }

    GWArithmetic::GWArithmetic(GWState& state) : _state(state)
    {
        _careful = new CarefulGWArithmetic(state);
    }

    GWArithmetic::~GWArithmetic()
    {
        if (_careful != nullptr && (GWArithmetic*)_careful != this)
        {
            delete _careful;
            _careful = nullptr;
        }
    }

    void GWArithmetic::alloc(GWNum& a)
    {
        a._gwnum = gwalloc(gwdata());
    }

    void GWArithmetic::free(GWNum& a)
    {
        gwfree(gwdata(), a._gwnum);
        a._gwnum = nullptr;
    }

    void GWArithmetic::copy(const GWNum& a, GWNum& res)
    {
        if (&a == &res)
            return;
        if (*res == nullptr)
            res.arithmetic().alloc(res);
        gwcopy(gwdata(), *a, *res);
    }

    void GWArithmetic::move(GWNum&& a, GWNum& res)
    {
        if (&a == &res)
            return;
        if (res._gwnum != nullptr)
            free(res);
        res._gwnum = a._gwnum;
        a._gwnum = nullptr;
    }

    void GWArithmetic::init(int32_t a, GWNum& res)
    {
        dbltogw(gwdata(), a, *res);
    }

    void GWArithmetic::init(const std::string& a, GWNum& res)
    {
        (popg() = a).to_GWNum(res);
    }

    int GWArithmetic::cmp(const GWNum& a, const GWNum& b)
    {
        if (&a == &b)
            return 0;
        return _state.giants->cmp(popg() = a, popg() = b);
    }

    int GWArithmetic::cmp(const GWNum& a, int32_t b)
    {
        if (b >= 0)
            return _state.giants->cmp(popg() = a, b);
        else
            return _state.giants->cmp(popg() = a, N() + b);
    }

    int GWArithmetic::cmp(const GWNum& a, const Giant& b)
    {
        return _state.giants->cmp(popg() = a, b);
    }

    void GWArithmetic::add(GWNum& a, GWNum& b, GWNum& res)
    {
        add(a, b, res, GWADD_DELAY_NORMALIZE);
    }

    void GWArithmetic::add(GWNum& a, int32_t b, GWNum& res)
    {
        if (*res != *a)
            copy(a, res);
        gwsmalladd(gwdata(), b, *res);
    }

    void GWArithmetic::add(GWNum& a, GWNum& b, GWNum& res, int options)
    {
        gwadd3o(gwdata(), *a, *b, *res, options);
    }

    void GWArithmetic::sub(GWNum& a, GWNum& b, GWNum& res)
    {
        sub(a, b, res, GWADD_DELAY_NORMALIZE);
    }

    void GWArithmetic::sub(GWNum& a, int32_t b, GWNum& res)
    {
        if (*res != *a)
            copy(a, res);
        gwsmalladd(gwdata(), -b, *res);
    }

    void GWArithmetic::sub(GWNum& a, GWNum& b, GWNum& res, int options)
    {
        gwsub3o(gwdata(), *a, *b, *res, options);
    }

    void GWArithmetic::neg(GWNum& a, GWNum& res)
    {
        if (*res != *a)
            copy(a, res);
        gwsmallmul(gwdata(), -1, *res);
    }

    void GWArithmetic::mul(GWNum& a, GWNum& b, GWNum& res)
    {
        mul(a, b, res, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
    }

    void GWArithmetic::mul(GWNum& a, int32_t b, GWNum& res)
    {
        if (*res != *a)
            copy(a, res);
        gwsmallmul(gwdata(), b, *res);
    }

    void GWArithmetic::mul(GWNum& a, GWNum& b, GWNum& res, int options)
    {
        gwmul3(gwdata(), *a, *b, *res, options);
    }

    void GWArithmetic::addsub(GWNum& a, GWNum& b, GWNum& res1, GWNum& res2, int options)
    {
        gwaddsub4o(gwdata(), *a, *b, *res1, *res2, options);
    }

    void GWArithmetic::addmul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options)
    {
        gwaddmul4(gwdata(), *a, *b, *c, *res, options);
    }

    void GWArithmetic::submul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options)
    {
        gwsubmul4(gwdata(), *a, *b, *c, *res, options);
    }

    void GWArithmetic::muladd(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options)
    {
        gwmuladd4(gwdata(), *a, *b, *c, *res, options);
    }

    void GWArithmetic::mulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options)
    {
        gwmulsub4(gwdata(), *a, *b, *c, *res, options);
    }

    void GWArithmetic::mulmuladd(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options)
    {
        gwmulmuladd5(gwdata(), *a, *b, *c, *d, *res, options);
    }

    void GWArithmetic::mulmulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options)
    {
        gwmulmulsub5(gwdata(), *a, *b, *c, *d, *res, options);
    }

    void GWArithmetic::fft(GWNum& a, GWNum& res)
    {
        gwfft(gwdata(), *a, *res);
    }

    void GWArithmetic::unfft(GWNum& a, GWNum& res)
    {
        gwunfft(gwdata(), *a, *res);
    }

    void GWArithmetic::div(GWNum& a, GWNum& b, GWNum& res)
    {
        GWNum inv_b = b;
        inv_b.inv();
        mul(a, inv_b, res);
    }

    void GWArithmetic::div(GWNum& a, int32_t b, GWNum& res)
    {
        GWNum inv_b(a.arithmetic());
        (popg() = b).inv(N()).to_GWNum(inv_b);
        mul(a, inv_b, res);
    }

    void GWArithmetic::gcd(GWNum& a, GWNum& b, GWNum& res)
    {
        (popg() = a).gcd(popg() = b).to_GWNum(res);
    }

    void GWArithmetic::inv(GWNum& a, GWNum& n, GWNum& res)
    {
        (popg() = a).inv(popg() = n).to_GWNum(res);
    }

    void GWArithmetic::inv(GWNum& a, GWNum& res)
    {
        (popg() = a).inv(N()).to_GWNum(res);
    }

    void GWArithmetic::mod(GWNum& a, GWNum& b, GWNum& res)
    {
        ((popg() = a)%(popg() = b)).to_GWNum(res);
    }

    const GWNumWrapper GWArithmetic::wrap(gwnum a)
    {
        return GWNumWrapper(*this, a);
    }

    std::string GWNum::to_string() const
    {
        return (arithmetic().popg() = *this).to_string();
    }

    void CarefulGWArithmetic::add(GWNum& a, int32_t b, GWNum& res)
    {
        GWNum tmp(*this);
        tmp = b;
        unfft(a, a);
        gwadd3o(gwdata(), *a, *tmp, *res, GWADD_FORCE_NORMALIZE);
    }

    void CarefulGWArithmetic::sub(GWNum& a, int32_t b, GWNum& res)
    {
        GWNum tmp(*this);
        tmp = b;
        unfft(a, a);
        gwsub3o(gwdata(), *a, *tmp, *res, GWADD_FORCE_NORMALIZE);
    }

    void CarefulGWArithmetic::add(GWNum& a, GWNum& b, GWNum& res, int options)
    {
        unfft(a, a);
        unfft(b, b);
        gwadd3o(gwdata(), *a, *b, *res, GWADD_FORCE_NORMALIZE);
    }

    void CarefulGWArithmetic::sub(GWNum& a, GWNum& b, GWNum& res, int options)
    {
        unfft(a, a);
        unfft(b, b);
        gwsub3o(gwdata(), *a, *b, *res, GWADD_FORCE_NORMALIZE);
    }

    void CarefulGWArithmetic::mul(GWNum& a, GWNum& b, GWNum& res, int options)
    {
        gwmul3_carefully(gwdata(), *a, *b, *res, options & (GWMUL_ADDINCONST | GWMUL_MULBYCONST));
    }

    void CarefulGWArithmetic::addsub(GWNum& a, GWNum& b, GWNum& res1, GWNum& res2, int options)
    {
        unfft(a, a);
        unfft(b, b);
        gwaddsub4o(gwdata(), *a, *b, *res1, *res2, GWADD_FORCE_NORMALIZE);
    }

    void CarefulGWArithmetic::addmul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options)
    {
        if (*c == *res)
        {
            GWNum tmp = c;
            add(a, b, res);
            mul(res, tmp, res, options);
        }
        else
        {
            add(a, b, res);
            mul(res, c, res, options);
        }
    }

    void CarefulGWArithmetic::submul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options)
    {
        if (*c == *res)
        {
            GWNum tmp = c;
            sub(a, b, res);
            mul(res, tmp, res, options);
        }
        else
        {
            sub(a, b, res);
            mul(res, c, res, options);
        }
    }

    void CarefulGWArithmetic::muladd(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options)
    {
        if (*c == *res || (options & GWMUL_MULBYCONST))
        {
            GWNum tmp(*this);
            unfft(c, tmp);
            if (options & GWMUL_MULBYCONST)
                gwsmallmul(gwdata(), mulbyconst(), *tmp);
            mul(a, b, res, options);
            add(res, tmp, res);
        }
        else
        {
            mul(a, b, res, options);
            add(res, c, res);
        }
    }

    void CarefulGWArithmetic::mulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options)
    {
        if (*c == *res || (options & GWMUL_MULBYCONST))
        {
            GWNum tmp(*this);
            unfft(c, tmp);
            if (options & GWMUL_MULBYCONST)
                gwsmallmul(gwdata(), mulbyconst(), *tmp);
            mul(a, b, res, options);
            sub(res, tmp, res);
        }
        else
        {
            mul(a, b, res, options);
            sub(res, c, res);
        }
    }

    void CarefulGWArithmetic::mulmuladd(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options)
    {
        GWNum tmp(*this);
        mul(c, d, tmp, options & (~GWMUL_ADDINCONST));
        mul(a, b, res, options);
        add(res, tmp, res);
    }

    void CarefulGWArithmetic::mulmulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options)
    {
        GWNum tmp(*this);
        mul(c, d, tmp, options & (~GWMUL_ADDINCONST));
        mul(a, b, res, options);
        sub(res, tmp, res);
    }

    void ReliableGWArithmetic::reset()
    {
        _op = 0;
        _suspect_ops.clear();
        _restart_flag = false;
        _failure_flag = false;
    }

    void ReliableGWArithmetic::restart(int op)
    {
        _op = op;
        _restart_flag = false;
    }

    void ReliableGWArithmetic::mul(GWNum& a, GWNum& b, GWNum& res, int options)
    {
        bool suspect = _suspect_ops.count(_op) != 0;
        if (!suspect)
        {
            // Try a normal operation.
            gwerror_checking(gwdata(), true);
            gwmul3(gwdata(), *a, *b, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                // If roundoff exceeds maximum, mark the operation as suspect.
                gw_clear_maxerr(gwdata());
                _suspect_ops.insert(_op);
                // If either source is the same as destination, we can't redo the operation.
                // Need to restart from a checkpoint.
                if (*res == *a || *res == *b)
                    _restart_flag = true;
                else
                    suspect = true;
            }
        }
        if (suspect)
        {
            // If operation is marked as suspect, make a copy of the source arguments.
            GWNum s1 = a;
            GWNum s2 = b;
            // Retry it with original arguments in case it's a hardware problem.
            gwerror_checking(gwdata(), true);
            gwmul3(gwdata(), *a, *b, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                // Not a hardware problem, retry carefully.
                gw_clear_maxerr(gwdata());
                gwerror_checking(gwdata(), true);
                gwmul3_carefully(gwdata(), *s1, (*a == *b) ? *s1 : *s2, *res, options & (~GWMUL_PRESERVE_S1) & (~GWMUL_PRESERVE_S2));
                gwerror_checking(gwdata(), false);
                if (gw_get_maxerr(gwdata()) > _max_roundoff)
                {
                    // All is lost, need larger FFT.
                    gw_clear_maxerr(gwdata());
                    _failure_flag = true;
                }
            }
            else // It's a hardware problem after all.
                _suspect_ops.erase(_op);
        }
        // Increase operation counter.
        _op++;
    }

    void ReliableGWArithmetic::addmul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options)
    {
        bool suspect = _suspect_ops.count(_op) != 0;
        if (!suspect)
        {
            gwerror_checking(gwdata(), true);
            gwaddmul4(gwdata(), *a, *b, *c, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                gw_clear_maxerr(gwdata());
                _suspect_ops.insert(_op);
                if (*res == *a || *res == *b || *res == *c)
                    _restart_flag = true;
                else
                    suspect = true;
            }
        }
        if (suspect)
        {
            GWNum s1 = a;
            GWNum s2 = b;
            GWNum s3 = c;
            gwerror_checking(gwdata(), true);
            gwaddmul4(gwdata(), *a, *b, *c, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                gw_clear_maxerr(gwdata());
                gwerror_checking(gwdata(), true);
                gwunfft(gwdata(), *s1, *s1);
                gwunfft(gwdata(), *s2, *s2);
                gwadd3o(gwdata(), *s1, *s2, *s1, GWADD_FORCE_NORMALIZE);
                gwmul3_carefully(gwdata(), *s1, *s3, *res, options & (~GWMUL_PRESERVE_S1) & (~GWMUL_PRESERVE_S2));
                gwerror_checking(gwdata(), false);
                if (gw_get_maxerr(gwdata()) > _max_roundoff)
                {
                    gw_clear_maxerr(gwdata());
                    _failure_flag = true;
                }
            }
            else
                _suspect_ops.erase(_op);
        }
        _op++;
    }

    void ReliableGWArithmetic::submul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options)
    {
        bool suspect = _suspect_ops.count(_op) != 0;
        if (!suspect)
        {
            gwerror_checking(gwdata(), true);
            gwsubmul4(gwdata(), *a, *b, *c, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                gw_clear_maxerr(gwdata());
                _suspect_ops.insert(_op);
                if (*res == *a || *res == *b || *res == *c)
                    _restart_flag = true;
                else
                    suspect = true;
            }
        }
        if (suspect)
        {
            GWNum s1 = a;
            GWNum s2 = b;
            GWNum s3 = c;
            gwerror_checking(gwdata(), true);
            gwsubmul4(gwdata(), *a, *b, *c, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                gw_clear_maxerr(gwdata());
                gwerror_checking(gwdata(), true);
                gwunfft(gwdata(), *s1, *s1);
                gwunfft(gwdata(), *s2, *s2);
                gwsub3o(gwdata(), *s1, *s2, *s1, GWADD_FORCE_NORMALIZE);
                gwmul3_carefully(gwdata(), *s1, *s3, *res, options & (~GWMUL_PRESERVE_S1) & (~GWMUL_PRESERVE_S2));
                gwerror_checking(gwdata(), false);
                if (gw_get_maxerr(gwdata()) > _max_roundoff)
                {
                    gw_clear_maxerr(gwdata());
                    _failure_flag = true;
                }
            }
            else
                _suspect_ops.erase(_op);
        }
        _op++;
    }

    void ReliableGWArithmetic::muladd(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options)
    {
        bool suspect = _suspect_ops.count(_op) != 0;
        if (!suspect)
        {
            gwerror_checking(gwdata(), true);
            gwmuladd4(gwdata(), *a, *b, *c, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                gw_clear_maxerr(gwdata());
                _suspect_ops.insert(_op);
                if (*res == *a || *res == *b || *res == *c)
                    _restart_flag = true;
                else
                    suspect = true;
            }
        }
        if (suspect)
        {
            GWNum s1 = a;
            GWNum s2 = b;
            GWNum s3 = c;
            gwerror_checking(gwdata(), true);
            gwmuladd4(gwdata(), *a, *b, *c, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                gw_clear_maxerr(gwdata());
                gwerror_checking(gwdata(), true);
                gwmul3_carefully(gwdata(), *s1, (*a == *b) ? *s1 : *s2, *res, options & (~GWMUL_PRESERVE_S1) & (~GWMUL_PRESERVE_S2) & (~GWMUL_STARTNEXTFFT));
                gwunfft(gwdata(), *s3, *s3);
                if (options & GWMUL_MULBYCONST)
                    gwsmallmul(gwdata(), mulbyconst(), *s3);
                gwadd3o(gwdata(), *res, *s3, *res, GWADD_FORCE_NORMALIZE);
                gwerror_checking(gwdata(), false);
                if (gw_get_maxerr(gwdata()) > _max_roundoff)
                {
                    gw_clear_maxerr(gwdata());
                    _failure_flag = true;
                }
            }
            else
                _suspect_ops.erase(_op);
        }
        _op++;
    }

    void ReliableGWArithmetic::mulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options)
    {
        bool suspect = _suspect_ops.count(_op) != 0;
        if (!suspect)
        {
            gwerror_checking(gwdata(), true);
            gwmulsub4(gwdata(), *a, *b, *c, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                gw_clear_maxerr(gwdata());
                _suspect_ops.insert(_op);
                if (*res == *a || *res == *b || *res == *c)
                    _restart_flag = true;
                else
                    suspect = true;
            }
        }
        if (suspect)
        {
            GWNum s1 = a;
            GWNum s2 = b;
            GWNum s3 = c;
            gwerror_checking(gwdata(), true);
            gwmulsub4(gwdata(), *a, *b, *c, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                gw_clear_maxerr(gwdata());
                gwerror_checking(gwdata(), true);
                gwmul3_carefully(gwdata(), *s1, (*a == *b) ? *s1 : *s2, *res, options & (~GWMUL_PRESERVE_S1) & (~GWMUL_PRESERVE_S2) & (~GWMUL_STARTNEXTFFT));
                gwunfft(gwdata(), *s3, *s3);
                if (options & GWMUL_MULBYCONST)
                    gwsmallmul(gwdata(), mulbyconst(), *s3);
                gwsub3o(gwdata(), *res, *s3, *res, GWADD_FORCE_NORMALIZE);
                gwerror_checking(gwdata(), false);
                if (gw_get_maxerr(gwdata()) > _max_roundoff)
                {
                    gw_clear_maxerr(gwdata());
                    _failure_flag = true;
                }
            }
            else
                _suspect_ops.erase(_op);
        }
        _op++;
    }

    void ReliableGWArithmetic::mulmuladd(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options)
    {
        bool suspect = _suspect_ops.count(_op) != 0;
        if (!suspect)
        {
            gwerror_checking(gwdata(), true);
            gwmulmuladd5(gwdata(), *a, *b, *c, *d, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                gw_clear_maxerr(gwdata());
                _suspect_ops.insert(_op);
                if (*res == *a || *res == *b || *res == *c || *res == *d)
                    _restart_flag = true;
                else
                    suspect = true;
            }
        }
        if (suspect)
        {
            GWNum s1 = a;
            GWNum s2 = b;
            GWNum s3 = c;
            GWNum s4 = d;
            gwerror_checking(gwdata(), true);
            gwmulmuladd5(gwdata(), *a, *b, *c, *d, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                gw_clear_maxerr(gwdata());
                gwerror_checking(gwdata(), true);
                gwmul3_carefully(gwdata(), *s1, (*a == *b) ? *s1 : *s2, *res, options & (~GWMUL_PRESERVE_S1) & (~GWMUL_PRESERVE_S2) & (~GWMUL_STARTNEXTFFT));
                gwmul3_carefully(gwdata(), *s3, (*c == *d) ? *s3 : *s4, *s3, options & (~GWMUL_PRESERVE_S1) & (~GWMUL_PRESERVE_S2) & (~GWMUL_STARTNEXTFFT) & (~GWMUL_ADDINCONST));
                gwadd3o(gwdata(), *res, *s3, *res, GWADD_FORCE_NORMALIZE);
                gwerror_checking(gwdata(), false);
                if (gw_get_maxerr(gwdata()) > _max_roundoff)
                {
                    gw_clear_maxerr(gwdata());
                    _failure_flag = true;
                }
            }
            else
                _suspect_ops.erase(_op);
        }
        _op++;
    }

    void ReliableGWArithmetic::mulmulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options)
    {
        bool suspect = _suspect_ops.count(_op) != 0;
        if (!suspect)
        {
            gwerror_checking(gwdata(), true);
            gwmulmulsub5(gwdata(), *a, *b, *c, *d, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                gw_clear_maxerr(gwdata());
                _suspect_ops.insert(_op);
                if (*res == *a || *res == *b || *res == *c || *res == *d)
                    _restart_flag = true;
                else
                    suspect = true;
            }
        }
        if (suspect)
        {
            GWNum s1 = a;
            GWNum s2 = b;
            GWNum s3 = c;
            GWNum s4 = d;
            gwerror_checking(gwdata(), true);
            gwmulmulsub5(gwdata(), *a, *b, *c, *d, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                gw_clear_maxerr(gwdata());
                gwerror_checking(gwdata(), true);
                gwmul3_carefully(gwdata(), *s1, (*a == *b) ? *s1 : *s2, *res, options & (~GWMUL_PRESERVE_S1) & (~GWMUL_PRESERVE_S2) & (~GWMUL_STARTNEXTFFT));
                gwmul3_carefully(gwdata(), *s3, (*c == *d) ? *s3 : *s4, *s3, options & (~GWMUL_PRESERVE_S1) & (~GWMUL_PRESERVE_S2) & (~GWMUL_STARTNEXTFFT) & (~GWMUL_ADDINCONST));
                gwsub3o(gwdata(), *res, *s3, *res, GWADD_FORCE_NORMALIZE);
                gwerror_checking(gwdata(), false);
                if (gw_get_maxerr(gwdata()) > _max_roundoff)
                {
                    gw_clear_maxerr(gwdata());
                    _failure_flag = true;
                }
            }
            else
                _suspect_ops.erase(_op);
        }
        _op++;
    }

    SerializedGWNum& SerializedGWNum::operator = (const GWNum& a)
    {
        int len;
        if (gwserialize(a.arithmetic().gwdata(), *a, NULL, 0, &len) != GWERROR_MALLOC)
            throw ArithmeticException();
        _data.resize(-len);
        if (gwserialize(a.arithmetic().gwdata(), *a, _data.data(), (int)_data.size(), &len) != 0 || _data.size() != len)
            throw ArithmeticException();
        return *this;
    }

    void SerializedGWNum::to_GWNum(GWNum& a) const
    {
        if (empty())
            a = 0;
        else
            gwdeserialize(a.arithmetic().gwdata(), _data.data(), (int)_data.size(), *a);
    }
}

int gwconvert(
    gwhandle *gwdata_s,	/* Handle initialized by gwsetup */
    gwhandle *gwdata_d,	/* Handle initialized by gwsetup */
    gwnum	s,
    gwnum	d)
{
    int	err_code;
    unsigned long limit;

    ASSERTG(gwdata_s->k == gwdata_d->k && gwdata_s->b == gwdata_d->b && gwdata_s->c == gwdata_d->c);

    /* Make sure data is not FFTed.  Caller should really try to avoid this scenario. */

    if (FFT_state(s) != NOT_FFTed) gwunfft(gwdata_s, s, s);

    /* Initialize the header */

    unnorms(d) = 0.0f;						/* Unnormalized adds count */
    *(uint32_t *)((char*)d - 28) = 0;					/* Has-been-pre-ffted flag */
    *(double *)((char*)d - 16) = 0.0;
    *(double *)((char*)d - 24) = 0.0;


    /* If this is a general-purpose mod, then only convert the needed words */
    /* which will be less than half the FFT length.  If this is a zero padded */
    /* FFT, then only convert a little more than half of the FFT data words. */
    /* For a DWT, convert all the FFT data. */

    if (gwdata_s->GENERAL_MOD) limit = gwdata_s->GW_GEN_MOD_MAX + 3;
    else if (gwdata_s->ZERO_PADDED_FFT) limit = gwdata_s->FFTLEN / 2 + 4;
    else limit = gwdata_s->FFTLEN;

    /* GENERAL_MOD has some strange cases we must handle.  In particular the */
    /* last fft word translated can be 2^bits and the next word could be -1, */
    /* this must be translated into zero, zero. */

    if (gwdata_s->GENERAL_MOD) {
        long	val, prev_val;
        while (limit < gwdata_s->FFTLEN) {
            err_code = get_fft_value(gwdata_s, s, limit, &val);
            if (err_code) return (err_code);
            if (val == -1 || val == 0) break;
            limit++;
            ASSERTG(limit <= gwdata_s->FFTLEN / 2 + 2);
            if (limit > gwdata_s->FFTLEN / 2 + 2) return (GWERROR_INTERNAL + 9);
        }
        while (limit > 1) {		/* Find top word */
            err_code = get_fft_value(gwdata_s, s, limit-1, &prev_val);
            if (err_code) return (err_code);
            if (val != prev_val || val < -1 || val > 0) break;
            limit--;
        }
        limit++;
    }

    /* If base is 2 we can simply copy the bits out of each FFT word */

    if (gwdata_s->b == 2) {
        long val;
        int64_t value;
        unsigned long i_s;
        unsigned long i_d, limit_d;
        int	bits_in_value, bits, bits1, bits2;
        long mask1, mask2, mask1i, mask2i;

        // Figure out how many FFT words we will need to set
        limit_d = gwdata_d->FFTLEN;
        if (gwdata_d->ZERO_PADDED_FFT && limit_d > gwdata_d->FFTLEN / 2 + 4) limit_d = gwdata_d->FFTLEN / 2 + 4;

        bits1 = gwdata_d->NUM_B_PER_SMALL_WORD;
        bits2 = bits1 + 1;
        mask1 = (1L << bits1) - 1;
        mask2 = (1L << bits2) - 1;
        mask1i = ~mask1;
        mask2i = ~mask2;

        /* Collect bits until we have all of them */

        value = 0;
        bits_in_value = 0;
        i_d = 0;

        for (i_s = 0; i_s < limit; i_s++) {
            err_code = get_fft_value(gwdata_s, s, i_s, &val);
            if (err_code) return (err_code);
            bits = gwdata_s->NUM_B_PER_SMALL_WORD;
            if (is_big_word(gwdata_s, i_s)) bits++;
            value += (int64_t)val << bits_in_value;
            bits_in_value += bits;

            for (; i_d < limit_d; i_d++) {
                if (i_d == limit_d - 1) {
                    if (i_s < limit - 1) break;
                    val = (long)value;
                }
                else {
                    int	big_word;
                    big_word = is_big_word(gwdata_d, i_d);
                    bits = big_word ? bits2 : bits1;
                    if (i_s < limit - 1 && bits > bits_in_value) break;
                    if (value >= 0)
                        val = (long)value & (big_word ? mask2 : mask1);
                    else {
                        val = (long)value | (big_word ? mask2i : mask1i);
                        value -= val;
                    }
                }
                set_fft_value(gwdata_d, d, i_d, val);
                value >>= bits;
                bits_in_value -= bits;
            }
        }

        /* Clear the upper words */

        for (; i_d < gwdata_d->FFTLEN; i_d++) {
            set_fft_value(gwdata_d, d, i_d, 0);
        }
    }

    /* Otherwise (base is not 2) we must do a radix conversion */

    else {
        //err_code = nonbase2_gwtogiant(gwdata, gg, v);
        //if (err_code) return (err_code);
        return -1;
    }

    /* Return success */

    gwdata_s->read_count += 1;
    gwdata_d->write_count += 1;
    return (0);
}

#define GWSERIALIZE_HEADER_SIZE 6
#define GWSERIALIZE_FLAG_IRRATIONAL 1
#define GWSERIALIZE_FLAG_GENERALMOD 2
#define GWSERIALIZE_FLAG_MMGW_MOD 4

int gwserialize(
    gwhandle *gwdata,  /* Handle initialized by gwsetup */
    gwnum gg,          /* Source gwnum */
    uint32_t *array,   /* Array to contain the serialized value */
    int arraylen,	   /* Maximum size of the array */
    int *arrayused)    /* Size of the serialized value */
{
    int	err_code;
    unsigned long limit;

/* Make sure data is not FFTed.  Caller should really try to avoid this scenario. */

	if (FFT_state (gg) != NOT_FFTed) gwunfft (gwdata, gg, gg);

/* Set result to zero in case of error.  If caller does not check the returned error code */
/* a result value of zero is less likely to cause problems/crashes. */

    *arrayused = 0;

/* If this is a general-purpose mod, then only convert the needed words */
/* which will be less than half the FFT length.  If this is a zero padded */
/* FFT, then only convert a little more than half of the FFT data words. */
/* For a DWT, convert all the FFT data. */

	if (gwdata->GENERAL_MOD) limit = gwdata->GW_GEN_MOD_MAX + 3;
    //else if (gwdata->GENERAL_MMGW_MOD) limit = 2*gwdata->FFTLEN;
    else if (gwdata->ZERO_PADDED_FFT) limit = gwdata->FFTLEN / 2 + 4;
	else limit = gwdata->FFTLEN;

/* GENERAL_MOD has some strange cases we must handle.  In particular the */
/* last fft word translated can be 2^bits and the next word could be -1, */
/* this must be translated into zero, zero. */

	if (gwdata->GENERAL_MOD) {
		long	val, prev_val;
		while (limit < gwdata->FFTLEN) {
			err_code = get_fft_value (gwdata, gg, limit, &val);
			if (err_code) return (err_code);
			if (val == -1 || val == 0) break;
			limit++;
			ASSERTG (limit <= gwdata->FFTLEN / 2 + 2);
			if (limit > gwdata->FFTLEN / 2 + 2) return (GWERROR_INTERNAL + 9);
		}
		while (limit > 1) {		/* Find top word */
			err_code = get_fft_value (gwdata, gg, limit-1, &prev_val);
			if (err_code) return (err_code);
			if (val != prev_val || val < -1 || val > 0) break;
			limit--;
		}
		limit++;
	}

    if ((int)limit + GWSERIALIZE_HEADER_SIZE > arraylen)
    {
        *arrayused = -((int)limit + GWSERIALIZE_HEADER_SIZE);
        return GWERROR_MALLOC;
    }

    array[0] = 0;
    array[1] = 0;
    array[2] = 0;
    array[3] = 0;
    array[4] = 0;
    array[5] = 0;

    int32_t val;
    unsigned long i;
    gwiter iter;
    for (i = 0, gwiter_init_zero(gwdata, &iter, gg); i < limit; i++, gwiter_next(&iter))
    {
        /*if (gwdata->GENERAL_MMGW_MOD && i == gwdata->FFTLEN)
            gwiter_init_zero(gwdata->negacyclic_gwdata, &iter, negacyclic_gwnum(gwdata, gg));*/
        err_code = gwiter_get_fft_value(&iter, &val);
        if (err_code) return (err_code);
        if (!gwdata->RATIONAL_FFT)
        {
            val <<= 1;
            if (gwiter_is_big_word(&iter))
                val |= 1;
        }
        array[i + GWSERIALIZE_HEADER_SIZE] = (uint32_t)val;
    }

    uint8_t* header = (uint8_t*)array;
    if (!gwdata->RATIONAL_FFT)
        header[0] |= GWSERIALIZE_FLAG_IRRATIONAL;
    if (gwdata->GENERAL_MOD)
        header[0] |= GWSERIALIZE_FLAG_GENERALMOD;
    if (gwdata->GENERAL_MMGW_MOD)
        header[0] |= GWSERIALIZE_FLAG_MMGW_MOD;
    header[1] = (uint8_t)gwdata->NUM_B_PER_SMALL_WORD;
    array[1] = (uint32_t)gwdata->b;
    array[2] = (uint32_t)gwdata->FFTLEN;
    array[3] = (uint32_t)(gwdata->GENERAL_MMGW_MOD ? gwdata->n : limit);
    *(uint64_t*)(array + 4) = (uint64_t)gwdata->k;
    *arrayused = (int)limit + GWSERIALIZE_HEADER_SIZE;

/* Return success */

	gwdata->read_count += 1;
	return (0);
}

void gwdeserialize(
    gwhandle *gwdata,	    /* Handle initialized by gwsetup */
    const uint32_t *array,  /* Array containing the binary value */
    int arraylen,	        /* Length of the array */
    gwnum g)	    	    /* Destination gwnum */
{
    unsigned long i, limit;
    int32_t val;
    gwiter iter;
    uint8_t* header = (uint8_t*)array;

    ASSERTG(arraylen >= GWSERIALIZE_HEADER_SIZE && array[2] != 0);

    limit = (header[0] & GWSERIALIZE_FLAG_MMGW_MOD) ? array[2] : array[3];

    if (!gwdata->RATIONAL_FFT == ((header[0] & GWSERIALIZE_FLAG_IRRATIONAL) != 0) &&
        gwdata->GENERAL_MOD == ((header[0] & GWSERIALIZE_FLAG_GENERALMOD) != 0) &&
        gwdata->GENERAL_MMGW_MOD == ((header[0] & GWSERIALIZE_FLAG_MMGW_MOD) != 0) &&
        gwdata->NUM_B_PER_SMALL_WORD == header[1] &&
        gwdata->b == array[1] &&
        gwdata->FFTLEN == array[2] &&
        (!gwdata->GENERAL_MMGW_MOD || gwdata->n == array[3]))
    {
        for (i = 0, gwiter_init_write_only(gwdata, &iter, g); i < limit; i++, gwiter_next(&iter))
        {
            /*if (gwdata->GENERAL_MMGW_MOD && i == gwdata->FFTLEN)
                gwiter_init_write_only(gwdata->negacyclic_gwdata, &iter, negacyclic_gwnum(gwdata, g));*/
            val = (int32_t)array[i + GWSERIALIZE_HEADER_SIZE];
            if ((header[0] & GWSERIALIZE_FLAG_IRRATIONAL))
            {
                ASSERTG(((val & 1) != 0) == gwiter_is_big_word(&iter));
                val >>= 1;
            }
            gwiter_set_fft_value(&iter, val);
        }

/* Clear the upper words */

		for ( ; i < gwdata->FFTLEN; i++, gwiter_next(&iter))
            gwiter_set_fft_value(&iter, 0);
    }

    else if (gwdata->b == array[1] && gwdata->FFTLEN < array[2] &&
             (gwdata->GENERAL_MOD != 0) == ((header[0] & GWSERIALIZE_FLAG_GENERALMOD) != 0) &&
             (gwdata->GENERAL_MMGW_MOD != 0) == ((header[0] & GWSERIALIZE_FLAG_MMGW_MOD) != 0))
        throw arithmetic::ArithmeticException("Can't deserialize to smaller transform.");

    else if (gwdata->b == array[1])
    {
        std::vector<int64_t> pow;
        int64_t value;
        int	bs_in_value, bs;
        unsigned long i_d, limit_d;
        bool more_data;

        pow.push_back(1);
        pow.push_back(gwdata->b);
        value = 0;
        bs_in_value = 0;
        i_d = 0;

        // Figure out how many FFT words we will need to set
        
        if (gwdata->GENERAL_MOD) limit_d = gwdata->GW_GEN_MOD_MAX + 3;
        //else if (gwdata->GENERAL_MMGW_MOD) limit_d = 2*gwdata->FFTLEN;
        else if (gwdata->ZERO_PADDED_FFT) limit_d = gwdata->FFTLEN / 2 + 4;
        else limit_d = gwdata->FFTLEN;

        /* Collect bs until we have all of them */

        gwiter_init_write_only(gwdata, &iter, g);
        for (i = 0; i < limit; i++)
        {
            bs = header[1];
            val = (int32_t)array[i + GWSERIALIZE_HEADER_SIZE];
            if ((header[0] & GWSERIALIZE_FLAG_IRRATIONAL))
            {
                if ((val & 1))
                    bs++;
                val >>= 1;
            }
            while (pow.size() <= bs_in_value)
                pow.push_back(pow.back()*gwdata->b);
            value += val*pow[bs_in_value];
            bs_in_value += bs;
            more_data = /*(!gwdata->GENERAL_MMGW_MOD || i != array[2] - 1) &&*/ i != limit - 1;

            for (; i_d < limit_d; i_d++, gwiter_next(&iter)) {
                /*if (gwdata->GENERAL_MMGW_MOD && i_d == gwdata->FFTLEN - 1)
                {
                    if (more_data)
                        break;
                    val = (int32_t)value;
                    value = 0;
                    bs_in_value = 0;
                    gwiter_set_fft_value(&iter, val);
                    i_d++;
                    gwiter_init_write_only(gwdata->negacyclic_gwdata, &iter, negacyclic_gwnum(gwdata, g));
                    break;
                }
                else*/ if (i_d == limit_d - 1) {
                    if (more_data) break;
                    val = (int32_t)value;
                }
                else {
                    bs = gwdata->NUM_B_PER_SMALL_WORD;
                    if (gwiter_is_big_word(&iter))
                        bs++;
                    if (more_data && bs > bs_in_value) break;
                    while (pow.size() <= bs)
                        pow.push_back(pow.back()*gwdata->b);
                    val = (int32_t)(value%pow[bs]);
                    value /= pow[bs];
                    if (val > pow[bs]/2)
                    {
                        val -= (int32_t)pow[bs];
                        value++;
                    }
                    else if (val < -pow[bs]/2)
                    {
                        val += (int32_t)pow[bs];
                        value--;
                    }
                    bs_in_value -= bs;
                }
                gwiter_set_fft_value(&iter, val);
            }
        }

        /* Clear the upper words */

        for (; i_d < gwdata->FFTLEN; i_d++, gwiter_next(&iter)) {
            gwiter_set_fft_value(&iter, 0);
        }

        if (gwdata->GENERAL_MOD && ((header[0] & GWSERIALIZE_FLAG_GENERALMOD) == 0) && *(uint64_t*)(array + 4) > 1)
        {
            if (*(uint64_t*)(array + 4) < GWSMALLMUL_MAX)
                gwsmallmul(gwdata, (double)*(uint64_t*)(array + 4), g);
            else
            {
                if (convert_factor == NULL)
                {
                    convert_factor = gwalloc(gwdata);
                    binarytogw(gwdata, array + 4, 2, convert_factor);
                }
                gwmul3_carefully(gwdata, g, convert_factor, g, 0);
            }
        }

        if (!gwdata->GENERAL_MOD && ((header[0] & GWSERIALIZE_FLAG_GENERALMOD) != 0) && gwdata->k > 1.0)
        {
            if (convert_factor == NULL)
            {
                giant Nalloc = NULL;
                giant N = gwdata->GW_MODULUS;
                if (N == NULL)
                {
                    Nalloc = popg(&gwdata->gdata, ((unsigned long)gwdata->bit_length >> 5) + 2);
                    N = Nalloc;
                    ultog(gwdata->b, N);
                    power(N, gwdata->n);
                    dblmulg(gwdata->k, N);
                    iaddg(gwdata->c, N);
                }
                giant invK = popg(&gwdata->gdata, ((unsigned long)gwdata->bit_length >> 5) + 2);
                dbltog(gwdata->k, invK);
                invg(N, invK);

                convert_factor = gwalloc(gwdata);
                gianttogw(gwdata, invK, convert_factor);

                if (Nalloc != NULL)
                    pushg(&gwdata->gdata, 1);
                pushg(&gwdata->gdata, 1);
            }

            gwmul3_carefully(gwdata, g, convert_factor, g, 0);
        }


        if (gwdata->GENERAL_MMGW_MOD && ((header[0] & GWSERIALIZE_FLAG_MMGW_MOD) == 0))
        {
            gwmul3_carefully(gwdata, g, gwdata->R2_4, g, 0);
            if (*(uint64_t*)(array + 4) > 1)
            {
                if (*(uint64_t*)(array + 4) < GWSMALLMUL_MAX)
                    gwsmallmul(gwdata, (double)*(uint64_t*)(array + 4), g);
                else
                {
                    if (convert_factor == NULL)
                    {
                        convert_factor = gwalloc(gwdata);
                        binarytogw(gwdata, array + 4, 2, convert_factor);
                    }
                    gwmul3_carefully(gwdata, g, convert_factor, g, 0);
                }
            }
        }

        if (((header[0] & GWSERIALIZE_FLAG_MMGW_MOD) != 0) &&
            (!gwdata->GENERAL_MMGW_MOD || (gwdata->GENERAL_MMGW_MOD && gwdata->n != array[3])))
        {
            if (convert_factor == NULL)
            {
                giant Nalloc = NULL;
                giant N = gwdata->GW_MODULUS;
                if (N == NULL)
                {
                    Nalloc = popg(&gwdata->gdata, ((unsigned long)gwdata->bit_length >> 5) + 2);
                    N = Nalloc;
                    ultog(gwdata->b, N);
                    power(N, gwdata->n);
                    dblmulg(gwdata->k, N);
                    iaddg(gwdata->c, N);
                }
                giant Rold = allocgiant((array[3] >> 5) + 4);
                itog(1, Rold);
                gshiftleft(array[3], Rold);
                sladdg(-1, Rold);
                dblmulg(gwdata->k, Rold);
                invg(N, Rold);

                convert_factor = gwalloc(gwdata);
                if (!gwdata->GENERAL_MMGW_MOD)
                {
                    dblmulg(2, Rold);
                    modg(N, Rold);
                    gianttogw(gwdata, Rold, convert_factor);
                }
                if (gwdata->GENERAL_MMGW_MOD && gwdata->n != array[3])
                {
                    giant R = allocgiant((((array[3] > gwdata->n ? array[3] : gwdata->n) + gwdata->n) >> 5) + 2);
                    itog(1, R);
                    gshiftleft(gwdata->n, R);
                    sladdg(-1, R);
                    squareg(R);
                    if (R->n[0] & 1)
                        addg(N, R);
                    gshiftright(1, R);
                    modg(N, R);
                    mulg(Rold, R);
                    modg(N, R);
                    gianttogw(gwdata->negacyclic_gwdata, R, convert_factor);
                    free(R);
                }

                if (Nalloc != NULL)
                    pushg(&gwdata->gdata, 1);
                free(Rold);
            }

            gwmul3_carefully(gwdata, g, convert_factor, g, 0);
        }
    }

    else // Radix conversion is necessary
        throw arithmetic::ArithmeticException("Can't deserialize.");

    /* Clear various flags, update counts */

    unnorms(g) = 0.0f;		/* Set unnormalized add counter */
    FFT_state(g) = NOT_FFTed;	/* Clear has been completely/partially FFTed flag */
    gwdata->write_count += 1;
}
