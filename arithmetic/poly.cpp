#define GDEBUG
#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include "gwnum.h"
#include "poly.h"
#include "exception.h"

namespace arithmetic
{
    GWNum Poly::eval(GWNum& x)
    {
        GWNum res(pm().gw());
        if (size() == 0 && monic())
            res = 1;
        else if (size() == 0)
            res = 0;
        else
            pm().gw().unfft((GWNum&)at(0), res);
        std::unique_ptr<GWNum> t;
        std::unique_ptr<GWNum> p;
        int i, d;
        d = degree();
        for (i = 1; i <= d; i++)
        {
            if (i == 2)
            {
                p.reset(new GWNum(pm().gw()));
                pm().gw().mul(x, x, *p, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
            }
            if (i > 2)
                pm().gw().mul(x, *p, *p, GWMUL_STARTNEXTFFT);
            if (i < size() && _poly[i] != nullptr)
            {
                GWNumWrapper a = at(i);
                if (!t)
                    t.reset(new GWNum(pm().gw()));
                if (i == 1)
                    pm().gw().mul(x, a, *t, GWMUL_FFT_S1 | GWMUL_PRESERVE_S2);
                else
                    pm().gw().mul(*p, a, *t, GWMUL_PRESERVE_S2);
                pm().gw().add(res, *t, res, i < d ? GWADD_DELAY_NORMALIZE : GWADD_FORCE_NORMALIZE);
            }
        }
        if (monic() && d > 1)
            pm().gw().add(res, *p, res, GWADD_FORCE_NORMALIZE);
        if (monic() && d == 1)
            pm().gw().add(res, x, res, GWADD_FORCE_NORMALIZE);
        return res;
    }

    Poly Poly::reciprocal(int precision, int options)
    {
        Poly res(pm(), precision, monic());
        pm().reciprocal(*this, res, options);
        return res;
    }

    int PolyMult::L3_CACHE_MB = 2500;

    int PolyMult::max_polymult_output(GWState& state)
    {
        int max_output;
        for (max_output = 2; max_output < (1 << 30) && gw_passes_safety_margin(state.gwdata(), polymult_safety_margin(max_output, max_output)); max_output <<= 1);
        return max_output;
    }

    PolyMult::PolyMult(GWArithmetic& gw, int max_threads) : _gw(gw)
    {
        GWASSERT(gw.state().polymult);
        _max_output = max_polymult_output(gw.state());
        polymult_init(pmdata(), gw.gwdata());
        polymult_set_max_num_threads(pmdata(), max_threads);
        pmdata()->L3_CACHE_SIZE = L3_CACHE_MB;
    }

    PolyMult::~PolyMult()
    {
        polymult_done(pmdata());
    }

    void PolyMult::set_threads(int threads)
    {
        polymult_set_num_threads(pmdata(), threads);
    }

    void PolyMult::alloc(Poly& a, int size)
    {
        if (a._freeable)
        {
            for (size_t i = size; i < a.size(); i++)
                gwfree(gw().gwdata(), a._poly[i]);
            a._poly.resize(size);
            for (auto it = a._poly.begin(); it != a._poly.end(); it++)
                if (*it == nullptr)
                    *it = gwalloc(gw().gwdata());
        }
        else
        {
            if (a._poly.size() < size)
                GWASSERT(0);
            a._poly.resize(size);
        }
        if (a._cache != nullptr)
            gwfree_array(gw().gwdata(), a._cache);
        a._cache = nullptr;
    }

    void PolyMult::free(Poly& a)
    {
        if (a._freeable)
            for (auto it = a._poly.begin(); it != a._poly.end(); it++)
                gwfree(gw().gwdata(), *it);
        a._poly.clear();
        if (a._cache != nullptr)
            gwfree_array(gw().gwdata(), a._cache);
        a._cache = nullptr;
        a._monic = false;
        a._freeable = true;
    }

    void PolyMult::copy(const Poly& a, Poly& res)
    {
        GWASSERT(!a.preprocessed());
        res.pm().alloc(res, a.size());
        for (size_t i = 0; i < a.size(); i++)
            gwcopy(gw().gwdata(), a._poly[i], res._poly[i]);
        res._monic = a._monic;
    }

    void PolyMult::move(Poly&& a, Poly& res)
    {
        if (!res.empty())
            res.pm().free(res);
        res._poly = std::move(a._poly);
        res._monic = a._monic;
        res._freeable = a._freeable;
        res._cache = a._cache;
        res._cache_size = a._cache_size;
        a._cache = nullptr;
        a._monic = false;
        a._freeable = true;
    }

    void PolyMult::fft(const Poly& a, Poly& res)
    {
        GWASSERT(!a.preprocessed());
        res.pm().alloc(res, a.size());
        for (size_t i = 0; i < a.size(); i++)
            gwfft(gw().gwdata(), a._poly[i], res._poly[i]);
        res._monic = a._monic;
    }

    void PolyMult::init(bool monic, Poly& res)
    {
        if (!res.empty())
            res.pm().free(res);
        res._monic = monic;
    }

    void PolyMult::init(const GWNum& a, bool monic, Poly& res)
    {
        res.pm().alloc(res, 1);
        gwcopy(gw().gwdata(), *a, res._poly[0]);
        res._monic = monic;
    }

    void PolyMult::init(GWNum&& a, bool monic, Poly& res)
    {
        if (!res.empty())
            res.pm().free(res);
        res._poly.push_back(a._gwnum);
        a._gwnum = nullptr;
        res._monic = monic;
        res._freeable = dynamic_cast<GWNumWrapper*>(&a) == nullptr;
    }

    void PolyMult::init(gwnum* data, size_t size, bool freeable, bool monic, Poly& res)
    {
        if (!res.empty())
            res.pm().free(res);
        res._poly.resize(size);
        memcpy(res._poly.data(), data, sizeof(gwnum)*size);
        res._monic = monic;
        res._freeable = freeable;
    }

    void PolyMult::mul(Poly& a, Poly& b, Poly& res, int options)
    {
        if (a.degree() < b.degree() || a.size() < b.size())
        {
            mul(b, a, res, options);
            return;
        }
        if (b.degree() < 0)
        {
            res.pm().free(res);
            return;
        }
        if (b.size() == 0 && !a.preprocessed()) // monic
        {
            if (&a != &res)
                copy(a, res);
            return;
        }

        int sa = a.size();
        int sb = b.size();
        res.pm().alloc(res, sa + sb - (a.monic() || b.monic() ? 0 : 1));

        if (sb == 1 && !a.preprocessed() && !b.preprocessed())
        {
            if (a.monic() && b.monic())
                gwadd3o(gw().gwdata(), a._poly[sa - 1], b._poly[0], res._poly[sa], GWADD_FORCE_NORMALIZE);
            else if (a.monic())
                gwcopy(gw().gwdata(), b._poly[0], res._poly[sa]);
            else if (b.monic())
                gwcopy(gw().gwdata(), a._poly[sa - 1], res._poly[sa]);
            for (int i = sa - 1; i >= 0; i--)
                if (b.monic() && i > 0)
                    gwmuladd4(gw().gwdata(), a._poly[i], b._poly[0], a._poly[i - 1], res._poly[i], (options & POLYMULT_STARTNEXTFFT ? GWMUL_STARTNEXTFFT : 0) | GWMUL_FFT_S2 | GWMUL_FFT_S3);
                else
                    gwmul3(gw().gwdata(), a._poly[i], b._poly[0], res._poly[i], (options & POLYMULT_STARTNEXTFFT ? GWMUL_STARTNEXTFFT : 0) | (i > 0 ? GWMUL_FFT_S2 : 0));
            res._monic = a.monic() && b.monic();
            return;
        }
        
        polymult(pmdata(), a.data(), sa, b.data(), sb, res.data(), res.size(), options | (a.monic() ? POLYMULT_INVEC1_MONIC : 0) | (b.monic() ? POLYMULT_INVEC2_MONIC : 0)/* | (pmdata()->num_threads > 1 ? POLYMULT_NO_UNFFT : 0)*/);
        /*if (pmdata()->num_threads > 1)
        {
            if ((options & POLYMULT_NEXTFFT))
                poly_unfft_fft_coefficients(pmdata(), res.data(), res.size());
            else
                poly_unfft_coefficients(pmdata(), res.data(), res.size());
        }*/

        res._monic = a.monic() && b.monic();
    }

    void PolyMult::mul(Poly&& a, Poly&& b, Poly& res, int options)
    {
        GWASSERT(!a._freeable || &a.pm().gw() == &res.pm().gw());
        GWASSERT(!b._freeable || &b.pm().gw() == &res.pm().gw());
        GWASSERT(!a.preprocessed());
        GWASSERT(!b.preprocessed());
        GWASSERT(a._freeable == b._freeable);
        if (a.degree() < b.degree() || a.size() < b.size())
        {
            mul(std::move(b), std::move(a), res, options);
            return;
        }
        if (b.degree() < 0)
        {
            res.pm().free(res);
            return;
        }
        if (&b == &res)
        {
            Poly tmp(*this, 0, false);
            move(std::move(b), tmp);
            mul(std::move(a), std::move(tmp), res, options);
            return;
        }

        int sa = a.size();
        int sb = b.size();
        if (&a != &res)
            move(std::move(a), res);
        else
            res.pm().alloc(res, sa); // free cache

        if (sb == 0) // monic
            return;

        res._poly.resize(sa + sb - (res.monic() || b.monic() ? 0 : 1));
        for (int i = 0; sa + i < res.size(); i++)
            res._poly[sa + i] = b._poly[i];

        if (sb == 1)
        {
            std::unique_ptr<GWNum> tmp(b.monic() ? new GWNum(gw()) : new GWNumWrapper(gw(), b._poly[0]));
            if (b.monic())
                gwfft(gw().gwdata(), b._poly[0], **tmp);
            if (res.monic() && b.monic())
                gwadd3o(gw().gwdata(), res._poly[sa - 1], b._poly[0], res._poly[sa], GWADD_FORCE_NORMALIZE);
            else if (b.monic())
                gwcopy(gw().gwdata(), res._poly[sa - 1], res._poly[sa]);
            for (int i = sa - 1; i >= 0; i--)
                if (b.monic() && i > 0)
                    gwmuladd4(gw().gwdata(), res._poly[i], **tmp, res._poly[i - 1], res._poly[i], (options & POLYMULT_STARTNEXTFFT ? GWMUL_STARTNEXTFFT : 0) | GWMUL_FFT_S2 | GWMUL_FFT_S3);
                else
                    gwmul3(gw().gwdata(), res._poly[i], **tmp, res._poly[i], (options & POLYMULT_STARTNEXTFFT ? GWMUL_STARTNEXTFFT : 0) | (i > 0 ? GWMUL_FFT_S2 : 0));
        }
        else
        {
            polymult(pmdata(), res._poly.data(), sa, b._poly.data(), sb, res.data(), res.size(), options | (res.monic() ? POLYMULT_INVEC1_MONIC : 0) | (b.monic() ? POLYMULT_INVEC2_MONIC : 0)/* | (pmdata()->num_threads > 1 ? POLYMULT_NO_UNFFT : 0)*/);
            /*if (pmdata()->num_threads > 1)
            {
                if ((options & POLYMULT_NEXTFFT))
                    poly_unfft_fft_coefficients(pmdata(), res.data(), res.size());
                else
                    poly_unfft_coefficients(pmdata(), res.data(), res.size());
            }*/
        }

        if (!res.monic() && !b.monic() && b._freeable)
            gwfree(gw().gwdata(), b._poly[sb - 1]);
        b._poly.clear();
        res._monic = res.monic() && b.monic();
        b._monic = false;
        b._freeable = true;
    }

    void PolyMult::mul_split(Poly& a, Poly& b, Poly& res_lo, Poly& res_hi, int size, int options)
    {
        if (a.degree() < b.degree() || a.size() < b.size())
        {
            mul_split(b, a, res_lo, res_hi, size, options);
            return;
        }
        int sa = a.size();
        int sb = b.size();
        int sr = sa + sb - (a.monic() || b.monic() ? 0 : 1);
        int slo = (sr <= size ? sr : size);
        int shi = (sr <= size ? 0 : sr - size);

        Poly tmp_lo(res_lo.pm());
        if (slo <= res_lo._poly.size())
            tmp_lo._poly.insert(tmp_lo._poly.begin(), res_lo._poly.begin(), res_lo._poly.begin() + slo);
        else
        {
            tmp_lo._poly = res_lo._poly;
            tmp_lo._freeable = res_lo._freeable;
            tmp_lo.pm().alloc(tmp_lo, slo);
        }
        Poly tmp_hi(res_hi.pm());
        if (shi <= res_hi._poly.size())
            tmp_hi._poly.insert(tmp_hi._poly.begin(), res_hi._poly.begin(), res_hi._poly.begin() + shi);
        else
        {
            tmp_hi._poly = res_hi._poly;
            tmp_hi._freeable = res_hi._freeable;
            tmp_hi.pm().alloc(tmp_hi, shi);
        }
        Poly tmp(*this);
        tmp._poly = tmp_lo._poly;
        tmp._poly.insert(tmp._poly.end(), tmp_hi._poly.begin(), tmp_hi._poly.end());

        GWASSERT(!(sb == 1 && !a.preprocessed() && !b.preprocessed()) || &a != &res_hi);
        if (sb == 1 && !a.preprocessed() && !b.preprocessed() && shi > 0 && b.data()[0] == tmp_hi._poly.data()[0])
        {
            Poly tmp_b(b);
            mul(a, tmp_b, tmp, options);
        }
        else
            mul(a, b, tmp, options);

        tmp._poly.clear();
        if (slo < res_lo._poly.size())
            res_lo.pm().alloc(res_lo, slo);
        if (shi < res_hi._poly.size())
            res_hi.pm().alloc(res_hi, shi);
        res_lo._poly = std::move(tmp_lo._poly);
        res_hi._poly = std::move(tmp_hi._poly);
        res_lo._monic = tmp._monic && (sr < size);
        res_hi._monic = tmp._monic && (sr >= size);
    }

    void PolyMult::mul_split(Poly&& a, Poly&& b, Poly& res_lo, Poly& res_hi, int size, int options)
    {
        int sa = a.size();
        int sb = b.size();
        mul(std::move(a), std::move(b), res_lo, options);
        res_hi.pm().free(res_hi);

        int sr = sa + sb - (a.monic() || b.monic() ? 0 : 1);
        if (sr >= size)
        {
            res_hi._poly.insert(res_hi._poly.begin(), res_lo._poly.begin() + size, res_lo._poly.end());
            res_lo._poly.erase(res_lo._poly.begin() + size, res_lo._poly.end());
            res_hi._monic = res_lo._monic;
            res_lo._monic = false;
            res_hi._freeable = res_lo._freeable;
        }
    }

    void PolyMult::preprocess(Poly& a, Poly& res, int size, int options)
    {
        if (a.size() > 1 && a._cache == nullptr)
        {
            res._cache = polymult_preprocess(pmdata(), a._poly.data(), a._poly.size(), size, size, options | POLYMULT_CIRCULAR | (a.monic() ? POLYMULT_INVEC1_MONIC : 0));
            res._cache_size = a._poly.size();
            if (res._freeable)
                for (auto it = res._poly.begin(); it != res._poly.end(); it++)
                    gwfree(gw().gwdata(), *it);
            res._poly.clear();
            res._monic = a.monic();
            res._freeable = true;
        }
        else
            copy(a, res);
    }

    void PolyMult::preprocess_and_mul(Poly& a, Poly& b, Poly& res, int size, int options)
    {
        GWASSERT(&a.pm().gw() == &res.pm().gw());
        GWASSERT(&b.pm().gw() == &res.pm().gw());
        GWASSERT(&a != &res);
        GWASSERT(&b != &res);
        GWASSERT(!a.preprocessed());
        GWASSERT(!b.preprocessed());
        GWASSERT(a._freeable == b._freeable);
        if (a.degree() < b.degree() || a.size() < b.size())
        {
            preprocess_and_mul(b, a, res, size, options);
            return;
        }
        if (b.degree() < 0)
        {
            res.pm().free(res);
            return;
        }

        int sa = a.size();
        int sb = b.size();
        if (sa > 1)
        {
            move(std::move(a), res);
            a._cache = polymult_preprocess(pmdata(), res._poly.data(), sa, size, size, POLYMULT_CIRCULAR | (options & POLYMULT_PRE_FFT ? POLYMULT_PRE_FFT : 0) | (res.monic() ? POLYMULT_INVEC1_MONIC : 0));
            a._cache_size = sa;
        }
        else
            copy(a, res);

        if (sb == 0) // monic
            return;
        if (sb == 1)
        {
            mul(res, b, res, options);
            return;
        }
        b._cache = polymult_preprocess(pmdata(), b._poly.data(), sb, size, size, POLYMULT_CIRCULAR | (options & POLYMULT_PRE_FFT ? POLYMULT_PRE_FFT : 0) | (b.monic() ? POLYMULT_INVEC1_MONIC : 0));
        b._cache_size = sb;

        res._poly.resize(sa + sb - (res.monic() || b.monic() ? 0 : 1));
        for (int i = 0; sa + i < res.size(); i++)
            res._poly[sa + i] = b._poly[i];

        polymult(pmdata(), a._cache, sa, b._cache, sb, res.data(), res.size(), options | (res.monic() ? POLYMULT_INVEC1_MONIC : 0) | (b.monic() ? POLYMULT_INVEC2_MONIC : 0)/* | (pmdata()->num_threads > 1 ? POLYMULT_NO_UNFFT : 0)*/);
        /*if (pmdata()->num_threads > 1)
        {
            if ((options & POLYMULT_NEXTFFT))
                poly_unfft_fft_coefficients(pmdata(), res.data(), res.size());
            else
                poly_unfft_coefficients(pmdata(), res.data(), res.size());
        }*/

        if (!res.monic() && !b.monic() && b._freeable)
            gwfree(gw().gwdata(), b._poly[sb - 1]);
        b._poly.clear();
        res._freeable = b._freeable;
        b._freeable = true;
        res._monic = res.monic() && b.monic();
        b._monic = false;
    }

    void PolyMult::poly_seize(Poly& a, Poly& res, Poly& to_free, int size)
    {
        res._poly.reserve(size);
        for (int i = 0; i < a._poly.size(); i++)
            if (a._poly[i] != nullptr)
            {
                if (res._poly.size() < size)
                    res._poly.push_back(a._poly[i]);
                else if (a._freeable)
                    to_free._poly.push_back(a._poly[i]);
            }
        res.pm().alloc(res, size);
    }

    void PolyMult::mul_range(Poly& a, Poly& b, Poly& res, int offset, int count, int options)
    {
        if (a.degree() < b.degree() || a.size() < b.size())
        {
            mul_range(b, a, res, offset, count, options);
            return;
        }
        if (b.degree() < 0 || a.degree() + b.degree() < offset)
        {
            res.pm().free(res);
            return;
        }
        Poly to_free(res.pm());
        if (b.size() == 0 && !a.preprocessed()) // monic
        {
            if (&a != &res)
            {
                res.pm().alloc(res, a.size() < offset + count ? a.size() - offset : count);
                for (int i = 0; i < res.size(); i++)
                    gwcopy(gw().gwdata(), a._poly[i + offset], res._poly[i]);
                res._monic = a.monic() && a.size() < offset + count;
            }
            else
            {
                res >>= offset;
                if (res.size() >= count)
                {
                    res._monic = false;
                    if (res._freeable)
                        to_free._poly.insert(to_free._poly.end(), res._poly.begin() + count, res._poly.end());
                    res._poly.resize(count);
                    GWASSERT(res.pm().gw().gwdata() == gw().gwdata());
                }
            }
            return;
        }

        int sa = a.size();
        int sb = b.size();
        int full = sa + sb - (a.monic() || b.monic() ? 0 : 1);
        int circular = offset + count;
        if (circular < full - offset)
            circular = full - offset;
        if (circular >= pmdata()->FFT_BREAK)
            circular = polymult_fft_size(circular);
        int size = full < offset + count ? full - offset : count;
        Poly tmp(res.pm());
        tmp._freeable = res._freeable;
        poly_seize(res, tmp, to_free, size);

        if (sb == 1 && !a.preprocessed() && !b.preprocessed())
        {
            GWASSERT(&a != &res || count == 1);
            if (sa >= offset && sa < offset + count)
            {
                if (a.monic() && b.monic())
                    gwadd3o(gw().gwdata(), a._poly[sa - 1], b._poly[0], tmp._poly[sa - offset], GWADD_FORCE_NORMALIZE);
                else if (a.monic())
                    gwcopy(gw().gwdata(), b._poly[0], tmp._poly[sa - offset]);
                else if (b.monic())
                    gwcopy(gw().gwdata(), a._poly[sa - 1], tmp._poly[sa - offset]);
            }
            for (int i = (sa < offset + count ? sa : offset + count) - 1; i >= offset; i--)
                if (b.monic() && i > 0)
                    gwmuladd4(gw().gwdata(), a._poly[i], b._poly[0], a._poly[i - 1], tmp._poly[i - offset], (options & POLYMULT_STARTNEXTFFT ? GWMUL_STARTNEXTFFT : 0) | GWMUL_FFT_S2 | GWMUL_FFT_S3);
                else
                    gwmul3(gw().gwdata(), a._poly[i], b._poly[0], tmp._poly[i - offset], (options & POLYMULT_STARTNEXTFFT ? GWMUL_STARTNEXTFFT : 0) | (i > 0 ? GWMUL_FFT_S2 : 0));
        }
        else
        {
            polymult2(pmdata(), a.data(), sa, b.data(), sb, tmp.data(), tmp.size(), nullptr, (full > circular ? circular : 0), offset, options | POLYMULT_MULMID | (full > circular ? POLYMULT_CIRCULAR : 0) | (a.monic() ? POLYMULT_INVEC1_MONIC : 0) | (b.monic() ? POLYMULT_INVEC2_MONIC : 0)/* | (pmdata()->num_threads > 1 ? POLYMULT_NO_UNFFT : 0)*/);
            /*if (pmdata()->num_threads > 1)
            {
                if ((options & POLYMULT_NEXTFFT))
                    poly_unfft_fft_coefficients(pmdata(), tmp.data(), tmp.size());
                else
                    poly_unfft_coefficients(pmdata(), tmp.data(), tmp.size());
            }*/
        }

        res._monic = a.monic() && b.monic() && full < offset + count;
        res._poly = std::move(tmp._poly);
    }

    void PolyMult::mul_twohalf(Poly& a, Poly& b, Poly& c, Poly& res1, Poly& res2, int half, int options)
    {
        if (b.degree() < 0 || a.degree() + b.degree() < half)
        {
            res1.pm().free(res1);
            mul_range(a, c, res2, half, half, options);
            return;
        }
        if (c.degree() < 0 || a.degree() + c.degree() < half)
        {
            res2.pm().free(res2);
            mul_range(a, b, res1, half, half, options);
            return;
        }
        if (b.size() == 0 || c.size() == 0) // monic
        {
            Poly tmp(a.pm(), a.size() - half, a.size() < 2*half && a._monic);
            for (int i = 0; i < tmp.size(); i++)
                gwcopy(gw().gwdata(), a._poly[i + half], tmp._poly[i]);
            if (b.size() == 0 && c.size() == 0)
            {
                copy(tmp, res1);
                copy(tmp, res2);
            }
            else if (b.size() == 0)
            {
                mul_range(a, c, res2, half, half, options);
                copy(tmp, res1);
            }
            else if (c.size() == 0)
            {
                mul_range(a, b, res1, half, half, options);
                copy(tmp, res2);
            }
            return;
        }

        int sa = a.size();
        int sb = b.size();
        int sc = c.size();
        int full1 = sa + sb - (a.monic() || b.monic() ? 0 : 1);
        int full2 = sa + sc - (a.monic() || c.monic() ? 0 : 1);
        int size1 = full1 < 2*half ? full1 - half : half;
        int size2 = full2 < 2*half ? full2 - half : half;

        if (sb == 1)
        {
            Poly tmp(res1.pm());
            mul_range(a, b, tmp, half, half, options);
            mul_range(a, c, res2, half, half, options);
            copy(tmp, res1);
            return;
        }
        if (sc == 1)
        {
            Poly tmp(res2.pm());
            mul_range(a, c, tmp, half, half, options);
            mul_range(a, b, res1, half, half, options);
            copy(tmp, res2);
            return;
        }

        res1.pm().alloc(res1, size1);
        res2.pm().alloc(res2, size2);

#ifdef NO_POLYMULT_SEVERAL
        a._cache = polymult_preprocess(pmdata(), a.data(), sa, 2*half, 2*half, POLYMULT_CIRCULAR | POLYMULT_PRE_FFT | (a.monic() ? POLYMULT_INVEC1_MONIC : 0));
        polymult2(pmdata(), a._cache, sa, b.data(), sb, res1.data(), res1.size(), nullptr, (full1 > 2*half ? 2*half : 0), 0, options | POLYMULT_MULHI | (full1 > 2*half ? POLYMULT_CIRCULAR : 0) | (a.monic() ? POLYMULT_INVEC1_MONIC : 0) | (b.monic() ? POLYMULT_INVEC2_MONIC : 0)/* | (pmdata()->num_threads > 1 ? POLYMULT_NO_UNFFT : 0)*/);
        polymult2(pmdata(), a._cache, sa, c.data(), sc, res2.data(), res2.size(), nullptr, (full2 > 2*half ? 2*half : 0), 0, options | POLYMULT_MULHI | (full2 > 2*half ? POLYMULT_CIRCULAR : 0) | (a.monic() ? POLYMULT_INVEC1_MONIC : 0) | (c.monic() ? POLYMULT_INVEC2_MONIC : 0)/* | (pmdata()->num_threads > 1 ? POLYMULT_NO_UNFFT : 0)*/);
        gwfree_array(gw().gwdata(), a._cache);
        a._cache = nullptr;
#endif
        polymult_arg args[2];
        polymult_arg& arg_b = args[0];
        polymult_arg& arg_c = args[1];
        arg_b.invec2 = b.data();
        arg_b.invec2_size = sb;
        arg_b.outvec = res1.data();
        arg_b.outvec_size = res1.size();
        arg_b.fmavec = nullptr;
        arg_b.circular_size = (full1 > 2*half ? 2*half : 0);
        arg_b.first_mulmid = 0;
        arg_b.options = (full1 > 2*half ? POLYMULT_CIRCULAR : 0) | (b.monic() ? POLYMULT_INVEC2_MONIC : 0);
        arg_c.invec2 = c.data();
        arg_c.invec2_size = sc;
        arg_c.outvec = res2.data();
        arg_c.outvec_size = res2.size();
        arg_c.fmavec = nullptr;
        arg_c.circular_size = (full2 > 2*half ? 2*half : 0);
        arg_c.first_mulmid = 0;
        arg_c.options = (full2 > 2*half ? POLYMULT_CIRCULAR : 0) | (c.monic() ? POLYMULT_INVEC2_MONIC : 0);
        polymult_several(pmdata(), a.data(), sa, args, 2, options | POLYMULT_MULHI | (a.monic() ? POLYMULT_INVEC1_MONIC : 0)/* | (pmdata()->num_threads > 1 ? POLYMULT_NO_UNFFT : 0)*/);
        /*if (pmdata()->num_threads > 1)
        {
            if ((options & POLYMULT_NEXTFFT))
                poly_unfft_fft_coefficients(pmdata(), res1.data(), res1.size());
            else
                poly_unfft_coefficients(pmdata(), res1.data(), res1.size());
            if ((options & POLYMULT_NEXTFFT))
                poly_unfft_fft_coefficients(pmdata(), res2.data(), res2.size());
            else
                poly_unfft_coefficients(pmdata(), res2.data(), res2.size());
        }*/

        res1._monic = a.monic() && b.monic() && full1 < 2*half;
        res2._monic = a.monic() && c.monic() && full2 < 2*half;
    }

    void PolyMult::mul_twohalf(Poly&& a, Poly& b, Poly& c, Poly& res1, Poly& res2, int half, int options)
    {
        GWASSERT(&a.pm().gw() == &res1.pm().gw());
        GWASSERT(&a.pm().gw() == &res2.pm().gw());
        //GWASSERT(b.preprocessed() || b.size() <= 1);
        //GWASSERT(c.preprocessed() || c.size() <= 1);
        if (b.degree() < 0 || a.degree() + b.degree() < half)
        {
            res1.pm().free(res1);
            res2 = std::move(a);
            mul_range(res2, c, res2, half, half, options);
            return;
        }
        if (c.degree() < 0 || a.degree() + c.degree() < half)
        {
            res2.pm().free(res2);
            res1 = std::move(a);
            mul_range(res1, b, res1, half, half, options);
            return;
        }
        if (b.size() == 0 && c.size() == 0) // monic
        {
            move(std::move(a), res1);
            res2.pm().free(res2);
            res2._poly.insert(res2._poly.begin(), res1._poly.begin() + half, res1._poly.end());
            res1._poly.resize(half);
            res2._monic = res1._monic;
            res2._freeable = res1._freeable;
            copy(res2, res1);
            return;
        }
        if (b.size() == 0) // monic
        {
            res1.pm().alloc(res1, a.size() - half);
            for (int i = 0; i < res1.size(); i++)
                gwcopy(gw().gwdata(), a._poly[i + half], res1._poly[i]);
            res1._monic = a.monic() && a.size() < 2*half;
            res2 = std::move(a);
            mul_range(res2, c, res2, half, half, options);
            return;
        }
        if (c.size() == 0) // monic
        {
            res2.pm().alloc(res2, a.size() - half);
            for (int i = 0; i < res2.size(); i++)
                gwcopy(gw().gwdata(), a._poly[i + half], res2._poly[i]);
            res2._monic = a.monic() && a.size() < 2*half;
            res1 = std::move(a);
            mul_range(res1, b, res1, half, half, options);
            return;
        }

        int sa = a.size();
        int sb = b.size();
        int sc = c.size();
        int full1 = sa + sb - (a.monic() || b.monic() ? 0 : 1);
        int full2 = sa + sc - (a.monic() || c.monic() ? 0 : 1);
        int size1 = full1 < 2*half ? full1 - half : half;
        int size2 = full2 < 2*half ? full2 - half : half;

        if (sb == 1)
        {
            mul_range(a, b, res1, half, half, options);
            res2 = std::move(a);
            mul_range(res2, c, res2, half, half, options);
            return;
        }
        if (sc == 1)
        {
            mul_range(a, c, res2, half, half, options);
            res1 = std::move(a);
            mul_range(res1, b, res1, half, half, options);
            return;
        }

        res1.pm().free(res1);
        res2.pm().free(res2);
        res1._freeable = a._freeable;
        res2._freeable = a._freeable;
        res1._poly = a._poly;
        if (sa > size1)
        {
            res2._poly.insert(res2._poly.begin(), res1._poly.begin() + size1, res1._poly.end());
            res1._poly.resize(size1);
        }
        if (sa > size1 + size2)
            res2._poly.resize(size2);
        res1.pm().alloc(res1, size1);
        res2.pm().alloc(res2, size2);

#ifdef NO_POLYMULT_SEVERAL
        a._cache = polymult_preprocess(pmdata(), a.data(), sa, 2*half, 2*half, POLYMULT_CIRCULAR | POLYMULT_PRE_FFT | (a.monic() ? POLYMULT_INVEC1_MONIC : 0));
        polymult2(pmdata(), a._cache, sa, b.data(), sb, res1.data(), res1.size(), nullptr, (full1 > 2*half ? 2*half : 0), 0, options | POLYMULT_MULHI | (full1 > 2*half ? POLYMULT_CIRCULAR : 0) | (a.monic() ? POLYMULT_INVEC1_MONIC : 0) | (b.monic() ? POLYMULT_INVEC2_MONIC : 0)/* | (pmdata()->num_threads > 1 ? POLYMULT_NO_UNFFT : 0)*/);
        polymult2(pmdata(), a._cache, sa, c.data(), sc, res2.data(), res2.size(), nullptr, (full2 > 2*half ? 2*half : 0), 0, options | POLYMULT_MULHI | (full2 > 2*half ? POLYMULT_CIRCULAR : 0) | (a.monic() ? POLYMULT_INVEC1_MONIC : 0) | (c.monic() ? POLYMULT_INVEC2_MONIC : 0)/* | (pmdata()->num_threads > 1 ? POLYMULT_NO_UNFFT : 0)*/);
        gwfree_array(gw().gwdata(), a._cache);
        a._cache = nullptr;
#else
        polymult_arg args[2];
        polymult_arg& arg_b = args[0];
        polymult_arg& arg_c = args[1];
        arg_b.invec2 = b.data();
        arg_b.invec2_size = sb;
        arg_b.outvec = res1.data();
        arg_b.outvec_size = res1.size();
        arg_b.fmavec = nullptr;
        arg_b.circular_size = (full1 > 2*half ? 2*half : 0);
        arg_b.first_mulmid = 0;
        arg_b.options = (full1 > 2*half ? POLYMULT_CIRCULAR : 0) | (b.monic() ? POLYMULT_INVEC2_MONIC : 0);
        arg_c.invec2 = c.data();
        arg_c.invec2_size = sc;
        arg_c.outvec = res2.data();
        arg_c.outvec_size = res2.size();
        arg_c.fmavec = nullptr;
        arg_c.circular_size = (full2 > 2*half ? 2*half : 0);
        arg_c.first_mulmid = 0;
        arg_c.options = (full2 > 2*half ? POLYMULT_CIRCULAR : 0) | (c.monic() ? POLYMULT_INVEC2_MONIC : 0);
        polymult_several(pmdata(), a.data(), sa, args, 2, options | POLYMULT_MULHI | (a.monic() ? POLYMULT_INVEC1_MONIC : 0)/* | (pmdata()->num_threads > 1 ? POLYMULT_NO_UNFFT : 0)*/);
        /*if (pmdata()->num_threads > 1)
        {
            if ((options & POLYMULT_NEXTFFT))
                poly_unfft_fft_coefficients(pmdata(), res1.data(), res1.size());
            else
                poly_unfft_coefficients(pmdata(), res1.data(), res1.size());
            if ((options & POLYMULT_NEXTFFT))
                poly_unfft_fft_coefficients(pmdata(), res2.data(), res2.size());
            else
                poly_unfft_coefficients(pmdata(), res2.data(), res2.size());
        }*/

#endif

        if (a._freeable)
            for (int i = size1 + size2; i < a.size(); i++)
                gwfree(gw().gwdata(), a._poly[i]);
        a._poly.clear();
        res1._monic = a.monic() && b.monic() && full1 < 2*half;
        res2._monic = a.monic() && c.monic() && full2 < 2*half;
        a._monic = false;
    }

    void PolyMult::fma_range(Poly& a, Poly& b, Poly& fma, Poly& res, int offset, int count, int options)
    {
        GWASSERT(!fma.preprocessed());
        if (a.degree() < b.degree() || a.size() < b.size())
        {
            fma_range(b, a, fma, res, offset, count, options);
            return;
        }
        Poly to_free(res.pm());
        if (b.degree() < 0 || a.degree() + b.degree() < offset)
        {
            if (fma.degree() < offset)
            {
                res.pm().free(res);
                return;
            }
            GWASSERT(0);
        }
        if (b.size() == 0 && !a.preprocessed()) // monic
        {
            Poly tmp(res.pm());
            tmp._freeable = res._freeable;
            tmp._monic = (a.monic() && fma.degree() < a.size()) ||
                (fma.monic() && a.degree() < fma.size() && ((options & POLYMULT_FMADD) || (options & POLYMULT_FNMADD)));
            int fma_full = (a.degree() > fma.degree() ? a.degree() : fma.degree()) + (tmp.monic() ? 0 : 1);
            poly_seize(res, tmp, to_free, fma_full < offset + count ? fma_full - offset : count);

            int i;
            for (i = 0; i < tmp.size() && i + offset < a.size() && i + offset < fma.size(); i++)
                if (options & POLYMULT_FMADD)
                    gwadd3o(gw().gwdata(), a._poly[i + offset], fma._poly[i + offset], tmp._poly[i], GWADD_FORCE_NORMALIZE);
                else if (options & POLYMULT_FMSUB)
                    gwsub3o(gw().gwdata(), a._poly[i + offset], fma._poly[i + offset], tmp._poly[i], GWADD_FORCE_NORMALIZE);
                else if (options & POLYMULT_FNMADD)
                    gwsub3o(gw().gwdata(), fma._poly[i + offset], a._poly[i + offset], tmp._poly[i], GWADD_FORCE_NORMALIZE);
            for (; i < tmp.size() && i + offset < a.size(); i++)
                if ((options & POLYMULT_FMADD) || (options & POLYMULT_FMSUB))
                    gwcopy(gw().gwdata(), a._poly[i + offset], tmp._poly[i]);
                else if (options & POLYMULT_FNMADD)
                {
                    dbltogw(gw().gwdata(), 0, tmp._poly[i]);
                    gwsub3o(gw().gwdata(), tmp._poly[i], a._poly[i + offset], tmp._poly[i], GWADD_GUARANTEED_OK);
                }
            for (; i < tmp.size() && i + offset < fma.size(); i++)
                if ((options & POLYMULT_FMADD) || (options & POLYMULT_FNMADD))
                    gwcopy(gw().gwdata(), fma._poly[i + offset], tmp._poly[i]);
                else if (options & POLYMULT_FMSUB)
                {
                    dbltogw(gw().gwdata(), 0, tmp._poly[i]);
                    gwsub3o(gw().gwdata(), tmp._poly[i], fma._poly[i + offset], tmp._poly[i], GWADD_GUARANTEED_OK);
                }
            if (!tmp.monic() && fma_full <= offset + count && a.size() < fma_full && fma.size() < fma_full)
                dbltogw(gw().gwdata(), 0, tmp._poly[a.size() - offset]);
            if (a.monic() && a.size() < offset + count && a.size() < fma_full)
                gwsmalladd(gw().gwdata(), (options & POLYMULT_FNMADD) ? -1 : 1, tmp._poly[a.size() - offset]);
            if (fma.monic() && fma.size() < offset + count && fma.size() < fma_full)
                gwsmalladd(gw().gwdata(), (options & POLYMULT_FMSUB) ? -1 : 1, tmp._poly[fma.size() - offset]);

            res._monic = tmp.monic() && fma_full < offset + count;
            res._poly = std::move(tmp._poly);
            return;
        }

        int sa = a.size();
        int sb = b.size();
        int full = sa + sb - (a.monic() || b.monic() ? 0 : 1);
        int circular = offset + count;
        if (circular < full - offset)
            circular = full - offset;
        if (circular >= pmdata()->FFT_BREAK)
            circular = polymult_fft_size(circular);
        int size = full < offset + count ? full - offset : count;
        Poly tmp(res.pm());
        tmp._freeable = res._freeable;
        int degree = full - (a.monic() && b.monic() ? 0 : 1);
        tmp._monic = (a.monic() && b.monic() && fma.degree() < full) ||
            (fma.monic() && degree < fma.size() && ((options & POLYMULT_FMADD) || (options & POLYMULT_FNMADD)));
        int fma_full = (degree > fma.degree() ? degree : fma.degree()) + (tmp.monic() ? 0 : 1);
        poly_seize(res, tmp, to_free, fma_full < offset + count ? fma_full - offset : count);

        int padding = offset + count - fma.size();
        if (padding > 0)
            fma._poly.insert(fma._poly.end(), padding, nullptr);
        polymult2(pmdata(), a.data(), sa, b.data(), sb, tmp.data(), size, fma.data() + offset, (full > circular ? circular : 0), offset, options | POLYMULT_MULMID | (full > circular ? POLYMULT_CIRCULAR : 0) | (a.monic() ? POLYMULT_INVEC1_MONIC : 0) | (b.monic() ? POLYMULT_INVEC2_MONIC : 0)/* | (pmdata()->num_threads > 1 ? POLYMULT_NO_UNFFT : 0)*/);
        /*if (pmdata()->num_threads > 1)
        {
            if ((options & POLYMULT_NEXTFFT))
                poly_unfft_fft_coefficients(pmdata(), tmp.data(), size);
            else
                poly_unfft_coefficients(pmdata(), tmp.data(), size);
        }*/
        if (padding > 0)
            fma._poly.erase(fma._poly.end() - padding, fma._poly.end());

        for (int i = size; i < tmp.size() && i + offset < fma.size(); i++)
            if ((options & POLYMULT_FMADD) || (options & POLYMULT_FNMADD))
                gwcopy(gw().gwdata(), fma._poly[i + offset], tmp._poly[i]);
            else if (options & POLYMULT_FMSUB)
            {
                dbltogw(gw().gwdata(), 0, tmp._poly[i]);
                gwsub3o(gw().gwdata(), tmp._poly[i], fma._poly[i + offset], tmp._poly[i], GWADD_GUARANTEED_OK);
            }
        if (!tmp.monic() && fma_full <= offset + count && full < fma_full && fma.size() < fma_full)
            dbltogw(gw().gwdata(), 0, tmp._poly[full - offset]);
        if (a.monic() && b.monic() && full < offset + count && full < fma_full)
            gwsmalladd(gw().gwdata(), (options & POLYMULT_FNMADD) ? -1 : 1, tmp._poly[full - offset]);
        if (fma.monic() && fma.size() < offset + count && fma.size() < fma_full)
            gwsmalladd(gw().gwdata(), (options & POLYMULT_FMSUB) ? -1 : 1, tmp._poly[fma.size() - offset]);

        res._monic = tmp.monic() && fma_full < offset + count;
        res._poly = std::move(tmp._poly);
    }

    void PolyMult::reciprocal(Poly& a, Poly& res, int options)
    {
        GWASSERT(&a.pm().gw() == &res.pm().gw());
        int i, j;
        int d;
        for (d = 2; d < res.size() + (a.monic() ? 1 : 0); d <<= 1);

        gwnum* g = res.data() + res.size() - 1;
        std::vector<gwnum> f(d);
        gwarray tmp = gwalloc_array(gw().gwdata(), d/2);

        int sa = a.size();
        if (a.monic())
        {
            dbltogw(gw().gwdata(), 0, *g);
            gwsub3o(gw().gwdata(), *g, a._poly[sa - 1], *g, GWADD_GUARANTEED_OK);
            for (i = 2; i < d; g -= i, i <<= 1)
            {
                for (j = 0; j < 2*i - 1; j++)
                    f[j] = (sa - 2*i + 1 + j) >= 0 ? a._poly[sa - 2*i + 1 + j] : nullptr;
                polymult2(pmdata(), g, i - 1, f.data(), 2*i - 1, tmp, i, nullptr, 2*i, i - 1, POLYMULT_CIRCULAR | POLYMULT_MULMID | POLYMULT_INVEC1_MONIC | POLYMULT_NEXTFFT/* | (pmdata()->num_threads > 1 && pmdata()->num_threads*5 < i ? POLYMULT_NO_UNFFT : 0)*/);
                /*if (pmdata()->num_threads > 1 && pmdata()->num_threads*5 < i)
                    poly_unfft_fft_coefficients(pmdata(), tmp, i);*/
                j = res.size() - 2*i + 1;
                polymult2(pmdata(), g, i - 1, tmp, i, res.data() + (j >= 0 ? j : 0), i + (j < 0 ? j : 0), nullptr, 0, i - 1 - (j < 0 ? j : 0), POLYMULT_MULMID | POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_NEGATE | (2*i < d ? POLYMULT_NEXTFFT : options)/* | (pmdata()->num_threads > 1 && pmdata()->num_threads*5 < i ? POLYMULT_NO_UNFFT : 0)*/);
                /*if (pmdata()->num_threads > 1 && pmdata()->num_threads*5 < i)
                {
                    if (2*i < d || (options & POLYMULT_NEXTFFT))
                        poly_unfft_fft_coefficients(pmdata(), res.data() + (j >= 0 ? j : 0), i + (j < 0 ? j : 0));
                    else
                        poly_unfft_coefficients(pmdata(), res.data() + (j >= 0 ? j : 0), i + (j < 0 ? j : 0));
                }*/
            }
        }
        else
        {
            gw().inv((GWNum&)gw().wrap(a._poly[sa - 1]), (GWNum&)gw().wrap(*g));
            for (i = 1; i < d; g -= i, i <<= 1)
            {
                for (j = 0; j < 2*i - 1; j++)
                    f[j] = (sa - 2*i + j) >= 0 ? a._poly[sa - 2*i + j] : nullptr;
                polymult2(pmdata(), g, i, f.data(), 2*i - 1, tmp, i, nullptr, 2*i, i - 1, POLYMULT_CIRCULAR | POLYMULT_MULMID | POLYMULT_NEXTFFT/* | (pmdata()->num_threads > 1 ? POLYMULT_NO_UNFFT : 0)*/);
                /*if (pmdata()->num_threads > 1)
                    poly_unfft_fft_coefficients(pmdata(), tmp, i);*/
                j = res.size() - 2*i;
                polymult2(pmdata(), g, i, tmp, i, res.data() + (j >= 0 ? j : 0), i + (j < 0 ? j : 0), nullptr, 0, i - 1 - (j < 0 ? j : 0), POLYMULT_MULMID | POLYMULT_INVEC2_NEGATE | (2*i < d ? POLYMULT_NEXTFFT : options)/* | (pmdata()->num_threads > 1 ? POLYMULT_NO_UNFFT : 0)*/);
                /*if (pmdata()->num_threads > 1)
                {
                    if (2*i < d || (options & POLYMULT_NEXTFFT))
                        poly_unfft_fft_coefficients(pmdata(), res.data() + (j >= 0 ? j : 0), i + (j < 0 ? j : 0));
                    else
                        poly_unfft_coefficients(pmdata(), res.data() + (j >= 0 ? j : 0), i + (j < 0 ? j : 0));
                }*/
            }
        }

        gwfree_array(gw().gwdata(), tmp);
        res._monic = a.monic();
    }

    void PolyMult::shiftleft(Poly& a, int b, Poly& res)
    {
        GWASSERT(&a.pm().gw() == &res.pm().gw());
        if (&a != &res)
            copy(a, res);
        res._poly.insert(res._poly.begin(), b, nullptr);
    }

    void PolyMult::shiftright(Poly& a, int b, Poly& res)
    {
        GWASSERT(&a.pm().gw() == &res.pm().gw());
        if (b > a.size())
        {
            res.pm().free(res);
            return;
        }
        if (&a == &res)
        {
            if (res._freeable)
                for (int i = 0; i < b && i < res.size(); i++)
                    gwfree(gw().gwdata(), res._poly[i]);
            res._poly.erase(res._poly.begin(), res._poly.begin() + b);
        }
        else
        {
            res.pm().alloc(res, (int)a.size() - b);
            for (int i = 0; i < res.size(); i++)
                gwcopy(gw().gwdata(), a._poly[b + i], res._poly[i]);
            res._monic = a._monic;
        }
    }

    void PolyMult::convert(const Poly& a, PolyMult& pm_res, Poly& res)
    {
        res._monic = a.monic();
        res.pm().alloc(res, a.size());
        for (int i = 0; i < a.size(); i++)
            if (gwconvert(gw().gwdata(), pm_res.gw().gwdata(), a._poly[i], res._poly[i]) != 0)
                throw InvalidFFTDataException();
    }

    void PolyMult::insert(GWNum&& a, Poly& res, size_t pos)
    {
        GWASSERT(res._freeable);
        res._poly.insert(res._poly.begin() + pos, *a);
        a._gwnum = nullptr;
    }

    GWNum PolyMult::remove(Poly& a, size_t pos)
    {
        GWASSERT(a._freeable);
        gwnum res = a._poly.at(pos);
        a._poly.erase(a._poly.begin() + pos);
        return GWNum(gw(), res);
    }
}
