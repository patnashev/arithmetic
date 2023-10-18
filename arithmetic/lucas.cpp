
#include <stdlib.h>
#include "gwnum.h"
#include "lucas.h"
#include "exception.h"

namespace arithmetic
{
    void LucasVArithmetic::copy(const LucasV& a, LucasV& res)
    {
        res.V() = a.V();
        res._parity = a._parity;
    }

    void LucasVArithmetic::move(LucasV&& a, LucasV& res)
    {
        res._V = std::move(a._V);
        res._parity = a._parity;
    }

    void LucasVArithmetic::init(LucasV& res)
    {
        res.V() = 2;
        res._parity = false;
    }

    void LucasVArithmetic::init(const GWNum& a, LucasV& res)
    {
        res.V() = a;
        res._parity = true;
    }

    void LucasVArithmetic::init(const GWNum& a, bool parity, LucasV& res)
    {
        res.V() = a;
        res._parity = parity;
    }

    void LucasVArithmetic::add(LucasV& a, LucasV& b, LucasV& a_minus_b, LucasV& res)
    {
        if (negativeQ() && b.parity())
            gw().muladd(a.V(), b.V(), a_minus_b.V(), res.V(), GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_FFT_S3 | GWMUL_STARTNEXTFFT);
        else
            gw().mulsub(a.V(), b.V(), a_minus_b.V(), res.V(), GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_FFT_S3 | GWMUL_STARTNEXTFFT);
        res._parity = (a.parity() != b.parity());
    }

    void LucasVArithmetic::dbl(LucasV& a, LucasV& res)
    {
        gw().setaddin(negativeQ() && a.parity() ? 2 : -2);
        gw().mul(a.V(), a.V(), res.V(), GWMUL_FFT_S1 | GWMUL_ADDINCONST | GWMUL_STARTNEXTFFT);
        res._parity = false;
    }

    void LucasVArithmetic::optimize(LucasV& a)
    {
        gwfft_for_fma(gw().gwdata(), *a.V(), *a.V());
    }

    void LucasUVArithmetic::copy(const LucasUV& a, LucasUV& res)
    {
        res.U() = a.U();
        res.V() = a.V();
        if (a._DU)
        {
            if (!res._DU)
                res._DU.reset(new GWNum(gw()));
            *res._DU = *a._DU;
        }
        else
            res._DU.reset();
        res._parity = a._parity;
    }

    void LucasUVArithmetic::move(LucasUV&& a, LucasUV& res)
    {
        res._U = std::move(a._U);
        res._V = std::move(a._V);
        res._DU = std::move(a._DU);
        res._parity = a._parity;
    }

    void LucasUVArithmetic::init(LucasUV& res)
    {
        res.U() = 0;
        res.V() = 1;
        res._DU.reset();
        res._parity = false;
    }

    void LucasUVArithmetic::init(GWNum& P, LucasUV& res)
    {
        if (!_half)
        {
            _half.reset(new GWNum(gw()));
            *_half = ((gw().N() + 1) >> 1);
        }
        res.U() = *_half;
        gw().mul(P, *_half, res.V(), GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
        res._DU.reset();
        res._parity = true;
    }

    void LucasUVArithmetic::init(GWNum& U, GWNum& V, LucasUV& res)
    {
        init(U, V, true, false, res);
    }

    void LucasUVArithmetic::init(GWNum& U, GWNum& V, bool parity, bool halved, LucasUV& res)
    {
        if (halved)
        {
            res.U() = U;
            res.V() = V;
        }
        else
        {
            if (!_half)
            {
                _half.reset(new GWNum(gw()));
                *_half = ((gw().N() + 1) >> 1);
            }
            gw().mul(U, *_half, res.U(), GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
            gw().mul(V, *_half, res.V(), GWMUL_STARTNEXTFFT);
        }
        res._DU.reset();
        res._parity = parity;
    }

    void LucasUVArithmetic::init_small(int index, LucasUV& res)
    {
        if (abs(index) >= _UV_small.size())
            throw ArithmeticException();

        int U = _UV_small[abs(index)].U;
        int V = _UV_small[abs(index)].V;
        int DU = _UV_small[abs(index)].DU;
        if (index < 0)
        {
            if (negativeQ() && (index & 1))
                V = -V;
            else
            {
                U = -U;
                DU = -DU;
            }
        }
        if (U & 1)
        {
            res.U() = U > 0 ? ((gw().N() + 1) >> 1) : ((gw().N() - 1) >> 1);
            gwsmallmul(gw().gwdata(), abs(U), *res.U());
        }
        else
            res.U() = U/2;
        if (V & 1)
        {
            res.V() = V > 0 ? ((gw().N() + 1) >> 1) : ((gw().N() - 1) >> 1);
            gwsmallmul(gw().gwdata(), abs(V), *res.V());
        }
        else
            res.V() = V/2;
        if (!res._DU)
            res._DU.reset(new GWNum(gw()));
        if (DU & 1)
        {
            *res._DU = DU > 0 ? ((gw().N() + 1) >> 1) : ((gw().N() - 1) >> 1);
            gwsmallmul(gw().gwdata(), abs(DU), **res._DU);
        }
        else
            *res._DU = DU/2;
        res._parity = (index & 1);
    }

    void LucasUVArithmetic::add(LucasUV& a, LucasUV& b, LucasUV& res)
    {
        add(a, b, res, GWMUL_STARTNEXTFFT);
    }

    void LucasUVArithmetic::sub(LucasUV& a, LucasUV& b, LucasUV& res)
    {
        add(a, b, res, GWMUL_STARTNEXTFFT | LUCASADD_NEGATIVE);
    }

    void LucasUVArithmetic::neg(LucasUV& a, LucasUV& res)
    {
        if (&res != &a)
            copy(a, res);
        if (negativeQ() && a.parity())
            gw().neg(res.V(), res.V());
        else
        {
            gw().neg(res.U(), res.U());
            if (res._DU)
                gw().neg(*res._DU, *res._DU);
        }
    }

    void LucasUVArithmetic::dbl(LucasUV& a, LucasUV& res)
    {
        dbl(a, res, GWMUL_STARTNEXTFFT);
    }

    void LucasUVArithmetic::add(LucasUV& a, LucasUV& b, LucasUV& res, int options)
    {
        GWNum* aU = a._U.get();
        if (aU == res._U.get())
        {
            if (!_tmp)
                _tmp.reset(new GWNum(gw()));
            std::swap(_tmp, res._U);
        }
        optimize(b);
        int optionsU = (options & LUCASADD_OPTIMIZE) ? (!_D ? (options & ~GWMUL_STARTNEXTFFT) : (options | GWMUL_STARTNEXTFFT)) : options;

        if (options & LUCASADD_NEGATIVE)
        {
            if (negativeQ() && b.parity())
            {
                gw().mulmulsub(a.V(), b.U(), *aU, b.V(), res.U(), GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_FFT_S3 | GWMUL_FFT_S4 | optionsU);
                if (options & LUCASADD_OPTIMIZE)
                    optimize(res);
                gw().mulmulsub(*aU, *b._DU, a.V(), b.V(), res.V(), options);
            }
            else
            {
                gw().mulmulsub(*aU, b.V(), a.V(), b.U(), res.U(), GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_FFT_S3 | GWMUL_FFT_S4 | optionsU);
                if (options & LUCASADD_OPTIMIZE)
                    optimize(res);
                gw().mulmulsub(a.V(), b.V(), *aU, *b._DU, res.V(), options);
            }
        }
        else
        {
            gw().mulmuladd(*aU, b.V(), a.V(), b.U(), res.U(), GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_FFT_S3 | GWMUL_FFT_S4 | optionsU);
            if (options & LUCASADD_OPTIMIZE)
                optimize(res);
            gw().mulmuladd(a.V(), b.V(), *aU, *b._DU, res.V(), options);
        }

        if (!(options & LUCASADD_OPTIMIZE))
            res._DU.reset();
        res._parity = (a.parity() != b.parity());
    }

    void LucasUVArithmetic::dbl(LucasUV& a, LucasUV& res, int options)
    {
        int optionsU = (options & LUCASADD_OPTIMIZE) ? (!_D ? (options & ~GWMUL_STARTNEXTFFT) : (options | GWMUL_STARTNEXTFFT)) : options;
        gw().setmulbyconst(2);
        gw().mul(a.V(), a.U(), res.U(), GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_MULBYCONST | optionsU);
        if (options & LUCASADD_OPTIMIZE)
            optimize(res);
        else
            res._DU.reset();

        gw().setpostaddin(negativeQ() && a.parity() ? 1 : -1);
        gw().mul(a.V(), a.V(), res.V(), GWMUL_MULBYCONST | GWMUL_ADDINCONST | options);
        //gw().mul(a.V(), a.V(), res.V(), GWMUL_MULBYCONST);
        //gwsmalladd(gw().gwdata(), negativeQ() && a.parity() ? 1 : -1, *res.V());
        /*if (!_half_fma)
        {
            if (!_half)
            {
                _half.reset(new GWNum(gw()));
                *_half = ((gw().N() + 1) >> 1);
            }
            _half_fma.reset(new GWNum(gw()));
            gwfft_for_fma(gw().gwdata(), **_half, **_half_fma);
        }
        gw().mulsub(a.V(), a.V(), *_half_fma, res.V(), GWMUL_MULBYCONST | options);*/
     
        res._parity = false;
    }

    void LucasUVArithmetic::dbl_add_small(LucasUV& a, int index, LucasUV& res, int options)
    {
        if (abs(index) >= _UV_small.size())
            throw ArithmeticException();

        gw().mul(a.V(), a.U(), res.U(), GWMUL_FFT_S1 | GWMUL_FFT_S2);
        if (!_half_fma)
        {
            if (!_half)
            {
                _half.reset(new GWNum(gw()));
                *_half = ((gw().N() + 1) >> 1);
            }
            _half_fma.reset(new GWNum(gw()));
            gwfft_for_fma(gw().gwdata(), **_half, **_half_fma);
        }
        if (negativeQ() && a.parity())
            gw().muladd(a.V(), a.V(), dynamic_cast<CarefulGWArithmetic*>(_gw) != nullptr ? *_half : *_half_fma, res.V(), 0);
        else
            gw().mulsub(a.V(), a.V(), dynamic_cast<CarefulGWArithmetic*>(_gw) != nullptr ? *_half : *_half_fma, res.V(), 0);

        if (!_tmp)
            _tmp.reset(new GWNum(gw()));
        *_tmp = res.U();
        if (!_tmp2)
            _tmp2.reset(new GWNum(gw()));
        *_tmp2 = res.V();
        gwsmallmul(gw().gwdata(), _UV_small[abs(index)].V, *res.U());
        gwsmallmul(gw().gwdata(), _UV_small[abs(index)].U, **_tmp2);
        gwsmallmul(gw().gwdata(), _UV_small[abs(index)].V, *res.V());
        gwsmallmul(gw().gwdata(), _UV_small[abs(index)].DU, **_tmp);
        if (index < 0)
        {
            if (negativeQ() && (index & 1))
            {
                gw().sub(*_tmp2, res.U(), res.U(), GWADD_FORCE_NORMALIZE);
                if (options & LUCASADD_OPTIMIZE)
                    optimize(res);
                gw().sub(*_tmp, res.V(), res.V(), GWADD_FORCE_NORMALIZE);
            }
            else
            {
                gw().sub(res.U(), *_tmp2, res.U(), GWADD_FORCE_NORMALIZE);
                if (options & LUCASADD_OPTIMIZE)
                    optimize(res);
                gw().sub(res.V(), *_tmp, res.V(), GWADD_FORCE_NORMALIZE);
            }

        }
        else
        {
            gw().add(res.U(), *_tmp2, res.U(), GWADD_FORCE_NORMALIZE);
            if (options & LUCASADD_OPTIMIZE)
                optimize(res);
            gw().add(res.V(), *_tmp, res.V(), GWADD_FORCE_NORMALIZE);
        }

        if (!(options & LUCASADD_OPTIMIZE))
            res._DU.reset();
        res._parity = (index & 1);
    }

    void LucasUVArithmetic::optimize(LucasUV& a)
    {
        if (!a._DU)
        {
            a._DU.reset(new GWNum(gw()));
            if (!_D)
            {
                gw().unfft(a.U(), *a._DU);
                gwsmallmul(gw().gwdata(), _UV_small[1].DU, **a._DU);
            }
            else
                gw().mul(*_D, a.U(), *a._DU, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
        }
    }

    void LucasUVArithmetic::mul(LucasUV& a, Giant& b, LucasUV& res)
    {
        int len = b.bitlen();
        int W;
        for (W = 2; W < 16 && /*3 + 3*(1 << (W - 1)) <= maxSize &&*/ (3 << (W - 2)) + len/0.69*(2 + 2/(W + 1.0)) > (3 << (W - 1)) + len/0.69*(2 + 2/(W + 2.0)); W++);
        std::vector<int16_t> naf_w;
        get_NAF_W(W, b, naf_w);
        mul(a, W, naf_w, res);
    }

    void LucasUVArithmetic::mul(LucasUV& a, int W, std::vector<int16_t>& naf_w, LucasUV& res)
    {
        int i, j;

        std::vector<std::unique_ptr<LucasUV>> u;
        for (i = 0; i < (1 << (W - 2)); i++)
            u.emplace_back(new LucasUV(*this));

        // Dictionary
        copy(a, *u[0]);
        optimize(*u[0]);
        if (W > 2)
        {
            dbl(a, res, GWMUL_STARTNEXTFFT);
            for (i = 1; i < (1 << (W - 2)); i++)
                add(res, *u[i - 1], *u[i], LUCASADD_OPTIMIZE | GWMUL_STARTNEXTFFT);
        }

        // Signed window
        copy(*u[naf_w.back()/2], res);
        for (i = (int)naf_w.size() - 2; i >= 0; i--)
        {
            if (naf_w[i] != 0)
            {
                for (j = 1; j < W; j++)
                    dbl(res, res, GWMUL_STARTNEXTFFT);
                dbl(res, res, GWMUL_STARTNEXTFFT);
                add(res, *u[abs(naf_w[i])/2], res, (i > 0 ? GWMUL_STARTNEXTFFT : 0) | (naf_w[i] < 0 ? LUCASADD_NEGATIVE : 0));
            }
            else
                dbl(res, res, (i > 0 ? GWMUL_STARTNEXTFFT : 0));
        }
        _tmp.reset();
    }
}
