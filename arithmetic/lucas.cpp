
#include <stdlib.h>
#include "gwnum.h"
#include "lucas.h"

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
        res._parity = (a.parity() && !b.parity()) || (!a.parity() && b.parity());
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
}