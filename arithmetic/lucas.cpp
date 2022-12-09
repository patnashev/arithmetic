
#include <stdlib.h>
#include "gwnum.h"
#include "lucas.h"

namespace arithmetic
{
    void LucasVArithmetic::copy(const LucasV& a, LucasV& res)
    {
        res.V() = a.V();
    }

    void LucasVArithmetic::move(LucasV&& a, LucasV& res)
    {
        res._V = std::move(a._V);
    }

    void LucasVArithmetic::init(LucasV& res)
    {
        res.V() = 2;
    }

    void LucasVArithmetic::init(const GWNum& a, LucasV& res)
    {
        res.V() = a;
    }

    void LucasVArithmetic::add(LucasV& a, LucasV& b, LucasV& a_minus_b, LucasV& res)
    {
        gw().mulsub(a.V(), b.V(), a_minus_b.V(), res.V(), GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_FFT_S3 | GWMUL_STARTNEXTFFT);
    }

    void LucasVArithmetic::dbl(LucasV& a, LucasV& res)
    {
        gw().setaddin(-2);
        gw().mul(a.V(), a.V(), res.V(), GWMUL_FFT_S1 | GWMUL_ADDINCONST | GWMUL_STARTNEXTFFT);
        gw().setaddin(0);
    }

    void LucasVArithmetic::optimize(LucasV& a)
    {
        gwfft_for_fma(gw().gwdata(), *a.V(), *a.V());
    }
}