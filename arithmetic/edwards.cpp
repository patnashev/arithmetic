
#include <stdlib.h>
#include <immintrin.h>
#include "cpuid.h"
#include "gwnum.h"
#include "edwards.h"
#include "exception.h"

namespace arithmetic
{
    void EdwardsArithmetic::copy(const EdPoint& a, EdPoint& res)
    {
        if (a.X)
        {
            if (!res.X)
                res.X.reset(new GWNum(gw()));
            *res.X = *a.X;
        }
        else
            res.X.reset();
        if (a.Y)
        {
            if (!res.Y)
                res.Y.reset(new GWNum(gw()));
            *res.Y = *a.Y;
        }
        else
            res.Y.reset();
        if (a.Z)
        {
            if (!res.Z)
                res.Z.reset(new GWNum(gw()));
            *res.Z = *a.Z;
        }
        else
            res.Z.reset();
        if (a.T)
        {
            if (!res.T)
                res.T.reset(new GWNum(gw()));
            *res.T = *a.T;
        }
        else
            res.T.reset();
    }

    void EdwardsArithmetic::move(EdPoint&& a, EdPoint& res)
    {
        res.X = std::move(a.X);
        res.Y = std::move(a.Y);
        res.Z = std::move(a.Z);
        res.T = std::move(a.T);
    }

    void EdwardsArithmetic::init(EdPoint& res)
    {
        res.X.reset();
        res.Y.reset();
        res.Z.reset();
        res.T.reset();
    }

    void EdwardsArithmetic::init(const GWNum& X, const GWNum& Y, EdPoint& res)
    {
        res.X.reset(new GWNum(gw()));
        *res.X = X;
        res.Y.reset(new GWNum(gw()));
        *res.Y = Y;
    }

    void EdwardsArithmetic::init(const GWNum& X, const GWNum& Y, const GWNum& Z, const GWNum& T, EdPoint& res)
    {
        res.X.reset(new GWNum(gw()));
        *res.X = X;
        res.Y.reset(new GWNum(gw()));
        *res.Y = Y;
        res.Z.reset(new GWNum(gw()));
        *res.Z = Z;
        res.T.reset(new GWNum(gw()));
        *res.T = T;
    }

    int EdwardsArithmetic::cmp(const EdPoint& a, const EdPoint& b)
    {
        int ret;
        GWNum tmp(gw());
        GWNum tmp2(gw());

        if (a.X)
            tmp = *a.X;
        else
            tmp = 1;
        if (b.Z)
            tmp *= *b.Z;
        if (b.X)
            tmp2 = *b.X;
        else
            tmp2 = 1;
        if (a.Z)
            tmp2 *= *a.Z;
        ret = gw().cmp(tmp, tmp2);
        if (ret != 0)
            return ret;
        if (a.Y)
            tmp = *a.Y;
        else
            tmp = 1;
        if (b.Z)
            tmp *= *b.Z;
        if (b.Y)
            tmp2 = *b.Y;
        else
            tmp2 = 1;
        if (a.Z)
            tmp2 *= *a.Z;
        return gw().cmp(tmp, tmp2);
    }

    void EdwardsArithmetic::add(EdPoint& a, EdPoint& b, EdPoint& res)
    {
        add(a, b, res, GWMUL_STARTNEXTFFT);
    }

    void EdwardsArithmetic::sub(EdPoint& a, EdPoint& b, EdPoint& res)
    {
        add(a, b, res, GWMUL_STARTNEXTFFT | EDADD_NEGATIVE);
    }

    void EdwardsArithmetic::neg(EdPoint& a, EdPoint& res)
    {
        if (&res != &a)
            copy(a, res);
        gw().neg(*res.X, *res.X);
        if (res.T)
            gw().neg(*res.T, *res.T);
    }

    void EdwardsArithmetic::dbl(EdPoint& a, EdPoint& res)
    {
        dbl(a, res, GWMUL_STARTNEXTFFT);
    }

    // https://eprint.iacr.org/2021/1061
    void EdwardsArithmetic::add(EdPoint& a, EdPoint& b, EdPoint& res, int options)
    {
        GWASSERT(&b != &res);
        bool safe11 = mul_safe(gw().gwdata(), 1, 1);

        if (!_tmp && (!(options & ED_PROJECTIVE) || &a == &res))
            _tmp.reset(new GWNum(gw()));
        if (!a.T)
        {
            if (a.Z)
                throw ArithmeticException("Projective coordinates are not suitable for add, call extend().");
            a.T.reset(new GWNum(gw()));
            gw().mul(*a.X, *a.Y, *a.T, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT_IF(b.Z || safe11 || (options & ED_PROJECTIVE)));
        }
        if (!b.T)
        {
            if (b.Z)
                throw ArithmeticException("Projective coordinates are not suitable for add, call extend().");
            b.T.reset(new GWNum(gw()));
            gw().mul(*b.X, *b.Y, *b.T, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT_IF(a.Z || safe11 || (options & ED_PROJECTIVE)));
        }
        if (!res.X)
            res.X.reset(new GWNum(gw()));
        if (!res.Y)
            res.Y.reset(new GWNum(gw()));
        if (!res.T)
            res.T.reset(new GWNum(gw()));

        if (a.Z)
        {
            if (!res.Z)
                res.Z.reset(new GWNum(gw()));
            if (options & EDADD_NEGATIVE)
                gw().setmulbyconst(-1);
            gw().mul(*b.T, *a.Z, *res.Z, GWMUL_FFT_S1 | ((options & EDADD_NEGATIVE) ? GWMUL_MULBYCONST : 0) | GWMUL_STARTNEXTFFT_IF(safe11 || (options & ED_PROJECTIVE))); // C = Z1 * T2
        }
        else
        {
            if (!res.Z)
                res.Z.reset(new GWNum(gw()));
            if ((safe11 || (options & ED_PROJECTIVE)) && !(options & EDADD_NEGATIVE))
                gw().copy(*b.T, *res.Z);
            else
                gw().unfft(*b.T, *res.Z);
            if (options & EDADD_NEGATIVE)
                gw().neg(*res.Z, *res.Z);
        }
        if (b.Z)
            gw().mul(*b.Z, *a.T, *res.T, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT_IF(safe11 || (options & ED_PROJECTIVE))); // D = T1 * Z2
        else if (&a != &res)
            gw().copy(*a.T, *res.T);
        gw().addsub(*res.T, *res.Z, *res.T, *res.Z, GWADD_DELAYNORM_IF(safe11 || (options & ED_PROJECTIVE))); // E = D + C, H = D - C
        if (options & EDADD_NEGATIVE)
        {
            gw().mulmuladd(*a.X, *b.Y, *a.Y, *b.X, (&a == &res) ? *_tmp : *res.X, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_FFT_S3 | GWMUL_FFT_S4 | GWMUL_STARTNEXTFFT); // F = X1*Y2 - Y1*[-]X2
            gw().mulmulsub(*a.Y, *b.Y, *a.X, *b.X, *res.Y, GWMUL_STARTNEXTFFT); // G = Y1*Y2 + X1*[-]X2
        }
        else
        {
            gw().mulmulsub(*a.X, *b.Y, *a.Y, *b.X, (&a == &res) ? *_tmp : *res.X, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_FFT_S3 | GWMUL_FFT_S4 | GWMUL_STARTNEXTFFT); // F = X1*Y2 - Y1*X2
            gw().mulmuladd(*a.Y, *b.Y, *a.X, *b.X, *res.Y, GWMUL_STARTNEXTFFT); // G = Y1*Y2 + X1*X2
        }
        if (&a == &res)
            swap(*res.X, *_tmp);
        if (!(options & ED_PROJECTIVE))
            gw().mul(*res.Z, *res.T, *_tmp, GWMUL_FFT_S1 | GWMUL_FFT_S2 | options); // T3 = E * H
        gw().mul(*res.Y, *res.Z, *res.Z, options); // Y3 = G * H
        gw().mul(*res.X, *res.T, *res.T, options); // X3 = E * F
        gw().mul(*res.X, *res.Y, *res.Y, options); // Z3 = F * G
        swap(*res.Y, *res.Z);
        swap(*res.X, *res.T);
        if (!(options & ED_PROJECTIVE))
            swap(*res.T, *_tmp);
        else
            res.T.reset();
    }

    // https://eprint.iacr.org/2021/1061
    void EdwardsArithmetic::dbl(EdPoint& a, EdPoint& res, int options)
    {
        bool safe21 = mul_safe(gw().gwdata(), 2, 1);
        bool safe11 = mul_safe(gw().gwdata(), 1, 1);
        bool save_write = safe11 && gw().gwdata()->PASS2_SIZE;

        if (!_tmp && (!(options & ED_PROJECTIVE) || !res.T))
            _tmp.reset(new GWNum(gw()));
        if (!res.X)
            res.X.reset(new GWNum(gw()));
        if (!res.Y)
            res.Y.reset(new GWNum(gw()));
        if (!res.T)
        {
            if (!(options & ED_PROJECTIVE))
                res.T.reset(new GWNum(gw()));
            else
                res.T.reset(_tmp.release());
        }

        gw().setmulbyconst(2);
        gw().mul(*a.X, *a.Y, *res.T, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT); // E = 2X1Y1
        if (a.Z)
        {
            if (!res.Z)
                res.Z.reset(new GWNum(gw()));
            gw().square(*a.Z, *res.Z, GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT_IF(safe21 || !safe11 || save_write)); // C = 2Z1^2
        }
        else
        {
            if (!res.Z)
                res.Z.reset(new GWNum(gw()));
            gw().init(2, *res.Z);
        }
        if (save_write)
        {
            gw().mulmuladd(*a.Y, *a.Y, *a.X, *a.X, *res.Y, GWMUL_STARTNEXTFFT); // G = Y1^2 + X1^2 
            gw().square(*a.X, *res.X, GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT); // G - H = 2X1^2
            if (!(options & ED_PROJECTIVE))
                gw().submul(*res.Y, *res.X, *res.T, *_tmp, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_FFT_S3 | ((options & GWMUL_STARTNEXTFFT) && (options & EDDBL_FOR_EXT_NORM_ADD) ? GWMUL_STARTNEXTFFT_IF(safe11) : options)); // H = G - (G - H), T3 = E * H
            gw().submul(*res.Z, *res.Y, *res.T, *res.T, options); // F = C - G, X3 = F * E
            gw().submul(*res.Z, *res.Y, *res.Y, *res.Z, options); // F = C - G, Z3 = F * G
            gw().submul(*res.Y, *res.X, *res.Y, *res.Y, options); // H = G - (G - H), Y3 = H * G
        }
        else
        {
            gw().square(*a.X, *res.X, GWMUL_STARTNEXTFFT_IF(safe21)); // A = X1^2
            gw().square(*a.Y, *res.Y, GWMUL_STARTNEXTFFT_IF(safe21)); // B = Y1^2
            gw().addsub(*res.Y, *res.X, *res.Y, *res.X, GWADD_DELAYNORM_IF(safe21 || safe11)); // G = B + A, H = B - A
            if (!(options & ED_PROJECTIVE))
                gw().mul(*res.X, *res.T, *_tmp, GWMUL_FFT_S1 | GWMUL_FFT_S2 | ((options & GWMUL_STARTNEXTFFT) && (options & EDDBL_FOR_EXT_NORM_ADD) ? GWMUL_STARTNEXTFFT_IF(safe11) : options)); // T3 = E * H
            gw().sub(*res.Z, *res.Y, *res.Z, GWADD_DELAYNORM_IF(safe21 || !safe11)); // F = C - G
            gw().mul(*res.Z, *res.T, *res.T, options); // X3 = F * E
            gw().mul(*res.Y, *res.Z, *res.Z, options); // Z3 = G * F
            gw().mul(*res.X, *res.Y, *res.Y, options); // Y3 = G * H
        }
        swap(*res.X, *res.T);
        if (!(options & ED_PROJECTIVE))
            swap(*res.T, *_tmp);
        else
            _tmp.reset(res.T.release());
    }

    void EdwardsArithmetic::mul(EdPoint& a, Giant& b, EdPoint& res)
    {
        int len = b.bitlen();
        int W;
        for (W = 2; W < 16 && /*3 + 3*(1 << (W - 1)) <= maxSize &&*/ (14 << (W - 2)) + len/0.69*(7 + 7/(W + 1.0)) > (14 << (W - 1)) + len/0.69*(7 + 7/(W + 2.0)); W++);
        std::vector<int16_t> naf_w;
        get_NAF_W(W, b, naf_w);
        mul(a, W, naf_w, res);
    }

    void EdwardsArithmetic::mul(EdPoint& a, int W, std::vector<int16_t>& naf_w, EdPoint& res)
    {
        int i, j;

        std::vector<std::unique_ptr<EdPoint>> u;
        for (i = 0; i < (1 << (W - 2)); i++)
            u.emplace_back(new EdPoint(*this));

        // Dictionary
        copy(a, *u[0]);
        if (W > 2)
        {
            dbl(a, res, GWMUL_STARTNEXTFFT);
            for (i = 1; i < (1 << (W - 2)); i++)
                add(*u[i - 1], res, *u[i], GWMUL_STARTNEXTFFT);
        }
        if (naf_w.size() > 100)
            normalize(u.begin(), u.end(), 0);

        // Signed window
        copy(*u[naf_w.back()/2], res);
        for (i = (int)naf_w.size() - 2; i >= 0; i--)
        {
            if (naf_w[i] != 0)
            {
                for (j = 1; j < W; j++)
                    dbl(res, res, GWMUL_STARTNEXTFFT | ED_PROJECTIVE);
                dbl(res, res, GWMUL_STARTNEXTFFT | (i > 0 ? 0 : EDDBL_FOR_EXT_NORM_ADD));
                add(res, *u[abs(naf_w[i])/2], res, (i > 0 ? GWMUL_STARTNEXTFFT | ED_PROJECTIVE : 0) | (naf_w[i] < 0 ? EDADD_NEGATIVE : 0));
            }
            else
                dbl(res, res, (i > 0 ? GWMUL_STARTNEXTFFT | ED_PROJECTIVE : 0));
        }
        _tmp.reset();
    }

    void EdwardsArithmetic::normalize(EdPoint& a, int options)
    {
        std::vector<EdPoint*> tmp;
        tmp.push_back(&a);
        normalize(tmp.begin(), tmp.end(), options);
    }

    template <typename Iter>
    void EdwardsArithmetic::normalize(Iter begin, Iter end, int options)
    {
        Iter it;
        Iter first = end;
        Iter last = end;
        for (it = begin; it != end; it++)
        {
            if (!(*it))
                continue;
            if (first == end)
                first = it;
            last = it;
            if (!(*it)->Z)
            {
                (*it)->Z.reset(new GWNum(gw()));
                *(*it)->Z = 1;
            }
            if (!(*it)->T)
                (*it)->T.reset(new GWNum(gw()));
        }
        if (first == end)
            return;
        std::swap((*first)->T, (*first)->Z);
        Iter prev = first;
        for ((it = first)++; it != end; it++)
            if (*it)
            {
                gw().mul(*(*prev)->T, *(*it)->Z, *(*it)->T, it != last ? GWMUL_STARTNEXTFFT : 0);
                prev = it;
            }
        try
        {
            gw().inv(*(*last)->T, *(*last)->T);
        }
        catch (const ArithmeticException&)
        {
            std::swap((*first)->T, (*first)->Z);
            throw;
        }
        for ((it = last)++, prev = last; it != first;)
        {
            it = prev;
            if (it != first)
            {
                for ((prev = it)--; !(*prev); prev--);
                gw().mul(*(*it)->T, *(*prev)->T, *(*prev)->T, GWMUL_STARTNEXTFFT);
                std::swap((*it)->T, (*prev)->T);
                gw().mul(*(*it)->Z, *(*prev)->T, *(*prev)->T, GWMUL_STARTNEXTFFT);
            }
            if ((*it)->X)
                gw().mul(*(*it)->T, *(*it)->X, *(*it)->X, options);
            if ((*it)->Y)
                gw().mul(*(*it)->T, *(*it)->Y, *(*it)->Y, options);
            (*it)->Z.reset();
            if ((*it)->X && (*it)->Y && !(options & ED_PROJECTIVE))
                gw().mul(*(*it)->X, *(*it)->Y, *(*it)->T, GWMUL_FFT_S1 | GWMUL_FFT_S2 | options);
            else
                (*it)->T.reset();
        }
    }

    template void EdwardsArithmetic::normalize(std::vector<std::unique_ptr<EdPoint>>::iterator begin, std::vector<std::unique_ptr<EdPoint>>::iterator end, int options);
    template void EdwardsArithmetic::normalize(std::vector<EdPoint*>::iterator begin, std::vector<EdPoint*>::iterator end, int options);

    GWNum EdwardsArithmetic::jinvariant(GWNum& ed_d)
    {
        // Returns j-invariant = 16*(1 + 14*d + d^2)^3/(d*(1 - d)^4)
        GWNum tmp = ((ed_d + 14)*ed_d + 1);
        return 16*square(tmp)*tmp/(ed_d*square(square(1 - ed_d)));
    }

    EdPoint EdwardsArithmetic::gen_curve(int seed, GWNum* ed_d)
    {
        GWArithmetic& gw = this->gw().carefully();

        int i, len;
        Giant tmp;
        GWNum s(gw), t(gw);
        s = 12;
        t = 40;
        tmp = seed;
        len = tmp.bitlen() - 1;
        for (i = 1; i <= len; i++)
        {
            GWNum lambda = (3*square(s) - 8)/(2*t);
            GWNum snext = lambda*lambda - s - s;
            t = lambda*(std::move(s) - snext) - t;
            s = std::move(snext);
            if (tmp.bit(len - i))
            {
                lambda = (t - 40)/(s - 12);
                snext = lambda*lambda - 12 - s;
                t = lambda*(std::move(s) - snext) - t;
                s = std::move(snext);
            }
        }

        // Asserting T^2 = S^3 - 8S - 32
        //GWASSERT(square(t) == (square(s) - 8)*s - 32);

        // SqrtD = (8A^2 - 1) * (8A^2 + 8A + 1) / (8A^2 + 4A + 1)^2
        // B = A * 2(4A + 1) / (8A^2 - 1)
        GWNum alpha = 1/((t + 25)/(s - 9) + 1);
        GWNum alpha8alpha(gw);
        gw.setmulbyconst(8);
        gw.square(alpha, alpha8alpha, GWMUL_MULBYCONST);
        alpha += alpha;
        GWNum sqrt_d = alpha8alpha - 1;
        alpha8alpha += alpha;
        GWNum beta = alpha8alpha/sqrt_d;
        alpha8alpha += 1;
        alpha8alpha += alpha;
        alpha += alpha;
        sqrt_d *= alpha8alpha + std::move(alpha);
        sqrt_d /= square(std::move(alpha8alpha));
        if (ed_d != nullptr)
            *ed_d = square(sqrt_d);

        GWNum x8 = 2*beta - 1;
        GWNum isdx8 = inv(x8*sqrt_d);
        GWNum x = x8*(4*beta - 3)/(6*beta - 5);
        GWNum y = x8*(t*(t + 50) - 104 - square(s)*(2*s - 27))/((t - 2 + 3*s)*(t + 16 + s));

        Giant gx, gy, gx8, gisdx8;
        (gx = x) %= gw.N();
        (gy = y) %= gw.N();
        (gx8 = x8) %= gw.N();
        (gisdx8 = isdx8) %= gw.N();
        Giant n1 = gw.N() - 1;
        Giant nx8 = gw.N() - gx8;
        Giant nisdx8 = gw.N() - gisdx8;

        // Checking torsion points
        ArithmeticException e("Torsion point.");
        if (gx == 0 && gy == 1)
            throw e;
        if (gx == 0 && gy == n1)
            throw e;
        if (gx == 1 && gy == 0)
            throw e;
        if (gx == n1 && gy == 0)
            throw e;
        if (gx == gx8 && gy == gx8)
            throw e;
        if (gx == nx8 && gy == gx8)
            throw e;
        if (gx == gx8 && gy == nx8)
            throw e;
        if (gx == nx8 && gy == nx8)
            throw e;
        if (gx == gisdx8 && gy == gisdx8)
            throw e;
        if (gx == nisdx8 && gy == gisdx8)
            throw e;
        if (gx == gisdx8 && gy == nisdx8)
            throw e;
        if (gx == nisdx8 && gy == nisdx8)
            throw e;

        // Asserting X^2 + Y^2 = 1 + d * X^2 * Y^2
        //GWASSERT(square(x) + square(y) == 1 + square(sqrt_d)*square(x)*square(y));

        EdPoint p(*this);
        p.X.reset(new GWNum(std::move(x)));
        p.Y.reset(new GWNum(std::move(y)));
        return p;
    }

    EdPoint EdwardsArithmetic::from_small(int32_t xa, int32_t xb, int32_t ya, int32_t yb, GWNum* ed_d)
    {
        Giant gxa, gxb, gya, gyb;
        gxa = xa;
        gxb = xb;
        gya = ya;
        gyb = yb;

        EdPoint p(*this);
        p.X.reset(new GWNum(gw()));
        p.Y.reset(new GWNum(gw()));
        p.Z.reset(new GWNum(gw()));
        p.T.reset(new GWNum(gw()));
        (gxa*gyb).to_GWNum(*p.X);
        (gya*gxb).to_GWNum(*p.Y);
        (gxb*gyb).to_GWNum(*p.Z);
        (gxa*gya).to_GWNum(*p.T);
        //GWASSERT(*p.X*(*p.Y) == *p.Z*(*p.T));

        if (ed_d != nullptr)
        {
            gxa.square();
            gxb.square();
            gya.square();
            gyb.square();
            (gxa*gyb + gya*gxb - gxb*gyb).to_GWNum(*ed_d);
            GWNum t(gw());
            inv(gxa*gya, gw().N()).to_GWNum(t);
            *ed_d *= t;
            //GWASSERT(square(*p.X) + square(*p.Y) == square(*p.Z) + *ed_d*square(*p.T));
        }

        return p;
    }

    bool EdwardsArithmetic::on_curve(EdPoint& a, GWNum& ed_d)
    {
        GWArithmetic& gw = this->gw().carefully();
        GWNum ed_d_a(gw), ed_d_b(gw);
        d_ratio(a, ed_d_a, ed_d_b);
        gw.mul(ed_d_b, ed_d, ed_d_b);
        gw.sub(ed_d_a, ed_d_b, ed_d_a);
        return (gw.popg() = ed_d_a)%gw.N() == 0;
    }

    void EdwardsArithmetic::d_ratio(EdPoint& a, GWNum& ed_d_a, GWNum& ed_d_b)
    {
        GWArithmetic& gw = this->gw().carefully();
        GWNum tmp(gw);
        gw.square(*a.X, ed_d_a, 0);
        gw.square(*a.Y, tmp, 0);
        if (a.T)
            gw.square(*a.T, ed_d_b, 0);
        else
            gw.mul(ed_d_a, tmp, ed_d_b);
        gw.add(ed_d_a, tmp, ed_d_a);
        if (a.Z)
        {
            gw.square(*a.Z, tmp, 0);
            gw.sub(ed_d_a, tmp, ed_d_a);
            if (!a.T)
                gw.mul(ed_d_a, tmp, ed_d_a);
        }
        else
            gw.sub(ed_d_a, 1, ed_d_a);
    }

    void EdPoint::serialize(Giant& X, Giant& Y, Giant& Z, Giant& T)
    {
        if (this->X)
            X = *this->X;
        else
            X = 0;
        if (this->Y)
            Y = *this->Y;
        else
            Y = 0;
        if (this->Z)
            Z = *this->Z;
        else
            Z = 1;
        if (this->T)
            T = *this->T;
        else
            T = 0;
    }

    void EdPoint::deserialize(const Giant& X, const Giant& Y, const Giant& Z, const Giant& T)
    {
        if (!X.empty() && X != 0)
        {
            this->X.reset(new GWNum(arithmetic().gw()));
            *this->X = X;
        }
        else
            this->X.reset();
        if (!Y.empty() && Y != 0)
        {
            this->Y.reset(new GWNum(arithmetic().gw()));
            *this->Y = Y;
        }
        else
            this->Y.reset();
        if (!Z.empty() && Z != 1)
        {
            this->Z.reset(new GWNum(arithmetic().gw()));
            *this->Z = Z;
        }
        else
            this->Z.reset();
        if (!T.empty() && T != 0)
        {
            this->T.reset(new GWNum(arithmetic().gw()));
            *this->T = T;
        }
        else
            this->T.reset();
    }

#ifdef NESTED_EDWARDS
    extern "C" unsigned long cache_line_offset(
        gwhandle *gwdata,	/* Handle initialized by gwsetup */
        unsigned long i);

    void NestedEdwardsArithmetic::set_gw(GWArithmetic& gw)
    {
        _gw = &gw;
        _offsets.clear();

        gwhandle* gwdata = gw.gwdata();
        if (gwdata->EXTRA_BITS < gwdata->fft_max_bits_per_word || !gwdata->ALL_COMPLEX_FFT || gwdata->k != 1.0)
            return;
        for (int i = 0; ; i++)
        {
            int index = i;

            // If we've read all the gwnum data, return FALSE
            if (index * 8 >= (int)gwdata->FFTLEN) break;

            {				// Process complex values where imaginary part is in the high half of the cache line
        //BUG - handle various FFT paddings - addr_offset may not be right as the asm code is allowed to store data in the gaps.
        // It might be easier/faster to just FFT the gap data.
        // Or create new gwnum cache_line_offset that returns cache lines in order skipping the gaps
                _offsets.push_back(cache_line_offset(gwdata, index));
            }
        }
    }

    void NestedEdwardsArithmetic::dbl(EdPoint& a, EdPoint& res, int options)
    {
        if (_offsets.empty() || !(options & ED_PROJECTIVE) || !a.Z)
        {
            EdwardsArithmetic::dbl(a, res, options);
            return;
        }

        gw().fft(*a.X, *a.X);
        gw().fft(*a.Y, *a.Y);
        gw().fft(*a.Z, *a.Z);
        if (!res.X)
            res.X.reset(new GWNum(gw()));
        if (!res.Y)
            res.Y.reset(new GWNum(gw()));
        if (!res.Z)
            res.Z.reset(new GWNum(gw()));
        for (int i = 0; i < _offsets.size(); i++)
        {
            int real_offset = _offsets[i];
            int imag_offset = _offsets[i] + 32;
            __m256d real_x = *(__m256d*)(((char*)**a.X) + real_offset);
            __m256d imag_x = *(__m256d*)(((char*)**a.X) + imag_offset);
            __m256d real_y = *(__m256d*)(((char*)**a.Y) + real_offset);
            __m256d imag_y = *(__m256d*)(((char*)**a.Y) + imag_offset);
            __m256d real_z = *(__m256d*)(((char*)**a.Z) + real_offset);
            __m256d imag_z = *(__m256d*)(((char*)**a.Z) + imag_offset);
            __m256d real_tmp, imag_tmp;
            // 2*X1*Y1
            __m256d real_t = _mm256_sub_pd(_mm256_mul_pd(real_x, real_y), _mm256_mul_pd(imag_x, imag_y));
            __m256d imag_t = _mm256_add_pd(_mm256_mul_pd(real_x, imag_y), _mm256_mul_pd(imag_x, real_y));
            real_t = _mm256_add_pd(real_t, real_t);
            imag_t = _mm256_add_pd(imag_t, imag_t);
            // X1^2
            real_tmp = _mm256_sub_pd(_mm256_mul_pd(real_x, real_x), _mm256_mul_pd(imag_x, imag_x));
            imag_x = _mm256_mul_pd(real_x, imag_x);
            imag_x = _mm256_add_pd(imag_x, imag_x);
            // Y1^2
            real_x = _mm256_sub_pd(_mm256_mul_pd(real_y, real_y), _mm256_mul_pd(imag_y, imag_y));
            imag_tmp = _mm256_mul_pd(real_y, imag_y);
            imag_tmp = _mm256_add_pd(imag_tmp, imag_tmp);
            // Y1^2 +/- X1^2
            real_y = _mm256_add_pd(real_x, real_tmp);
            real_x = _mm256_sub_pd(real_x, real_tmp);
            imag_y = _mm256_add_pd(imag_tmp, imag_x);
            imag_x = _mm256_sub_pd(imag_tmp, imag_x);
            // 2*Z1^2
            real_tmp = _mm256_sub_pd(_mm256_mul_pd(real_z, real_z), _mm256_mul_pd(imag_z, imag_z));
            imag_z = _mm256_mul_pd(real_z, imag_z);
            imag_z = _mm256_add_pd(imag_z, imag_z);
            imag_z = _mm256_add_pd(imag_z, imag_z);
            real_z = _mm256_add_pd(real_tmp, real_tmp);
            // 2*Z1^2 - (Y1^2 + X1^2)
            real_z = _mm256_sub_pd(real_z, real_y);
            imag_z = _mm256_sub_pd(imag_z, imag_y);

            // 2*X1*Y1 * (2*Z1^2 - Y1^2 - X1^2)
            *(__m256d*)(((char*)**res.X) + real_offset) = _mm256_sub_pd(_mm256_mul_pd(real_t, real_z), _mm256_mul_pd(imag_t, imag_z));
            *(__m256d*)(((char*)**res.X) + imag_offset) = _mm256_add_pd(_mm256_mul_pd(real_t, imag_z), _mm256_mul_pd(imag_t, real_z));
            // (Y1^2 - X1^2) * (Y1^2 + X1^2)
            *(__m256d*)(((char*)**res.Y) + real_offset) = _mm256_sub_pd(_mm256_mul_pd(real_x, real_y), _mm256_mul_pd(imag_x, imag_y));
            *(__m256d*)(((char*)**res.Y) + imag_offset) = _mm256_add_pd(_mm256_mul_pd(real_x, imag_y), _mm256_mul_pd(imag_x, real_y));
            // (Y1^2 + X1^2) * (2*Z1^2 - Y1^2 - X1^2)
            *(__m256d*)(((char*)**res.Z) + real_offset) = _mm256_sub_pd(_mm256_mul_pd(real_y, real_z), _mm256_mul_pd(imag_y, imag_z));
            *(__m256d*)(((char*)**res.Z) + imag_offset) = _mm256_add_pd(_mm256_mul_pd(real_y, imag_z), _mm256_mul_pd(imag_y, real_z));
        }
        FFT_state(**res.X) = FULLY_FFTed;
        FFT_state(**res.Y) = FULLY_FFTed;
        FFT_state(**res.Z) = FULLY_FFTed;
        gwunfft2(gw().gwdata(), **res.X, **res.X, options);
        gwunfft2(gw().gwdata(), **res.Y, **res.Y, options);
        gwunfft2(gw().gwdata(), **res.Z, **res.Z, options);
        gw().gwdata()->fft_count += 8;
    }
#endif
}
