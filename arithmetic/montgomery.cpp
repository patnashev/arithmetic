
#include <deque>
#include <stdlib.h>
#include "gwnum.h"
#include "montgomery.h"
#include "exception.h"

namespace arithmetic
{
    void MontgomeryArithmetic::copy(const EdY& a, EdY& res)
    {
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
        if (a.ZpY)
        {
            if (!res.ZpY)
                res.ZpY.reset(new GWNum(gw()));
            *res.ZpY = *a.ZpY;
        }
        else
            res.ZpY.reset();
        if (a.ZmY)
        {
            if (!res.ZmY)
                res.ZmY.reset(new GWNum(gw()));
            *res.ZmY = *a.ZmY;
        }
        else
            res.ZmY.reset();
    }

    void MontgomeryArithmetic::move(EdY&& a, EdY& res)
    {
        res.Y = std::move(a.Y);
        res.Z = std::move(a.Z);
        res.ZpY = std::move(a.ZpY);
        res.ZmY = std::move(a.ZmY);
    }

    void MontgomeryArithmetic::init(EdY& a)
    {
        a.Y.reset(new GWNum(gw()));
        *a.Y = 1;
        a.Z.reset();
        a.ZpY.reset();
        a.ZmY.reset();
    }

    void MontgomeryArithmetic::init(const EdPoint& a, EdY& res)
    {
        if (!res.Y)
            res.Y.reset(new GWNum(gw()));
        *res.Y = *a.Y;
        if (a.Z)
        {
            res.Z.reset(new GWNum(gw()));
            *res.Z = *a.Z;
        }
        else
            res.Z.reset();
        res.ZpY.reset();
        res.ZmY.reset();
    }

    int MontgomeryArithmetic::cmp(const EdY& a, const EdY& b)
    {
        if (!a.Y && !b.Y)
            return 0;
        if (!a.Y && b.Y)
            return -1;
        if (a.Y && !b.Y)
            return 1;
        if (!a.Z && !b.Z)
            return gw().cmp(*a.Y, *b.Y);
        if (!a.Z && b.Z)
            return gw().cmp(*a.Y*(*b.Z), *b.Y);
        if (a.Z && !b.Z)
            return gw().cmp(*a.Y, *b.Y*(*a.Z));
        return gw().cmp(*a.Y*(*b.Z), *b.Y*(*a.Z));
    }

    void MontgomeryArithmetic::add(EdY& a, EdY& b, EdY& a_minus_b, EdY& res)
    {
        bool safe1 = square_safe(gw().gwdata(), 1);
        bool safe11 = mul_safe(gw().gwdata(), 1, 1);

        if (!a_minus_b.Y)
        {
            dbl(a, res);
            return;
        }
        optimize(a_minus_b);
        bool normalized_a = !a.Z;
        bool normalized_b = !b.Z;
        if (!res.Y)
            res.Y.reset(new GWNum(gw()));
        if (!res.Z)
            res.Z.reset(new GWNum(gw()));
        if (!res.ZpY)
            res.ZpY.reset(new GWNum(gw()));
        if (!res.ZmY)
            res.ZmY.reset(new GWNum(gw()));
        std::unique_ptr<GWNum> tmp;
        if (&a_minus_b == &res || (!normalized_b && (&b == &res)))
            tmp.reset(new GWNum(gw()));

        if (!normalized_b && (&b == &res))
            std::swap(tmp, res.Z);
        if (!normalized_a)
            gw().mul(*a.Z, *b.Y, *res.Z, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT_IF(safe1));
        else
            gw().copy(*b.Y, *res.Z);
        if (!normalized_b)
            gw().mul(*a.Y, (&b == &res) ? *tmp : *b.Z, *res.Y, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT_IF(safe1));
        else if (&a != &res)
            gw().copy(*a.Y, *res.Y);
        gw().addsub(*res.Y, *res.Z, *res.Z, *res.Y, GWADD_DELAYNORM_IF(safe1));
        gw().square(*res.Z, *res.Z, GWMUL_STARTNEXTFFT);
        gw().square(*res.Y, *res.Y, GWMUL_STARTNEXTFFT);
        gw().mulmuladd(*a_minus_b.ZmY, *res.Z, *a_minus_b.ZpY, *res.Y, (&a_minus_b == &res) ? *tmp : *res.ZpY, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_FFT_S3 | GWMUL_FFT_S4 | GWMUL_STARTNEXTFFT_IF(safe11));
        gw().mulmulsub(*a_minus_b.ZmY, *res.Z, *a_minus_b.ZpY, *res.Y, *res.ZmY, GWMUL_STARTNEXTFFT_IF(safe11 && ((!normalized_a && !normalized_b) || safe1)));
        if (&a_minus_b == &res)
            std::swap(res.ZpY, tmp);
        if (safe11)
            gw().fft(*res.ZpY, *res.ZpY);
        if (safe11 && ((!normalized_a && !normalized_b) || safe1))
            gw().fft(*res.ZmY, *res.ZmY);
        gw().copy(*res.ZpY, *res.Z);
        gw().copy(*res.ZmY, *res.Y);
        gw().addsub(*res.ZpY, *res.ZmY, *res.ZpY, *res.ZmY, GWADD_DELAYNORM_IF(safe11));
    }

    void MontgomeryArithmetic::dbl(EdY& a, EdY& res)
    {
        bool safe1 = square_safe(gw().gwdata(), 1);
        bool safe11 = mul_safe(gw().gwdata(), 1, 1);

        bool normalized = !a.Z;
        if (!res.Y)
            res.Y.reset(new GWNum(gw()));
        if (!res.Z)
            res.Z.reset(new GWNum(gw()));
        if (!res.ZpY)
            res.ZpY.reset(new GWNum(gw()));
        if (!res.ZmY)
            res.ZmY.reset(new GWNum(gw()));

        // d = 1 - 1/Ad4
        // t1 = zz*(yy - d*yy)
        // t2 = (zz - yy)*(zz - d*yy)
        // y_2 = t1 - t2
        // z_2 = t1 + t2

        gw().square(*a.Y, *res.Y, GWMUL_STARTNEXTFFT_IF(safe11)); // yy
        if (!normalized)
            gw().square(*a.Z, *res.Z, GWMUL_STARTNEXTFFT_IF(safe11)); // zz
        else
            gw().init(1, *res.Z);
        gw().sub(*res.Z, *res.Y, *res.ZmY, GWADD_DELAYNORM_IF(safe11)); // zz - yy
        gw().mul(ed_d(), *res.Y, *res.ZpY, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT); // d*yy
        gw().submul(*res.Z, *res.ZpY, *res.ZmY, *res.ZmY, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT_IF(safe1)); // t2
        gw().submul(*res.Y, *res.ZpY, *res.Z, *res.Z, GWMUL_STARTNEXTFFT_IF(safe1)); // t1
        gw().copy(*res.ZmY, *res.Y);
        gw().copy(*res.Z, *res.ZpY);
        gw().addsub(*res.Z, *res.Y, *res.Z, *res.Y, GWADD_DELAYNORM_IF(safe1));
    }

    void MontgomeryArithmetic::optimize(EdY& a)
    {
        bool safe11 = mul_safe(gw().gwdata(), 1, 1);

        //normalize(a);
        if (!a.ZpY)
        {
            a.ZpY.reset(new GWNum(gw()));
            a.ZmY.reset(new GWNum(gw()));
            if (a.Z)
                *a.ZpY = *a.Z;
            else
                *a.ZpY = 1;
            *a.ZmY = *a.Y;
            gw().addsub(*a.ZpY, *a.ZmY, *a.ZpY, *a.ZmY, GWADD_DELAYNORM_IF(safe11));
        }
    }

    void MontgomeryArithmetic::normalize(EdY& a)
    {
        std::vector<EdY*> tmp;
        tmp.push_back(&a);
        normalize(tmp.begin(), tmp.end());
    }

    template <typename Iter>
    void MontgomeryArithmetic::normalize(Iter begin, Iter end)
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
            if (!(*it)->ZpY)
                (*it)->ZpY.reset(new GWNum(gw()));
            (*it)->ZmY.reset();
        }
        if (first == end)
            return;
        std::swap((*first)->ZpY, (*first)->Z);
        Iter prev = first;
        for ((it = first)++; it != end; it++)
            if (*it)
            {
                gw().mul(*(*prev)->ZpY, *(*it)->Z, *(*it)->ZpY, it != last ? GWMUL_STARTNEXTFFT : 0);
                prev = it;
            }
        try
        {
            gw().inv(*(*last)->ZpY, *(*last)->ZpY);
        }
        catch (const ArithmeticException&)
        {
            std::swap((*first)->ZpY, (*first)->Z);
            throw;
        }
        for ((it = last)++, prev = last; it != first;)
        {
            it = prev;
            if (it != first)
            {
                for ((prev = it)--; !(*prev); prev--);
                gw().mul(*(*it)->ZpY, *(*prev)->ZpY, *(*prev)->ZpY, GWMUL_STARTNEXTFFT);
                std::swap((*it)->ZpY, (*prev)->ZpY);
                gw().mul(*(*it)->Z, *(*prev)->ZpY, *(*prev)->ZpY, GWMUL_STARTNEXTFFT);
            }
            if ((*it)->Y)
                gw().mul(*(*it)->ZpY, *(*it)->Y, *(*it)->Y, GWMUL_STARTNEXTFFT);
            (*it)->Z.reset();
            (*it)->ZpY.reset();
        }
    }

    template void MontgomeryArithmetic::normalize(std::vector<std::unique_ptr<EdY>>::iterator begin, std::vector<std::unique_ptr<EdY>>::iterator end);
    template void MontgomeryArithmetic::normalize(std::vector<EdY*>::iterator begin, std::vector<EdY*>::iterator end);
    template void MontgomeryArithmetic::normalize(std::deque<std::unique_ptr<EdY>>::iterator begin, std::deque<std::unique_ptr<EdY>>::iterator end);

    void EdY::serialize(Giant& Y, Giant& Z)
    {
        if (this->Y)
            Y = *this->Y;
        else
            Y = 0;
        if (this->Z)
            Z = *this->Z;
        else
            Z = 1;
    }

    void EdY::deserialize(const Giant& Y, const Giant& Z)
    {
        if (Y != 0)
        {
            this->Y.reset(new GWNum(arithmetic().gw()));
            *this->Y = Y;
        }
        else
            this->Y.reset();
        if (Z != 1)
        {
            this->Z.reset(new GWNum(arithmetic().gw()));
            *this->Z = Z;
        }
        else
            this->Z.reset();
        ZpY.reset();
        ZmY.reset();
    }
}
