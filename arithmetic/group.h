#pragma once

#include <memory>
#include <vector>

#include "arithmetic.h"

namespace arithmetic
{
    void get_NAF_W(int W, Giant& a, std::vector<int16_t>& res);

    template<class Element>
    class GroupArithmetic
    {
    public:
        virtual void copy(const Element& a, Element& res) = 0;
        virtual void move(Element&& a, Element& res) = 0;
        virtual void init(Element& a) = 0;
        virtual void add(Element& a, Element& b, Element& res) = 0;
        virtual void sub(Element& a, Element& b, Element& res) = 0;
        virtual void neg(Element& a, Element& res) = 0;
        virtual void dbl(Element& a, Element& res) = 0;

        virtual void mul(Element& a, int W, std::vector<int16_t>& naf_w, Element& res)
        {
            int i, j;

            std::vector<std::unique_ptr<Element>> u;
            for (i = 0; i < (1 << (W - 2)); i++)
                u.emplace_back(new Element(a.arithmetic()));

            // Dictionary
            copy(a, *u[0]);
            if (W > 2)
            {
                dbl(a, res);
                for (i = 1; i < (1 << (W - 2)); i++)
                    add(*u[i - 1], res, *u[i]);
            }

            // Signed window
            copy(*u[naf_w.back()/2], res);
            for (i = (int)naf_w.size() - 2; i >= 0; i--)
            {
                if (naf_w[i] != 0)
                {
                    for (j = 0; j < W; j++)
                        dbl(res, res);
                    if (naf_w[i] > 0)
                        add(res, *u[naf_w[i]/2], res);
                    else
                        sub(res, *u[-naf_w[i]/2], res);
                }
                else
                    dbl(res, res);
            }
        }
    };

    template<class Arithmetic, class Element>
    class GroupElement
    {
    public:
        GroupElement(Arithmetic& arithmetic) : _arithmetic(arithmetic) { }

        Arithmetic& arithmetic() const { return _arithmetic; }

        Element& operator += (Element& a)
        {
            arithmetic().add((Element&)*this, a, (Element&)*this);
            return (Element&)*this;
        }

        friend Element operator + (Element& a, Element& b)
        {
            Element res(a.arithmetic());
            res.arithmetic().add(a, b, res);
            return res;
        }

        friend Element operator + (Element&& a, Element& b)
        {
            Element res(std::move(a));
            res += ((void*)&a == (void*)&b) ? res : b;
            return res;
        }

        friend Element operator + (Element& a, Element&& b)
        {
            Element res(std::move(b));
            res += ((void*)&a == (void*)&b) ? res : a;
            return res;
        }

        friend Element operator + (Element&& a, Element&& b)
        {
            return std::move(a) + b;
        }

        Element& operator -= (Element& a)
        {
            arithmetic().sub((Element&)*this, a, (Element&)*this);
            return (Element&)*this;
        }

        friend Element operator - (Element& a, Element& b)
        {
            Element res(a.arithmetic());
            res.arithmetic().sub(a, b, res);
            return res;
        }

        friend Element operator - (Element&& a, Element& b)
        {
            Element res(std::move(a));
            res -= ((void*)&a == (void*)&b) ? res : b;
            return res;
        }

        friend Element operator - (Element& a, Element&& b)
        {
            Element res(std::move(b));
            res -= ((void*)&a == (void*)&b) ? res : a;
            res.arithmetic().neg(res, res);
            return res;
        }

        friend Element operator - (Element&& a, Element&& b)
        {
            return std::move(a) - b;
        }

        friend Element operator - (Element& a)
        {
            Element res(a.arithmetic());
            res.arithmetic().neg(a, res);
            return res;
        }

        friend Element operator - (Element&& a)
        {
            Element res(std::move(a));
            res.arithmetic().neg(res, res);
            return res;
        }

        Element& operator *= (Giant& a)
        {
            arithmetic().mul((Element&)*this, a, (Element&)*this);
            return (Element&)*this;
        }

        friend Element operator * (Element& a, Giant& b)
        {
            Element res(a.arithmetic());
            res.arithmetic().mul(a, b, res);
            return res;
        }

        friend Element operator * (Giant& a, Element& b)
        {
            Element res(b.arithmetic());
            res.arithmetic().mul(b, a, res);
            return res;
        }

        friend Element operator * (Element&& a, Giant& b)
        {
            Element res(std::move(a));
            res *= b;
            return res;
        }

        friend Element operator * (Giant& a, Element&& b)
        {
            Element res(std::move(b));
            res *= a;
            return res;
        }

        friend void swap(Element& a, Element& b)
        {
            Element t(std::move(a));
            a = std::move(b);
            b = std::move(t);
        }

    private:
        Arithmetic& _arithmetic;
    };

    int get_DAC_S_d(int e, int start, int end, int *maxlen);
    extern const size_t precomputed_DAC_S_d_len;
    extern const int precomputed_DAC_S_d[];

    template<class Element>
    class DifferentialGroupArithmetic
    {
    public:
        virtual void copy(const Element& a, Element& res) = 0;
        virtual void move(Element&& a, Element& res) = 0;
        virtual void init(Element& a) = 0;
        virtual void add(Element& a, Element& b, Element& a_minus_b, Element& res) = 0;
        virtual void dbl(Element& a, Element& res) = 0;

        virtual void mul(Element& a, int32_t b, Element& res)
        {
            if (b == 0)
                init(res);
            else if (b == 1)
                res = a;
            else if (b == 2)
                dbl(a, res);
            else if (b == 3 && &a != &res)
            {
                dbl(a, res);
                add(a, res, a, res);
            }
            else if (b == 3)
            {
                Element tmp = res;
                dbl(tmp, res);
                add(tmp, res, tmp, res);
            }
            else
            {
                int len;
                for (len = 1; b >= (1 << len); len++);
                int i = b;
                if (len > 2)
                    i >>= len - 2;
                len -= 2;
                int j;
                for (j = 1; (b & j) == 0; j <<= 1);
                while ((1 << len) < j) { i >>= 1; len++; }

                Element tmp = a;
                Element res2(res.arithmetic());
                if (i == 1)
                {
                    if (&a != &res)
                        res = tmp;
                    if ((1 << len) > j)
                        dbl(res, res2);
                }
                if (i == 2)
                {
                    dbl(tmp, res);
                    if ((1 << len) > j)
                        add(tmp, res, tmp, res2);
                }
                if (i == 3)
                {
                    dbl(tmp, res2);
                    add(tmp, res2, tmp, res);
                    if ((1 << len) > j)
                       dbl(res2, res2);
                }

                if ((1 << len) > j*2)
                    optimize(tmp);
                for (i = (1 << (len - 1)); i >= j; i >>= 1)
                {
                    if (b & i)
                    {
                        add(res, res2, tmp, res);
                        if (i > j)
                            dbl(res2, res2);
                    }
                    else
                    {
                        if (i > j)
                            add(res, res2, tmp, res2);
                        dbl(res, res);
                    }
                }
                for (; i; i >>= 1)
                    dbl(res, res);
            }
        }

        virtual void mul(Element& a, int32_t b, Element& res1, Element& res2)
        {
            if (b == 0)
            {
                res2 = a;
                init(res1);
                return;
            }
            int len;
            for (len = 1; b >= (1 << len); len++);
            int i = b;
            if (len > 2)
                i >>= len - 2;
            len -= 2;
            Element tmp = a;
            if (i == 1)
            {
                res1 = tmp;
                dbl(res1, res2);
            }
            if (i == 2)
            {
                dbl(tmp, res1);
                add(tmp, res1, tmp, res2);
            }
            if (i == 3)
            {
                dbl(tmp, res2);
                add(tmp, res2, tmp, res1);
                dbl(res2, res2);
            }
            if (len <= 0)
                return;

            if (len > 1)
                optimize(tmp);
            for (i = (1 << (len - 1)); i; i >>= 1)
            {
                if (b & i)
                {
                    add(res1, res2, tmp, res1);
                    dbl(res2, res2);
                }
                else
                {
                    add(res1, res2, tmp, res2);
                    dbl(res1, res1);
                }
            }
        }

        // https://eprint.iacr.org/2017/293.pdf
        virtual void mul(Element& a, int32_t prime, size_t index, Element& res)
        {
            if (prime < 14)
            {
                mul(a, prime, res);
                return;
            }
            std::vector<int> chain;
            int len = 60;
            int e = prime;
            int d = 1;
            if (index >= 0 && index < precomputed_DAC_S_d_len)
                d = precomputed_DAC_S_d[index];
            else
                d = get_DAC_S_d(e, (int)(e/1.618) - 100, (int)(e/1.618) + 100, &len);
            while (d != 0)
            {
                if (e/2 < d)
                {
                    d = e - d;
                    chain.push_back(0);
                }
                else if ((e & 1) == 0 && d < e/4)
                {
                    chain.push_back(2);
                    e = e/2;
                }
                else
                {
                    chain.push_back(1);
                    e = e - d;
                }
            }
            d = 1;
            e = 2;
            int ed = 1;
            bool ed_init = false;
            Element Td = a;
            Element Ted(a.arithmetic());
            Element Te = std::move(res);
            dbl(Td, Te);
            for (int i = (int)chain.size() - 3; i >= 0; i--)
            {
                if (ed != e - d)
                    printf("error");
                if (chain[i] == 2)
                {
                    ed = ed + e;
                    e = 2*e;
                    add(ed_init ? Ted : Td, Te, Td, Ted);
                    dbl(Te, Te);
                    ed_init = true;
                }
                else if (chain[i] == 1)
                {
                    ed = e;
                    e = e + d;
                    add(Te, Td, ed_init ? Ted : Td, Ted);
                    swap(Te, Ted);
                    ed_init = true;
                }
                else
                {
                    ed = d;
                    d = e - d;
                    if (ed_init)
                        swap(Td, Ted);
                }
            }
            if (e != prime)
                printf("error");
            res = std::move(Te);
        }

        virtual void mul(Element& a, Giant& b, Element& res1, Element& res2)
        {
            Element tmp = a;
            optimize(tmp);
            res1 = a;
            dbl(res1, res2);

            int len = b.bitlen() - 1;
            for (int i = 1; i <= len; i++)
            {
                if (b.bit(len - i))
                {
                    add(res1, res2, tmp, res1);
                    dbl(res2, res2);
                }
                else
                {
                    add(res1, res2, tmp, res2);
                    dbl(res1, res1);
                }
            }
        }

        virtual void optimize(Element& a) = 0;

    private:
    };

    template<class Arithmetic, class Element>
    class DifferentialGroupElement
    {
    public:
        DifferentialGroupElement(Arithmetic& arithmetic) : _arithmetic(arithmetic) { }

        Arithmetic& arithmetic() const { return _arithmetic; }

        Element& operator *= (int32_t a)
        {
            arithmetic().mul((Element&)*this, a, (Element&)*this);
            return (Element&)*this;
        }

        friend Element operator * (Element& a, int32_t b)
        {
            Element res(a.arithmetic());
            res.arithmetic().mul(a, b, res);
            return res;
        }

        friend Element operator * (int32_t a, Element& b)
        {
            Element res(b.arithmetic());
            res.arithmetic().mul(b, a, res);
            return res;
        }

        friend Element operator * (Element&& a, int32_t b)
        {
            Element res(std::move(a));
            res *= b;
            return res;
        }

        friend Element operator * (int32_t a, Element&& b)
        {
            Element res(std::move(b));
            res *= a;
            return res;
        }

        friend void swap(Element& a, Element& b)
        {
            Element t(std::move(a));
            a = std::move(b);
            b = std::move(t);
        }

    private:
        Arithmetic& _arithmetic;
    };
}
