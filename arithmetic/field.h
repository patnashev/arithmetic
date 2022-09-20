#pragma once

#include <string>

namespace arithmetic
{
    template<class Element>
    class FieldArithmetic
    {
    public:
        virtual void alloc(Element& a) = 0;
        virtual void free(Element& a) = 0;
        virtual void copy(const Element& a, Element& res) = 0;
        virtual void move(Element&& a, Element& res) = 0;
        virtual void init(int32_t a, Element& res) = 0;
        virtual void init(const std::string& a, Element& res) = 0;
        virtual int cmp(const Element& a, const Element& b) = 0;
        virtual int cmp(const Element& a, int32_t b) = 0;
        virtual void add(Element& a, Element& b, Element& res) = 0;
        virtual void add(Element& a, int32_t b, Element& res) = 0;
        virtual void sub(Element& a, Element& b, Element& res) = 0;
        virtual void sub(Element& a, int32_t b, Element& res) = 0;
        virtual void neg(Element& a, Element& res) = 0;
        virtual void mul(Element& a, Element& b, Element& res) = 0;
        virtual void mul(Element& a, int32_t b, Element& res) = 0;
        virtual void div(Element& a, Element& b, Element& res) = 0;
        virtual void div(Element& a, int32_t b, Element& res) = 0;
        virtual void gcd(Element& a, Element& b, Element& res) = 0;
        virtual void inv(Element& a, Element& n, Element& res) = 0;
        virtual void mod(Element& a, Element& b, Element& res) = 0;
    };

    template<class Arithmetic, class Element>
    class FieldElement
    {
    public:
        FieldElement(Arithmetic& arithmetic) : _arithmetic(arithmetic) { }

        Arithmetic& arithmetic() const { return _arithmetic; }

        virtual std::string to_string() const = 0;

        Element& operator = (int32_t a)
        {
            arithmetic().init(a, (Element&)*this);
            return (Element&)*this;
        }

        Element& operator = (const std::string& a)
        {
            arithmetic().init(a, (Element&)*this);
            return (Element&)*this;
        }

        friend bool operator == (const Element& a, const Element& b)
        {
            return a.arithmetic().cmp(a, b) == 0;
        }

        friend bool operator != (const Element& a, const Element& b)
        {
            return a.arithmetic().cmp(a, b) != 0;
        }

        friend bool operator > (const Element& a, const Element& b)
        {
            return a.arithmetic().cmp(a, b) > 0;
        }

        friend bool operator < (const Element& a, const Element& b)
        {
            return a.arithmetic().cmp(a, b) < 0;
        }

        friend bool operator >= (const Element& a, const Element& b)
        {
            return a.arithmetic().cmp(a, b) >= 0;
        }

        friend bool operator <= (const Element& a, const Element& b)
        {
            return a.arithmetic().cmp(a, b) <= 0;
        }

        friend bool operator == (const Element& a, int b)
        {
            return a.arithmetic().cmp(a, b) == 0;
        }

        friend bool operator != (const Element& a, int b)
        {
            return a.arithmetic().cmp(a, b) != 0;
        }

        friend bool operator > (const Element& a, int b)
        {
            return a.arithmetic().cmp(a, b) > 0;
        }

        friend bool operator < (const Element& a, int b)
        {
            return a.arithmetic().cmp(a, b) < 0;
        }

        friend bool operator >= (const Element& a, int b)
        {
            return a.arithmetic().cmp(a, b) >= 0;
        }

        friend bool operator <= (const Element& a, int b)
        {
            return a.arithmetic().cmp(a, b) <= 0;
        }

        friend bool operator == (int a, const Element& b)
        {
            return b.arithmetic().cmp(b, a) == 0;
        }

        friend bool operator != (int a, const Element& b)
        {
            return b.arithmetic().cmp(b, a) != 0;
        }

        friend bool operator > (int a, const Element& b)
        {
            return b.arithmetic().cmp(b, a) < 0;
        }

        friend bool operator < (int a, const Element& b)
        {
            return b.arithmetic().cmp(b, a) > 0;
        }

        friend bool operator >= (int a, const Element& b)
        {
            return b.arithmetic().cmp(b, a) <= 0;
        }

        friend bool operator <= (int a, const Element& b)
        {
            return b.arithmetic().cmp(b, a) >= 0;
        }

        Element& add(Element& a)
        {
            arithmetic().add((Element&)*this, a, (Element&)*this);
            return (Element&)*this;
        }

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

        Element& operator += (int32_t a)
        {
            arithmetic().add((Element&)*this, a, (Element&)*this);
            return (Element&)*this;
        }

        friend Element operator + (Element& a, int32_t b)
        {
            Element res(a.arithmetic());
            res.arithmetic().add(a, b, res);
            return res;
        }

        friend Element operator + (int32_t a, Element& b)
        {
            Element res(b.arithmetic());
            res.arithmetic().add(b, a, res);
            return res;
        }

        friend Element operator + (Element&& a, int32_t b)
        {
            Element res(std::move(a));
            res += b;
            return res;
        }

        friend Element operator + (int32_t a, Element&& b)
        {
            Element res(std::move(b));
            res += a;
            return res;
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

        Element& operator -= (int32_t a)
        {
            arithmetic().sub((Element&)*this, a, (Element&)*this);
            return (Element&)*this;
        }

        friend Element operator - (Element& a, int32_t b)
        {
            Element res(a.arithmetic());
            res.arithmetic().sub(a, b, res);
            return res;
        }

        friend Element operator - (int32_t a, Element& b)
        {
            Element res(b.arithmetic());
            res.arithmetic().init(a, res);
            res.arithmetic().sub(res, b, res);
            return res;
        }

        friend Element operator - (Element&& a, int32_t b)
        {
            Element res(std::move(a));
            res -= b;
            return res;
        }

        friend Element operator - (int32_t a, Element&& b)
        {
            Element res(std::move(b));
            res -= a;
            res.arithmetic().neg(res, res);
            return res;
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

        Element& square()
        {
            arithmetic().mul((Element&)*this, (Element&)*this, (Element&)*this);
            return (Element&)*this;
        }

        friend Element square(Element& a)
        {
            Element res(a.arithmetic());
            res.arithmetic().mul(a, a, res);
            return res;
        }

        friend Element square(Element&& a)
        {
            Element res(std::move(a));
            res.square();
            return res;
        }

        Element& operator *= (Element& a)
        {
            arithmetic().mul((Element&)*this, a, (Element&)*this);
            return (Element&)*this;
        }

        Element& operator *= (Element&& a)
        {
            *this *= a;
            return (Element&)*this;
        }

        friend Element operator * (Element& a, Element& b)
        {
            Element res(a.arithmetic());
            res.arithmetic().mul(a, b, res);
            return res;
        }

        friend Element operator * (Element&& a, Element& b)
        {
            Element res(std::move(a));
            res *= ((void*)&a == (void*)&b) ? res : b;
            return res;
        }

        friend Element operator * (Element& a, Element&& b)
        {
            Element res(std::move(b));
            res *= ((void*)&a == (void*)&b) ? res : a;
            return res;
        }

        friend Element operator * (Element&& a, Element&& b)
        {
            return std::move(a) * b;
        }

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

        Element& operator /= (Element& a)
        {
            arithmetic().div((Element&)*this, a, (Element&)*this);
            return (Element&)*this;
        }

        Element& operator /= (int32_t a)
        {
            arithmetic().div((Element&)*this, a, (Element&)*this);
            return (Element&)*this;
        }

        friend Element operator / (Element& a, int32_t b)
        {
            Element res(a.arithmetic());
            res.arithmetic().div(a, b, res);
            return res;
        }

        friend Element operator / (Element&& a, int32_t b)
        {
            Element res(std::move(a));
            res /= b;
            return res;
        }

        Element& operator %= (Element& a)
        {
            arithmetic().mod((Element&)*this, a, (Element&)*this);
            return (Element&)*this;
        }

        friend Element operator % (Element& a, Element& b)
        {
            Element res(a.arithmetic());
            res.arithmetic().mod(a, b, res);
            return res;
        }

        friend Element operator % (Element&& a, Element& b)
        {
            Element res(std::move(a));
            res %= ((void*)&a == (void*)&b) ? res : b;
            return res;
        }

        friend Element operator % (Element& a, Element&& b)
        {
            return a % b;
        }

        friend Element operator % (Element&& a, Element&& b)
        {
            return std::move(a) % b;
        }

        Element& gcd(Element& a)
        {
            arithmetic().gcd((Element&)*this, a, (Element&)*this);
            return (Element&)*this;
        }

        friend Element gcd(Element& a, Element& b)
        {
            Element res(a.arithmetic());
            res.arithmetic().gcd(a, b, res);
            return res;
        }

        friend Element gcd(Element&& a, Element& b)
        {
            Element res(std::move(a));
            res.gcd(((void*)&a == (void*)&b) ? res : b);
            return res;
        }

        friend Element gcd(Element& a, Element&& b)
        {
            Element res(std::move(b));
            res.gcd(((void*)&a == (void*)&b) ? res : a);
            return res;
        }

        friend Element gcd(Element&& a, Element&& b)
        {
            return gcd(std::move(a), b);
        }

        Element& inv(Element& n)
        {
            arithmetic().inv((Element&)*this, n, (Element&)*this);
            return (Element&)*this;
        }

        friend Element inv(Element& a, Element& n)
        {
            Element res(a.arithmetic());
            res.arithmetic().inv(a, n, res);
            return res;
        }

        friend Element inv(Element&& a, Element& b)
        {
            Element res(std::move(a));
            res.inv(((void*)&a == (void*)&b) ? res : b);
            return res;
        }

        friend Element inv(Element& a, Element&& b)
        {
            return inv(a, b);
        }

        friend Element inv(Element&& a, Element&& b)
        {
            return inv(std::move(a), b);
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
