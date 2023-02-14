#pragma once

#include "group.h"

namespace arithmetic
{
    class LucasV;

    class LucasVArithmetic : public DifferentialGroupArithmetic<LucasV>
    {
        friend class LucasV;
    public:
        using Element = LucasV;

    public:
        LucasVArithmetic(bool negativeQ = false) : _gw(nullptr), _negativeQ(negativeQ) { }
        LucasVArithmetic(GWArithmetic& gw, bool negativeQ = false) : _gw(&gw), _negativeQ(negativeQ) { }
        virtual ~LucasVArithmetic() { }

        virtual void copy(const LucasV& a, LucasV& res) override;
        virtual void move(LucasV&& a, LucasV& res) override;
        virtual void init(LucasV& res) override;
        virtual void init(const GWNum& P, LucasV& res);
        virtual void init(const GWNum& P, bool parity, LucasV& res);
        virtual void add(LucasV& a, LucasV& b, LucasV& a_minus_b, LucasV& res) override;
        virtual void dbl(LucasV& a, LucasV& res) override;
        virtual void optimize(LucasV& a) override;

        GWArithmetic& gw() { return *_gw; }
        void set_gw(GWArithmetic& gw) { _gw = &gw; }
        bool negativeQ() { return _negativeQ; }

    private:
        GWArithmetic* _gw;
        bool _negativeQ;
    };

    class LucasV : public DifferentialGroupElement<LucasVArithmetic, LucasV>
    {
        friend class LucasVArithmetic;
    public:
        using Arithmetic = LucasVArithmetic;

    public:
        LucasV(LucasVArithmetic& arithmetic) : DifferentialGroupElement<LucasVArithmetic, LucasV>(arithmetic), _V(new GWNum(arithmetic.gw()))
        {
            arithmetic.init(*this);
        }
        template<class T>
        LucasV(LucasVArithmetic& arithmetic, T&& P, bool parity = true) : DifferentialGroupElement<LucasVArithmetic, LucasV>(arithmetic), _V(new GWNum(arithmetic.gw()))
        {
            *_V = std::forward<T>(P);
            _parity = parity;
        }
        LucasV(LucasVArithmetic& arithmetic, gwnum a) : DifferentialGroupElement<LucasVArithmetic, LucasV>(arithmetic), _V(new GWNumWrapper(arithmetic.gw(), a))
        {
        }
        ~LucasV()
        {
        }
        LucasV(const LucasV& a) : DifferentialGroupElement<LucasVArithmetic, LucasV>(a.arithmetic()), _V(new GWNum(a.arithmetic().gw()))
        {
            arithmetic().init(a.V(), a.parity(), *this);
        }
        LucasV(LucasV&& a) noexcept : DifferentialGroupElement<LucasVArithmetic, LucasV>(a.arithmetic()), _V(std::move(a._V)), _parity(a._parity)
        {
        }

        LucasV& operator = (const LucasV& a)
        {
            arithmetic().copy(a, *this);
            return *this;
        }
        LucasV& operator = (LucasV&& a) noexcept
        {
            arithmetic().move(std::move(a), *this);
            return *this;
        }

        LucasV& operator = (const GWNum& P)
        {
            arithmetic().init(P, *this);
            return *this;
        }
        friend bool operator == (const LucasV& a, const LucasV& b)
        {
            return a.V() == b.V();
        }
        friend bool operator != (const LucasV& a, const LucasV& b)
        {
            return a.V() != b.V();
        }

        GWNum& V() { return *_V; }
        const GWNum& V() const { return *_V; }
        bool parity() const { return _parity; }

    private:
        std::unique_ptr<GWNum> _V;
        bool _parity;
    };
}
