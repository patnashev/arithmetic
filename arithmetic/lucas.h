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
        LucasVArithmetic() : _gw(nullptr) { }
        LucasVArithmetic(GWArithmetic& gw) : _gw(&gw) { }
        virtual ~LucasVArithmetic() { }

        virtual void copy(const LucasV& a, LucasV& res) override;
        virtual void move(LucasV&& a, LucasV& res) override;
        virtual void init(LucasV& res) override;
        virtual void init(const GWNum& P, LucasV& res);
        virtual void add(LucasV& a, LucasV& b, LucasV& a_minus_b, LucasV& res) override;
        virtual void dbl(LucasV& a, LucasV& res) override;
        virtual void optimize(LucasV& a) override;

        GWArithmetic& gw() { return *_gw; }
        void set_gw(GWArithmetic& gw) { _gw = &gw; }

    private:
        GWArithmetic* _gw;
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
        LucasV(LucasVArithmetic& arithmetic, const GWNum& P) : DifferentialGroupElement<LucasVArithmetic, LucasV>(arithmetic), _V(new GWNum(arithmetic.gw()))
        {
            arithmetic.init(P, *this);
        }
        LucasV(LucasVArithmetic& arithmetic, gwnum a) : DifferentialGroupElement<LucasVArithmetic, LucasV>(arithmetic), _V(new GWNumWrapper(arithmetic.gw(), a))
        {
        }
        ~LucasV()
        {
        }
        LucasV(const LucasV& a) : DifferentialGroupElement<LucasVArithmetic, LucasV>(a.arithmetic()), _V(new GWNum(a.arithmetic().gw()))
        {
            arithmetic().init(a.V(), *this);
        }
        LucasV(LucasV&& a) noexcept : DifferentialGroupElement<LucasVArithmetic, LucasV>(a.arithmetic()), _V(std::move(a._V))
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

    private:
        std::unique_ptr<GWNum> _V;
    };
}
