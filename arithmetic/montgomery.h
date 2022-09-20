#pragma once

#include "group.h"
#include "edwards.h"

namespace arithmetic
{
    class EdY;

    class MontgomeryArithmetic : public DifferentialGroupArithmetic<EdY>
    {
        friend class EdY;
    public:
        using Element = EdY;

    public:
        MontgomeryArithmetic(GWNum& ed_d) : _gw(nullptr), _ed_d(ed_d) { }
        MontgomeryArithmetic(GWArithmetic& gw, GWNum& ed_d) : _gw(&gw), _ed_d(ed_d) { }
        virtual ~MontgomeryArithmetic() { }

        virtual void copy(const EdY& a, EdY& res) override;
        virtual void move(EdY&& a, EdY& res) override;
        virtual void init(EdY& res) override;
        virtual void init(const EdPoint& a, EdY& res);
        virtual int cmp(const EdY& a, const EdY& b);
        virtual void add(EdY& a, EdY& b, EdY& a_minus_b, EdY& res) override;
        virtual void dbl(EdY& a, EdY& res) override;
        virtual void optimize(EdY& a) override;

        virtual void normalize(EdY& a);
        template <typename Iter>
        void normalize(Iter begin, Iter end);

        GWArithmetic& gw() { return *_gw; }
        void set_gw(GWArithmetic& gw) { _gw = &gw; }
        GWNum& ed_d() { return _ed_d; }

    private:
        GWArithmetic* _gw;
        GWNum& _ed_d;
    };

    class EdY : public DifferentialGroupElement<MontgomeryArithmetic, EdY>
    {
        friend class MontgomeryArithmetic;
    public:
        using Arithmetic = MontgomeryArithmetic;

    public:
        EdY(MontgomeryArithmetic& arithmetic) : DifferentialGroupElement<MontgomeryArithmetic, EdY>(arithmetic)
        {
            //arithmetic.init(*this);
        }
        EdY(MontgomeryArithmetic& arithmetic, const EdPoint& a) : DifferentialGroupElement<MontgomeryArithmetic, EdY>(arithmetic)
        {
            arithmetic.init(a, *this);
        }
        ~EdY()
        {
        }
        EdY(const EdY& a) : DifferentialGroupElement<MontgomeryArithmetic, EdY>(a.arithmetic())
        {
            arithmetic().copy(a, *this);
        }
        EdY(EdY&& a) noexcept : DifferentialGroupElement<MontgomeryArithmetic, EdY>(a.arithmetic())
        {
            arithmetic().move(std::move(a), *this);
        }

        EdY& operator = (const EdY& a)
        {
            arithmetic().copy(a, *this);
            return *this;
        }
        EdY& operator = (EdY&& a) noexcept
        {
            arithmetic().move(std::move(a), *this);
            return *this;
        }

        EdY& operator = (const EdPoint& a)
        {
            arithmetic().init(a, *this);
            return *this;
        }
        friend bool operator == (const EdY& a, const EdY& b)
        {
            return a.arithmetic().cmp(a, b) == 0;
        }
        friend bool operator != (const EdY& a, const EdY& b)
        {
            return a.arithmetic().cmp(a, b) != 0;
        }
        EdY& normalize()
        {
            arithmetic().normalize(*this);
            return *this;
        }

        void serialize(Giant& Y, Giant& Z);
        void deserialize(const Giant& Y, const Giant& Z);

    public:
        std::unique_ptr<GWNum> Y;
        std::unique_ptr<GWNum> Z;

        std::unique_ptr<GWNum> ZpY;
        std::unique_ptr<GWNum> ZmY;
    };

}
