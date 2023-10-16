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
        bool& parity() { return _parity; }
        const bool& parity() const { return _parity; }

    private:
        std::unique_ptr<GWNum> _V;
        bool _parity;
    };

    class LucasUV;

    class LucasUVArithmetic : public GroupArithmetic<LucasUV>
    {
        friend class LucasUV;
    public:
        using Element = LucasUV;
        static const int LUCASADD_OPTIMIZE = 0x100000;
        static const int LUCASADD_NEGATIVE = 0x200000;

    public:
        LucasUVArithmetic(GWArithmetic& gw, int P, bool negativeQ = false) : _gw(&gw), _negativeQ(negativeQ)
        {
            if (P < (GWSMALLMUL_MAX - 4)/P)
            {
                int D = P*P - (negativeQ ? -4 : 4);
                _UV_small.emplace_back(0, 2, 0);
                _UV_small.emplace_back(1, P, D);
                while (true)
                {
                    int U = (_UV_small.back().V + P*_UV_small.back().U)/2;
                    int V = (P*_UV_small.back().V + D*_UV_small.back().U)/2;
                    if (U > GWSMALLMUL_MAX/D)
                        break;
                    _UV_small.emplace_back(U, V, D*U);
                }
            }
            else
            {
                Giant tmp;
                tmp = P;
                tmp *= P;
                tmp -= (negativeQ ? -4 : 4);
                _D.reset(new GWNum(gw));
                *_D = tmp;
            }
            
        }
        template<class T>
        LucasUVArithmetic(GWArithmetic& gw, T&& D, bool negativeQ = false) : _gw(&gw), _negativeQ(negativeQ), _D(new GWNum(gw))
        {
            *_D = std::forward<T>(D);
        }
        virtual ~LucasUVArithmetic() { }

        virtual void copy(const LucasUV& a, LucasUV& res) override;
        virtual void move(LucasUV&& a, LucasUV& res) override;
        virtual void init(LucasUV& res) override;
        virtual void init(GWNum& P, LucasUV& res);
        virtual void init(GWNum& U, GWNum& V, LucasUV& res);
        virtual void init(GWNum& U, GWNum& V, bool parity, bool halved, LucasUV& res);
        virtual void init_small(int index, LucasUV& res);
        virtual void add(LucasUV& a, LucasUV& b, LucasUV& res) override;
        virtual void sub(LucasUV& a, LucasUV& b, LucasUV& res) override;
        virtual void add(LucasUV& a, LucasUV& b, LucasUV& res, int options);
        virtual void neg(LucasUV& a, LucasUV& res) override;
        virtual void dbl(LucasUV& a, LucasUV& res) override;
        virtual void dbl(LucasUV& a, LucasUV& res, int options);
        virtual void dbl_add_small(LucasUV& a, int index, LucasUV& res, int options);
        virtual void mul(LucasUV& a, Giant& b, LucasUV& res);
        virtual void mul(LucasUV& a, int W, std::vector<int16_t>& naf_w, LucasUV& res) override;
        virtual void optimize(LucasUV& a);

        GWArithmetic& gw() { return *_gw; }
        void set_gw(GWArithmetic& gw) { _gw = &gw; }
        bool negativeQ() { return _negativeQ; }
        int max_small() { if (_D) return 0; return (int)_UV_small.size() - 1; }

    private:
        GWArithmetic* _gw;
        std::unique_ptr<GWNum> _D;
        std::unique_ptr<GWNum> _half;
        std::unique_ptr<GWNum> _half_fma;
        std::unique_ptr<GWNum> _tmp;
        std::unique_ptr<GWNum> _tmp2;

        bool _negativeQ;
        struct LucasUV_small
        {
            int U, V, DU;
            LucasUV_small(int U, int V, int DU) { this->U = U; this->V = V; this->DU = DU; }
        };
        std::vector<LucasUV_small> _UV_small;
    };

    class LucasUV : public GroupElement<LucasUVArithmetic, LucasUV>
    {
        friend class LucasUVArithmetic;
    public:
        using Arithmetic = LucasUVArithmetic;

    public:
        LucasUV(LucasUVArithmetic& arithmetic) : GroupElement<LucasUVArithmetic, LucasUV>(arithmetic), _U(new GWNum(arithmetic.gw())), _V(new GWNum(arithmetic.gw()))
        {
            //arithmetic.init(*this);
        }
        template<class T>
        LucasUV(LucasUVArithmetic& arithmetic, T&& P) : GroupElement<LucasUVArithmetic, LucasUV>(arithmetic), _U(new GWNum(arithmetic.gw())), _V(new GWNum(arithmetic.gw()))
        {
            V() = std::forward<T>(P);
            arithmetic.init(V(), *this);
        }
        template<class TU, class TV>
        LucasUV(LucasUVArithmetic& arithmetic, TU&& U, TV&& V, bool parity = true, bool halved = true) : GroupElement<LucasUVArithmetic, LucasUV>(arithmetic), _U(new GWNum(arithmetic.gw())), _V(new GWNum(arithmetic.gw()))
        {
            U() = std::forward<TU>(U);
            V() = std::forward<TV>(V);
            if (!halved)
                arithmetic.init(U(), V(), *this);
        }
        ~LucasUV()
        {
        }
        LucasUV(const LucasUV& a) : GroupElement<LucasUVArithmetic, LucasUV>(a.arithmetic()), _U(new GWNum(a.arithmetic().gw())), _V(new GWNum(a.arithmetic().gw()))
        {
            arithmetic().copy(a, *this);
        }
        LucasUV(LucasUV&& a) noexcept : GroupElement<LucasUVArithmetic, LucasUV>(a.arithmetic())
        {
            arithmetic().move(std::move(a), *this);
        }

        LucasUV& operator = (const LucasUV& a)
        {
            arithmetic().copy(a, *this);
            return *this;
        }
        LucasUV& operator = (LucasUV&& a) noexcept
        {
            arithmetic().move(std::move(a), *this);
            return *this;
        }

        friend bool operator == (const LucasUV& a, const LucasUV& b)
        {
            return a.U() == b.U() && a.V() == b.V();
        }
        friend bool operator != (const LucasUV& a, const LucasUV& b)
        {
            return a.U() != b.U() || a.V() != b.V();
        }

        GWNum& U() { return *_U; }
        const GWNum& U() const { return *_U; }
        GWNum& V() { return *_V; }
        const GWNum& V() const { return *_V; }
        bool parity() const { return _parity; }

    private:
        std::unique_ptr<GWNum> _U;
        std::unique_ptr<GWNum> _V;
        std::unique_ptr<GWNum> _DU;
        bool _parity;
    };
}
