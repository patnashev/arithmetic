#pragma once

#include <memory>
#include <vector>
#include "giant.h"
#include "arithmetic.h"

class File;

class InputNum
{
public:
    static const int ZERO = 0;
    static const int GENERIC = 1;
    static const int KBNC = 2;
    static const int FACTORIAL = 3;
    static const int PRIMORIAL = 4;

    static const int ALGEBRAIC_SIMPLE = 0;
    static const int ALGEBRAIC_CYCLOTOMIC = 1; //     X^3 -+ 1 = (X^2 +- X + 1)(X -+ 1)
    static const int ALGEBRAIC_QUAD = 2;       //  1/4 X^4 + 1 = (1/2 X^2 + X + 1)(1/2 X^2 - X + 1)
    static const int ALGEBRAIC_HEX = 3;        // 1/27 x^6 + 1 = (1/3 X^2 + X + 1)(1/3 X^2 - X + 1)(1/3 X^2 + 1)

public:
    InputNum() { _gk = 0; _gb = 0; _gd = 1; }
    InputNum(int k, int b, int n, int c) { init(k, b, n, c); }
    template<class T>
    InputNum(T&& b) { _type = GENERIC; _gk = 1; _gb = std::forward<T>(b); _n = 0; _gd = 1; _c = 0; process(); }

    template<class TK, class TB>
    void init(TK&& k, TB&& b, int n, int c) { _type = KBNC; _gk = std::forward<TK>(k); _gb = std::forward<TB>(b); _n = n; _gd = 1; _c = c; _custom_k.clear(); _custom_b.clear(); _custom_d.clear(); _algebraic_type = ALGEBRAIC_SIMPLE; _algebraic_k = 0; process(); }
    template<class T>
    void init(T&& b) { _type = GENERIC; _gk = 1; _gb = std::forward<T>(b); _n = 0; _gd = 1; _c = 0; _custom_k.clear(); _custom_b.clear(); _custom_d.clear(); _algebraic_type = ALGEBRAIC_SIMPLE; _algebraic_k = 0; process(); }
    bool read(File& file);
    void write(File& file);
    bool parse(const std::string& s, bool c_required = true);
    void setup(arithmetic::GWState& state);
    void print_info();

    //deprecated
    static uint64_t parse_numeral(const std::string& s);

    bool empty() const { return _gb == 0; }
    int type() { return _type; }
    uint64_t k() { return _gk.size() == 2 && _gk.bitlen() < 52 ? *(uint64_t*)(_gk.data()) : _gk.size() == 1 ? *(_gk.data()) : 0; }
    uint32_t b() { return _gb.size() == 1 ? *(_gb.data()) : 0; }
    uint32_t n() { return _type == GENERIC ? 1 : _n; }
    uint32_t d() { return _gd.size() == 1 ? *(_gd.data()) : 0; }
    int32_t c() { return _c; }
    arithmetic::Giant& gk() { return _gk; }
    arithmetic::Giant& gb() { return _gb; }
    arithmetic::Giant& gd() { return _gd; }
    int gfn() { return _gfn; }
    int multifactorial() { return _multifactorial; }
    int algebraic_type() { return _algebraic_type; }
    int algebraic_k() { return _algebraic_k; }
    arithmetic::Giant value() { return _type == GENERIC ? _gb : _type != KBNC ? _gk*_gb/_gd + _c : _gk*power(_gb, _n)/_gd + _c; }
    uint32_t mod(uint32_t modulus);
    uint32_t fingerprint() { return mod(3417905339UL); }

    int bitlen();

    bool is_half_factored();
    void add_factor(arithmetic::Giant& factor);
    std::vector<std::pair<arithmetic::Giant, int>>& b_factors() { return _b_factors; }
    arithmetic::Giant& b_cofactor() { return _b_cofactor; }
    std::vector<std::pair<arithmetic::Giant, int>>& factors() { return _factors; }
    arithmetic::Giant& cofactor() { return _cofactor; }
    std::vector<int> factorize_minus1(int depth);
    void factorize_f_p();

    const std::string& input_text() { return _input_text; }
    const std::string& display_text() { return _display_text; }

private:
    void process();
    std::string build_text(int max_len = -1);

private:
    int _type = ZERO;
    arithmetic::Giant _gk;
    arithmetic::Giant _gb;
    uint32_t _n = 0;
    arithmetic::Giant _gd;
    int32_t _c = 0;
    std::vector<std::pair<arithmetic::Giant, int>> _b_factors;
    arithmetic::Giant _b_cofactor;
    std::vector<std::pair<arithmetic::Giant, int>> _factors;
    arithmetic::Giant _cofactor;
    std::string _input_text;
    std::string _display_text;
    std::string _custom_k;
    std::string _custom_b;
    std::string _custom_d;
    int _gfn = 0;
    int _multifactorial = 0;
    int _algebraic_type = ALGEBRAIC_SIMPLE;
    int32_t _algebraic_k = 0;
};
