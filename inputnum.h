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

public:
    InputNum() { _gk = 0; _gb = 0; }
    InputNum(int k, int b, int n, int c) { init(k, b, n, c); }
    template<class T>
    InputNum(T&& b) { _type = GENERIC; _gk = 1; _gb = std::forward<T>(b); _n = 0; _c = 0; process(); }

    template<class TK, class TB>
    void init(TK&& k, TB&& b, int n, int c) { _type = KBNC; _gk = std::forward<TK>(k); _gb = std::forward<TB>(b); _n = n; _c = c; _custom_k.clear(); _custom_b.clear(); _cyclotomic = 0; _cyclotomic_k = 0; _hex_k = 0; process(); }
    template<class T>
    void init(T&& b) { _type = GENERIC; _gk = 1; _gb = std::forward<T>(b); _n = 0; _c = 0; _custom_k.clear(); _custom_b.clear(); _cyclotomic = 0; _cyclotomic_k = 0; _hex_k = 0; process(); }
    bool read(File& file);
    void write(File& file);
    bool parse(const std::string& s);
    void setup(arithmetic::GWState& state);
    bool is_base2();
    void to_base2(InputNum& k, InputNum& base2);
    bool is_factorized_half();
    void add_factor(arithmetic::Giant& factor);
    void print_info();

    static uint64_t parse_numeral(const std::string& s);

    bool empty() const { return _gb == 0; }
    int type() { return _type; }
    uint32_t k() { return _gk.size() == 1 ? *(_gk.data()) : 0; }
    uint32_t b() { return _gb.size() == 1 ? *(_gb.data()) : 0; }
    uint32_t n() { return _type == KBNC ? _n : 1; }
    int32_t c() { return _c; }
    arithmetic::Giant& gk() { return _gk; }
    arithmetic::Giant& gb() { return _gb; }
    int gfn() { return _gfn; }
    int cyclotomic() { return _cyclotomic; }
    arithmetic::Giant value() { return _type == GENERIC ? _gb : _type != KBNC ? _gk*_gb + _c : _gk*power(_gb, _n) + _c; }
    uint32_t fingerprint();
    int bitlen();

    std::vector<std::pair<arithmetic::Giant, int>>& b_factors() { return _b_factors; }
    arithmetic::Giant& b_cofactor() { return _b_cofactor; }
    std::vector<int> factorize_minus1(int depth);

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
    int32_t _c = 0;
    std::vector<std::pair<arithmetic::Giant, int>> _b_factors;
    arithmetic::Giant _b_cofactor;
    std::vector<std::pair<arithmetic::Giant, int>> _factors;
    arithmetic::Giant _cofactor;
    std::string _input_text;
    std::string _display_text;
    std::string _custom_k;
    std::string _custom_b;
    int _gfn = 0;
    int _cyclotomic = 0;
    uint32_t _cyclotomic_k = 0;
    uint32_t _hex_k = 0;
};
