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

    template<class TK, class TB>
    void init(TK k, TB b, int n, int c) { _type = KBNC; _gk = k; _gb = b; _n = n; _c = c; process(); }
    bool read(File& file);
    void write(File& file);
    bool parse(const std::string& s);
    void setup(arithmetic::GWState& state);
    bool is_base2();
    void to_base2(InputNum& k, InputNum& base2);
    bool is_factorized_half();
    void add_factor(arithmetic::Giant& factor);

    static uint64_t parse_numeral(const std::string& s);

    bool empty() const { return _gb == 0; }
    uint32_t k() { return _gk.size() == 1 ? *(_gk.data()) : 0; }
    uint32_t b() { return _gb.size() == 1 ? *(_gb.data()) : 0; }
    uint32_t n() { return _type == KBNC ? _n : 1; }
    int32_t c() { return _c; }
    arithmetic::Giant& gk() { return _gk; }
    arithmetic::Giant& gb() { return _gb; }
    int gfn() { return _gfn; }
    arithmetic::Giant value() { return _type == GENERIC ? _gb : _type != KBNC ? _gk*_gb + _c : _gk*power(_gb, _n) + _c; }
    uint32_t fingerprint();
    int bitlen();

    std::vector<std::pair<arithmetic::Giant, int>>& b_factors() { return _b_factors; }
    std::unique_ptr<arithmetic::Giant>& b_cofactor() { return _b_cofactor; }
    std::vector<int> factorize_minus1(int depth);

    const std::string& input_text() { return _input_text; }
    const std::string& display_text() { return _display_text; }

private:
    void add_factor(uint32_t factor, int power = 1);
    void add_factor(const arithmetic::Giant& factor, int power = 1);
    void process();
    std::string build_text(int max_len = -1);

private:
    int _type = ZERO;
    arithmetic::Giant _gk;
    arithmetic::Giant _gb;
    uint32_t _n = 0;
    int32_t _c = 0;
    std::vector<std::pair<arithmetic::Giant, int>> _b_factors;
    std::unique_ptr<arithmetic::Giant> _b_cofactor;
    std::string _input_text;
    std::string _display_text;
    int _gfn = 0;
};
