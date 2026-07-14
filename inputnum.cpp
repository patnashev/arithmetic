
#include <climits>
#include <cctype>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <stdlib.h>
#include <map>
#include "gwnum.h"
#include "cpuid.h"
#include "inputnum.h"
#include "file.h"
#include "edwards.h"
#include "integer.h"
#include "exception.h"

using namespace arithmetic;

bool InputNum::read(File& file)
{
    std::unique_ptr<Reader> reader(file.get_reader());
    if (!reader)
        return false;
    if (reader->type() != 0)
        return false;
    if (!reader->read(_gb))
        return false;
    _type = GENERIC;
    _gk = 1;
    _n = 0;
    _c = 0;
    _input_text = file.filename();
    _display_text = file.filename();
    return true;
}

void InputNum::write(File& file)
{
    std::unique_ptr<Writer> writer(file.get_writer(0, 0));
    if (_type != KBNC)
        writer->write(_gk*_gb + _c);
    else
        writer->write(_gk*power(_gb, _n) + _c);
    file.commit_writer(*writer);
}

void add_factor(std::vector<std::pair<arithmetic::Giant, int>>& factors, const Giant& factor, int power)
{
    auto it = factors.begin();
    for (; it != factors.end() && it->first != factor; it++);
    if (it != factors.end())
    {
        it->second += power;
        if (it->second == 0)
            factors.erase(it);
    }
    else
        factors.emplace_back(factor, power);
}

void add_factor(std::vector<std::pair<arithmetic::Giant, int>>& factors, uint32_t factor, int power)
{
    Giant tmp;
    tmp = factor;
    add_factor(factors, tmp, power);
}

void InputNum::add_factor(Giant& factor)
{
    if (!_b_cofactor.empty() && _b_cofactor%factor == 0)
    {
        ::add_factor(_b_factors, factor, 1);
        if (std::abs(_c) == 1)
            ::add_factor(_factors, factor, _n);
        _b_cofactor /= factor;
        if (_b_cofactor == 1)
            _b_cofactor.arithmetic().free(_b_cofactor);
        if (!_cofactor.empty())
        {
            _cofactor /= power(factor, _n);
            if (_cofactor == 1)
                _cofactor.arithmetic().free(_cofactor);
        }
    }
    else if (!_cofactor.empty() && _cofactor%factor == 0)
    {
        ::add_factor(_factors, factor, 1);
        _cofactor /= factor;
        if (_cofactor == 1)
            _cofactor.arithmetic().free(_cofactor);
    }
}

void factorize(Giant& N, std::vector<std::pair<arithmetic::Giant, int>>& factors, Giant& cofactor, std::function<bool(Giant&, uint32_t)> is_factor = nullptr, uint32_t s = 10)
{
    uint32_t i, j;
    int power;
    Giant tmp = N;
    for (i = 0; !tmp.bit(i); i++);
    if (i > 0)
    {
        add_factor(factors, 2, i);
        tmp >>= i;
    }
    std::vector<bool> bitmap;
    if (tmp > 1)
    {
        bitmap.resize((size_t)1 << (s - 1), false);
        std::vector<std::pair<int, int>> smallprimes;
        for (i = 1; i < bitmap.size(); i++)
            if (!bitmap[i])
            {
                smallprimes.emplace_back(i*2 + 1, (i*2 + 1)*(i*2 + 1)/2);
                if (i < ((size_t)1 << (s/2 - 1)))
                    for (; smallprimes.back().second < bitmap.size(); smallprimes.back().second += smallprimes.back().first)
                        bitmap[smallprimes.back().second] = true;
                if ((is_factor != nullptr && is_factor(tmp, i*2 + 1)) || (is_factor == nullptr && tmp%(i*2 + 1) == 0))
                {
                    for (power = 1, tmp /= i*2 + 1; tmp%(i*2 + 1) == 0; power++, tmp /= i*2 + 1);
                    add_factor(factors, i*2 + 1, power);
                }
            }
        if (tmp > 1 && tmp < (1 << (2*s)))
        {
            add_factor(factors, tmp, 1);
            tmp = 1;
        }
        for (j = 0; j < s && tmp > 1; j += 5)
        {
            bitmap.resize((size_t)1 << (s - 1 + (j + 5 < s ? j + 5 : s)), false);
            for (auto it = smallprimes.begin(); it != smallprimes.end(); it++)
                for (; it->second < bitmap.size(); it->second += it->first)
                    bitmap[it->second] = true;
            for (i = 1 << (s - 1 + j); i < bitmap.size(); i++)
                if (!bitmap[i] && ((is_factor != nullptr && is_factor(tmp, i*2 + 1)) || (is_factor == nullptr && tmp%(i*2 + 1) == 0)))
                {
                    for (power = 1, tmp /= i*2 + 1; tmp%(i*2 + 1) == 0; power++, tmp /= i*2 + 1);
                    add_factor(factors, i*2 + 1, power);
                }
            if (tmp > 1 && (uint32_t)tmp.bitlen() <= 2*(s + (j + 5 < s ? j + 5 : s)))
            {
                add_factor(factors, tmp, 1);
                tmp = 1;
            }
        }
    }
    if (tmp > 1)
    {
        if (cofactor.empty())
            cofactor = std::move(tmp);
        else
            cofactor *= tmp;
    }

/*    if (_b_cofactor)
    {
        int len;
        GWState gwstate;
        gwstate.setup(*_b_cofactor);
        GWArithmetic gw(gwstate);
        {
            GWNum y(gw);
            y = 3;
            gw.setmulbyconst(3);
            len = _b_cofactor->bitlen() - 1;
            for (i = 1; i <= len; i++)
                gw.mul(y, y, y, _b_cofactor->bit(len - i) ? GWMUL_MULBYCONST : 0);
            if (y == 3)
            {
                add_factor(*_b_cofactor, 1);
                _b_cofactor.reset();
            }
        }
    }

    if (_b_cofactor)
    {
        int len = _b_cofactor->bitlen()/10;
        if (len > 2*s)
            len = 2*s;
        tmp.arithmetic().alloc(tmp, (1 << (len - 5))/0.69);
        Giant tmp2(tmp.arithmetic(), tmp.capacity() < 8192 ? tmp.capacity() : 8192);
        tmp = 0;
        tmp2 = 1 << len;
        for (i = 1; i < bitmap.size() && i < (1 << len); i++)
            if (!bitmap[i])
            {
                uint32_t p = i*2 + 1;
                j = p;
                if (p <= (1 << (len/2 + 1)))
                {
                    uint32_t pp = (1 << len)/p;
                    while (j <= pp)
                        j *= p;
                }
                tmp2 *= j;
                if (tmp2.size() > 8190)
                {
                    if (tmp == 0)
                        tmp = tmp2;
                    else
                        tmp *= tmp2;
                    tmp2 = 1;
                }
            }
        if (tmp == 0)
            tmp = std::move(tmp2);
        else
            tmp *= tmp2;
        len = tmp.bitlen();
        int W;
        for (W = 2; W < 16 && (14 << (W - 2)) + len/0.69*(7 + 7/(W + 1.0)) >(14 << (W - 1)) + len/0.69*(7 + 7/(W + 2.0)); W++);
        std::vector<int16_t> naf_w;
        get_NAF_W(W, tmp, naf_w);
        GWState gwstate;
        gwstate.setup(*_b_cofactor);
        GWArithmetic gw(gwstate);
        for (j = 0; true; j++)
        {
            double timer = getHighResTimer();
            {
                EdwardsArithmetic ed(gw);
                GWNum ed_d(gw);
                EdPoint P = ed.gen_curve(*(int *)&timer, &ed_d);
                ed.mul(P, W, naf_w, P);
                if (!ed.on_curve(P, ed_d))
                {
                    gwstate.done();
                    gwstate.next_fft_count++;
                    gwstate.setup(*_b_cofactor);
                }
                tmp = gcd(*P.X, *_b_cofactor);
            }
            if (tmp != 1)
            {
                GWState rgwstate;
                rgwstate.setup(tmp);
                GWArithmetic rgw(rgwstate);
                GWNum x(rgw);
                x = 3;
                rgw.setmulbyconst(3);
                len = tmp.bitlen() - 1;
                for (i = 1; i <= len; i++)
                    rgw.mul(x, x, x, tmp.bit(len - i) ? GWMUL_MULBYCONST : 0);
                if (x == 3)
                {
                    add_factor(tmp, 1);
                    *_b_cofactor /= tmp;
                    gwstate.done();
                    gwstate.setup(*_b_cofactor);
                    GWNum y(gw);
                    y = 3;
                    gw.setmulbyconst(3);
                    len = _b_cofactor->bitlen() - 1;
                    for (i = 1; i <= len; i++)
                        gw.mul(y, y, y, _b_cofactor->bit(len - i) ? GWMUL_MULBYCONST : 0);
                    if (y == 3)
                    {
                        add_factor(*_b_cofactor, 1);
                        _b_cofactor.reset();
                        break;
                    }
                }
            }
        }
    }*/
}

class Expr
{
public:
    Expr()
    {
        value = 0;
    }

    bool evaluate(bool factorize = true)
    {
        if (str_func.empty())
        {
            if (args.size() != 1 || str_args(0).empty())
                return error(0, "value expected");
            value = str_args(0);
            if (factorize)
                ::factorize(value, factors, cofactor);
        }
        else if (str_func == "p")
        {
            if (args.size() != 1 || str_args(0).empty())
                return error(0, "prime index expected");
            std::string::const_iterator it;
            for (it = str_args(0).begin(); it != str_args(0).end() && std::isdigit(*it); it++);
            if (it != str_args(0).end())
            {
                InputNum recursive;
                InputNum::ParseResult res = recursive.parse(str_args(0), false);
                if (!res)
                    return error(args[0].first + res.pos, res.message);
                value = recursive.value();
            }
            else
                value = str_args(0);
            if (value > 105097565)
                return error(args[0].first, "prime index too big");
            if (value == 0)
                value = 1;
            else
            {
                value = *PrimeIterator::get(value.data()[0] - 1);
                if (factorize)
                    ::add_factor(factors, value, 1);
            }
        }
        else if (str_func == "()")
        {
            if (args.size() != 1 || str_args(0).empty())
                return error(0, "expression expected");
            InputNum recursive;
            InputNum::ParseResult res = recursive.parse(str_args(0), false);
            if (!res)
                return error(args[0].first + res.pos, res.message);
            value = recursive.value();
            str_value = "(" + recursive.input_text() + ")";
            if (factorize)
            {
                if (recursive.c() == 0)
                {
                    factors = std::move(recursive.factors());
                    cofactor = std::move(recursive.cofactor());
                }
                else
                    ::factorize(value, factors, cofactor);
            }
        }
        else
            return error(0, "unknown function");
        if (value.empty() || value == 0)
            return error(0, "zero value");

        return true;
    }

    void merge_factors(std::vector<std::pair<arithmetic::Giant, int>>& factors_, Giant& cofactor_, int power)
    {
        for (auto& factor : factors)
            ::add_factor(factors_, factor.first, factor.second*power);
        if (!cofactor.empty())
        {
            cofactor.power(abs(power));
            if (power < 0)
            {
                if (cofactor_.empty() || cofactor_%cofactor != 0)
                    throw ArithmeticException();
                cofactor_ /= cofactor;
            }
            else if (!cofactor_.empty())
                cofactor_ *= cofactor;
            else
                cofactor_ = cofactor;
        }
    }

    std::string& str_args(int index) { return args[index].second; }

private:
    bool error(int pos, const char* msg) { error_pos = pos; error_msg = msg; return false; }
    bool error(int pos, std::string& msg) { error_pos = pos; error_msg = std::move(msg); return false; }

public:
    std::string str_func;
    std::vector<std::pair<int, std::string>> args;
    
    Giant value;
    std::string str_value;

    std::vector<std::pair<arithmetic::Giant, int>> factors;
    Giant cofactor;

    int error_pos = -1;
    std::string error_msg;
};

class Token
{
public:
    Token() { }

    bool iterate(const std::string& str, std::string::const_iterator& it)
    {
        pos = (int)(it - str.begin());
        if (it == str.end())
            return false;
        if (*it == '+' || *it == '-' || *it == '*' || *it == '/' || *it == '^' || *it == '!' || *it == '#')
        {
            oper = *it;
            it++;
            return true;
        }

        std::string::const_iterator it_pos = it;
        std::string::const_iterator it_s;
        if (std::isdigit(*it))
        {
            expr.reset(new Expr());
            for (it_s = it; it != str.end() && std::isdigit(*it); it++);
            expr->args.emplace_back((int)(it_s - it_pos), std::string(it_s, it));
            return true;
        }
        if (std::isalpha(*it))
        {
            expr.reset(new Expr());
            for (it_s = it; it != str.end() && std::isalpha(*it); it++);
            expr->str_func = std::string(it_s, it);
            if (it == str.end())
                return true;
            if (std::isdigit(*it))
            {
                for (it_s = it; it != str.end() && std::isdigit(*it); it++);
                expr->args.emplace_back((int)(it_s - it_pos), std::string(it_s, it));
                return true;
            }
        }
        if (*it == '(')
        {
            if (!expr)
            {
                expr.reset(new Expr());
                expr->str_func = "()";
            }
            while (it != str.end() && *it != ')')
            {
                it++;
                int recursion = 0;
                for (it_s = it; it != str.end() && (recursion > 0 || (*it != ')' && *it != ',')); it++)
                {
                    if (*it == '(')
                        recursion++;
                    if (*it == ')')
                        recursion--;
                }
                expr->args.emplace_back((int)(it_s - it_pos), std::string(it_s, it));
            }
            if (it != str.end() && *it == ')')
                it++;
            return true;
        }

        return false;
    }

public:
    int pos = 0;
    char oper = 0;
    std::unique_ptr<Expr> expr;
};

Giant factorial(Giant& factor)
{
    uint32_t n = factor.data()[0];
    uint32_t multifactorial = factor.data()[1];
    Giant res, tmp;
    res = 1;
    tmp = 1;
    uint32_t i = n%multifactorial;
    while (i < 2)
        i += multifactorial;
    for (; i <= n; i += multifactorial)
    {
        tmp *= i;
        if (tmp.size() > 8192 || i == n)
        {
            res *= tmp;
            tmp = 1;
        }
    }
    return res;
}

Giant primorial(Giant& factor)
{
    int n = (int)factor.data()[0];
    Giant res, tmp;
    res = 1;
    tmp = 1;
    int last = 1;
    for (auto it = PrimeIterator::get(); *it <= n; it++)
    {
        last = *it;
        tmp *= last;
        if (tmp.size() > 8192)
        {
            res *= tmp;
            tmp = 1;
        }
    }
    if (tmp != 1)
        res *= tmp;
    factor.data()[0] = last;
    return res;
}

InputNum::ParseResult InputNum::parse(const std::string& s, bool c_required)
{
    Giant tmp;

    int type = KBNC;
    Giant gk, gb, gd;
    uint32_t n = 0;
    int64_t c = 0;
    std::string custom_k;
    std::string custom_b;
    std::string custom_d;
    uint32_t multifactorial = 0;
    int algebraic_type = ALGEBRAIC_SIMPLE;
    int32_t algebraic_k = 0;
    std::vector<std::pair<arithmetic::Giant, int>> factors;
    Giant cofactor;
    std::vector<std::pair<arithmetic::Giant, int>> b_factors;
    Giant b_cofactor;

    std::vector<Token> tokens;
    std::string::const_iterator its;
    for (its = s.begin(); its != s.end() && std::isspace(*its); its++);
    if (its != s.end() && *its == '\"')
        its++;
    while (its != s.end())
    {
        Token token;
        if (!token.iterate(s, its))
            break;
        tokens.push_back(std::move(token));
    }
    if (its != s.end() && *its == '\"')
        its++;
    for (; its != s.end() && std::isspace(*its); its++);
    if (its != s.end())
        return InputNum::ParseResult(false, (int)(its - s.begin()), "excess symbols");

    auto it_expr = tokens.begin();
    if (it_expr == tokens.end())
        return InputNum::ParseResult(false, 0, "no value");
    auto it_oper = it_expr + 1;

    if (it_expr->expr && it_expr->expr->str_func == "()")
    {
        bool bf = false;
        auto itF = it_oper;
        while (itF != tokens.end() && itF->oper == '/')
        {
            itF++;
            if (itF == tokens.end() || !itF->expr)
            {
                bf = false;
                break;
            }
            bf = true;
            itF++;
        }
        if (itF != tokens.end())
            bf = false;
        if (bf)
        {
            InputNum recursive;
            InputNum::ParseResult res = recursive.parse(it_expr->expr->str_args(0));
            if (!res)
                return InputNum::ParseResult(false, it_expr->pos + it_expr->expr->args[0].first + res.pos, res.message);
            it_expr++;

            Giant gf;
            std::string custom_f;
            while (it_expr != tokens.end())
            {
                it_expr++;
                if (!it_expr->expr->evaluate(false))
                    return InputNum::ParseResult(false, it_expr->pos + it_expr->expr->error_pos, it_expr->expr->error_msg);
                if (gf.empty())
                {
                    gf = std::move(it_expr->expr->value);
                    custom_f = std::move(it_expr->expr->str_value);
                }
                else
                {
                    if (!custom_f.empty() || !it_expr->expr->str_value.empty())
                    {
                        if (custom_f.empty())
                            custom_f = gf.to_string();
                        if (it_expr->expr->str_value.empty())
                            it_expr->expr->str_value = it_expr->expr->value.to_string();
                        custom_f += "/";
                        custom_f += it_expr->expr->str_value;
                    }
                    gf *= it_expr->expr->value;
                }
                it_expr++;
            }
            if (recursive.value()%gf != 0)
                return InputNum::ParseResult(false, (int)s.size(), "not divisible");

            _type = recursive._type;
            _gk = std::move(recursive._gk);
            _gb = std::move(recursive._gb);
            _n = recursive._n;
            _c = recursive._c;
            _custom_k = std::move(recursive._custom_k);
            _custom_b = std::move(recursive._custom_b);
            _custom_d = std::move(recursive._custom_d);
            _multifactorial = recursive._multifactorial;
            _algebraic_type = recursive._algebraic_type;
            _algebraic_k = recursive._algebraic_k;
            _factors = std::move(recursive._factors);
            _cofactor = std::move(recursive._cofactor);
            _b_factors = std::move(recursive._b_factors);
            _b_cofactor = std::move(recursive._b_cofactor);

            _input_text = "(" + recursive._input_text + ")/" + (!custom_f.empty() ? custom_f : gf.to_string());
            if (recursive._gf != 1 && (!custom_f.empty() || !recursive._custom_f.empty()))
            {
                if (custom_f.empty())
                    custom_f = gf.to_string();
                if (recursive._custom_f.empty())
                    recursive._custom_f = recursive._gf.to_string();
                custom_f = recursive._custom_f + "/" + custom_f;
            }
            if (_type == GENERIC)
            {
                if (!_custom_b.empty())
                    _custom_b = "(" + _custom_b + ")/" + (!custom_f.empty() ? custom_f : gf.to_string());
                custom_f.clear();
                _gb /= gf;
                gf = 1;
            }
            _gf = std::move(gf)*recursive._gf;
            _custom_f = std::move(custom_f);
            _display_text = build_text(30);

            return InputNum::ParseResult(true);
        }
    }

    if (it_expr->expr && !it_expr->expr->str_func.empty() && it_expr->expr->str_func != "p" && (it_oper == tokens.end()))
    {
        if (it_expr->expr->str_func == "Phi") // (X +- 1)*X + 1
        {
            if (it_expr->expr->args.size() != 2)
                return InputNum::ParseResult(false, it_expr->pos, "two arguments expected");
            int cyclotomic = stoi(it_expr->expr->str_args(0));
            if (cyclotomic != 3 && cyclotomic != 6)
                return InputNum::ParseResult(false, it_expr->pos + it_expr->expr->args[0].first, "only 3 and 6 are supported");
            if (!it_expr->expr->str_args(1).empty() && it_expr->expr->str_args(1)[0] == '-')
            {
                cyclotomic = (cyclotomic == 3 ? 6 : 3);
                it_expr->expr->str_args(1) = it_expr->expr->str_args(1).substr(1);
            }
            InputNum recursive;
            InputNum::ParseResult res = recursive.parse(it_expr->expr->str_args(1) + "+1");
            if (!res)
                return InputNum::ParseResult(false, it_expr->pos + it_expr->expr->args[1].first + res.pos, res.message);
            type = recursive.type();
            gk = recursive.value();
            custom_k = "(" + recursive.input_text() + ")";
            if (cyclotomic == 6)
            {
                gk -= 2;
                custom_k[custom_k.size() - 3] = '-';
            }
            factors = std::move(recursive.factors());
            cofactor = std::move(recursive.cofactor());
            b_factors = std::move(recursive.b_factors());
            b_cofactor = std::move(recursive.b_cofactor());
            ::factorize(gk, factors, cofactor);
            if (recursive.k() != 1)
            {
                gk *= recursive.gk();
                custom_k += "*" + recursive.gk().to_string();
            }
            gb = recursive.gb();
            custom_b = recursive._custom_b;
            n = recursive._n;
            multifactorial = recursive._multifactorial;
            custom_d = recursive._custom_d;
            c = 1;
            if (recursive.k() > 0 && recursive.k() < 1000)
            {
                algebraic_type = ALGEBRAIC_CYCLOTOMIC;
                algebraic_k = (cyclotomic == 6 ? -1 : 1)*(int32_t)recursive.k();
            }
        }
        else if (it_expr->expr->str_func == "Quad" || // (1/2 X +- 1)*X + 1
                 it_expr->expr->str_func == "Hex")    // (1/3 X +- 1)*X + 1
        {
            if (it_expr->expr->args.size() != 1)
                return InputNum::ParseResult(false, it_expr->pos, "one argument expected");
            bool neg = false;
            if (!it_expr->expr->str_args(0).empty() && it_expr->expr->str_args(0)[0] == '-')
            {
                neg = true;
                it_expr->expr->str_args(0) = it_expr->expr->str_args(0).substr(1);
            }
            int divisor = (it_expr->expr->str_func == "Quad" ? 2 : 3);
            InputNum recursive;
            InputNum::ParseResult res = recursive.parse(it_expr->expr->str_args(0) + "+1");
            if (!res)
                return InputNum::ParseResult(false, it_expr->pos + it_expr->expr->args[0].first + res.pos, res.message);
            type = recursive.type();
            gk = recursive.value();
            gk -= 1;
            if (gk%divisor != 0)
                return InputNum::ParseResult(false, it_expr->pos + it_expr->expr->args[0].first, "should be divisible by " + std::to_string(divisor));
            gk /= divisor;
            recursive._c = (neg ? -1 : 1);
            gk += recursive._c;
            factors = std::move(recursive.factors());
            cofactor = std::move(recursive.cofactor());
            b_factors = std::move(recursive.b_factors());
            b_cofactor = std::move(recursive.b_cofactor());
            ::factorize(gk, factors, cofactor);
            if (recursive.k() != 1)
            {
                gk *= recursive.gk();
                custom_k = "*" + recursive.gk().to_string();
            }
            gb = recursive.gb();
            custom_b = recursive._custom_b;
            n = recursive._n;
            multifactorial = recursive._multifactorial;
            custom_d = recursive._custom_d;
            c = 1;
            if (it_expr->expr->str_func == "Quad" && recursive.k() > 0 && recursive.k() < 100)
            {
                algebraic_type = ALGEBRAIC_QUAD;
                algebraic_k = (int32_t)(recursive.k()*recursive._c);
            }
            if (it_expr->expr->str_func == "Hex" && recursive.k() > 0 && recursive.k() < 32)
            {
                algebraic_type = ALGEBRAIC_HEX;
                algebraic_k = (int32_t)(recursive.k()*recursive._c);
            }
            if (type == KBNC)
            {
                if (recursive.gk()%divisor == 0)
                    recursive.gk() /= divisor;
                else
                {
                    recursive.gk() *= recursive.gb()/divisor;
                    recursive._n--;
                }
                custom_k = "(" + recursive.build_text() + ")" + custom_k;
            }
            else if (type == FACTORIAL || type == PRIMORIAL)
            {
                std::string st = recursive.build_text();
                custom_k = "(" + st.substr(0, st.size() - 2) + "/" + std::to_string(divisor) + st.substr(st.size() - 2) + ")" + custom_k;
            }
        }
        else
            return InputNum::ParseResult(false, it_expr->pos, "unknown function");
        it_expr++;
    }
    else
    {
        gk = 1;
        while (it_oper != tokens.end() && it_oper->oper == '*' && it_expr->expr)
        {
            if (!it_expr->expr->evaluate())
                return InputNum::ParseResult(false, it_expr->pos + it_expr->expr->error_pos, it_expr->expr->error_msg);
            if (gk == 1)
            {
                gk = std::move(it_expr->expr->value);
                custom_k = std::move(it_expr->expr->str_value);
                it_expr->expr->merge_factors(factors, cofactor, 1);
            }
            else
            {
                if (!custom_k.empty() || !it_expr->expr->str_value.empty())
                {
                    if (custom_k.empty())
                        custom_k = gk.to_string();
                    if (it_expr->expr->str_value.empty())
                        it_expr->expr->str_value = it_expr->expr->value.to_string();
                    custom_k += "*";
                    custom_k += it_expr->expr->str_value;
                }
                gk *= it_expr->expr->value;
                it_expr->expr->merge_factors(factors, cofactor, 1);
            }

            it_oper++;
            it_expr = it_oper;
            if (it_oper != tokens.end())
                it_oper++;
        }

        if (it_expr == tokens.end())
            return InputNum::ParseResult(false, (int)s.size(), "unexpected end");
        if (!it_expr->expr)
            return InputNum::ParseResult(false, it_expr->pos, "unexpected operator");
        if (it_oper != tokens.end() && it_oper->expr)
            return InputNum::ParseResult(false, it_oper->pos, "unexpected value");

        if (it_oper == tokens.end() || it_oper->oper == '^' || it_oper->oper == '+' || it_oper->oper == '-' || it_oper->oper == '/')
        {
            if (!it_expr->expr->evaluate())
                return InputNum::ParseResult(false, it_expr->pos + it_expr->expr->error_pos, it_expr->expr->error_msg);
            gb = std::move(it_expr->expr->value);
            custom_b = std::move(it_expr->expr->str_value);
            it_expr->expr->merge_factors(b_factors, b_cofactor, 1);
            Expr* b_expr = it_expr->expr.get();
            it_expr++;
            if (it_expr != tokens.end())
                it_expr++;

            if (it_oper != tokens.end() && it_oper->oper == '^')
            {
                if (it_expr == tokens.end())
                    return InputNum::ParseResult(false, (int)s.size(), "unexpected end");
                if (!it_expr->expr)
                    return InputNum::ParseResult(false, it_expr->pos, "unexpected operator");
                if (!it_expr->expr->evaluate(false))
                    return InputNum::ParseResult(false, it_expr->pos + it_expr->expr->error_pos, it_expr->expr->error_msg);
                if (it_expr->expr->value.size() > 1)
                    return InputNum::ParseResult(false, it_expr->pos, "exponent too big");
                n = it_expr->expr->value.data()[0];

                it_expr++;
                it_oper = it_expr;
                if (it_expr != tokens.end())
                    it_expr++;
            }
            else
                n = 1;
            b_expr->merge_factors(factors, cofactor, (int)n);
        }
        else if (it_oper->oper == '!')
        {
            type = FACTORIAL;
            if (!it_expr->expr->evaluate(false))
                return InputNum::ParseResult(false, it_expr->pos + it_expr->expr->error_pos, it_expr->expr->error_msg);
            if (it_expr->expr->value.size() > 1)
                return InputNum::ParseResult(false, it_expr->pos, "factorial too big");
            n = it_expr->expr->value.data()[0];

            it_oper++;
            for (multifactorial = 1; it_oper != tokens.end() && it_oper->oper == '!'; it_oper++, multifactorial++);
            it_expr = it_oper;
            if (multifactorial == 1 && it_expr != tokens.end() && it_expr->expr)
            {
                if (!it_expr->expr->evaluate(false))
                    return InputNum::ParseResult(false, it_expr->pos + it_expr->expr->error_pos, it_expr->expr->error_msg);
                if (it_expr->expr->value > n)
                    return InputNum::ParseResult(false, it_expr->pos, "multifactorial step too big");
                multifactorial = it_expr->expr->value.data()[0];

                it_expr++;
                it_oper = it_expr;
            }
            if (it_expr != tokens.end())
                it_expr++;

            Giant factor;
            factor = ((uint64_t)multifactorial << 32) + n;
            factor.arithmetic().neg(factor, factor);
            gb = factorial(factor);
            ::add_factor(factors, factor, 1);
        }
        else if (it_oper->oper == '#')
        {
            type = PRIMORIAL;
            if (!it_expr->expr->evaluate(false))
                return InputNum::ParseResult(false, it_expr->pos + it_expr->expr->error_pos, it_expr->expr->error_msg);
            if (it_expr->expr->value.bitlen() > 31)
                return InputNum::ParseResult(false, it_expr->pos, "primorial too big");

            Giant factor;
            factor.arithmetic().neg(it_expr->expr->value, factor);
            gb = primorial(factor);
            n = factor.data()[0];
            ::add_factor(factors, factor, 1);

            it_oper++;
            it_expr = it_oper;
            if (it_expr != tokens.end())
                it_expr++;
        }
        else
            return InputNum::ParseResult(false, it_oper->pos, "unexpected symbol");

        if (it_oper != tokens.end())
        {
            if (it_oper->expr)
                return InputNum::ParseResult(false, it_oper->pos, "unexpected value");
            if (it_expr == tokens.end())
                return InputNum::ParseResult(false, (int)s.size(), "unexpected end");
            if (!it_expr->expr)
                return InputNum::ParseResult(false, it_expr->pos, "unexpected operator");
        }

        gd = 1;
        while (it_expr != tokens.end() && it_expr->expr && it_oper->oper == '/')
        {
            if (!it_expr->expr->evaluate())
                return InputNum::ParseResult(false, it_expr->pos + it_expr->expr->error_pos, it_expr->expr->error_msg);
            if (gd == 1)
            {
                gd = std::move(it_expr->expr->value);
                custom_d = std::move(it_expr->expr->str_value);
                it_expr->expr->merge_factors(factors, cofactor, -1);
            }
            else
            {
                if (!custom_d.empty() || !it_expr->expr->str_value.empty())
                {
                    if (custom_d.empty())
                        custom_d = gd.to_string();
                    if (it_expr->expr->str_value.empty())
                        it_expr->expr->str_value = it_expr->expr->value.to_string();
                    custom_d += "/";
                    custom_d += it_expr->expr->str_value;
                }
                gd *= it_expr->expr->value;
                it_expr->expr->merge_factors(factors, cofactor, -1);
            }

            it_expr++;
            it_oper = it_expr;
            if (it_expr != tokens.end())
                it_expr++;
        }
        if (gd > 1 && gk*(type == KBNC ? power(gb, n) : gb)%gd != 0)
            return InputNum::ParseResult(false, it_oper != tokens.end() ? it_oper->pos : (int)s.size(), "not divisible");

        if (it_oper != tokens.end())
        {
            if (it_oper->expr)
                return InputNum::ParseResult(false, it_oper->pos, "unexpected value");
            if (it_expr == tokens.end())
                return InputNum::ParseResult(false, (int)s.size(), "unexpected end");
            if (!it_expr->expr)
                return InputNum::ParseResult(false, it_expr->pos, "unexpected operator");
        }

        if (it_expr != tokens.end() && (it_oper->oper == '+' || it_oper->oper == '-'))
        {
            bool minus = (it_oper->oper == '-');
            if (!it_expr->expr->evaluate(false))
                return InputNum::ParseResult(false, it_expr->pos + it_expr->expr->error_pos, it_expr->expr->error_msg);
            if (it_expr->expr->value.bitlen() > 63)
                return InputNum::ParseResult(false, it_expr->pos, "C too big");
            if (it_expr->expr->value.size() == 1)
                c = it_expr->expr->value.data()[0];
            if (it_expr->expr->value.size() == 2)
                c = *(int64_t*)it_expr->expr->value.data();
            if (minus)
                c = -c;
            it_expr++;
            it_oper = it_expr;
        }
        else if (it_oper == tokens.end() && type == KBNC && gk == 1 && n == 1 && gd == 1)
        {
            type = GENERIC;
            n = 0;
        }
        else if (c_required)
            return InputNum::ParseResult(false, it_oper != tokens.end() ? it_oper->pos : (int)s.size(), "C expected");
    }
    if (it_oper != tokens.end())
        return InputNum::ParseResult(false, it_oper->pos, it_oper->expr ? "unexpected value" : "unexpected operator");
    if (it_expr != tokens.end())
        return InputNum::ParseResult(false, it_expr->pos, it_expr->expr ? "unexpected value" : "unexpected operator");

    _type = type;
    _gk = std::move(gk);
    _gb = std::move(gb);
    _n = n;
    _c = c;
    _gf = 1;
    _custom_k = std::move(custom_k);
    _custom_b = std::move(custom_b);
    _custom_d = std::move(custom_d);
    _custom_f.clear();
    _multifactorial = multifactorial;
    _algebraic_type = algebraic_type;
    _algebraic_k = algebraic_k;
    _factors = std::move(factors);
    _cofactor = std::move(cofactor);
    _b_factors = std::move(b_factors);
    _b_cofactor = std::move(b_cofactor);
    std::sort(_factors.begin(), _factors.end(), [](std::pair<arithmetic::Giant, int>& a, std::pair<arithmetic::Giant, int>& b) { return (a.second > 0 && b.second < 0) || a.first < b.first; });

    process(gd, true);
    return InputNum::ParseResult(true);
}

void InputNum::process()
{
    arithmetic::Giant gd;
    gd = 1;
    process(gd, false);
}

void InputNum::process(arithmetic::Giant& gd, bool factored)
{
    if (!factored)
    {
        _factors.clear();
        if (!_cofactor.empty())
            _cofactor.arithmetic().free(_cofactor);
        _b_factors.clear();
        if (!_b_cofactor.empty())
            _b_cofactor.arithmetic().free(_b_cofactor);

        if (_type != GENERIC)
        {
            factorize(_gk, _factors, _cofactor);
            if (_type == KBNC)
            {
                factorize(_gb, _b_factors, _b_cofactor);
                for (auto& factor : _b_factors)
                    ::add_factor(_factors, factor.first, factor.second*_n);
            }
            if (_type == FACTORIAL)
            {
                Giant factorial;
                factorial = ((uint64_t)_multifactorial << 32) + _n;
                factorial.arithmetic().neg(factorial, factorial);
                ::add_factor(_factors, factorial, 1);
            }
            if (_type == PRIMORIAL)
            {
                Giant primorial;
                primorial = _n;
                primorial.arithmetic().neg(primorial, primorial);
                ::add_factor(_factors, primorial, 1);
            }

            std::sort(_factors.begin(), _factors.end(), [](std::pair<arithmetic::Giant, int>& a, std::pair<arithmetic::Giant, int>& b) { return a.first < b.first; });
            if (!_b_cofactor.empty() && !_cofactor.empty())
                _cofactor *= power(_b_cofactor, _n);
            else if (!_b_cofactor.empty())
                _cofactor = power(_b_cofactor, _n);
            if (gd > 1)
            {
                std::vector<std::pair<arithmetic::Giant, int>> factors;
                arithmetic::Giant cofactor;
                factorize(gd, factors, cofactor);
                for (auto& factor : factors)
                    ::add_factor(_factors, factor.first, -factor.second);
                if (!cofactor.empty())
                {
                    if (_cofactor.empty() || _cofactor%cofactor != 0)
                        throw ArithmeticException();
                    _cofactor /= cofactor;
                }
            }
        }
    }

    if (_type == KBNC && _c == 1 && gd == 1 && _algebraic_type == ALGEBRAIC_SIMPLE)
    {
        if (_gk == 1 && _n > 1 && (_n & (_n - 1)) == 0)
            for (_gfn = 1; (1UL << _gfn) < _n; _gfn++);
        if ((abs(_gk.bitlen() - log2(_gb)*_n) < _n/30 + 10))
        {
            Giant tmp = _gb;
            tmp.power(_n);
            Giant cyclo_tmp = tmp - _gk;
            if (cyclo_tmp == 1)
            {
                _algebraic_type = ALGEBRAIC_CYCLOTOMIC;
                _algebraic_k = -1;
            }
            else if (cyclo_tmp == -1)
            {
                _algebraic_type = ALGEBRAIC_CYCLOTOMIC;
                _algebraic_k = 1;
            }
            else if (_gb%2 == 0)
            {
                tmp /= 2;
                tmp -= _gk;
                if (tmp == 1)
                {
                    _algebraic_type = ALGEBRAIC_QUAD;
                    _algebraic_k = -1;
                }
                if (tmp == -1)
                {
                    _algebraic_type = ALGEBRAIC_QUAD;
                    _algebraic_k = 1;
                }
            }
            else if (_gb%3 == 0)
            {
                tmp /= 3;
                tmp -= _gk;
                if (tmp == 1)
                {
                    _algebraic_type = ALGEBRAIC_HEX;
                    _algebraic_k = -1;
                }
                if (tmp == -1)
                {
                    _algebraic_type = ALGEBRAIC_HEX;
                    _algebraic_k = 1;
                }
            }
        }
    }

    if (gd > 1 && _custom_d.empty())
        _custom_d = gd.to_string();
    _input_text = build_text();

    if (_type == KBNC && _gb != 1)
    {
        if (_b_cofactor.empty())
        {
            int a = _b_factors[0].second;
            for (auto it = ++_b_factors.begin(); it != _b_factors.end(); it++)
                a = gcd(a, it->second);
            if (a > 1)
            {
                for (auto it = _b_factors.begin(); it != _b_factors.end(); it++)
                    it->second /= a;
                _n *= a;
                _gb = 1;
                for (auto it = _b_factors.begin(); it != _b_factors.end(); it++)
                    _gb *= power(it->first, it->second);
                _custom_b.clear();
            }
        }

        while (_gk%_gb == 0 && _algebraic_type == ALGEBRAIC_SIMPLE)
        {
            _gk /= _gb;
            _n++;
            _custom_k.clear();
        }
    }

    if (gd > 1)
    {
        Giant k_d = gcd(_gk, gd);
        if (k_d > 1)
        {
            _gk /= k_d;
            gd /= k_d;
        }
        uint32_t n = 0;
        if (gd > 1 && _type == KBNC)
        {
            Giant tmp = gcd(gd, _gb);
            if (tmp > 1)
                n = (uint32_t)(log2(gd)/log2(tmp)) + 10;
            if (n > _n - 1)
                n = _n - 1;
            if (n > 0)
            {
                tmp = power(_gb, n);
                if (tmp%gd == 0)
                {
                    tmp /= gd;
                    while (tmp%_gb == 0)
                    {
                        tmp /= _gb;
                        n--;
                    }
                    _gk *= tmp;
                    _n -= n;
                    gd = 1;
                    if (_gk.size() <= 2)
                    {
                        _custom_k.clear();
                        _custom_d.clear();
                    }
                    else
                    {
                        if (!_custom_k.empty() && (tmp != 1 || k_d != 1))
                        {
                            if (tmp != 1)
                                _custom_k += "*(" + (!_custom_b.empty() ? _custom_b : _gb.to_string()) + (n > 1 ? "^" + std::to_string(n) : "") + ")";
                            else if (k_d != 1)
                                _custom_d = k_d.to_string();
                        }
                        else
                            _custom_d.clear();
                    }
                }
            }
        }
        if (gd > 1)
        {
            GWASSERT(_gb%gd == 0);
            _gb /= gd;
        }
        if (n == 0)
        {
            if (_custom_k.empty() && _custom_b.empty())
                _custom_d.clear();
            if (!_custom_k.empty() && _custom_b.empty())
            {
                if (k_d == 1)
                    _custom_d.clear();
                if (k_d != 1 && gd != 1)
                    _custom_d = k_d.to_string();
            }
            if (_custom_k.empty() && !_custom_b.empty())
            {
                if (gd == 1)
                    _custom_d.clear();
                if (k_d != 1 && gd != 1)
                    _custom_d = gd.to_string();
            }
        }
    }

    _display_text = build_text(30);
}

std::string InputNum::build_text(int max_len)
{
    std::string res;
    res.reserve(32);
    if (_gf != 1)
        res.append(1, '(');
    if (_gk != 1 || !_custom_k.empty())
    {
        std::string sk = !_custom_k.empty() ? _custom_k : _gk.to_string();
        if (sk.size() > max_len && max_len > 0)
        {
            res.append(sk, 0, max_len/2);
            res.append(3, '.');
            res.append(sk, sk.size() - max_len/2, max_len/2);
        }
        else
            res.append(sk);
        if (_gb != 1)
            res.append(1, '*');
    }
    if (_type == FACTORIAL)
    {
        res.append(std::to_string(_n));
        if (_multifactorial <= 3)
            res.append(_multifactorial, '!');
        else
            res.append(1, '!').append(std::to_string(_multifactorial));
    }
    else if (_type == PRIMORIAL)
    {
        res.append(std::to_string(_n));
        res.append(1, '#');
    }
    else
    {
        std::string sb = !_custom_b.empty() ? _custom_b : _gb.to_string();
        if (sb.size() > max_len && max_len > 0)
        {
            res.append(sb, 0, max_len/2);
            res.append(3, '.');
            res.append(sb, sb.size() - max_len/2, max_len/2);
        }
        else
            res.append(sb);
        if (_type == KBNC && _n != 1)
        {
            res.append(1, '^');
            res.append(std::to_string(_n));
        }
    }
    if (!_custom_d.empty())
    {
        res.append(1, '/');
        std::string sd = _custom_d;
        if (sd.size() > max_len && max_len > 0)
        {
            res.append(sd, 0, max_len/2);
            res.append(3, '.');
            res.append(sd, sd.size() - max_len/2, max_len/2);
        }
        else
            res.append(sd);
    }
    if (_c != 0)
    {
        if (_c > 0)
            res.append(1, '+');
        if (_c < 0)
            res.append(1, '-');
        res.append(std::to_string(std::abs(_c)));
    }
    if (_gf != 1)
    {
        res.append(1, ')');
        res.append(1, '/');
        std::string sf = !_custom_f.empty() ? _custom_f : _gf.to_string();
        if (sf.size() > max_len && max_len > 0)
        {
            res.append(sf, 0, max_len/2);
            res.append(3, '.');
            res.append(sf, sf.size() - max_len/2, max_len/2);
        }
        else
            res.append(sf);
    }

    return res;
}

void InputNum::setup(GWState& state)
{
    if (_type == GENERIC)
    {
        state.setup(_gb);
    }
    else if (_type != KBNC)
    {
        state.setup((_gk*_gb + _c)/_gf);
    }
    else if (_algebraic_type != ALGEBRAIC_SIMPLE && b() != 0 && !state.force_mod_type)
    {
        if (_algebraic_type == ALGEBRAIC_CYCLOTOMIC)
        {
            uint32_t k = (uint32_t)abs(_algebraic_k);
            state.known_factors = (k*power(_gb, _n) + (_algebraic_k < 0 ? 1 : -1))*_gf;
            state.setup(k*k*k, b(), 3*_n, _algebraic_k < 0 ? 1 : -1);
        }
        if (_algebraic_type == ALGEBRAIC_QUAD)
        {
            Giant x = (-_algebraic_k)*power(_gb, _n);
            x <<= 1;
            state.known_factors = (value() + x)*_gf;
            uint64_t k_squared = (uint64_t)(_algebraic_k*_algebraic_k);
            if (k_squared%2 == 0)
                state.setup(k_squared*k_squared/4, b(), 4*_n, 1);
            else if (b()%4 == 0)
                state.setup(k_squared*k_squared*b()/4, b(), 4*_n - 1, 1);
            else
                state.setup(k_squared*k_squared*b()*b()/4, b(), 4*_n - 2, 1);
        }
        if (_algebraic_type == ALGEBRAIC_HEX)
        {
            Giant x = (-_algebraic_k)*power(_gb, _n);
            Giant val = value() + x;
            state.known_factors = (val + std::move(x))*val*_gf;
            uint64_t k_squared = (uint64_t)(_algebraic_k*_algebraic_k);
            if (k_squared%3 == 0)
                state.setup(k_squared*k_squared*k_squared/27, b(), 6*_n, 1);
            else
                state.setup(k_squared*k_squared*k_squared*b()*b()*b()/27, b(), 6*_n - 3, 1);
        }
        if (*state.N%3417905339UL != fingerprint())
            throw ArithmeticException();
        if (state.gwdata()->GENERAL_MOD || state.gwdata()->GENERAL_MMGW_MOD)
        {
            state.done();
            state.known_factors = 1;
            _algebraic_type = ALGEBRAIC_SIMPLE;
            _algebraic_k = 0;
            setup(state);
        }
        return;
    }
    else if (k() != 0 && b() != 0 && std::abs(_c) < (1ULL << 30)  && (_gf == 1 || !state.force_mod_type))
    {
        state.known_factors = _gf;
        state.setup(k(), b(), _n, _c);
        if (_gf != 1)
        {
            if (state.gwdata()->GENERAL_MOD || state.gwdata()->GENERAL_MMGW_MOD)
            {
                state.done();
                state.known_factors = 1;
                state.setup((_gk*power(_gb, _n) + _c)/_gf);
            }
            else
            {
                if (*state.N%3417905339UL != fingerprint())
                    throw ArithmeticException();
                return;
            }
        }
    }
    else
    {
        state.setup((_gk*power(_gb, _n) + _c)/_gf);
    }
    if (state.fingerprint != fingerprint())
        throw ArithmeticException();
}

int InputNum::bitlen()
{
    if (_type == GENERIC)
        return _gb.bitlen();
    int res;
    if (_type != KBNC)
        res = _gk.bitlen() + _gb.bitlen() - _gf.bitlen();
    else if (b() == 2)
        res = _gk.bitlen() + _n + (_c > 0 || _gk != 1 ? 1 : 0) - _gf.bitlen();
    else
        res = (int)std::ceil(log2(_gk) + log2(_gb)*_n - log2(_gf));
    if (res <= 64)
        res = value().bitlen();
    return res;
}

uint64_t InputNum::parse_numeral(const std::string& s)
{
    size_t pos = 0;
    uint64_t val = std::stoull(s, &pos);
    for (auto it = s.begin() + pos; it != s.end(); it++)
        switch (*it)
        {
        case 'k':
        case 'K':
            val *= 1000;
            break;
        case 'M':
            val *= 1000000;
            break;
        case 'G':
            val *= 1000000000;
            break;
        case 'T':
            val *= 1000000000000ULL;
            break;
        case 'P':
            val *= 1000000000000000ULL;
            break;
        }
    return val;
}

uint32_t InputNum::mod(uint32_t modulus)
{
    if (_type == GENERIC)
        return _gb%modulus;

    uint32_t f_inv = 1;
    if (_gf != 1)
    {
        f_inv = inv(_gf%modulus, modulus);
        if (!f_inv)
            return value()%modulus;
    }

    uint64_t result = 1;
    if (_type != KBNC)
    {
        result = _gb%modulus;
        if (_gk > 1)
            result *= _gk%modulus;
    }
    else
    {
        if (!_cofactor.empty())
            result = _cofactor%modulus;
        for (auto& factor : _factors)
        {
            uint32_t b;
            if (factor.first < modulus)
                b = factor.first.data()[0];
            else
                b = factor.first%modulus;
            uint64_t mult = 1;
            for (uint32_t i = 1; i <= (uint32_t)factor.second; i <<= 1)
            {
                if (i == 1)
                    mult = b;
                else
                {
                    mult *= mult;
                    if (mult >= (1ULL << 32))
                        mult %= modulus;
                }
                if ((factor.second & i) != 0)
                {
                    result *= mult;
                    if (result >= (1ULL << 32))
                        result %= modulus;
                }
            }
        }
    }

    if (_c >= 0)
        result += _c;
    else
    {
        uint32_t c = (uint32_t)((-_c)%modulus);
        if (result < c)
            result += (uint64_t)modulus << 31;
        result -= c;
    }

    if (f_inv != 1)
    {
        result %= modulus;
        result *= f_inv;
    }
    return (uint32_t)(result%modulus);
}

bool InputNum::is_half_factored()
{
    if (std::abs(c()) != 1)
        return false;
    if (_cofactor.empty() || _cofactor.bitlen()*2 + 10 < bitlen())
        return true;
    if (_cofactor.bitlen()*2 > bitlen() + 10)
        return false;
    Giant tmp;
    tmp = 1;
    for (auto& factor : _factors)
        tmp *= power(factor.first, factor.second);
    if (_c == 1)
        return tmp > _cofactor;
    tmp -= 1;
    return tmp > _cofactor;
}

std::vector<int> InputNum::factorize_minus1(int depth)
{
    uint32_t s = (depth + 1)/2;
    if (s%2 == 1)
        s++;
    Giant minus1 = value() - 1;
    std::vector<std::pair<arithmetic::Giant, int>> factors;
    Giant cofactor;
    factorize(minus1, factors, cofactor, [&](Giant& x, uint32_t p) { return _gf%p != 0 ? mod(p) == 1 : x%p == 0; }, s);

    std::vector<int> divisors;
    for (auto& factor : factors)
    {
        if (factor.first > INT_MAX)
            continue;
        int base = (int)(factor.first.data()[0]);
        int divisor = base;
        for (int i = 1; i < factor.second && divisor < INT_MAX/base; i++, divisor *= base);
        divisors.push_back(divisor);
    }

    return divisors;
}

void InputNum::expand_factors()
{
    bool f_p_exist = false;
    for (auto& f : _factors)
        if (f.first < 0)
            f_p_exist = true;
    if (!f_p_exist)
        return;

    std::map<uint32_t, int> factors;
    for (auto& f : _factors)
        if (f.first > 0 && f.first.size() == 1)
            factors[*f.first.data()] = f.second;

    for (auto& f : _factors)
        if (f.first < 0)
        {
            if (f.first.size() == 2)
            {
                uint32_t n = f.first.data()[0];
                uint32_t multifactorial = f.first.data()[1];
                uint32_t i = n%multifactorial;
                while (i < 2)
                    i += multifactorial;
                for (; i <= n; i += multifactorial)
                {
                    uint32_t j = i;
                    for (auto it = PrimeIterator::get(); (*it)*(*it) <= (int)j; it++)
                        if (j%(*it) == 0)
                        {
                            int& power = factors[*it];
                            do
                            {
                                power += f.second;
                                j /= *it;
                            } while (j%(*it) == 0);
                        }
                    if (j != 1)
                        factors[j] += f.second;
                }
            }
            if (f.first.size() == 1)
            {
                uint32_t n = f.first.data()[0];
                for (auto it = PrimeIterator::get(); *it <= (int)n; it++)
                    factors[*it] += f.second;
            }
        }

    Giant tmp;
    std::vector<std::pair<Giant, int>> new_factors;
    for (auto& f : factors)
        if (f.second != 0)
        {
            tmp = f.first;
            new_factors.emplace_back(std::move(tmp), f.second);
        }
    for (auto& f : _factors)
        if (f.first > 0 && f.first.size() > 1)
            new_factors.emplace_back(std::move(f.first), f.second);
    _factors = std::move(new_factors);
}

std::vector<int> InputNum::factorize_small()
{
    bool trial_test = (bitlen() <= 40);
    Giant N;
    if (trial_test)
        N = value();
    else
        N = ((std::move(N) = 2) << 61) - 1;

    std::vector<std::pair<arithmetic::Giant, int>> factors;
    arithmetic::Giant cofactor;
    if (trial_test)
        factorize(N, factors, cofactor);
    else
    {
        if (mod(2) == 0)
            ::add_factor(factors, 2, 1);
        factorize(N, factors, cofactor, [&](Giant& x, uint32_t p) { bool res = (mod(p) == 0); if (res) x *= p; return res; });
    }

    std::vector<int> res;
    for (auto& factor : factors)
    {
        if (factor.first.bitlen() < 32)
            res.push_back((int)factor.first.data()[0]);
        else
            res.push_back(0);
    }
    return res;
}

void InputNum::print_info()
{
    if (_type == ZERO)
        return;

    Giant N = value();
    Giant tmp;
    std::string st;

    std::cout << display_text() << ", " << N.bitlen() << " bits, " << N.digits() << " digits." << std::endl;

    bool mod_available = true;
    for (auto& factor : _factors)
        if (factor.first < 0)
            mod_available = false;

    std::vector<std::pair<arithmetic::Giant, int>> factors;
    arithmetic::Giant cofactor;
    factorize(N, factors, cofactor, [&](Giant& x, uint32_t p) { return mod_available && _gf%p != 0 ? mod(p) == 0 : x%p == 0; });
    if (factors.empty())
        std::cout << "No small factors." << std::endl;
    else if (factors[0].first == N)
        std::cout << N.to_string() << " is prime." << std::endl;
    else
    {
        std::cout << "Small factors:";
        for (auto& factor : factors)
        {
            std::cout << " " << factor.first.to_string();
            if (factor.second > 1)
                std::cout << "^" << factor.second;
        }
        std::cout << std::endl;
    }

    if (type() == KBNC && std::abs(_c) == 1 && (k() == 1 || _cofactor.empty() || (!_b_cofactor.empty() && _cofactor == power(_b_cofactor, _n))))
    {
        uint32_t n = _n;
        if (k() != 1)
        {
            n = gcd(_n, (uint32_t)_factors[0].second);
            for (auto it = ++_factors.begin(); it != _factors.end(); it++)
                n = gcd(n, (uint32_t)it->second);
        }

        uint32_t l = 0;
        if (n > 1 && _c == 1)
        {
            for (l = 1; !(n & l); l <<= 1);
            if (l == n)
                l = 0;
        }
        if (n > 1 && _c == -1 && (k() != 1 || b() != 2))
            l = 1;
        if (n > 1 && _c == -1 && k() == 1 && b() == 2)
        {
            PrimeIterator it = PrimeIterator::get();
            int sqrt_n = (int)std::sqrt(n);
            for (; n%(*it) != 0 && *it < sqrt_n; it++);
            if (n%(*it) == 0)
                l = *it;
        }

        if (l != 0)
        {
            st = "";
            for (auto& factor : _factors)
            {
                if (!st.empty())
                    st += "*";
                st += factor.first.to_string();
                if (factor.second/n*l > 1)
                {
                    st += "^";
                    st += std::to_string(factor.second/n*l);
                }
            }
            if (!_b_cofactor.empty())
            {
                if (!st.empty())
                    st += "*";
                st += _b_cofactor.to_string();
                if (_n/n*l > 1)
                {
                    st += "^";
                    st += std::to_string(_n/n*l);
                }
            }
            if (_c > 0)
                st += "+1";
            else
                st += "-1";
            std::cout << "Algebraic factor: " << st << std::endl;
        }
        else if (_c == 1 && _cofactor.empty() && _factors[0].first == 2 && _factors[0].second%4 == 2)
        {
            l = 4;
            for (auto it = ++_factors.begin(); it != _factors.end(); it++)
                if (it->second%4 != 0)
                    l = 0;
            if (l != 0)
            {
                st = "2" + (_factors[0].second/2 > 1 ? "^" + std::to_string(_factors[0].second/2) : "");
                for (auto it = ++_factors.begin(); it != _factors.end(); it++)
                    st += "*" + it->first.to_string() +  "^" + std::to_string(it->second/2);
                st += "+2" + (_factors[0].second/4 + 1 > 1 ? "^" + std::to_string(_factors[0].second/4 + 1) : "");
                for (auto it = ++_factors.begin(); it != _factors.end(); it++)
                    st += "*" + it->first.to_string() + (it->second/4 > 1 ? "^" + std::to_string(it->second/4) : "");
                std::cout << "Algebraic factor: " << st << "+1" << std::endl;
            }
        }
        else if (_c == 1 && _cofactor.empty() && ((_factors[0].first == 3 && _factors[0].second%6 == 3) || (_factors.size() > 1 && _factors[1].first == 3 && _factors[1].second%6 == 3)))
        {
            l = 6;
            for (auto it = _factors.begin(); it != _factors.end(); it++)
                if (it->first != 3 && it->second%6 != 0)
                    l = 0;
            if (l != 0)
            {
                st = "";
                for (auto it = ++_factors.begin(); it != _factors.end(); it++)
                    st += (!st.empty() ? "*" : "") + it->first.to_string() + (it->second/3 > 1 ? "^" + std::to_string(it->second/3) : "");
                std::cout << "Algebraic factor: " << st << "+1" << std::endl;
            }
        }
    }
    if (type() == KBNC && _algebraic_type == ALGEBRAIC_CYCLOTOMIC && abs(_algebraic_k) == 1)
    {
        uint32_t d = _n;
        while (d%3 == 0)
            d /= 3;
        if (_algebraic_k < 0)
            while (!(d & 1))
                d >>= 1;

        if (d > 1)
        {
            d = _n/d;
            std::string sb = !_custom_b.empty() ? _custom_b : _gb.to_string();
            if (sb.size() > 30)
            {
                st = "";
                st.append(sb, 0, 30/2);
                st.append(3, '.');
                st.append(sb, sb.size() - 30/2, 30/2);
                sb = st;
            }
            else
                st = sb;
            st += "^" + std::to_string(2*d) + (_algebraic_k > 0 ? "+" : "-") + sb + (d > 1 ? "^" + std::to_string(d) : "") + "+1";
            std::cout << "Algebraic factor: " << st << std::endl;
        }
    }

    factors.clear();
    if (!cofactor.empty())
        cofactor.arithmetic().free(cofactor);
    if (c() == 1)
    {
        factors = _factors;
        cofactor = _cofactor;
    }
    else
    {
        N -= 1;
        factorize(N, factors, cofactor, [&](Giant& x, uint32_t p) { return mod_available && _gf%p != 0 ? mod(p) == 1 : x%p == 0; });
        N += 1;
    }
    if (factors.empty())
        std::cout << "No small factors of N-1." << std::endl;
    else
    {
        tmp = 1;
        st.clear();
        for (auto& factor : factors)
        {
            if (!st.empty())
                st += " ";
            if (factor.second < 0)
                st += "/";
            Giant value;
            if (factor.first < 0 && factor.first.size() == 2)
            {
                st += std::to_string(factor.first.data()[0]);
                int multifactorial = factor.first.data()[1];
                if (multifactorial <= 3)
                    st.append(multifactorial, '!');
                else
                    st.append(1, '!').append(std::to_string(multifactorial));
                if (!cofactor.empty())
                    value = factorial(factor.first);
            }
            else if (factor.first < 0 && factor.first.size() == 1)
            {
                st += std::to_string(factor.first.data()[0]) + "#";
                if (!cofactor.empty())
                    value = primorial(factor.first);
            }
            else
            {
                st += factor.first.to_string();
                if (!cofactor.empty())
                    value = factor.first;
            }
            if (!cofactor.empty())
            {
                value.power(abs(factor.second));
                if (factor.second > 0)
                    tmp *= value;
                else
                    tmp /= value;
            }
            if (abs(factor.second) != 1)
            {
                st += "^";
                st += std::to_string(abs(factor.second));
            }
        }
        std::cout << "Factors of N-1";
        if (cofactor.empty())
            std::cout << " (fully factored)";
        else if (square(tmp) > N)
            std::cout << " (more than a half (" << std::fixed << std::setprecision(1) << (log2(tmp)/log2(N)*100) << "%) is factored)";
        else if (2*square(tmp)*tmp > N)
            std::cout << " (more than a third (" << std::fixed << std::setprecision(1) << (log2(tmp)/log2(N)*100) << "%) is factored)";
        else
            std::cout << " (" << std::fixed << std::setprecision(3) << (log2(tmp)/log2(N)*100) << "% is factored)";
        std::cout << ": " << st << std::endl;
    }

    factors.clear();
    if (!cofactor.empty())
        cofactor.arithmetic().free(cofactor);
    if (c() == -1)
    {
        factors = _factors;
        cofactor = _cofactor;
    }
    else
    {
        N += 1;
        factorize(N, factors, cofactor, [&](Giant& x, uint32_t p) { return mod_available && _gf%p != 0 ? mod(p) == p - 1 : x%p == 0; });
        N -= 1;
    }
    if (factors.empty())
        std::cout << "No small factors of N+1." << std::endl;
    else
    {
        tmp = 1;
        st.clear();
        for (auto& factor : factors)
        {
            if (!st.empty())
                st += " ";
            if (factor.second < 0)
                st += "/";
            Giant value;
            if (factor.first < 0 && factor.first.size() == 2)
            {
                st += std::to_string(factor.first.data()[0]);
                int multifactorial = factor.first.data()[1];
                if (multifactorial <= 3)
                    st.append(multifactorial, '!');
                else
                    st.append(1, '!').append(std::to_string(multifactorial));
                if (!cofactor.empty())
                    value = factorial(factor.first);
            }
            else if (factor.first < 0 && factor.first.size() == 1)
            {
                st += std::to_string(factor.first.data()[0]) + "#";
                if (!cofactor.empty())
                    value = primorial(factor.first);
            }
            else
            {
                st += factor.first.to_string();
                if (!cofactor.empty())
                    value = factor.first;
            }
            if (!cofactor.empty())
            {
                value.power(abs(factor.second));
                if (factor.second > 0)
                    tmp *= value;
                else
                    tmp /= value;
            }
            if (abs(factor.second) != 1)
            {
                st += "^";
                st += std::to_string(abs(factor.second));
            }
        }
        std::cout << "Factors of N+1";
        if (cofactor.empty())
            std::cout << " (fully factored)";
        else if (square(tmp - 1) > N)
            std::cout << " (more than a half (" << std::fixed << std::setprecision(1) << (log2(tmp)/log2(N)*100) << "%) is factored)";
        else if (2*square(tmp)*(tmp - 1) > N)
            std::cout << " (more than a third (" << std::fixed << std::setprecision(1) << (log2(tmp)/log2(N)*100) << "%) is factored)";
        else
            std::cout << " (" << std::fixed << std::setprecision(3) << (log2(tmp)/log2(N)*100) << "% is factored)";
        std::cout << ": " << st << std::endl;
    }

    std::cout << "(2/N)=" << kronecker(2, N) << ", (3/N)=" << kronecker(3, N) << ", (5/N)=" << kronecker(5, N) << ", (7/N)=" << kronecker(7, N) << ", (11/N)=" << kronecker(11, N) << "." << std::endl;
}
