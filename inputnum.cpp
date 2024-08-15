
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

template<class It>
void iterate_digits(bool& prime, It& it, It& first, const It& last)
{
    prime = it != last && *it == 'p';
    if (prime)
        it++;
    for (first = it; it != last && std::isdigit(*it); it++);
}

template<class It>
void iterate_func(std::string& f, It& it, const It& last)
{
    f.clear();
    It it_s;
    for (it_s = it; it != last && std::isalpha(*it); it++);
    if (it != last && *it == '(')
        f = std::string(it_s, it);
    else
        it = it_s;
}

template<class It>
void iterate_args(std::string& f, std::vector<std::string>& args, It& it, const It& last)
{
    args.clear();
    if (f.empty())
        f = "()";
    it++;
    while (it != last && *it != ')')
    {
        int recursion = 0;
        It it_s = it;
        for (it_s = it; it != last && (recursion > 0 || (*it != ')' && *it != ',')); it++)
        {
            if (*it == '(')
                recursion++;
            if (*it == ')')
                recursion--;
        }
        args.emplace_back(it_s, it);
        if (*it == ',')
            it++;
    }
    if (it != last && *it == ')')
        it++;
}

int get_prime(int n)
{
    if (n == 0)
        return 1;
    PrimeIterator primes = PrimeIterator::get();
    for (int i = 1; i < n; i++, primes++);
    return *primes;
}

template<class It, class T, typename std::enable_if_t<std::is_integral<T>::value, bool> = true>
bool parse_digits(bool prime, const It& first, const It& last, T& value)
{
    if (last == first || last - first > (prime ? 8 : 9))
        return false;
    if (prime)
        value = get_prime(stoi(std::string(first, last)));
    else
        value = stoi(std::string(first, last));
    if (value == 0)
        return false;
    return true;
}

template<class It, class T, typename std::enable_if_t<!std::is_integral<T>::value, bool> = true>
bool parse_digits(bool prime, It& first, It& last, T& value)
{
    if (last == first || (prime && last - first > 8))
        return false;
    if (prime)
        value = get_prime(stoi(std::string(first, last)));
    else
        value = std::string(first, last);
    if (value == 0)
        return false;
    return true;
}

bool parse_func(const std::string& f, const std::vector<std::string>& args, Giant& value, std::string& custom_val)
{
    if (f == "()")
    {
        if (args.size() > 1)
            return false;
        InputNum recursive;
        if (!recursive.parse(args[0], false))
            return false;
        value = recursive.value();
        custom_val = "(" + recursive.input_text() + ")";
    }
    else
    {
        custom_val = f;
        for (auto it_a = args.begin(); it_a != args.end(); it_a++)
            custom_val += (it_a == args.begin() ? "(" : ",") + *it_a;
        custom_val += ")";
        return false; // not implemented
    }
    return true;
}

bool InputNum::parse(const std::string& s, bool c_required)
{
    std::string::const_iterator it, it_s;
    int type = KBNC;
    Giant gk, gb, gd;
    uint32_t n = 0;
    int32_t c = 0;
    Giant tmp;
    bool prime;
    std::string f;
    std::vector<std::string> args;
    std::string custom_k;
    std::string custom_b;
    std::string custom_d;
    int multifactorial = 0;
    int algebraic_type = ALGEBRAIC_SIMPLE;
    int32_t algebraic_k = 0;

    for (it = s.begin(); it != s.end() && std::isspace(*it); it++);
    if (it != s.end() && *it == '\"')
        it++;

    iterate_func(f, it, s.end());
    if (!f.empty() || (it != s.end() && *it == '('))
        iterate_args(f, args, it, s.end());
    else
        iterate_digits(prime, it, it_s, s.end());

    if (!f.empty() && (it == s.end() || *it == '\"'))
    {
        if (f == "Phi" && args.size() == 2)
        {
            int cyclotomic = stoi(args[0]);
            if (cyclotomic != 3 && cyclotomic != 6)
                return false;
            if (!args[1].empty() && args[1][0] == '-')
            {
                cyclotomic = (cyclotomic == 3 ? 6 : 3);
                args[1] = args[1].substr(1);
            }
            InputNum recursive;
            if (!recursive.parse(args[1] + "+1"))
                return false;
            type = recursive.type();
            gk = recursive.value();
            custom_k = "(" + recursive.input_text() + ")";
            if (cyclotomic == 6)
            {
                gk -= 2;
                custom_k[custom_k.size() - 3] = '-';
            }
            if (recursive.k() != 1)
            {
                gk *= recursive.gk();
                custom_k += "*" + recursive.gk().to_string();
            }
            gb = recursive.gb();
            custom_b = recursive._custom_b;
            n = recursive._n;
            multifactorial = recursive._multifactorial;
            gd = recursive.gd();
            custom_d = recursive._custom_d;
            c = 1;
            if (recursive.k() > 0 && recursive.k() < 1000 && recursive.d() == 1)
            {
                algebraic_type = ALGEBRAIC_CYCLOTOMIC;
                algebraic_k = (cyclotomic == 6 ? -1 : 1)*(int32_t)recursive.k();
            }
        }
        else if (f == "Quad" && args.size() == 1)
        {
            bool neg = false;
            if (!args[0].empty() && args[0][0] == '-')
            {
                neg = true;
                args[0] = args[0].substr(1);
            }
            InputNum recursive;
            if (!recursive.parse(args[0] + "+1"))
                return false;
            type = recursive.type();
            gk = recursive.value();
            gk -= 1;
            if (gk%2 != 0)
                return false;
            gk /= 2;
            recursive._c = (neg ? -1 : 1);
            gk += recursive._c;
            if (recursive.k() != 1)
            {
                gk *= recursive.gk();
                custom_k = "*" + recursive.gk().to_string();
            }
            gb = recursive.gb();
            custom_b = recursive._custom_b;
            n = recursive._n;
            multifactorial = recursive._multifactorial;
            gd = recursive.gd();
            custom_d = recursive._custom_d;
            c = 1;
            if (recursive.k() > 0 && recursive.k() < 100 && recursive.d() == 1)
            {
                algebraic_type = ALGEBRAIC_QUAD;
                algebraic_k = (int32_t)recursive.k()*recursive._c;
            }
            if (type == KBNC)
            {
                if (recursive.gk()%2 == 0)
                    recursive.gk() /= 2;
                else
                {
                    recursive.gk() *= recursive.gb()/2;
                    recursive._n--;
                }
                custom_k = "(" + recursive.build_text() + ")" + custom_k;
            }
            else if (type == FACTORIAL || type == PRIMORIAL)
            {
                std::string st = recursive.build_text();
                custom_k = "(" + st.substr(0, st.size() - 2) + "/2" + st.substr(st.size() - 2) + ")" + custom_k;
            }
        }
        else if (f == "Hex" && args.size() == 1)
        {
            bool neg = false;
            if (!args[0].empty() && args[0][0] == '-')
            {
                neg = true;
                args[0] = args[0].substr(1);
            }
            InputNum recursive;
            if (!recursive.parse(args[0] + "+1"))
                return false;
            type = recursive.type();
            gk = recursive.value();
            gk -= 1;
            if (gk%3 != 0)
                return false;
            gk /= 3;
            recursive._c = (neg ? -1 : 1);
            gk += recursive._c;
            if (recursive.k() != 1)
            {
                gk *= recursive.gk();
                custom_k = "*" + recursive.gk().to_string();
            }
            gb = recursive.gb();
            custom_b = recursive._custom_b;
            n = recursive._n;
            multifactorial = recursive._multifactorial;
            gd = recursive.gd();
            custom_d = recursive._custom_d;
            c = 1;
            if (recursive.k() > 0 && recursive.k() < 32 && recursive.d() == 1)
            {
                algebraic_type = ALGEBRAIC_HEX;
                algebraic_k = (int32_t)recursive.k()*recursive._c;
            }
            if (type == KBNC)
            {
                if (recursive.gk()%3 == 0)
                    recursive.gk() /= 3;
                else
                {
                    recursive.gk() *= recursive.gb()/3;
                    recursive._n--;
                }
                custom_k = "(" + recursive.build_text() + ")" + custom_k;
            }
            else if (type == FACTORIAL || type == PRIMORIAL)
            {
                std::string st = recursive.build_text();
                custom_k = "(" + st.substr(0, st.size() - 2) + "/3" + st.substr(st.size() - 2) + ")" + custom_k;
            }
            // Montgomery reduction is faster than k^6/27 * b^6n + 1
            //algebraic_type = ALGEBRAIC_SIMPLE;
            //algebraic_k = 0;
        }
        else
            return false;
    }
    else
    {
        gk = 1;
        while (it != s.end() && *it == '*')
        {
            tmp = 1;
            std::string tmp_k;
            if (!f.empty())
            {
                if (!parse_func(f, args, gk == 1 ? gk : tmp, gk == 1 ? custom_k : tmp_k))
                    return false;
            }
            else if (!parse_digits(prime, it_s, it, gk == 1 ? gk : tmp))
                return false;
            if (tmp != 1)
            {
                if (!custom_k.empty() || !tmp_k.empty())
                {
                    if (custom_k.empty())
                        custom_k = gk.to_string();
                    if (tmp_k.empty())
                        tmp_k = tmp.to_string();
                    custom_k += "*";
                    custom_k += tmp_k;
                }
                gk *= tmp;
            }

            it++;
            iterate_func(f, it, s.end());
            if (!f.empty() || (it != s.end() && *it == '('))
                iterate_args(f, args, it, s.end());
            else
                iterate_digits(prime, it, it_s, s.end());
        }

        if (it == s.end() || *it == '\"' || *it == '^' || *it == '+' || *it == '-' || *it == '/')
        {
            if (!f.empty())
            {
                if (!parse_func(f, args, gb, custom_b))
                    return false;
            }
            else if (!parse_digits(prime, it_s, it, gb))
                return false;

            if (it != s.end() && *it == '^')
            {
                iterate_digits(prime, ++it, it_s, s.end());
                if (!parse_digits(prime, it_s, it, n))
                    return false;
            }
            else
                n = 1;
        }
        else if (*it == '!')
        {
            type = FACTORIAL;
            if (!parse_digits(prime, it_s, it, n))
                return false;
            for (multifactorial = 1, it++; it != s.end() && *it == '!'; it++, multifactorial++);
            if (multifactorial == 1)
            {
                iterate_digits(prime, it, it_s, s.end());
                if (it != it_s && !parse_digits(prime, it_s, it, multifactorial))
                    return false;
            }
            gb = 1;
            tmp = 1;
            uint32_t i = n%multifactorial;
            while (i < 2)
                i += multifactorial;
            for (; i <= n; i += multifactorial)
            {
                tmp *= i;
                if (tmp.size() > 8192 || i == n)
                {
                    gb *= tmp;
                    tmp = 1;
                }
            }
        }
        else if (*it == '#')
        {
            type = PRIMORIAL;
            if (it == it_s || it - it_s > (prime ? 8 : 9))
                return false;
            n = stoi(std::string(it_s, it));
            it++;
            gb = 1;
            tmp = 1;
            int last = 1;
            PrimeIterator primes = PrimeIterator::get();
            for (uint32_t i = 0; prime ? i < n : *primes <= (int)n; i++, primes++)
            {
                last = *primes;
                tmp *= last;
                if (tmp.size() > 8192)
                {
                    gb *= tmp;
                    tmp = 1;
                }
            }
            n = last;
            if (tmp != 1)
                gb *= tmp;
        }
        else
            return false;

        gd = 1;
        while (it != s.end() && *it == '/')
        {
            it++;
            iterate_func(f, it, s.end());
            if (!f.empty() || (it != s.end() && *it == '('))
                iterate_args(f, args, it, s.end());
            else
                iterate_digits(prime, it, it_s, s.end());

            tmp = 1;
            std::string tmp_d;
            if (!f.empty())
            {
                if (!parse_func(f, args, gd == 1 ? gd : tmp, gd == 1 ? custom_d : tmp_d))
                    return false;
            }
            else if (!parse_digits(prime, it_s, it, gd == 1 ? gd : tmp))
                return false;
            if (tmp != 1)
            {
                if (!custom_d.empty() || !tmp_d.empty())
                {
                    if (custom_d.empty())
                        custom_d = gd.to_string();
                    if (tmp_d.empty())
                        tmp_d = tmp.to_string();
                    custom_d += "/";
                    custom_d += tmp_d;
                }
                gd *= tmp;
            }
        }
        if (gd > 1 && gk*(type == KBNC ? power(gb, n) : gb)%gd != 0)
            return false;

        if (it != s.end() && (*it == '+' || *it == '-'))
        {
            bool minus = (*it == '-');
            iterate_digits(prime, ++it, it_s, s.end());
            if (!parse_digits(prime, it_s, it, c))
                return false;
            if (minus)
                c = -c;
        }
        else if (type == KBNC && gk == 1 && n == 1)
        {
            type = GENERIC;
            n = 0;
            gb /= gd;
            if (!custom_b.empty() && custom_d.empty())
                custom_d = gd.to_string();
            gd = 1;
        }
        else if (c_required)
            return false;
    }
    if (it != s.end() && *it == '\"')
        it++;
    for (; it != s.end() && std::isspace(*it); it++);
    if (it != s.end())
        return false;

    _type = type;
    _gk = std::move(gk);
    _gb = std::move(gb);
    _n = n;
    _gd = std::move(gd);
    _c = c;
    _custom_k = std::move(custom_k);
    _custom_b = std::move(custom_b);
    _custom_d = std::move(custom_d);
    _multifactorial = multifactorial;
    _algebraic_type = algebraic_type;
    _algebraic_k = algebraic_k;
    process();
    return true;
}

void add_factor(std::vector<std::pair<arithmetic::Giant, int>>& factors, const Giant& factor, int power)
{
    auto it = factors.begin();
    for (; it != factors.end() && it->first != factor; it++);
    if (it != factors.end())
        it->second += power;
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
        if (abs(_c) == 1)
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
            if (tmp > 1 && (uint32_t)tmp.bitlen() <= 2*(s - 1 + (j + 5 < s ? j + 5 : s)))
            {
                add_factor(factors, tmp, 1);
                tmp = 1;
            }
        }
    }
    if (tmp > 1)
        cofactor = std::move(tmp);

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

void InputNum::process()
{
    _factors.clear();
    if (!_cofactor.empty())
        _cofactor.arithmetic().free(_cofactor);
    _b_factors.clear();
    if (!_b_cofactor.empty())
        _b_cofactor.arithmetic().free(_b_cofactor);

    if (_type != KBNC)
    {
        /*Giant tmp;
        _b_factors.reserve(_n);
        PrimeIterator primes = PrimeIterator::get();
        for (i = 0; i < _n; i++, primes++)
        {
            tmp = *primes;
            _b_factors.emplace_back(tmp, 1);
        }*/
    }
    else
    {
        factorize(_gb, _b_factors, _b_cofactor);
    }
    
    if (_type != GENERIC)
    {
        factorize(_gk, _factors, _cofactor);
        for (auto& factor : _b_factors)
            ::add_factor(_factors, factor.first, factor.second*_n);
        std::sort(_factors.begin(), _factors.end(), [](std::pair<arithmetic::Giant, int>& a, std::pair<arithmetic::Giant, int>& b) { return a.first < b.first; });
        if (!_b_cofactor.empty() && !_cofactor.empty())
            _cofactor *= power(_b_cofactor, _n);
        else if (!_b_cofactor.empty())
            _cofactor = power(_b_cofactor, _n);
        if (_gd > 1)
        {
            std::vector<std::pair<arithmetic::Giant, int>> factors;
            arithmetic::Giant cofactor;
            factorize(_gd, factors, cofactor);
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

    if (_type == KBNC && _c == 1 && _gd == 1 && _algebraic_type == ALGEBRAIC_SIMPLE)
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

    _input_text = build_text();

    if (_gd > 1)
    {
        Giant tmp = gcd(_gk, _gd);
        if (tmp > 1)
        {
            if (!_custom_k.empty() && _custom_d.empty())
                _custom_d = _gd.to_string();
            if (_custom_k.empty() && !_custom_d.empty())
                _custom_k = _gk.to_string();
            _gk /= tmp;
            _gd /= tmp;
        }
        if (_gd > 1 && _type == KBNC)
        {
            uint32_t n = (uint32_t)((32 + log2(_gd))/log2(_gb));
            if (n > _n - 1)
                n = _n - 1;
            if (n > 0)
            {
                tmp = power(_gb, n);
                if (tmp%_gd == 0)
                {
                    tmp /= _gd;
                    while (tmp%_gb == 0)
                    {
                        tmp /= _gb;
                        n--;
                    }
                    if (!_custom_k.empty())
                    {
                        if (_custom_d.empty())
                            _custom_k += "*" + tmp.to_string();
                        else
                            _custom_k += "*(" + (!_custom_b.empty() ? _custom_b : _gb.to_string()) + "^" + std::to_string(n) + ")";
                    }
                    else if (!_custom_d.empty())
                        _custom_d.clear();
                    _gk *= tmp;
                    _n -= n;
                    _gd = 1;
                }
            }
        }
        if (_gd > 1 && _type != KBNC)
        {
            _gb /= _gd;
            if (_custom_d.empty())
                _custom_d = _gd.to_string();
            _gd = 1;
        }
    }

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

    _display_text = build_text(30);
}

std::string InputNum::build_text(int max_len)
{
    std::string res;
    res.reserve(32);
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
    if (_gd > 1 || !_custom_d.empty())
    {
        res.append(1, '/');
        std::string sd = !_custom_d.empty() ? _custom_d : _gd.to_string();
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
        res.append(std::to_string(abs(_c)));
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
        state.setup(_gk*_gb/_gd + _c);
    }
    else if (_algebraic_type != ALGEBRAIC_SIMPLE && b() != 0 && !state.force_mod_type)
    {
        if (_algebraic_type == ALGEBRAIC_CYCLOTOMIC)
        {
            uint32_t k = (uint32_t)abs(_algebraic_k);
            state.known_factors = k*power(_gb, _n) + (_algebraic_k < 0 ? 1 : -1);
            state.setup(k*k*k, b(), 3*_n, _algebraic_k < 0 ? 1 : -1);
        }
        if (_algebraic_type == ALGEBRAIC_QUAD)
        {
            Giant x = (-_algebraic_k)*power(_gb, _n);
            x <<= 1;
            state.known_factors = value() + x;
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
            state.known_factors = (val + std::move(x))*val;
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
    else if (k() != 0 && b() != 0 && d() == 1)
    {
        state.setup(k(), b(), _n, _c);
    }
    else
    {
        state.setup(_gk*power(_gb, _n)/_gd + _c);
    }
    if (state.fingerprint != fingerprint())
        throw ArithmeticException();
}

int InputNum::bitlen()
{
    if (_type == GENERIC)
        return _gb.bitlen();
    else if (_type != KBNC)
        return _gk.bitlen() + _gb.bitlen() - _gd.bitlen();
    else if (b() == 2)
        return _gk.bitlen() + _n - _gd.bitlen() + 1;
    return (int)std::ceil(log2(_gk) + log2(_gb)*_n - log2(_gd));
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
        uint32_t c = (uint32_t)(-_c);
        if (result < c)
            result += (uint64_t)modulus << 31;
        result -= c;
    }
    return (uint32_t)(result%modulus);
}

bool InputNum::is_half_factored()
{
    if (abs(_c) != 1)
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
    factorize(minus1, factors, minus1, [&](Giant& x, uint32_t p) { return mod(p) == 1; }, s);

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

void InputNum::factorize_f_p()
{
    if (_type != FACTORIAL && _type != PRIMORIAL)
        return;
    std::map<uint32_t, int> factors;
    for (auto& f : _factors)
        if (f.first.size() == 1)
            factors[*f.first.data()] = f.second;

    if (_type == FACTORIAL)
    {
        uint32_t i = _n%_multifactorial;
        while (i < 2)
            i += _multifactorial;
        for (; i <= _n; i += _multifactorial)
        {
            uint32_t j = i;
            for (auto it = PrimeIterator::get(); (*it)*(*it) <= (int)j; it++)
                if (j%(*it) == 0)
                {
                    int& power = factors[*it];
                    do
                    {
                        power++;
                        j /= *it;
                    } while (j%(*it) == 0);
                }
            if (j != 1)
                factors[j]++;
        }
    }
    if (_type == PRIMORIAL)
    {
        for (auto it = PrimeIterator::get(); *it <= (int)_n; it++)
            factors[*it]++;
    }

    Giant tmp;
    std::vector<std::pair<Giant, int>> new_factors;
    for (auto& f : factors)
    {
        tmp = f.first;
        new_factors.emplace_back(std::move(tmp), f.second);
    }
    for (auto& f : _factors)
        if (f.first.size() > 1)
            new_factors.emplace_back(std::move(f.first), f.second);
    _factors = std::move(new_factors);
}

void InputNum::print_info()
{
    if (_type == ZERO)
        return;

    Giant N = value();
    Giant tmp;
    std::string st;

    std::cout << display_text() << ", " << N.bitlen() << " bits, " << N.digits() << " digits." << std::endl;

    std::vector<std::pair<arithmetic::Giant, int>> factors;
    arithmetic::Giant cofactor;
    factorize(N, factors, cofactor, [&](Giant& x, uint32_t p) { return mod(p) == 0; });
    if (factors.empty())
        std::cout << "No small factors." << std::endl;
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

    if (_type == KBNC && abs(_c) == 1 && d() == 1 && (k() == 1 || _cofactor.empty() || (!_b_cofactor.empty() && _cofactor == power(_b_cofactor, _n))))
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
    if (_type == KBNC && _algebraic_type == ALGEBRAIC_CYCLOTOMIC && abs(_algebraic_k) == 1)
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
    if (_c == 1)
    {
        factors = _factors;
        cofactor = _cofactor;
    }
    else
    {
        N -= 1;
        factorize(N, factors, cofactor, [&](Giant& x, uint32_t p) { return mod(p) == 1; });
        N += 1;
    }
    if (factors.empty() && (_c != 1 || (_type != FACTORIAL && _type != PRIMORIAL)))
        std::cout << "No small factors of N-1." << std::endl;
    else
    {
        tmp = 1;
        st.clear();
        if (_c == 1 && _type == FACTORIAL)
        {
            tmp = _gb;
            st = std::to_string(_n);
            if (_multifactorial <= 3)
                st.append(_multifactorial, '!');
            else
                st.append(1, '!').append(std::to_string(_multifactorial));
        }
        if (_c == 1 && _type == PRIMORIAL)
        {
            tmp = _gb;
            st = std::to_string(_n) + "#";
        }
        for (auto& factor : factors)
        {
            if (factor.second == 0)
                continue;
            if (!st.empty())
                st += " ";
            if (factor.second < 0)
                st += "/";
            st += factor.first.to_string();
            if (abs(factor.second) > 1)
            {
                st += "^";
                st += std::to_string(abs(factor.second));
            }
            if (!cofactor.empty())
            {
                if (factor.second > 0)
                    tmp *= power(factor.first, factor.second);
                else
                    tmp /= power(factor.first, -factor.second);
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
    if (_c == -1)
    {
        factors = _factors;
        cofactor = _cofactor;
    }
    else
    {
        N += 1;
        factorize(N, factors, cofactor, [&](Giant& x, uint32_t p) { return mod(p) == p - 1; });
        N -= 1;
    }
    if (factors.empty() && (_c != -1 || (_type != FACTORIAL && _type != PRIMORIAL)))
        std::cout << "No small factors of N+1." << std::endl;
    else
    {
        tmp = 1;
        st.clear();
        if (_c == -1 && _type == FACTORIAL)
        {
            tmp = _gb;
            st = std::to_string(_n);
            if (_multifactorial <= 3)
                st.append(_multifactorial, '!');
            else
                st.append(1, '!').append(std::to_string(_multifactorial));
        }
        if (_c == -1 && _type == PRIMORIAL)
        {
            tmp = _gb;
            st = std::to_string(_n) + "#";
        }
        for (auto& factor : factors)
        {
            if (factor.second == 0)
                continue;
            if (!st.empty())
                st += " ";
            if (factor.second < 0)
                st += "/";
            st += factor.first.to_string();
            if (abs(factor.second) > 1)
            {
                st += "^";
                st += std::to_string(abs(factor.second));
            }
            if (!cofactor.empty())
            {
                if (factor.second > 0)
                    tmp *= power(factor.first, factor.second);
                else
                    tmp /= power(factor.first, -factor.second);
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
