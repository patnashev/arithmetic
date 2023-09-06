
#include <climits>
#include <cctype>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <stdlib.h>
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
        if (!recursive.parse(args[0]))
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

bool InputNum::parse(const std::string& s)
{
    std::string::const_iterator it, it_s;
    int type = KBNC;
    Giant gk, gb;
    uint32_t n = 0;
    int32_t c = 0;
    Giant tmp;
    bool prime;
    std::string f;
    std::vector<std::string> args;
    std::string custom_k;
    std::string custom_b;
    int32_t cyclotomic_k = 0;
    int32_t hex_k = 0;

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
            n = recursive._n;
            c = 1;
            cyclotomic_k = (cyclotomic == 6 ? -1 : 1)*(recursive.k() > 0 && recursive.k() < 1000 ? recursive.k() : 0);
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
            n = recursive._n;
            c = 1;
            hex_k = (recursive.k() > 0 && recursive.k() < 32 ? recursive.k() : 0)*recursive._c;
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
        }
        else
            return false;
    }
    else
    {
        if (it != s.end() && *it == '*')
        {
            if (!f.empty())
            {
                if (!parse_func(f, args, gk, custom_k))
                    return false;
            }
            else if (!parse_digits(prime, it_s, it, gk))
                return false;

            it++;
            iterate_func(f, it, s.end());
            if (!f.empty() || (it != s.end() && *it == '('))
                iterate_args(f, args, it, s.end());
            else
                iterate_digits(prime, it, it_s, s.end());
        }
        else
            gk = 1;
        if (it == s.end() || *it == '^' || *it == '\"')
        {
            if (!f.empty())
            {
                if (!parse_func(f, args, gb, custom_b))
                    return false;
            }
            else if (!parse_digits(prime, it_s, it, gb))
                return false;
        }
        else if (*it == '!')
        {
            type = FACTORIAL;
            if (!parse_digits(prime, it_s, it, n))
                return false;
            gb = 1;
            tmp = 1;
            for (uint32_t i = 2; i <= n; i++)
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
        if (it != s.end())
        {
            if (*it == '^')
            {
                iterate_digits(prime, ++it, it_s, s.end());
                if (!parse_digits(prime, it_s, it, n))
                    return false;
            }
            else
                it++;
            bool minus = false;
            if (it == s.end())
                return false;
            else if (*it == '+')
                minus = false;
            else if (*it == '-')
                minus = true;
            else
                return false;
            iterate_digits(prime, ++it, it_s, s.end());
            if (!parse_digits(prime, it_s, it, c))
                return false;
            if (minus)
                c = -c;
        }
        else
        {
            n = 1;
            c = 0;
        }
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
    _c = c;
    _custom_k = custom_k;
    _custom_b = custom_b;
    _cyclotomic_k = cyclotomic_k;
    _hex_k = hex_k;
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
    
    if (abs(_c) == 1)
    {
        factorize(_gk, _factors, _cofactor);
        for (auto& factor : _b_factors)
            ::add_factor(_factors, factor.first, factor.second*_n);
        std::sort(_factors.begin(), _factors.end(), [](std::pair<arithmetic::Giant, int>& a, std::pair<arithmetic::Giant, int>& b) { return a.first < b.first; });
        if (!_b_cofactor.empty() && !_cofactor.empty())
            _cofactor *= power(_b_cofactor, _n);
        else if (!_b_cofactor.empty())
            _cofactor = power(_b_cofactor, _n);
    }

    if (_type == KBNC && _c == 1 && _cyclotomic_k == 0 && _hex_k == 0)
    {
        if (_gk == 1 && _n > 1 && (_n & (_n - 1)) == 0)
            for (_gfn = 1; (1UL << _gfn) < _n; _gfn++);
        if ((abs(_gk.bitlen() - log2(_gb)*_n) < _n/30 + 10))
        {
            Giant tmp = _gb;
            tmp.power(_n);
            Giant cyclo_tmp = tmp - _gk;
            if (cyclo_tmp == 1)
                _cyclotomic_k = -1;
            if (cyclo_tmp == -1)
                _cyclotomic_k = 1;
            if (_gb%3 == 0)
            {
                tmp /= 3;
                tmp -= _gk;
                if (tmp == 1)
                    _hex_k = -1;
                if (tmp == -1)
                    _hex_k = 1;
            }
        }
    }

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

        while (_gk%_gb == 0 && _cyclotomic_k == 0 && _hex_k == 0)
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
    if (_gk != 1)
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
        res.append(1, '!');
    }
    else if (_type == PRIMORIAL)
    {
        res.append(std::to_string(_n));
        res.append(1, '#');
    }
    else if (_gb != 1)
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
        state.setup(_gk*_gb + _c);
    }
    else if ((_cyclotomic_k != 0 || _hex_k != 0) && b() != 0 && !state.force_general_mod)
    {
        if (_cyclotomic_k != 0)
        {
            uint32_t k = (uint32_t)abs(_cyclotomic_k);
            state.known_factors = k*power(_gb, _n) + (_cyclotomic_k < 0 ? 1 : -1);
            state.setup(k*k*k, b(), 3*_n, _cyclotomic_k < 0 ? 1 : -1);
        }
        if (_hex_k != 0)
        {
            Giant x = (-_hex_k)*power(_gb, _n);
            Giant val = value() + x;
            state.known_factors = (val + std::move(x))*val;
            uint64_t hex_k = (uint64_t)(_hex_k*_hex_k);
            if (hex_k%3 == 0)
                state.setup(hex_k*hex_k*hex_k/27, b(), 6*_n, 1);
            else
                state.setup(hex_k*hex_k*hex_k*b()*b()*b()/27, b(), 6*_n - 3, 1);
        }
        if (*state.N%3417905339UL != fingerprint())
            throw ArithmeticException();
        if (state.gwdata()->GENERAL_MOD)
        {
            state.done();
            state.known_factors = 1;
            _cyclotomic_k = 0;
            _hex_k = 0;
            setup(state);
        }
        return;
    }
    else if (k() != 0 && b() != 0)
    {
        state.setup(k(), b(), _n, _c);
    }
    else
    {
        state.setup(_gk*power(_gb, _n) + _c);
    }
    if (state.fingerprint != fingerprint())
        throw ArithmeticException();
}

int InputNum::bitlen()
{
    if (_type == GENERIC)
        return _gb.bitlen();
    else if (_type != KBNC)
        return _gk.bitlen() + _gb.bitlen() - 1;
    else if (b() == 2)
        return _gk.bitlen() + _n;
    return (int)std::ceil(log2(_gk) + log2(_gb)*_n);
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
    uint32_t b;
    if (_gb < modulus)
        b = _gb.data()[0];
    else
        b = _gb%modulus;
    if (_type == GENERIC)
        return b;
    uint64_t result;
    if (_type == KBNC)
    {
        result = (_n & 1) ? b : 1;
        uint64_t mult = b;
        for (uint32_t i = 2; i <= _n; i <<= 1)
        {
            mult *= mult;
            if (mult >= (1ULL << 32))
                mult %= modulus;
            if ((_n & i) != 0)
            {
                result *= mult;
                if (result >= (1ULL << 32))
                    result %= modulus;
            }
        }
    }
    else
        result = b;
    if (_gk > 1)
    {
        if (_gk < modulus)
            result *= _gk.data()[0];
        else
            result *= _gk%modulus;
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
            st = std::to_string(_n) + "!";
        }
        if (_c == 1 && _type == PRIMORIAL)
        {
            tmp = _gb;
            st = std::to_string(_n) + "#";
        }
        for (auto& factor : factors)
        {
            if (!st.empty())
                st += " ";
            st += factor.first.to_string();
            if (factor.second > 1)
            {
                st += "^";
                st += std::to_string(factor.second);
            }
            if (!cofactor.empty())
                tmp *= power(factor.first, factor.second);
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
            st = std::to_string(_n) + "!";
        }
        if (_c == -1 && _type == PRIMORIAL)
        {
            tmp = _gb;
            st = std::to_string(_n) + "#";
        }
        for (auto& factor : factors)
        {
            if (!st.empty())
                st += " ";
            st += factor.first.to_string();
            if (factor.second > 1)
            {
                st += "^";
                st += std::to_string(factor.second);
            }
            if (!cofactor.empty())
                tmp *= power(factor.first, factor.second);
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
