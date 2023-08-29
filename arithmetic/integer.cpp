
#include <cmath>
#include "integer.h"

namespace arithmetic
{
    uint32_t gcd(uint32_t a, uint32_t b)
    {
        while (a)
        {
            std::swap(a, b);
            a = a%b;
        }
        return b;
    }

    bool is_prime(uint32_t a)
    {
        for (auto it = PrimeIterator::get(); (uint32_t)*it*(uint32_t)(*it) <= a; it++)
            if (a%(uint32_t)(*it) == 0)
                return false;
        return true;
    }

    uint32_t phi(uint32_t N)
    {
        uint32_t res = 1;
        for (auto it = PrimeIterator::get(); (uint32_t)*it*(uint32_t)(*it) <= N; it++)
            if (N%(uint32_t)(*it) == 0)
            {
                N /= (uint32_t)*it;
                res *= (uint32_t)*it - 1;
                while (N%(uint32_t)(*it) == 0)
                {
                    N /= (uint32_t)*it;
                    res *= (uint32_t)*it;
                }
            }
        if (N > 1)
            res *= N - 1;
        return res;
    }

    int jacobi(uint32_t a, uint32_t b)
    {
        /* Computes Jacobi (a, b) */

        uint32_t  jdvs, jdvd, jr;
        int resul; uint32_t s, t, u, v;
        jdvs = a;
        jdvd = b;
        resul = 1;
        while (jdvs)
            if (jdvs==1)				/* Finished ! */
                return(resul);
            else {
                v = jdvd;
                s = (v-1)>>1;  		/* (dvd-1)/2 */
                t = (v+1)>>1;  		/* (dvd+1)/2 */
                while (!(jdvs & 1)) {	/* While dvs is even */
                                /* resul *= (-1)**(dvd**2-1)/8; */
                    if (t & 1) {
                        if ((s>>1) & 1)
                            resul = -resul;
                    }
                    else
                        if ((t>>1) & 1)
                            resul = -resul;
                    jdvs >>= 1;		/* dvs /= 2; */
                }
                if (jdvs==1)			/* Finished ! */
                    return(resul);
                else {
                    u = (jdvs-1)>>1; 	/* (dvs-1)/2 */
                    if (s & u & 1)	/* resul *= (-1)**(dvd-1)*(dvs-1)/4; */
                        resul = -resul;
                    jr = jdvd%jdvs;
                    jdvd = jdvs;
                    jdvs = jr;
                }      			/* dvs != 1 */
            }					/* dvs != 1 */
        return(jdvd);				/* a and b are not coprime, */
                                /* so, return their gcd. */
    }

    int kronecker(uint32_t a, uint32_t b)
    {
        /* Computes Kronecker (a, b) = (a/b) */
        /* Note : returns gcd (a, b) if a and b are not coprime. */

        uint32_t  jdvs, jdvd, jr, jmul;
        int resul; uint32_t s, t, u, v, amod8;
        jdvs = a;
        jdvd = b;
        jmul = 1;
        amod8 = a & 7;
        resul = 1;
        while (!(jdvd&1)) {			/* if b is even... */
            jdvd >>= 1;
            if ((amod8 == 3) || (amod8 == 5))
                resul = -resul;
            else if (!(jdvs & 1)) {
                jdvs >>= 1;
                jmul <<= 1;
            }
        }		/* Now, jdvd is odd and jmul is the largest power of 2 dividing both a and b. */
        if (jdvd == 1 || jdvs == 1)
            return((jmul == 1) ? resul : jmul);	/* Finished ! */
        while (jdvs)
            if (jdvs==1)				/* Finished ! */
                return((jmul == 1) ? resul : jmul);
            else {
                v = jdvd;
                s = (v-1)>>1;  		/* (dvd-1)/2 */
                t = (v+1)>>1;  		/* (dvd+1)/2 */
                while (!(jdvs & 1)) {	/* While dvs is even */
                                /* resul *= (-1)**(dvd**2-1)/8; */
                    if (t & 1) {
                        if ((s>>1) & 1)
                            resul = -resul;
                    }
                    else
                        if ((t>>1) & 1)
                            resul = -resul;
                    jdvs >>= 1;		/* dvs /= 2; */
                }
                if (jdvs==1)			/* Finished ! */
                    return((jmul == 1) ? resul : jmul);
                else {
                    u = (jdvs-1)>>1; 	/* (dvs-1)/2 */
                    if (s & u & 1)	/* resul *= (-1)**(dvd-1)*(dvs-1)/4; */
                        resul = -resul;
                    jr = jdvd%jdvs;
                    jdvd = jdvs;
                    jdvs = jr;
                }      			/* dvs != 1 */
            }					/* dvs != 1 */
        return(jdvd*jmul);			/* a and b are not coprime, */
                                /* so, return their gcd. */
    }

    std::unique_ptr<PrimeList> PrimeList::_list65536;

    PrimeList::PrimeList(int max)
    {
        int i, j, k;
        std::vector<bool> bitmap(max/2, false);
        k = 0;
        for (i = 1; i < max/2; i++)
            if (!bitmap[i])
            {
                k++;
                if (i > 16384)
                    continue;
                for (j = (i*2 + 1)*(i*2 + 1)/2; j < bitmap.size(); j += i*2 + 1)
                    bitmap[j] = true;
            }
        primes.clear();
        primes.reserve(k + 1);
        primes.push_back(2);
        for (i = 1; i < max/2; i++)
            if (!bitmap[i])
                primes.push_back(i*2 + 1);
    }

    template<class T, class IT>
    void sieve_range_t(T start, T end, std::vector<T>& list, IT& it, const IT& it_end)
    {
        int i, p;
        list.clear();
#ifndef __clang__
        list.reserve((int)((std::expint(log(end)) - (start > 2 ? std::expint(log(start)) : 0))*1.1) + 100);
#endif
        if (start <= 2 && 2 < end)
            list.push_back(2);
        if (!(start & 1))
            start++;
        if (!(end & 1))
            end++;
        std::vector<bool> bitmap((end - start)/2, false);
        if (start == 1 && end > 1)
            bitmap[0] = true;
        if (*it == 2)
            it++;
        for (; it != it_end && (*it)*(T)(*it) < end; it++)
        {
            p = *it;
            i = start%p;
            if (i != 0)
                i = p - i;
            if (i & 1)
                i += p;
            i >>= 1;
            for (; i < bitmap.size(); i += p)
                bitmap[i] = true;
            if (start <= p)
                list.push_back(p);
        }
        for (i = 0; i < bitmap.size(); i++)
            if (!bitmap[i])
                list.push_back(start + i*2);
    }

    void PrimeList::sieve_range(int start, int end, std::vector<int>& list)
    {
        std::vector<int>::iterator it = primes.begin();
        sieve_range_t<int,std::vector<int>::iterator>(start, end, list, it, primes.size() < 4792 ? primes.end() : primes.begin() + 4792);
    }

    PrimeIterator PrimeList::begin()
    {
        return PrimeIterator(*this);
    }

    PrimeIterator& PrimeIterator::operator++()
    {
        return (*this) += 1;
    }

    void PrimeIterator::operator++(int)
    {
        (*this) += 1;
    }

    PrimeIterator& PrimeIterator::operator+=(int offset)
    {
        _cur += offset;
        if (_cur > 105097564)
            _cur = 105097564;
        if (_cur < _list.size())
            return *this;
        while (_cur - _range_pos >= _range.size())
        {
            int start;
            if (_range_pos == 0)
            {
                _range_pos = _list.size();
                start = _list[_list.size() - 1] + 2;
            }
            else
            {
                _range_pos += _range.size();
                start = _range[_range.size() - 1] + 2;
            }
            int range = _list[_list.size() - 1];
            if (start < (int)((1UL << 31) - range))
                _list.sieve_range(start, start + range, _range);
            else
            {
                _list.sieve_range(start, (int)((1UL << 31) - 1), _range);
                _range.push_back((int)((1UL << 31) - 1));
            }
        }
        return *this;
    }

    std::unique_ptr<PrimeIterator> _iter105097564;

    const PrimeIterator& PrimeIterator::max()
    {
        if (!_iter105097564)
        {
            _iter105097564.reset(new PrimeIterator(PrimeList::primes_16bit()));
            _iter105097564->_cur = 105097564;
            _iter105097564->_range_pos = 105097564;
            _iter105097564->_range.push_back((int)((1UL << 31) - 1));
        }
        return *_iter105097564;
    }

    void PrimeIterator::sieve_range(uint64_t start, uint64_t end, std::vector<uint64_t>& list)
    {
        sieve_range_t<uint64_t, PrimeIterator>(start, end, list, *this, max());
    }
}
