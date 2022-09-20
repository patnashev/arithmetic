#include <iostream>

#include "gwnum.h"
#include "arithmetic.h"
#include "lucas.h"
#include "edwards.h"
#include "montgomery.h"

using namespace arithmetic;

int main(int argc, char *argv[])
{
    int i;
    Giant tmp;
    int len;
    
    GWState gwstate;
    gwstate.setup(1, 2, 4096, 1);
    GWArithmetic gw(gwstate);
    //CarefulGWArithmetic gw(gwstate);
    GWNum s(gw.carefully()), t(gw.carefully());
    s = 12;
    t = 40;
    tmp = 7636607;
    len = tmp.bitlen() - 1;
    for (i = 1; i <= len; i++)
    {
        GWNum lambda = (3*square(s) - 8)/(2*t);
        GWNum snext = lambda*lambda - s - s;
        t = lambda*(std::move(s) - snext) - t;
        s = std::move(snext);
        if (tmp.bit(len - i))
        {
            lambda = (t - 40)/(s - 12);
            snext = lambda*lambda - 12 - s;
            t = lambda*(std::move(s) - snext) - t;
            s = std::move(snext);
        }
    }
    GWNum alpha = 1/((t + 25)/(s - 9) + 1);
    GWNum beta = 2*alpha*(4*alpha + 1)/(8*square(alpha) - 1);
    GWNum d = (2*square(2*beta - 1) - 1)/square(square((2*beta - 1)));
    GWNum X1 = (2*beta - 1)*(4*beta - 3)/(6*beta - 5);
    GWNum Y1 = (2*beta - 1)*(t*(t + 50) - 104 - square(s)*(2*s - 27))/((t - 2 + 3*s)*(t + 16 + s));
    
    EdwardsArithmetic ed(gw);
    EdPoint P1(ed, X1, Y1);
    EdPoint P = ed.gen_curve(7636607, nullptr);
    std::cout << ed.on_curve(P, d) << std::endl;
    P1.Z.reset(new GWNum(gw));
    *P1.X *= 2;
    *P1.Y *= 2;
    *P1.Z = 2;
    P1.extend();

    tmp = 1 << 5;
    tmp = std::move(tmp) * 2949403 * 542441 * 399577 * 381529 * 235489 * 149441 * 55061 * 5849 * 4159 * 2477 * 13 * 3;
    len = tmp.bitlen() - 1;
    for (i = 1; i <= len; i++)
    {
        if (tmp.bit(len - i))
        {
            ed.dbl(P, P, GWMUL_STARTNEXTFFT);
            ed.add(P, P1, P, GWMUL_STARTNEXTFFT | ed.ED_PROJECTIVE);
        }
        else
            ed.dbl(P, P, GWMUL_STARTNEXTFFT | ed.ED_PROJECTIVE);
    }

    std::vector<int16_t> naf_w;
    get_NAF_W(5, tmp, naf_w);
    ed.mul(P1, 5, naf_w, P1);

    std::cout << (P == P1) << std::endl;
    std::cout << (*P.X != *P1.X) << std::endl;
    P.normalize();
    P1.normalize();
    std::cout << (*P.X == *P1.X) << std::endl;

    std::cout << gcd(*P.X, gw.N()).to_string() << std::endl;

    P = ed.from_small(17, 19, 17, 33, &d);
    MontgomeryArithmetic mont(gw, d);
    EdY ma(mont, P);
    EdY mb = ma;
    mont.dbl(mb, mb);
    mont.add(mb, ma, ma, mb);
    tmp = 3;
    P *= tmp;
    P.normalize();
    mb.normalize();
    std::cout << (*P.Y == *mb.Y) << std::endl;

    GiantsArithmetic giants;
    Giant a(giants);
    a = 1;
    Giant b = a;
    a <<= 31;
    if (0 < a)
        a = "8";
    a = std::move(a) + a;
    a = std::move(a)*a;
    a += b;
    swap(a, b);

    a = -1 + (a + b) + (5 + a);
    b = 250;
    a = std::move(a)%b;

    std::cout << a.to_string() << std::endl;

    b = 5;
    a = 2;
    a *= b;
    a = -(10 - (a - b)*2 - 5 - a);
    a = -std::move(a);

    std::cout << a.to_string() << std::endl;

    a = "123456789012345678901234567890";
    b = "123456789012345678901234567890123456789012345678901234567890";
    a += b;

    std::cout << a.to_string() << std::endl;

    LucasVArithmetic lucas(gw.carefully());
    LucasV V(lucas), V1(lucas), V2(lucas);
    GWNum lucasP(lucas.gw());
    lucasP = 6;
    lucasP /= 5;
    V = lucasP;
    V1 = lucasP;
    for (i = 2; i < 100; i++)
    {
        //V *= i;
        //lucas.mul(V, i, V);
        //swap(V, V2);
        lucas.mul(V, i, -1, V);
        lucas.mul(V1, (tmp = i), V1, V2);
        if (V.V() != V1.V())
            printf("error");
    }

    GWState gwstateProth;
    gwstateProth.setup(224027, 2, 99763, 1);
    //gwstateProth.setup(227753, 2, 91397, 1);
    ReliableGWArithmetic rgw(gwstateProth);
    GWNum x(rgw);
    do
    {
        rgw.restart();
        x = 3;
        rgw.setmulbyconst(3);
        tmp = 224027;
        len = tmp.bitlen() - 1;
        for (i = 1; i <= len; i++)
            rgw.carefully().mul(x, x, x, tmp.bit(len - i) ? GWMUL_MULBYCONST : 0);
        for (i = 1; i <= 99762; i++)
            rgw.mul(x, x, x, GWMUL_STARTNEXTFFT);
    } while (rgw.restart_flag() && !rgw.failure_flag());

    tmp = x;
    if (tmp == rgw.N() - 1)
        std::cout << "prime\n";
    else
        printf("%016" PRIX64 "\n", *(uint64_t*)tmp.data());

    GWState gwstateP;
    gwstateP.setup(5, 2, 127, 1);
    GWArithmetic gwP(gwstateP);
    FFT fft(gwP, 8);
    Poly poly1(fft);
    poly1[3] = 1; poly1[2] = 2; poly1[1] = 3; poly1[0] = 4;
    fft.transform(poly1, poly1);
    Poly poly2(fft);
    poly2[3] = 5; poly2[2] = 6; poly2[1] = 7; poly2[0] = 8;
    fft.transform(poly2, poly2);
    for (i = 0; i < poly1.size(); i++)
        poly1[i] *= poly2[i];
    fft.inv_transform(poly1, poly1);
    for (int i = 7; i >= 0; i--)
        std::cout << poly1[i].to_string() << std::endl;
 
    return 0;
}
