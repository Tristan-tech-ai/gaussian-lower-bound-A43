"""check_hand_J2.py (11 Nov 2026): numbers behind the hand argument for the bound p > exp(-2) c* (paper A, Lemma 2(c)).

Written BEFORE the run (what the script must show; a reading, not a verdict script, no registration):
  (1) the short form J(2) = -5/2 + 2 sqrt(6) phi(0) + 7 Phibar(sqrt 6) - sqrt(6) phi(sqrt 6) agrees with the form of the paper,
      J = 2[-lambda ja + R jb - jc/2] + 2 Phibar(R), and with a direct quadrature of 2 int_0^inf phi(u) G_2(u) du;
  (2) J(2) <= -5/2 + 2 sqrt(3/pi) + exp(-3)/sqrt(12 pi) < -5/2 + 2 + 1/6 < 0 (the hand inequalities);
  (3) rho(2) = 2 e^2/(e + 1) < 4, and the constant (1/4)(4 rho)^(-1/6) lies above 2^(-8/3) and above 3/20;
  (4) the chord bound g(x0 + u) <= (1 - u/R) g(x0) + (u/R) g(x0 + R) + u(R - u)/2 holds on random Gaussian mixtures.
Usage: py check_hand_J2.py   (python-flint at 256 bits for (1) to (3), mpmath for the quadrature, floats for (4); no file is read or written)
"""
import math
import random

import mpmath as mp
from flint import arb, ctx

ctx.prec = 256
SQ2 = arb(2).sqrt()
PI = arb.pi()
E = arb(1).exp()


def phi(x):
    return (-(x * x) / 2).exp() / (2 * PI).sqrt()


def Phibar(x):
    return (1 - (x / SQ2).erf()) / 2


def J_paper(lam):
    R = (2 * (1 + lam)).sqrt()
    ja = arb(1) / 2 - Phibar(R)
    jb = phi(arb(0)) - phi(R)
    jc = ja - R * phi(R)
    return 2 * (-lam * ja + R * jb - jc / 2) + 2 * Phibar(R)


def J_short(lam):
    R = (2 * (1 + lam)).sqrt()
    return -(lam + arb(1) / 2) + 2 * R * phi(arb(0)) + (R * R + 1) * Phibar(R) - R * phi(R)


def show(name, x, digits=15):
    print("%-44s %s   (radius %s)" % (name, x.str(digits, radius=False), x.rad().str(3)))


def constants(lam):
    ell = (-lam).exp()
    L1 = (1 - ell) / lam
    rho = (E - 1) / L1
    c0max = (4 * rho) ** (arb(-1) / 3)
    return rho, c0max, c0max.sqrt() / 4


print("precision: %d bits" % ctx.prec)
print("(1) three forms of J(2)")
two = arb(2)
R6 = arb(6).sqrt()
jp = J_paper(two)
js = J_short(two)
jh = -arb(5) / 2 + 2 * R6 * phi(arb(0)) + 7 * Phibar(R6) - R6 * phi(R6)
show("J(2), form of the paper", jp)
show("J(2), short form", js)
show("J(2), -5/2 + 2 sqrt6 phi(0) + 7 Phibar - sqrt6 phi", jh)
mp.mp.dps = 40
Rm = mp.sqrt(6)
quad = 2 * mp.quad(lambda u: mp.npdf(u) * (-2 + Rm * u - u * u / 2), [0, Rm]) + 2 * mp.quad(lambda u: mp.npdf(u), [Rm, mp.inf])
print("%-44s %s" % ("J(2), direct quadrature (mpmath, 40 digits)", mp.nstr(quad, 15)))
ok1 = bool(abs(jp - js) < arb("1e-60")) and bool(abs(jp - jh) < arb("1e-60")) and abs(quad - mp.mpf(jp.mid().str(40, radius=False))) < mp.mpf("1e-30")
print("three forms agree: %s" % ok1)

print("(2) the hand inequalities")
mills = phi(R6) / R6
show("Phibar(sqrt 6)", Phibar(R6))
show("phi(sqrt 6)/sqrt 6 (bound on Phibar)", mills)
show("7 Phibar(sqrt 6) - sqrt 6 phi(sqrt 6)", 7 * Phibar(R6) - R6 * phi(R6))
show("exp(-3)/sqrt(12 pi)", (-arb(3)).exp() / (12 * PI).sqrt())
show("2 sqrt(3/pi)", 2 * (3 / PI).sqrt())
bound = -arb(5) / 2 + 2 * (3 / PI).sqrt() + (-arb(3)).exp() / (12 * PI).sqrt()
show("-5/2 + 2 sqrt(3/pi) + exp(-3)/sqrt(12 pi)", bound)
show("-5/2 + 2 + 1/6", -arb(5) / 2 + 2 + arb(1) / 6)
ok2 = (bool(Phibar(R6) < mills) and bool(jp < bound) and bool(2 * (3 / PI).sqrt() < 2)
       and bool((-arb(3)).exp() / (12 * PI).sqrt() < arb(1) / 6) and bool(bound < 0))
print("Phibar < phi/R, J(2) < bound, 2 sqrt(3/pi) < 2, exp(-3)/sqrt(12 pi) < 1/6, bound < 0: %s" % ok2)

print("(3) constants that follow from the depth lambda")
for name, lam in (("2", two), ("3/2", arb(3) / 2), ("1.1648", arb("1.1648"))):
    rho, c0max, kappa = constants(lam)
    show("J(%s)" % name, J_paper(lam))
    show("rho(%s) = (e - 1) lambda/(1 - exp(-lambda))" % name, rho)
    show("c0max(%s) = (4 rho)^(-1/3)" % name, c0max)
    show("constant(%s) = (1/4)(4 rho)^(-1/6)" % name, kappa)
rho2, c0max2, kappa2 = constants(two)
show("2 e^2/(e + 1)", 2 * E * E / (E + 1))
show("2^(-8/3)", arb(2) ** (arb(-8) / 3))
show("1 + sqrt 3", 1 + arb(3).sqrt())
ok3 = (bool(abs(rho2 - 2 * E * E / (E + 1)) < arb("1e-60")) and bool(rho2 < 4) and bool(E < 1 + arb(3).sqrt())
       and bool(kappa2 > arb(2) ** (arb(-8) / 3)) and bool(arb(2) ** (arb(-8) / 3) > arb(3) / 20))
print("rho(2) = 2e^2/(e+1) < 4, e < 1 + sqrt 3, constant(2) > 2^(-8/3) > 3/20: %s" % ok3)

print("(4) chord bound on random Gaussian mixtures (floats, seed 20261111)")
rng = random.Random(20261111)
RF = math.sqrt(6.0)


def logp(y, xs, ws):
    terms = [math.log(w) - 0.5 * (y - x) ** 2 for x, w in zip(xs, ws)]
    m = max(terms)
    return m + math.log(sum(math.exp(t - m) for t in terms)) - 0.5 * math.log(2 * math.pi)


worst = -1e9
cases = 0
for A in (0.5, 3.0, 20.0):
    for K in (1, 2, 5, 20):
        for _ in range(50):
            xs = [rng.uniform(-A, A) for _ in range(K)]
            raw = [rng.random() + 1e-3 for _ in range(K)]
            ws = [r / sum(raw) for r in raw]
            x0 = rng.uniform(-A, A)
            for sign in (1.0, -1.0):
                g0 = logp(x0, xs, ws)
                gR = logp(x0 + sign * RF, xs, ws)
                for i in range(41):
                    u = RF * i / 40
                    lhs = logp(x0 + sign * u, xs, ws)
                    rhs = (1 - u / RF) * g0 + (u / RF) * gR + 0.5 * u * (RF - u)
                    worst = max(worst, lhs - rhs)
                    cases += 1
print("points tested: %d; largest value of left side minus right side: %.3e" % (cases, worst))
ok4 = worst < 1e-9
print("chord bound holds at every tested point: %s" % ok4)
print("CHECK hand_J2: %s" % ("ALL FOUR STATEMENTS HOLD" if (ok1 and ok2 and ok3 and ok4) else "NOT ALL STATEMENTS HOLD"))
