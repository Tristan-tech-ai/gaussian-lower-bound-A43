"""sharp_depth_ball_arithmetic.py: the computation of Appendix C of the paper (the sharp depth), in ball arithmetic.

No step of the proof of Theorem 1 uses it; it gives the sharper constant 0.1660 in place of 3/20. It shows:
  (1) J(1.1648) is a ball that lies entirely below 0;
  (2) J(1.1647) is a ball that lies entirely above 0, so the root of J lies between 1.1647 and 1.1648;
  (3) with lambda1 = 1.1648, ell = exp(-lambda1), L1 = (1 - ell)/lambda1, rho = (e - 1)/L1, c0max = (4 rho)^(-1/3) and
      kappa = sqrt(c0max)/4 = (1/4)(4 rho)^(-1/6), the ball of kappa lies above 0.1660.
J(lambda) = 2[-lambda ja + R jb - jc/2] + 2 Phibar(R), R = sqrt(2(1 + lambda)), ja = Phi(R) - 1/2, jb = phi(0) - phi(R), jc = ja - R phi(R).
This is the form of J in Appendix C written in other terms; proofs/check_hand_J2.py checks that the two forms agree.
Usage: py sharp_depth_ball_arithmetic.py     (python-flint, 256 bits; reads and writes no file)
"""

from flint import arb, ctx

ctx.prec = 256
SQ2 = arb(2).sqrt()
PI = arb.pi()


def phi(x):
    return (-(x * x) / 2).exp() / (2 * PI).sqrt()


def Phi(x):
    return (1 + (x / SQ2).erf()) / 2


def Phibar(x):
    return (1 - (x / SQ2).erf()) / 2


def J(lam):
    R = (2 * (1 + lam)).sqrt()
    ja = Phi(R) - arb(1) / 2
    jb = phi(arb(0)) - phi(R)
    jc = ja - R * phi(R)
    return 2 * (-lam * ja + R * jb - jc / 2) + 2 * Phibar(R)


def show(name, x, digits=15):
    print("%-34s %s   (radius %s)" % (name, x.str(digits, radius=False), x.rad().str(3)))


print("precision: %d bits" % ctx.prec)
lam1 = arb("1.1648")
lam0 = arb("1.1647")
j1 = J(lam1)
j0 = J(lam0)
show("J(1.1648)", j1)
show("J(1.1647)", j0)
print("J(1.1648) < 0: %s" % bool(j1 < 0))
print("J(1.1647) > 0: %s" % bool(j0 > 0))
ell = (-lam1).exp()
L1 = (1 - ell) / lam1
e1 = arb(1).exp()
rho = (e1 - 1) / L1
c0max = (4 * rho) ** (arb(-1) / 3)
kappa = c0max.sqrt() / 4
show("ell = exp(-1.1648)", ell)
show("L(1, ell) = (1 - ell)/1.1648", L1)
show("L(1, e) = e - 1", e1 - 1)
show("rho = (e - 1)/L(1, ell)", rho)
show("c0max = (4 rho)^(-1/3)", c0max)
show("kappa = sqrt(c0max)/4", kappa)
print("kappa > 0.1660: %s" % bool(kappa > arb("0.1660")))
print(
    "delta_inf = 2 c0^(3/2) sqrt(rho) at c0 = c0max equals 1: %s"
    % (2 * c0max ** (arb(3) / 2) * rho.sqrt()).str(15, radius=False)
)
ok = bool(j1 < 0) and bool(j0 > 0) and bool(kappa > arb("0.1660"))
print("CHECK sharp_depth: %s" % ("ALL THREE STATEMENTS HOLD" if ok else "NOT ALL STATEMENTS HOLD"))
