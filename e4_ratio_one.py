#!/usr/bin/env python3
"""Certified enclosure showing C_* < 1 in the source's two-sided exponential model.

Label: CERTIFIED.  Method:
 1. eta is isolated by 200 steps of bisection on [9/10, 1] in exact rational
    arithmetic (fractions.Fraction), using the sign of
        P(s) = (s-1)(1+s)^2(4-s)^2 - 3(4-s)^2 + 8(1+s)^2,
    which equals (1+s)^2(4-s)^2 psi'(s) and so has the sign of psi'(s); psi' is
    strictly increasing on (-1,4), so eta is its only zero there.
 2. All further quantities (kappa_inf, a_+-, beta_+-, G_+-, C_*, R_*) are evaluated
    in mpmath.iv interval arithmetic (outward rounding at 80 digits) starting
    from the rational isolating interval for eta.
 3. Every printed endpoint is rounded OUTWARD (lower endpoints down, upper
    endpoints up) to 40 decimal digits, so the printed intervals are themselves
    valid enclosures.
 4. The four-point check of the rational factorisation of chi is a consistency
    test only; the factorisation itself is exact algebra (manuscript Sec. 10.2).  The structural sign results (C_* < 1 for spectrally
negative and C_* > 1 for spectrally positive processes) are proved in the
manuscript; this script only certifies the numerical value for the source model.

Model: psi(s) = -s + s^2/2 - 3s/(1+s) + 2s/(4-s), lambda_* = 2.
Usage: python3 e4_ratio_one.py [--out DIR]
"""
import argparse
import json
from fractions import Fraction as Fr
from pathlib import Path

from mpmath import iv, mp

ROOT = Path(__file__).resolve().parent


def dpsi(s):
    """Exact psi'(s) for rational s."""
    return -1 + s - Fr(3) / (1 + s) ** 2 + Fr(8) / (4 - s) ** 2


def P(s):
    return (s - 1) * (1 + s) ** 2 * (4 - s) ** 2 - 3 * (4 - s) ** 2 + 8 * (1 + s) ** 2


def bracket_eta(steps=200):
    """Bisection with exact rationals; psi' is strictly increasing on (-1,4)."""
    lo, hi = Fr(9, 10), Fr(1)
    assert dpsi(lo) < 0 < dpsi(hi)
    for _ in range(steps):
        mid = (lo + hi) / 2
        assert (P(mid) < 0) == (dpsi(mid) < 0)
        if P(mid) < 0:
            lo = mid
        else:
            hi = mid
    assert dpsi(lo) < 0 < dpsi(hi)
    return lo, hi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, default=ROOT / 'checks')
    out = ap.parse_args().out
    out.mkdir(parents=True, exist_ok=True)
    iv.dps = 80
    mp.dps = 80
    lo, hi = bracket_eta()
    eta = iv.mpf([iv.mpf(lo.numerator) / lo.denominator, iv.mpf(hi.numerator) / hi.denominator])
    psi = lambda s: -s + s * s / 2 - 3 * s / (1 + s) + 2 * s / (4 - s)
    kinf = -psi(eta)
    bp, bm = 4 - eta, 1 + eta
    disc = iv.sqrt((5 - 2 * eta) ** 2 + 32 * kinf / eta ** 2)
    rminus = (-(5 - 2 * eta) - disc) / 2
    rplus = (-(5 - 2 * eta) + disc) / 2
    a_p, a_m = -eta - rminus, eta + rplus
    Gp = lambda s: bp * (a_p - s) / (a_p * (bp - s))
    Gm = lambda s: bm * (a_m - s) / (a_m * (bm - s))
    C = Gp(2 - eta) / Gm(eta)
    R = kinf * Gp(-eta) / Gm(eta)
    # factorisation check: chi(s) - s^2(a_+-s)(a_-+s)/(2(beta_+-s)(beta_-+s)) encloses 0
    fact_ok = True
    for t in ('-0.5', '0.25', '0.75', '1.25'):
        x = iv.mpf(t)
        chi = psi(eta + x) - psi(eta)
        fac = x * x * (a_p - x) * (a_m + x) / (2 * (bp - x) * (bm + x))
        d = chi - fac
        fact_ok &= bool(d.a <= 0 <= d.b) and float(mp.mpf(d.b) - mp.mpf(d.a)) < 1e-50
    assert fact_ok
    assert a_p.a > bp.b and a_m.a > bm.b      # orderings a_+ > beta_+, a_- > beta_-
    assert C.b < 1                            # C_* < 1, certified
    assert float(mp.mpf(C.b) - mp.mpf(C.a)) < 1e-50
    def outw(x, n=40):
        q = mp.mpf(10) ** n
        return [mp.nstr(mp.floor(mp.mpf(x.a) * q) / q, n + 5), mp.nstr(mp.ceil(mp.mpf(x.b) * q) / q, n + 5)]
    s = lambda x, n=30: mp.nstr(mp.mpf(x.mid), n)
    res = {
        'label': 'certified: exact rational bracket of eta, outward-rounded interval arithmetic',
        'model': 'psi(s) = -s + s^2/2 - 3s/(1+s) + 2s/(4-s); lambda_* = 2',
        'eta_bisection_steps': 200,
        'eta': s(eta), 'kappa_inf': s(kinf), 'a_plus': s(a_p), 'a_minus': s(a_m),
        'beta_plus': s(bp), 'beta_minus': s(bm),
        'factorisation_enclosures_contain_zero': fact_ok,
        'C_star_outward_interval_40d': outw(C),
        'R_star_outward_interval_40d': outw(R),
        'eta_outward_interval_45d': outw(eta, 45),
        'kappa_inf_outward_interval_40d': outw(kinf),
        'a_plus_outward_interval_40d': outw(a_p), 'a_minus_outward_interval_40d': outw(a_m),
        'beta_plus_outward_interval_40d': outw(bp), 'beta_minus_outward_interval_40d': outw(bm),
        'G_plus_at_lambda_star_minus_eta_outward_40d': outw(Gp(2 - eta)),
        'G_minus_at_eta_outward_40d': outw(Gm(eta)),
        'G_plus_at_minus_eta_outward_40d': outw(Gp(-eta)),
        'bisection_polynomial': '(s-1)(1+s)^2(4-s)^2 - 3(4-s)^2 + 8(1+s)^2 on [9/10,1], 200 exact rational steps',
        'factorisation_check_kind': 'consistency test at s = -0.5, 0.25, 0.75, 1.25',
        'C_star_enclosure_width_below': '1e-50',
        'C_star_strictly_less_than_one': True,
        'G_plus_at_lambda_star_minus_eta': s(Gp(2 - eta)),
        'G_minus_at_eta': s(Gm(eta)),
        'brownian_reference': 'C_* = 1 exactly (no jumps: G_+ = G_- = 1)',
    }
    (out / 'e4_ratio_one.json').write_text(json.dumps(res, indent=2, sort_keys=True) + '\n')
    print('C_* in', res['C_star_outward_interval_40d'], '(< 1 certified)')


if __name__ == '__main__':
    main()
