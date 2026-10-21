"""Algebraic evaluation of a candidate QSD constant; no asymptotic theorem assumed.

Run with Python 3, standard library only. The two formulas arise from Wiener-Hopf
factors and formal boundary profiles. Agreement checks the algebra, not the
probabilistic passage from finite intervals to stationary overshoots.
"""

import argparse
from decimal import Decimal as D, localcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def calculate(out):
    lines = []
    with localcontext() as ctx:
        ctx.prec = 70

        def psi(s):
            return -s + s * s / 2 - 3 * s / (s + 1) + 2 * s / (4 - s)

        def dpsi(s):
            return -1 + s - 3 / (s + 1) ** 2 + 8 / (4 - s) ** 2

        lo, hi = D(0), D(2)
        for _ in range(240):
            mid = (lo + hi) / 2
            if dpsi(mid) > 0:
                hi = mid
            else:
                lo = mid
        eta = (lo + hi) / 2
        rho0 = -psi(eta)
        d, e, b = 4 - eta, 1 + eta, 2 - eta
        p = 5 - 2 * eta
        disc = (p * p + 32 * rho0 / (eta * eta)).sqrt()
        rminus, rplus = (-p - disc) / 2, (-p + disc) / 2
        A, B = -eta - rminus, eta + rplus

        def Gplus(s):
            return d * (A - s) / (A * (d - s))

        def Gminus(s):
            return e * (B - s) / (B * (e - s))

        ratio = Gplus(b) / Gminus(eta)
        ellminus, ellplus = 1 / d - 1 / A, 1 / e - 1 / B
        left = 1 / eta**2 + ellminus * (1 / eta - 1 / (eta + A))
        right = 1 / b**2 + ellplus * (1 / b - 1 / (b + B))
        boundary_ratio = right / left
        error = abs(ratio - boundary_ratio)
        factor_error = D(0)
        for s in (D('-0.5'), D('0.25'), D('0.75'), D('1.25')):
            kappa = psi(eta + s) - psi(eta)
            factored = s * s * (A - s) * (B + s) / (2 * (d - s) * (e + s))
            factor_error = max(factor_error, abs(kappa - factored))

        for name, val in (
            ('eta', eta), ('kappa_inf', rho0), ('a_plus', A), ('a_minus', B),
            ('Gplus(2-eta)', Gplus(b)), ('Gminus(eta)', Gminus(eta)),
            ('prefactor from overshoots', ratio),
            ('prefactor from formal boundary profiles', boundary_ratio),
            ('absolute difference', error), ('max factorization residual', factor_error),
        ):
            lines.append(f'{name}: {val}')
        assert error < D('1e-60')
        assert factor_error < D('1e-60')
    out.mkdir(parents=True, exist_ok=True)
    (out / 'prefactor_decimal.txt').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, default=ROOT / 'checks')
    calculate(ap.parse_args().out)
