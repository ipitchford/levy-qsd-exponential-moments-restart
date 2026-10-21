#!/usr/bin/env python3
"""Independent constant and finite-interval checks for levy_qsd_restart.tex.

This file does not import the manuscript's verification scripts.  The exact
source model is reduced independently to a four-mode boundary-value problem.
The numerical checks are consistency checks, not substitutes for the proofs.
"""

from decimal import Decimal, localcontext
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq


SCRIPT_DIR = Path(__file__).resolve().parent
OUT = SCRIPT_DIR/'checks'/'manuscript_constant_check.json'


ROUNDOFF_KEYS = ('residual', 'minimum_unnormalized_density', 'identity_error')


def clean(o, key=''):
    """Deterministic JSON/CSV values: floats to 10 significant digits;
    round-off-level residuals replaced by a threshold string."""
    if isinstance(o, dict):
        return {k: clean(v, k) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v, key) for v in o]
    if isinstance(o, float) or type(o).__name__.startswith('float'):
        o = float(o)
        if any(t in key for t in ROUNDOFF_KEYS):
            return '<1e-10' if abs(o) < 1e-10 else float(f'{o:.3g}')
        return float(f'{o:.10g}')
    return o


def decimal_constants():
    """Compute the constants at 85 digits, using three formulas for R*."""
    with localcontext() as ctx:
        ctx.prec = 85
        D = Decimal
        one, two, three, four = map(D, ('1', '2', '3', '4'))
        eta = D('.94')
        for _ in range(14):
            f = eta-one-three/(eta+one)**2+D(8)/(four-eta)**2
            fp = one+D(6)/(eta+one)**3+D(16)/(four-eta)**3
            eta -= f/fp
        psi = lambda s: -s+s*s/two-three*s/(one+s)+two*s/(four-s)
        c = -psi(eta)
        v = one+D(6)/(one+eta)**3+D(16)/(four-eta)**3
        d, e = four-eta, one+eta
        disc = ((D(5)-two*eta)**2+D(32)*c/eta**2).sqrt()
        rminus = (-(D(5)-two*eta)-disc)/two
        rplus = (-(D(5)-two*eta)+disc)/two
        A, F = -eta-rminus, eta+rplus
        # Paired ladder-clock normalisation c_+=1, c_-=1/2.
        phip = lambda q: q*(q+A)/(q+d)
        phim = lambda q: q*(q+F)/(two*(q+e))
        mp, mm = A/d, F/(two*e)
        Gp = lambda s: d*(A-s)/(A*(d-s))
        Gm = lambda s: e*(F-s)/(F*(e-s))
        b = two-eta
        C_overshoot = Gp(b)/Gm(eta)
        # pi^2 cancels from R(b)/L(eta).
        C_ladder = (mm/(b*phim(b)))/(mp/(eta*phip(eta)))
        R_overshoot = c*Gp(-eta)/Gm(eta)
        R_ladder = (mm/mp)*phip(eta)**2
        # pi^2 cancels from pi^2*v*G_+(-eta)/(4*L(eta)).
        R_endpoint = v*eta*phip(eta)*Gp(-eta)/(two*mp)
        # Closed expressions use d+eta=4, e-eta=1, and d-b=2.
        C_closed = F*d*(A+eta-two)/(two*A*e*(F-eta))
        R_closed = eta**2*(eta+A)**2*F*d/(D(32)*A*e)
        errs = {
            'psi_lambda_star': abs(psi(two)),
            'eta_stationarity': abs(eta-one-three/(eta+one)**2+D(8)/(four-eta)**2),
            'ladder_mean_product': abs(mp*mm-v/two),
            'C_two_formulas': abs(C_overshoot-C_ladder),
            'R_overshoot_vs_ladder': abs(R_overshoot-R_ladder),
            'R_endpoint_vs_ladder': abs(R_endpoint-R_ladder),
            'C_closed_vs_overshoot': abs(C_closed-C_overshoot),
            'R_closed_vs_ladder': abs(R_closed-R_ladder),
        }
        # Clear all denominators in the rational centred exponent.  The
        # coefficient vector is obtained directly from the original psi,
        # using polynomial arithmetic on Decimal coefficients.
        def padd(a, b):
            a, b = list(a), list(b)
            return [(a[i] if i<len(a) else D(0))+
                    (b[i] if i<len(b) else D(0)) for i in range(max(len(a),len(b)))]
        def pmul(a, b):
            out = [D(0)]*(len(a)+len(b)-1)
            for i, x in enumerate(a):
                for j, y in enumerate(b):
                    out[i+j] += x*y
            return out
        scale = lambda a, k: [k*x for x in a]
        z = [eta, one]
        # 2(d-s)(e+s)*(psi(eta+s)+c).
        denom = pmul([d,-one],[e,one])
        poly = padd(scale(pmul(denom, padd(padd(scale(z,-one),scale(pmul(z,z),one/two)),[c])),two),
                    padd(scale(pmul(z,[d,-one]),-D(6)),scale(pmul(z,[e,one]),D(4))))
        factored = [D(0),D(0),A*F,A-F,-one]
        errs['factorization_coefficient_max'] = max(abs(x-y) for x,y in zip(poly,factored))
        assert A>d and F>e
        assert max(errs.values()) < D('1e-77'), errs
        return {
            'eta': str(eta), 'c': str(c), 'v': str(v),
            'd':str(d), 'e':str(e), 'A':str(A), 'D':str(F),
            'C_star':str(C_overshoot), 'R_star':str(R_ladder),
            'G_plus_b':str(Gp(b)), 'G_minus_eta':str(Gm(eta)),
            'G_plus_minus_eta':str(Gp(-eta)),
            'max_decimal_error':str(max(errs.values())),
            'individual_errors':{k:str(z) for k,z in errs.items()},
        }


def finite_interval(L, eta, c, v):
    """Solve the exact differential-integral eigenproblem in a stable basis.

    q''/2+q'+(rho-5)q+8U+3V=0,
    U'=q-4U, V'=V-q,
    q(0)=q(L)=U(0)=V(L)=0.
    Its characteristic polynomial is
    r^4+5r^3+(2rho-8)r^2+(6rho-28)r-8rho.

    Work with h(y)=exp(eta*y)*q(y). Anchor the positive fast mode at L,
    the negative fast mode at zero, and use the real/imaginary slow modes.
    """
    d, e = 4-eta, 1+eta

    def modes(alpha):
        rho = c+alpha
        # Expand directly in centred s=r+eta. Using psi'(eta)=0 avoids
        # cancellation in the constant/linear coefficients near alpha=0.
        roots = np.roots([1,5-4*eta,2*alpha-v*d*e,
                          2*alpha*(d-e),-2*alpha*d*e])
        real = sorted([z.real for z in roots if abs(z.imag)<1e-9])
        slow = [z for z in roots if z.imag>1e-9]
        if len(real)!=2 or len(slow)!=1:
            raise RuntimeError(('unexpected roots',L,alpha,roots))
        sn, sp = real
        s = slow[0]
        en, ep, ec = np.exp(sn*L), np.exp(-sp*L), np.exp(s*L)
        M = np.array([
            [1, ep, 1, 0],
            [en, 1, ec.real, ec.imag],
            [1/(sn+d), ep/(sp+d), (1/(s+d)).real, (1/(s+d)).imag],
            [en/(e-sn), 1/(e-sp), (ec/(e-s)).real, (ec/(e-s)).imag],
        ])
        return sn, sp, s, M

    alpha_scale = np.pi**2*v/(2*L**2)
    grid = alpha_scale*np.geomspace(.12,2.0,180)
    dets = [np.linalg.det(modes(x)[3]) for x in grid]
    brackets = [(x,y) for x,y,fx,fy in zip(grid[:-1],grid[1:],dets[:-1],dets[1:]) if fx*fy<0]
    assert brackets, ('no eigenvalue bracket',L)
    alpha = brentq(lambda a:np.linalg.det(modes(a)[3]),*brackets[0],xtol=1e-16,rtol=1e-13)
    sn, sp, s, M = modes(alpha)
    _, sv, vh = np.linalg.svd(M)
    coeff = vh[-1]

    def basis(y):
        z = np.exp(s*y)
        return np.array([np.exp(sn*y),np.exp(sp*(y-L)),z.real,z.imag])

    if coeff@basis(L/2)<0:
        coeff = -coeff

    def left_integral(b):
        slow = -np.expm1((s-b)*L)/(b-s)
        terms = np.array([
            -np.expm1((sn-b)*L)/(b-sn),
            (np.exp(-b*L)-np.exp(-sp*L))/(sp-b),
            slow.real, slow.imag,
        ])
        return float(coeff@terms)

    def right_integral(b):
        slow = (np.exp(s*L)-np.exp(-b*L))/(s+b)
        terms = np.array([
            (np.exp(sn*L)-np.exp(-b*L))/(sn+b),
            -np.expm1(-(sp+b)*L)/(sp+b),
            slow.real,slow.imag,
        ])
        return float(coeff@terms)

    def derivative(y):
        slow = s*np.exp(s*y)
        terms = np.array([sn*np.exp(sn*y),sp*np.exp(sp*(y-L)),slow.real,slow.imag])
        return float(coeff@terms)

    Z = left_integral(eta)
    total = left_integral(0)
    C_L = right_integral(2-eta)/Z
    # Actual upper exit flux, combining Gaussian boundary flux and positive
    # exponential jump exit intensity.  This does not use the proposed R*.
    R_L = (-.5*derivative(L)+2*right_integral(d))/Z
    lower_rate = (.5*derivative(0)+3*left_integral(eta+1))/Z
    upper_rate = np.exp(-eta*L)*R_L
    total_rate_err = abs(lower_rate+upper_rate-(c+alpha))
    vals = np.array([coeff@basis(y) for y in np.linspace(0,L,2001)])
    boundary_error = float(np.linalg.norm(M@coeff,np.inf))
    assert Z>0 and total>0
    assert vals.min()>-2e-10
    assert boundary_error<2e-12
    assert total_rate_err<1e-9, (L, lower_rate, upper_rate, c+alpha, total_rate_err)
    return {
        'L':L, 'alpha_L2':float(alpha*L**2),
        'scaled_moment_C_L':C_L, 'scaled_restart_rate_R_L':R_L,
        'L2_left_mass':L**2*Z/total,
        'boundary_residual':boundary_error,
        'total_exit_rate_residual':total_rate_err,
        'minimum_unnormalized_density':float(vals.min()),
        'second_smallest_singular_value':float(sv[-2]),
    }


def brownian_checks():
    a = .8
    L = 180.0
    k = np.pi/L
    moment = lambda theta: (1+np.exp((theta-a)*L))/(1+np.exp(-a*L))*(a*a+k*k)/((theta-a)**2+k*k)
    critical = moment(a)/L**2
    expected_critical = 2*a*a/np.pi**2
    exact_critical = 2*(a*a+k*k)/(np.pi**2*(1+np.exp(-a*L)))
    assert abs(critical-exact_critical)<1e-15
    # Include a negative moment parameter and both positive regimes.
    checks = []
    for theta in (-.7,.2,1.3):
        observed = moment(theta)*np.exp(-max(theta-a,0)*L)
        expected = a*a/(theta-a)**2
        checks.append({'theta':theta,'scaled_moment':observed,'limit':expected})
        assert abs(observed/expected-1)<.002
    assert abs(critical/expected_critical-1)<.002
    return {'a':a,'L':L,'critical_scaled_moment':critical,
            'critical_limit_2a2_over_pi2':expected_critical,
            'critical_exact_identity_error':abs(critical-exact_critical),
            'other_regimes':checks}


def main():
    constants = decimal_constants()
    eta,c,v = (float(constants[k]) for k in ('eta','c','v'))
    intervals = [finite_interval(L,eta,c,v) for L in (12,20,40,80,160,320,640,1280)]
    C, R = float(constants['C_star']), float(constants['R_star'])
    for row in intervals:
        row['C_relative_error'] = row['scaled_moment_C_L']/C-1
        row['R_relative_error'] = row['scaled_restart_rate_R_L']/R-1
    assert abs(intervals[-1]['C_relative_error'])<5e-4
    assert abs(intervals[-1]['R_relative_error'])<5e-4
    ap = argparse.ArgumentParser(); ap.add_argument('--out', type=Path, default=OUT.parent)
    out = ap.parse_args().out; out.mkdir(parents=True, exist_ok=True)
    rename = {'c':'kappa_inf','d':'beta_plus','e':'beta_minus','A':'a_plus','D':'a_minus'}
    result = {'decimal_constants':{rename.get(k,k):v for k,v in constants.items()},'finite_interval':intervals,
              'brownian':brownian_checks(),
              'status':'all checks passed'}
    (out/'manuscript_constant_check.json').write_text(json.dumps(clean(result),indent=2,sort_keys=True)+'\n')
    print('C* =',constants['C_star'])
    print('R* =',constants['R_star'])
    print('L      L^2 alpha       C_L             R_L              flux residual')
    for x in intervals:
        print(f"{x['L']:3d}  {x['alpha_L2']:.11f}  {x['scaled_moment_C_L']:.11f}  {x['scaled_restart_rate_R_L']:.11f}  {x['total_exit_rate_residual']:.2e}")
    print('All checks passed; wrote',out/'manuscript_constant_check.json')


if __name__ == '__main__':
    main()
