#!/usr/bin/env python3
"""Independent finite checks and illustrations for the QSD manuscript.

The proof is analytic. Floating-point checks are diagnostics, not a proof of
the infinite-domain theorem. Exact rational identities use fractions.Fraction.
Run with Python 3, numpy, scipy and matplotlib. No network or randomness.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import os
from fractions import Fraction as F
from pathlib import Path

import numpy as np
import scipy
from scipy import linalg, optimize, special

ROOT = Path(__file__).resolve().parent


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


def qmat(n: int, original: bool = False) -> np.ndarray:
    """Killed lattice generator on {1,...,n}; endpoints are 0 and n+1.

    Centred rates: +2 at 1, +1 at 1, -1 at 3, total variance 8.
    Original rates after inverse tilt eta=log 2: +2 at 1/4,
    +1 at 1/2, -1 at 6. Diagonals retain all killing rates.
    """
    rates = {2: .25, 1: .5, -1: 6.} if original else {2: 1., 1: 1., -1: 3.}
    q = -sum(rates.values()) * np.eye(n)
    for z, rate in rates.items():
        for i in range(n):
            j = i + z
            if 0 <= j < n:
                q[i, j] += rate
    return q


def qmat_exact(n: int, original: bool = False) -> list[list[F]]:
    rates = {2: F(1,4), 1: F(1,2), -1: F(6)} if original else {2:F(1),1:F(1),-1:F(3)}
    q = [[F(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        q[i][i] = -sum(rates.values())
        for z,rate in rates.items():
            if 0<=i+z<n: q[i][i+z]+=rate
    return q


def mm(a,b):
    n=len(a)
    return [[sum(a[i][k]*b[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def charpoly_exact(a):
    """Faddeev-LeVerrier, with exact rational arithmetic."""
    n=len(a)
    b=[[F(i==j) for j in range(n)] for i in range(n)]
    coeffs=[F(1)]
    for k in range(1,n+1):
        b=mm(a,b)
        c=-sum(b[i][i] for i in range(n))/k
        coeffs.append(c)
        for i in range(n): b[i][i]+=c
    assert all(x==0 for row in b for x in row)
    return [str(c) for c in coeffs]


def brownian_log_moment(theta, eta, length):
    k=math.pi/length
    c=theta-eta
    return (np.logaddexp(0.,c*length)-np.logaddexp(0.,-eta*length)
            +math.log(eta*eta+k*k)-math.log(c*c+k*k))


def original_brownian_log_hit(length, budget, drift=1., sigma=1.):
    sd=sigma*math.sqrt(budget)
    a=special.log_ndtr((-length-drift*budget)/sd)
    b=-2*drift*length/sigma**2+special.log_ndtr((drift*budget-length)/sd)
    return float(np.logaddexp(a,b))


def source_example():
    def psi(s): return -s+s*s/2-3*s/(1+s)+2*s/(4-s)
    def dpsi(s): return -1+s-3/(1+s)**2+8/(4-s)**2
    def ddpsi(s): return 1+6/(1+s)**3+16/(4-s)**3
    eta=optimize.brentq(dpsi,0.,2.,xtol=1e-14)
    assert abs(psi(2.))<1e-14
    assert abs(dpsi(eta))<1e-12
    return dict(eta=eta,lambda_star=2.,moment_exponent=2.-eta,
                tilted_variance=ddpsi(eta),limiting_killing_rate=-psi(eta),
                eigenvalue_correction=math.pi**2*ddpsi(eta)/2,
                derivative_residual=dpsi(eta))


def run(out: Path, plot: bool):
    out.mkdir(parents=True,exist_ok=True)
    exact=[]
    for n in [2,3,4,5]:
        a=qmat_exact(n); b=qmat_exact(n,True)
        # K0 = H^{-1} K H + 7/4 I, H_ii=2^(i+1).
        for i in range(n):
            for j in range(n):
                transformed=b[i][j]*F(2)**(j-i)+(F(7,4) if i==j else 0)
                assert transformed==a[i][j]
        exact.append(dict(n=n,centred_matrix=[[str(x) for x in row] for row in a],
                          characteristic_polynomial=charpoly_exact(a)))
    eta=math.log(2.)
    ls=math.log((-3+math.sqrt(105))/2)
    lattice=[]
    for n in [15,31,63,127,255]:
        q=qmat(n)
        eig,vec=linalg.eig(q.T)
        k=int(np.argmax(eig.real))
        ev=eig[k]
        p=vec[:,k].real
        if p.sum()<0:p=-p
        assert abs(ev.imag)<1e-10 and p.min()>0
        p/=p.sum()
        rho=-ev.real
        length=n+1
        y=np.arange(1,n+1,dtype=float)
        logp=np.log(p)
        logden=special.logsumexp(logp-eta*y)
        logm=special.logsumexp(logp+(ls-eta)*y)-logden
        sine=np.sin(math.pi*y/length);sine/=sine.sum()
        residual=linalg.norm(q.T@p+rho*p,np.inf)
        assert residual<1e-11
        lattice.append(dict(n=n,length=length,scaled_killing_rate=length**2*rho,
                            limiting_scaled_rate=4*math.pi**2,sine_total_variation=.5*np.abs(p-sine).sum(),
                            log_moment_per_length=logm/length,limiting_moment_rate=ls-eta,
                            eigen_residual=residual))
    # These are convergence diagnostics; no finite test proves the limit.
    assert lattice[-1]['sine_total_variation']<lattice[0]['sine_total_variation']
    assert abs(lattice[-1]['scaled_killing_rate']-4*math.pi**2)<abs(lattice[0]['scaled_killing_rate']-4*math.pi**2)
    brownian=[]
    for drift,sigma in [(1.,1.),(.5,1.3),(2.,.7)]:
        a=drift/sigma**2
        for length in [1.,3.,10.,30.]:
            logm=brownian_log_moment(2*a,a,length)
            assert abs(logm-a*length)<5e-13
            rho=drift**2/(2*sigma**2)+sigma**2*math.pi**2/(2*length**2)
            assert abs(length**2*(rho-drift**2/(2*sigma**2))-sigma**2*math.pi**2/2)<2e-12
            brownian.append(dict(drift=drift,sigma=sigma,length=length,log_moment_residual=logm-a*length))
        # At theta=eta the numerator integral is 2L/pi; this tests the
        # critical factor of two independently of the lambda_star identity.
        length=1e5/a
        critical=math.exp(brownian_log_moment(a,a,length))/length**2
        predicted=2*a*a/math.pi**2
        assert abs(critical/predicted-1)<2e-9
    saturation=[]
    for length in [4.,8.,12.,16.,24.,32.]:
        budget=math.exp(1.5*length)
        rho=.5*(1+math.pi**2/length**2)
        hazard=rho*budget/(1+math.exp(length))
        log_pr=math.log(-math.expm1(-hazard))
        log_po=original_brownian_log_hit(length,budget)
        logratio=log_pr-log_po-2.5*length
        saturation.append(dict(length=length,restart_probability=math.exp(log_pr),
                               source_normalised_ratio=math.exp(logratio),
                               predicted_ratio=math.exp(-length/2),
                               budget=budget))
    assert saturation[-1]['source_normalised_ratio']<1e-6
    data=dict(status='PASS',scope='finite arithmetic checks and numerical diagnostics only; analytic proof is separate',
              exact_rational_lattice_checks=exact,brownian_identities=brownian,
              lattice_diagnostics=lattice,source_example=source_example(),saturation_counterexample=saturation)
    (out/'verification.json').write_text(json.dumps(clean(data),indent=2,sort_keys=True)+'\n')
    for name,rows in [('lattice_diagnostics',lattice),('saturation_counterexample',saturation)]:
        with (out/(name+'.csv')).open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(clean(rows))
    if plot: make_plots(out,lattice)
    print(json.dumps(dict(status=data['status'],source_example=data['source_example'],
                          last_lattice=lattice[-1],last_saturation=saturation[-1]),indent=2))


def make_plots(out,lattice):
    import tempfile
    os.environ.setdefault('MPLCONFIGDIR',tempfile.mkdtemp(prefix='mplcfg_'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,
                         'figure.facecolor':'#fbf8f2','axes.facecolor':'#fbf8f2','savefig.facecolor':'#fbf8f2'})
    plum='#73365b';ochre='#b97726';teal='#26746d';grey='#5e5964'
    fig,axs=plt.subplots(1,2,figsize=(12.8,5.6),layout='constrained')
    q=np.linspace(-.5,2.5,501)
    eta=1.
    for length,c in [(8,ochre),(32,teal),(128,plum)]:
        axs[0].plot(q,[brownian_log_moment(x,eta,length)/length for x in q],color=c,label=f'L = {length}')
    axs[0].plot(q,np.maximum(q-1,0),'--',color=grey,label='Exact limiting exponent')
    axs[0].set(xlabel=r'Moment parameter $\theta$',ylabel=r'$L^{-1}\log M_L(\theta)$',title='Moment growth has a threshold')
    axs[0].legend(frameon=False,fontsize=9)
    length=20
    b=np.linspace(.4,1.8,401)
    rho=.5*(1+math.pi**2/length**2)
    hazard=rho*np.exp(b*length)/(1+math.exp(length))
    prob=-np.expm1(-hazard)
    axs[1].plot(b,prob,color=plum,label='Exact restart probability')
    axs[1].plot(b,np.minimum(hazard,2.0),color=ochre,ls='--',label='Linear approximation (clipped at 2)')
    axs[1].axhline(1,color=grey,lw=.8)
    axs[1].axvline(1,color=teal,ls=':',lw=1)
    axs[1].set(xlabel=r'Budget exponent $b$, with $B=e^{bL}$',ylabel='Success probability',ylim=(0,1.55),title='The restart benefit saturates')
    axs[1].legend(frameon=False,fontsize=9,loc='upper left')
    fig.suptitle('Exact Brownian benchmark: drift −1, variance rate 1',fontsize=15)
    fig.savefig(out/'moment_and_saturation.png',dpi=180,metadata={'Software':None})
    fig.savefig(out/'moment_and_saturation.pdf',metadata={'CreationDate':None,'Producer':None,'Creator':None})
    plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(12.8,5.2),layout='constrained')
    for n,c in [(15,ochre),(63,teal),(255,plum)]:
        eig,vec=linalg.eig(qmat(n).T);p=vec[:,np.argmax(eig.real)].real
        if p.sum()<0:p=-p
        p/=p.sum();length=n+1
        axs[0].plot(np.arange(1,n+1)/length,length*p,color=c,label=f'L = {length}')
    u=np.linspace(0,1,401)
    axs[0].plot(u,math.pi/2*np.sin(math.pi*u),ls='--',color=grey,label='Sine density')
    axs[0].set(xlabel='Scaled position u',ylabel='Scaled quasi-stationary mass',title='A non-reversible lattice check')
    axs[0].legend(frameon=False,fontsize=9)
    axs[1].plot([x['length'] for x in lattice],[x['scaled_killing_rate'] for x in lattice],'o-',color=plum)
    axs[1].axhline(4*math.pi**2,color=teal,ls='--',label=r'$\pi^2 v/2=4\pi^2$')
    axs[1].set(xlabel='Interval length L',ylabel=r'$L^2\alpha_L$',title='The decay rate approaches its prediction',xscale='log')
    axs[1].legend(frameon=False,fontsize=9)
    fig.savefig(out/'lattice_scaling.png',dpi=180,metadata={'Software':None})
    fig.savefig(out/'lattice_scaling.pdf',metadata={'CreationDate':None,'Producer':None,'Creator':None})
    plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,default=ROOT/'checks');parser.add_argument('--no-plots',action='store_true')
    args=parser.parse_args();run(args.out,not args.no_plots)
