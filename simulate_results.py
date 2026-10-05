#!/usr/bin/env python3
"""Reproduce simulation results and figures for the Optical Review manuscript.

All reported quantities are generated from this Monte Carlo model. No experimental
measurements are used. The code uses fixed seeds and writes raw CSV files.
"""
from __future__ import annotations
import json, math, platform, sys
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy.stats import ncx2
from scipy.integrate import quad
from scipy.special import i0e
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'
FIG = ROOT / 'figures'
DATA.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

CONFIG = {
    'wavelength_m': 1550e-9,
    'link_distance_m': 1000.0,
    'waist_contribution_m': 0.010,
    'aperture_radius_m': 0.025,
    'ring_radius_m': 0.075,
    'pointing_mean_m': [0.030, -0.015],
    'pointing_sigma_major_m': 0.035,
    'pointing_sigma_minor_m': 0.015,
    'pointing_ellipse_angle_deg': 30.0,
    'sigma_ln': 0.38,
    'scintillation_corr_length_m': 0.055,
    'tail_alpha': 0.05,
    'blocked_attenuation_uniform': [0.0, 0.12],
    'micro_obscuration_probability': 0.01,
    'micro_obscuration_uniform': [0.0, 0.4],
    'optimization_bounds': {'x0_m': [-0.10, 0.10], 'y0_m': [-0.10, 0.10], 'theta_urad': [60.0, 350.0]},
    'test_samples_per_state': 60000,
    'repeat_trials': 10,
    'repeat_samples_per_state': 20000,
    'base_seed': 20260919,
}

DESIGNS = {
    'Centroid-fixed': dict(x0=0.0, y0=0.0, theta_urad=100.0),
    'Bias-cancelled': dict(x0=-0.030, y0=0.015, theta_urad=100.0),
    'Divergence-only robust': dict(x0=-0.030, y0=0.015, theta_urad=126.33756146),
    'Center-only robust': dict(x0=-0.03958562, y0=0.01155312, theta_urad=100.0),
    'Equal-state mean joint': dict(x0=-0.0295457768, y0=0.0149357509, theta_urad=60.0),
    'Proposed robust joint': dict(x0=-0.0390755073, y0=0.0119549665, theta_urad=105.707386),
}

STATE_BLOCK = {
    'none': (),
    'left': (3, 4),
    'upper_right': (1, 2),
    'lower': (5, 6),
}
STATES = list(STATE_BLOCK.keys())

# receiver centers: center + six vertices at 0,60,... deg
angles = np.deg2rad(np.arange(0, 360, 60.0))
AP_POS = np.vstack(([0.0, 0.0], np.column_stack((CONFIG['ring_radius_m']*np.cos(angles), CONFIG['ring_radius_m']*np.sin(angles)))))
Dmat = np.linalg.norm(AP_POS[:,None,:]-AP_POS[None,:,:], axis=2)
C = np.exp(- (Dmat / CONFIG['scintillation_corr_length_m'])**(5.0/3.0))
# numerical safety
C += 1e-12*np.eye(7)
Lcorr = np.linalg.cholesky(C)


def cov_pointing(sig_major=None, sig_minor=None, phi_deg=None):
    sig_major = CONFIG['pointing_sigma_major_m'] if sig_major is None else sig_major
    sig_minor = CONFIG['pointing_sigma_minor_m'] if sig_minor is None else sig_minor
    phi = math.radians(CONFIG['pointing_ellipse_angle_deg'] if phi_deg is None else phi_deg)
    R = np.array([[math.cos(phi), -math.sin(phi)],[math.sin(phi), math.cos(phi)]])
    return R @ np.diag([sig_major**2, sig_minor**2]) @ R.T


def beam_radius(theta_urad):
    theta = theta_urad * 1e-6
    return math.sqrt(CONFIG['waist_contribution_m']**2 + (theta*CONFIG['link_distance_m'])**2)


def aperture_fraction(d, w):
    """Exact encircled fraction for a decentered normalized 2D Gaussian."""
    d = np.asarray(d)
    return ncx2.cdf((2*CONFIG['aperture_radius_m']/w)**2, df=2, nc=(2*d/w)**2)


def lower_tail_mean(x, alpha=0.05):
    x = np.asarray(x)
    k = max(1, int(math.ceil(alpha*len(x))))
    # exact sample lower-tail conditional mean (average k smallest observations)
    return np.partition(x, k-1)[:k].mean()


def generate_randoms(n, seed, sig_major=None, sigma_ln=None):
    rng = np.random.default_rng(seed)
    Sigma = cov_pointing(sig_major=sig_major)
    p = rng.multivariate_normal(np.array(CONFIG['pointing_mean_m']), Sigma, size=n)
    z = rng.standard_normal((n,7)) @ Lcorr.T
    sln = CONFIG['sigma_ln'] if sigma_ln is None else sigma_ln
    h = np.exp(sln*z - 0.5*sln**2)
    u_block = rng.random((n,7))
    u_micro = rng.random((n,7))
    u_micro_att = rng.random((n,7))
    return p, h, u_block, u_micro, u_micro_att


def collected(design, state, rnd):
    p,h,u_block,u_micro,u_micro_att = rnd
    w = beam_radius(design['theta_urad'])
    center = p + np.array([design['x0'], design['y0']])
    d = np.linalg.norm(center[:,None,:]-AP_POS[None,:,:], axis=2)
    g = aperture_fraction(d,w)
    b = np.ones_like(g)
    for idx in STATE_BLOCK[state]:
        lo,hi = CONFIG['blocked_attenuation_uniform']
        b[:,idx] = lo + (hi-lo)*u_block[:,idx]
    # low-rate independent micro-obscuration stressor
    event = u_micro < CONFIG['micro_obscuration_probability']
    lo,hi = CONFIG['micro_obscuration_uniform']
    micro_att = lo + (hi-lo)*u_micro_att
    b = np.where(event, b*micro_att, b)
    return np.sum(b*h*g, axis=1)


def metrics(x):
    return {
        'mean_eta': float(np.mean(x)),
        'std_eta': float(np.std(x, ddof=1)),
        'lcvar5': float(lower_tail_mean(x, 0.05)),
        'p1': float(np.quantile(x,0.01)),
        'p5': float(np.quantile(x,0.05)),
        'p10': float(np.quantile(x,0.10)),
        'outage_lt_0p15': float(np.mean(x<0.15)),
        'outage_lt_0p20': float(np.mean(x<0.20)),
    }


def eval_design(design, n=60000, seed=20260919, sig_major=None, sigma_ln=None):
    rows=[]; arrays={}
    for si,state in enumerate(STATES):
        rnd = generate_randoms(n, seed+1009*si, sig_major=sig_major, sigma_ln=sigma_ln)
        x=collected(design,state,rnd)
        arrays[state]=x
        m=metrics(x); m['state']=state
        rows.append(m)
    return pd.DataFrame(rows), arrays


def quadrature_fraction(d,w):
    # Radial integration after angular integration; i0e avoids overflow.
    a=CONFIG['aperture_radius_m']
    def integrand(r):
        x=4*r*d/(w*w)
        # exp(-2(r^2+d^2)/w^2) I0(x) = exp(-2(r-d)^2/w^2) i0e(x)
        return (4*r/(w*w))*math.exp(-2*(r-d)*(r-d)/(w*w))*float(i0e(x))
    return quad(integrand,0,a,epsabs=1e-12,epsrel=1e-11,limit=200)[0]


def make_fig1():
    fig,ax=plt.subplots(figsize=(6.2,5.3))
    th=np.linspace(0,2*np.pi,300)
    for i,(x,y) in enumerate(AP_POS):
        ax.plot(x+CONFIG['aperture_radius_m']*np.cos(th), y+CONFIG['aperture_radius_m']*np.sin(th), lw=1.2)
        ax.text(x,y,str(i),ha='center',va='center',fontsize=8)
    # pointing covariance 2-sigma ellipse about biased mean with proposed static offset shown separately
    vals,vecs=np.linalg.eigh(cov_pointing())
    pts=(vecs @ (2*np.sqrt(vals)[:,None]*np.vstack((np.cos(th),np.sin(th))))).T + np.array(CONFIG['pointing_mean_m'])
    ax.plot(pts[:,0],pts[:,1],'--',lw=1.1,label='2-sigma pointing ellipse')
    d=DESIGNS['Proposed robust joint']; w=beam_radius(d['theta_urad'])
    bc=np.array(CONFIG['pointing_mean_m'])+np.array([d['x0'],d['y0']])
    ax.plot(bc[0]+w*np.cos(th),bc[1]+w*np.sin(th),':',lw=1.2,label='Proposed 1/e^2 footprint at mean pointing')
    ax.scatter([0],[0],marker='+',s=50,label='Array centroid')
    ax.scatter([bc[0]],[bc[1]],marker='x',s=45,label='Mean proposed footprint center')
    ax.set_xlabel('Receiver-plane x (m)'); ax.set_ylabel('Receiver-plane y (m)')
    ax.axis('equal'); ax.grid(alpha=.2); ax.legend(fontsize=7,loc='upper right')
    fig.tight_layout();
    for ext in ['png','pdf']:
        fig.savefig(FIG/f'Fig1_system_geometry.{ext}',dpi=600 if ext=='png' else None,bbox_inches='tight')
    plt.close(fig)


def make_ecdf_fig(arr_mean, arr_rob):
    fig,axs=plt.subplots(2,2,figsize=(6.8,5.6),sharex=True,sharey=True)
    for ax,state in zip(axs.ravel(),STATES):
        for label,arr in [('Equal-state mean joint',arr_mean[state]),('Proposed robust joint',arr_rob[state])]:
            xs=np.sort(arr); ys=np.arange(1,len(xs)+1)/len(xs)
            # emphasize lower tail but show overall CDF
            ax.plot(xs,ys,lw=1.0,label=label)
        ax.axhline(.05,ls='--',lw=.8)
        ax.set_title(state.replace('_',' '),fontsize=9)
        ax.grid(alpha=.2)
    axs[1,0].set_xlabel('Collected-power factor $\\eta$'); axs[1,1].set_xlabel('Collected-power factor $\\eta$')
    axs[0,0].set_ylabel('Empirical CDF'); axs[1,0].set_ylabel('Empirical CDF')
    axs[0,0].legend(fontsize=7)
    fig.tight_layout()
    for ext in ['png','pdf']:
        fig.savefig(FIG/f'Fig2_state_cdfs.{ext}',dpi=600 if ext=='png' else None,bbox_inches='tight')
    plt.close(fig)


def make_lcvar_fig(state_metrics):
    pv=state_metrics.pivot(index='state',columns='design',values='lcvar5').loc[STATES]
    fig,ax=plt.subplots(figsize=(7.0,4.2))
    x=np.arange(len(STATES)); width=0.13
    for j,col in enumerate(pv.columns):
        ax.bar(x+(j-(len(pv.columns)-1)/2)*width,pv[col].values,width,label=col)
    ax.set_xticks(x); ax.set_xticklabels([s.replace('_',' ') for s in STATES])
    ax.set_ylabel('Lower-tail conditional mean (5%)')
    ax.grid(axis='y',alpha=.2); ax.legend(fontsize=6,ncol=2)
    fig.tight_layout()
    for ext in ['png','pdf']:
        fig.savefig(FIG/f'Fig3_state_lcvar.{ext}',dpi=600 if ext=='png' else None,bbox_inches='tight')
    plt.close(fig)


def make_sensitivity_fig(df_jit, df_turb):
    fig,ax=plt.subplots(figsize=(6.2,4.2))
    for name,g in df_jit.groupby('design'):
        ax.plot(g['sigma_major_m']*1000,g['worst_lcvar5'],marker='o',label=name)
    ax.set_xlabel('Major-axis pointing jitter standard deviation (mm)'); ax.set_ylabel('Worst-state 5% lower-tail mean')
    ax.grid(alpha=.2); ax.legend(fontsize=7); fig.tight_layout()
    for ext in ['png','pdf']:
        fig.savefig(FIG/f'Fig4_pointing_sensitivity.{ext}',dpi=600 if ext=='png' else None,bbox_inches='tight')
    plt.close(fig)

    fig,ax=plt.subplots(figsize=(6.2,4.2))
    for name,g in df_turb.groupby('design'):
        ax.plot(g['sigma_ln'],g['worst_lcvar5'],marker='s',label=name)
    ax.set_xlabel('Lognormal standard-deviation parameter $\\sigma_{ln}$'); ax.set_ylabel('Worst-state 5% lower-tail mean')
    ax.grid(alpha=.2); ax.legend(fontsize=7); fig.tight_layout()
    for ext in ['png','pdf']:
        fig.savefig(FIG/f'Fig5_scintillation_sensitivity.{ext}',dpi=600 if ext=='png' else None,bbox_inches='tight')
    plt.close(fig)


def make_tradeoff_fig(agg):
    fig,ax=plt.subplots(figsize=(6.4,4.5))
    for _,r in agg.iterrows():
        ax.scatter(r['equal_state_mean_eta'],r['worst_state_lcvar5'],s=35)
        ax.annotate(r['design'],(r['equal_state_mean_eta'],r['worst_state_lcvar5']),xytext=(4,3),textcoords='offset points',fontsize=7)
    ax.set_xlabel('Equal-state mean collected-power factor'); ax.set_ylabel('Worst-state 5% lower-tail mean')
    ax.grid(alpha=.2); fig.tight_layout()
    for ext in ['png','pdf']:
        fig.savefig(FIG/f'Fig6_mean_tail_tradeoff.{ext}',dpi=600 if ext=='png' else None,bbox_inches='tight')
    plt.close(fig)


def main():
    with open(ROOT/'reproducibility'/'config.json','w') as f: json.dump(CONFIG,f,indent=2)
    with open(ROOT/'reproducibility'/'software_versions.txt','w') as f:
        f.write(f'Python {sys.version.split()[0]}\nNumPy {np.__version__}\nSciPy {scipy.__version__}\nPandas {pd.__version__}\nMatplotlib {matplotlib.__version__}\nPlatform {platform.platform()}\n')

    # Exact formula versus independent numerical quadrature
    pairs=[(0,.08),(.03,.08),(.07,.08),(.05,.12),(.10,.16)]
    qrows=[]
    for d,w in pairs:
        ex=float(aperture_fraction(d,w)); qu=quadrature_fraction(d,w)
        qrows.append(dict(displacement_m=d,beam_radius_m=w,noncentral_chi2=ex,direct_quadrature=qu,abs_difference=abs(ex-qu)))
    pd.DataFrame(qrows).to_csv(DATA/'aperture_quadrature_validation.csv',index=False)

    # Nominal independent test evaluation
    allrows=[]; arrays={}
    for j,(name,dsgn) in enumerate(DESIGNS.items()):
        df,arr=eval_design(dsgn,n=CONFIG['test_samples_per_state'],seed=CONFIG['base_seed']+20000*j)
        df.insert(0,'design',name); allrows.append(df); arrays[name]=arr
    state_df=pd.concat(allrows,ignore_index=True)
    state_df.to_csv(DATA/'state_metrics.csv',index=False)

    aggs=[]
    for name,g in state_df.groupby('design',sort=False):
        blocked=g[g.state!='none']
        aggs.append(dict(
            design=name,
            x0_m=DESIGNS[name]['x0'],y0_m=DESIGNS[name]['y0'],theta_urad=DESIGNS[name]['theta_urad'],
            beam_radius_m=beam_radius(DESIGNS[name]['theta_urad']),
            equal_state_mean_eta=g.mean_eta.mean(),
            worst_state_lcvar5=g.lcvar5.min(),
            average_state_lcvar5=g.lcvar5.mean(),
            worst_state_p5=g.p5.min(),
            blocked_lcvar_spread=blocked.lcvar5.max()-blocked.lcvar5.min()
        ))
    agg=pd.DataFrame(aggs)
    agg.to_csv(DATA/'design_metrics.csv',index=False)

    # Independent repeated trials with per-trial worst-state metric and equal-state mean
    rep=[]
    for trial in range(CONFIG['repeat_trials']):
        for j,(name,dsgn) in enumerate(DESIGNS.items()):
            df,_=eval_design(dsgn,n=CONFIG['repeat_samples_per_state'],seed=CONFIG['base_seed']+500000+trial*10000+j*137)
            rep.append(dict(trial=trial+1,design=name,worst_state_lcvar5=df.lcvar5.min(),equal_state_mean_eta=df.mean_eta.mean(),worst_state_p5=df.p5.min()))
    repdf=pd.DataFrame(rep); repdf.to_csv(DATA/'repeated_trials.csv',index=False)
    reps=[]
    for name,g in repdf.groupby('design',sort=False):
        n=len(g); tcrit=2.262157 # df=9, two-sided 95%
        for metric in ['worst_state_lcvar5','equal_state_mean_eta','worst_state_p5']:
            mean=g[metric].mean(); sd=g[metric].std(ddof=1); half=tcrit*sd/math.sqrt(n)
            reps.append(dict(design=name,metric=metric,mean=mean,sd=sd,ci95_low=mean-half,ci95_high=mean+half,n_trials=n))
    pd.DataFrame(reps).to_csv(DATA/'repeated_trial_summary.csv',index=False)

    # Sensitivity - fixed designs, independent CRN-style seeds for pairwise comparability
    jitrows=[]
    for si,sig in enumerate([.020,.030,.035,.040,.050,.060]):
        for name in ['Equal-state mean joint','Proposed robust joint']:
            df,_=eval_design(DESIGNS[name],n=30000,seed=CONFIG['base_seed']+700000+si*5000,sig_major=sig)
            jitrows.append(dict(sigma_major_m=sig,design=name,worst_lcvar5=df.lcvar5.min(),equal_state_mean_eta=df.mean_eta.mean()))
    jitdf=pd.DataFrame(jitrows); jitdf.to_csv(DATA/'sensitivity_jitter.csv',index=False)

    turbrows=[]
    for ti,sln in enumerate([.15,.25,.38,.50,.65]):
        for name in ['Equal-state mean joint','Proposed robust joint']:
            df,_=eval_design(DESIGNS[name],n=30000,seed=CONFIG['base_seed']+800000+ti*5000,sigma_ln=sln)
            turbrows.append(dict(sigma_ln=sln,design=name,worst_lcvar5=df.lcvar5.min(),equal_state_mean_eta=df.mean_eta.mean()))
    turbd=pd.DataFrame(turbrows); turbd.to_csv(DATA/'sensitivity_turbulence.csv',index=False)

    # MC convergence of the two key designs
    conv=[]
    for ni,n in enumerate([2000,5000,10000,20000,40000]):
        for name in ['Equal-state mean joint','Proposed robust joint']:
            df,_=eval_design(DESIGNS[name],n=n,seed=CONFIG['base_seed']+900000+ni*7000)
            conv.append(dict(samples_per_state=n,design=name,worst_lcvar5=df.lcvar5.min(),equal_state_mean_eta=df.mean_eta.mean()))
    convdf=pd.DataFrame(conv); convdf.to_csv(DATA/'monte_carlo_convergence.csv',index=False)

    # A concise manuscript comparison table with changes relative to selected references
    ref=agg.set_index('design')
    prop=ref.loc['Proposed robust joint']
    comps=[]
    for base in ['Centroid-fixed','Bias-cancelled','Divergence-only robust','Center-only robust','Equal-state mean joint']:
        b=ref.loc[base]
        comps.append(dict(
            baseline=base,
            baseline_worst_lcvar5=b.worst_state_lcvar5,
            proposed_worst_lcvar5=prop.worst_state_lcvar5,
            relative_change_worst_lcvar_percent=100*(prop.worst_state_lcvar5/b.worst_state_lcvar5-1),
            baseline_equal_state_mean=b.equal_state_mean_eta,
            proposed_equal_state_mean=prop.equal_state_mean_eta,
            relative_change_mean_percent=100*(prop.equal_state_mean_eta/b.equal_state_mean_eta-1)
        ))
    pd.DataFrame(comps).to_csv(DATA/'baseline_comparisons.csv',index=False)

    make_fig1()
    # Use same independent raw arrays from nominal evaluation for key comparison
    make_ecdf_fig(arrays['Equal-state mean joint'],arrays['Proposed robust joint'])
    make_lcvar_fig(state_df)
    make_sensitivity_fig(jitdf,turbd)
    make_tradeoff_fig(agg)

    print('\nAGGREGATE DESIGN METRICS')
    print(agg.to_string(index=False,float_format=lambda x:f'{x:.6f}'))
    print('\nREPEATED TRIAL SUMMARY')
    print(pd.DataFrame(reps).to_string(index=False,float_format=lambda x:f'{x:.6f}'))
    print('\nQUADRATURE VALIDATION')
    print(pd.DataFrame(qrows).to_string(index=False,float_format=lambda x:f'{x:.12g}'))

if __name__=='__main__':
    main()
