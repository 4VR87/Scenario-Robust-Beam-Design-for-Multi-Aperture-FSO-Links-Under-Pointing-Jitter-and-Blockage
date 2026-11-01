"""New evidence for reviewer comments on parameter assumptions and validation.

The independent solver is isolated in green_independent_solver.py. This driver
compares its results with the original analytic/Monte Carlo implementation.
All new results are simulations and retained designs are never reselected.
"""
from pathlib import Path
import sys,json,time,platform
import numpy as np
import pandas as pd
import scipy
from scipy.stats import t
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'repository'))
import simulate_results as original
import green_independent_solver as independent
OUT=ROOT/'green_revision_data';OUT.mkdir(exist_ok=True)
ALL_DESIGNS=json.loads((ROOT/'pink_revision_data/designs.json').read_text())
DESIGNS={key:ALL_DESIGNS[key] for key in ['S3_legacy_mean','S3_legacy_joint','S6_center_alpha05','S6_joint_alpha05']}
PATTERNS=independent.PAIRS

def primary_power(vector,channel):
    p,bh=channel
    distance=np.linalg.norm((p+np.array(vector[:2]))[:,None,:]-original.AP_POS[None,:,:],axis=2)
    return np.sum(bh*original.aperture_fraction(distance,original.beam_radius(vector[2])),axis=1)

def primary_channel(n,seed,pair):
    p,h,u,um,ua=original.generate_randoms(n,seed)
    b=np.ones_like(h)
    for j in pair:b[:,j]=.12*u[:,j]
    return p,h*np.where(um<.01,b*.4*ua,b)

def summarize(frame,groups,metrics):
    rows=[]
    for keys,g in frame.groupby(groups,sort=False):
        if not isinstance(keys,tuple):keys=(keys,)
        for metric in metrics:
            x=g[metric].to_numpy();mean=x.mean();sd=x.std(ddof=1)
            half=t.ppf(.975,len(x)-1)*sd/np.sqrt(len(x))
            rows.append(dict(zip(groups,keys))|dict(metric=metric,mean=mean,sd=sd,
                ci95_low=mean-half,ci95_high=mean+half,repetitions=len(x)))
    return pd.DataFrame(rows)

def validate():
    start=time.perf_counter()
    couplings={key:independent.QuadratureCoupling(vector[2]) for key,vector in DESIGNS.items()}
    checks=[]
    distances=np.linspace(0,.5499,276)+.000037
    for name,vector in DESIGNS.items():
        width=couplings[name].width
        low=independent.disk_integral(distances,width,16,64)
        high=independent.disk_integral(distances,width,64,256)
        mid=independent.disk_integral(distances,width,32,128)
        reference=original.aperture_fraction(distances,width)
        checks.append(dict(design=name,width_m=width,
            low_high_max_abs=float(np.max(np.abs(low-high))),
            mid_high_max_abs=float(np.max(np.abs(mid-high))),
            high_analytic_max_abs=float(np.max(np.abs(high-reference))),
            interpolated_analytic_max_abs=float(np.max(np.abs(couplings[name](distances)-reference)))))
    pd.DataFrame(checks).to_csv(OUT/'quadrature_checks.csv',index=False)
    paired=[]
    for si,pair in enumerate(PATTERNS):
        channel=independent.sobol_channel(8192,202610240500+1009*si,pair)
        for name,vector in DESIGNS.items():
            first=primary_power(vector,channel)
            second=independent.collected(vector,channel,couplings[name])
            paired.append(dict(design=name,state='none' if not pair else '-'.join(map(str,pair)),
                max_abs_power_error=float(np.max(np.abs(first-second))),
                mean_abs_power_error=float(np.mean(np.abs(first-second))),
                primary_L5=float(original.lower_tail_mean(first)),
                independent_L5=independent.lower_tail(second)))
    pd.DataFrame(paired).to_csv(OUT/'paired_power_checks.csv',index=False)
    rows=[];states=[]
    for n in [4096,16384,65536]:
        for trial in range(8):
            for method in ['Analytic_MC','Independent_RQMC']:
                by={name:[] for name in DESIGNS};means={name:[] for name in DESIGNS}
                for si,pair in enumerate(PATTERNS):
                    seed=(202610240600 if method=='Analytic_MC' else 202610240700)+1000000*trial+1009*si
                    channel=primary_channel(n,seed,pair) if method=='Analytic_MC' else independent.sobol_channel(n,seed,pair)
                    for name,vector in DESIGNS.items():
                        power=primary_power(vector,channel) if method=='Analytic_MC' else independent.collected(vector,channel,couplings[name])
                        tail=float(original.lower_tail_mean(power)) if method=='Analytic_MC' else independent.lower_tail(power)
                        mean=float(np.mean(power));by[name].append(tail);means[name].append(mean)
                        states.append(dict(n=n,trial=trial+1,method=method,design=name,
                            state='none' if not pair else '-'.join(map(str,pair)),L5=tail,mean_eta=mean))
                for name in DESIGNS:
                    rows.append(dict(n=n,trial=trial+1,method=method,design=name,
                        worst_L5=min(by[name]),mean_eta=float(np.mean(means[name]))))
            print('Independent validation',n,'trial',trial+1,'complete',flush=True)
        pd.DataFrame(rows).to_csv(OUT/'independent_validation_trials.csv',index=False)
    df=pd.DataFrame(rows)
    summarize(df,['n','method','design'],['worst_L5','mean_eta']).to_csv(OUT/'independent_validation_summary.csv',index=False)
    pd.DataFrame(states).to_csv(OUT/'independent_validation_state_trials.csv',index=False)
    meta=dict(geometry_max_abs_difference=float(np.max(np.abs(independent.POSITIONS-original.AP_POS))),
        pointing_covariance_max_abs_difference=float(np.max(np.abs(independent.POINTING_COV-original.cov_pointing()))),
        scintillation_covariance_max_abs_difference=float(np.max(np.abs(independent.LATENT_COV-original.C))),
        lookup_step_m=.0002,lookup_max_m=.55,radial_nodes=32,angular_nodes=128,
        independent_replications=8,sample_sizes=[4096,16384,65536],states=[list(p) for p in PATTERNS],
        elapsed_s=time.perf_counter()-start,
        interpretation='Independent numerical implementation of the same assumed channel; not experimental or wave-optics validation.')
    (OUT/'independent_validation_config.json').write_text(json.dumps(meta,indent=2))

def sensitivity():
    baseline=dict(sigma_ln=.38,corr_length_m=.055,block_max=.12,micro_probability=.01,micro_max=.4)
    cases=[('Nominal','nominal',baseline)]
    for field,values in [('sigma_ln',[.15,.65]),('corr_length_m',[.020,.120]),('block_max',[0.,.30]),
                         ('micro_probability',[0.,.05]),('micro_max',[.1,.8])]:
        for value in values:cases.append((f'{field}={value:g}',field,baseline|{field:value}))
    designs={name:ALL_DESIGNS[name] for name in ['S3_legacy_joint','S6_mean','S6_joint_alpha05']}
    separation=np.linalg.norm(original.AP_POS[:,None,:]-original.AP_POS[None,:,:],axis=2)
    roots={}
    for length in [.020,.055,.120]:
        cov=np.exp(-(separation/length)**(5/3))+1e-12*np.eye(7)
        roots[length]=np.linalg.cholesky(cov)
    rows=[]
    for trial in range(6):
        tails={(label,name):[] for label,_,_ in cases for name in designs}
        for si,pair in enumerate(PATTERNS):
            rng=np.random.default_rng(202610240900+1000000*trial+1009*si)
            n=20000
            pointing=rng.multivariate_normal([.030,-.015],original.cov_pointing(),size=n)
            standard=rng.standard_normal((n,7));u=rng.random((n,7));um=rng.random((n,7));ua=rng.random((n,7))
            geometric={}
            for name,v in designs.items():
                distance=np.linalg.norm((pointing+np.array(v[:2]))[:,None,:]-original.AP_POS[None,:,:],axis=2)
                geometric[name]=original.aperture_fraction(distance,original.beam_radius(v[2]))
            for label,family,params in cases:
                s=params['sigma_ln'];h=np.exp(s*(standard@roots[params['corr_length_m']].T)-.5*s*s)
                b=np.ones_like(h)
                for aperture in pair:b[:,aperture]=params['block_max']*u[:,aperture]
                b*=np.where(um<params['micro_probability'],params['micro_max']*ua,1.)
                for name in designs:
                    power=np.sum(geometric[name]*h*b,axis=1)
                    tails[label,name].append(float(original.lower_tail_mean(power)))
        for label,family,params in cases:
            for name in designs:
                rows.append(dict(trial=trial+1,case=label,family=family,design=name,worst_L5=min(tails[label,name]),**params))
        print('Parameter sensitivity trial',trial+1,'complete',flush=True)
    df=pd.DataFrame(rows);df.to_csv(OUT/'parameter_sensitivity_trials.csv',index=False)
    summarize(df,['case','family','design'],['worst_L5']).to_csv(OUT/'parameter_sensitivity_summary.csv',index=False)
    (OUT/'parameter_sensitivity_config.json').write_text(json.dumps(dict(n_per_state=20000,trials=6,
        base_seed=202610240900,trial_increment=1000000,state_increment=1009,
        scenarios=[dict(case=label,family=family,**params) for label,family,params in cases],
        policy='One parameter at a time, no retuning, common random numbers across cases and designs within trial/state.'),indent=2))

if __name__=='__main__':
    (OUT/'environment.json').write_text(json.dumps(dict(python=sys.version,numpy=np.__version__,
        scipy=scipy.__version__,pandas=pd.__version__,platform=platform.platform(),
        independent_module='green_independent_solver.py',primary_module='repository/simulate_results.py'),indent=2))
    validate();sensitivity()
    print('ALL GREEN REVIEW EXPERIMENTS COMPLETE',flush=True)
