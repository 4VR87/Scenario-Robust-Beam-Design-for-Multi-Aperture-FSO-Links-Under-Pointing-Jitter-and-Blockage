#!/usr/bin/env python3
"""Fixed-seed sample-average optimization for manuscript design parameters.

This script is intentionally separate because it is more computationally expensive than
simulate_results.py. It uses common random numbers across candidate designs.
"""
from pathlib import Path
import sys, math
import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution
sys.path.insert(0,str(Path(__file__).resolve().parent))
import simulate_results as sim

ROOT=Path(__file__).resolve().parents[1]
N=5000
SEED=20260919+12345
# Cache one common-random-number sample per scenario state.
RNDS={s:sim.generate_randoms(N,SEED+1009*i) for i,s in enumerate(sim.STATES)}

def dsgn(v): return {'x0':float(v[0]),'y0':float(v[1]),'theta_urad':float(v[2])}

def state_stats(v):
    dd=dsgn(v); l=[]; m=[]
    for s in sim.STATES:
        x=sim.collected(dd,s,RNDS[s]); l.append(sim.lower_tail_mean(x,.05)); m.append(float(np.mean(x)))
    return np.array(l),np.array(m)

def robust_obj(v): return -float(np.min(state_stats(v)[0]))
def mean_obj(v): return -float(np.mean(state_stats(v)[1]))

def run(fun,bounds,seed,maxiter=55,popsize=11):
    return differential_evolution(fun,bounds,seed=seed,maxiter=maxiter,popsize=popsize,tol=2e-4,polish=True,workers=1,updating='immediate')

rows=[]
# Full 3D designs
bounds=[(-.10,.10),(-.10,.10),(60.,350.)]
for name,fun,seed in [('Proposed robust joint',robust_obj,181),('Equal-state mean joint',mean_obj,182)]:
    r=run(fun,bounds,seed)
    lc,me=state_stats(r.x)
    rows.append(dict(design=name,x0_m=r.x[0],y0_m=r.x[1],theta_urad=r.x[2],sample_objective=-r.fun,worst_lcvar5=lc.min(),equal_state_mean_eta=me.mean(),nit=r.nit,nfev=r.nfev,success=r.success))
# center-only with theta=100

def center_obj(v): return robust_obj([v[0],v[1],100.0])
r=run(center_obj,[(-.10,.10),(-.10,.10)],183,maxiter=45,popsize=10)
v=[r.x[0],r.x[1],100.0]; lc,me=state_stats(v)
rows.append(dict(design='Center-only robust',x0_m=v[0],y0_m=v[1],theta_urad=v[2],sample_objective=-r.fun,worst_lcvar5=lc.min(),equal_state_mean_eta=me.mean(),nit=r.nit,nfev=r.nfev,success=r.success))
# divergence-only at bias cancellation

def div_obj(v): return robust_obj([-.03,.015,v[0]])
r=run(div_obj,[(60.,350.)],184,maxiter=35,popsize=10)
v=[-.03,.015,r.x[0]]; lc,me=state_stats(v)
rows.append(dict(design='Divergence-only robust',x0_m=v[0],y0_m=v[1],theta_urad=v[2],sample_objective=-r.fun,worst_lcvar5=lc.min(),equal_state_mean_eta=me.mean(),nit=r.nit,nfev=r.nfev,success=r.success))

out=pd.DataFrame(rows)
out.to_csv(ROOT/'data'/'optimization_results.csv',index=False)
print(out.to_string(index=False))
