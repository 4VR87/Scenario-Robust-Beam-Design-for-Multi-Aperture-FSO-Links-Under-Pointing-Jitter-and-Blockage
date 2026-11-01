"""Add the second reviewer's studies to the preserved manuscript OOXML."""
PDATA=ROOT/'pink_revision_data'
opt6=pd.read_csv(PDATA/'optimizations.csv')
sum6=pd.read_csv(PDATA/'independent_trial_summary.csv')
trial6=pd.read_csv(PDATA/'independent_trials.csv')
state6=pd.read_csv(PDATA/'independent_state_metrics.csv')
curve6=pd.read_csv(PDATA/'link_power_sweep.csv')
mix6=pd.read_csv(PDATA/'randomized_sector_holdout.csv').set_index('design')
quant6=pd.read_csv(PDATA/'quantized_control_trials.csv')
check6=json.loads((PDATA/'model_checks.json').read_text())
cmp6=json.loads((PDATA/'s6_joint_center_comparison.json').read_text())
link6=json.loads((PDATA/'link_config.json').read_text())
seed6=opt6[opt6.run.str.contains('S6_joint_seed')].sort_values('run')

def sm6(design,metric='worst_L5'):
    return sum6[(sum6.design==design)&(sum6.metric==metric)].iloc[0]
def value6(design,metric='worst_L5'):
    return float(sm6(design,metric)['mean'])
def ci6(design,metric='worst_L5',digits=5):
    r=sm6(design,metric)
    return f'{r["mean"]:.{digits}f} ({r.ci95_low:.{digits}f}–{r.ci95_high:.{digits}f})'
def upright(text):
    r=mr(text);mp=r.find(qn('m:rPr'))
    if mp is None:mp=OxmlElement('m:rPr');r.insert(0,mp)
    st=OxmlElement('m:sty');st.set(qn('m:val'),'p');mp.append(st)
    return r
def findp(prefix):
    found=[p for p in D.paragraphs if p.text.startswith(prefix)]
    if len(found)!=1:raise RuntimeError(f'Expected one paragraph {prefix}: {len(found)}')
    return found[0]

# Changed text is shaded light green in Word and in the PDF rendering.
def convert_to_green(doc):
    for h in list(doc.element.xpath('.//w:highlight')):
        rp=h.getparent();rp.remove(h)
        shade=rp.find(qn('w:shd'))
        if shade is None:shade=OxmlElement('w:shd');rp.append(shade)
        shade.set(qn('w:val'),'clear');shade.set(qn('w:color'),'auto');shade.set(qn('w:fill'),'DFF2D8')
    for shade in doc.element.xpath('.//w:shd[@w:fill="FFFF00"]'):
        shade.set(qn('w:fill'),'DFF2D8')

replace(P[0],'Risk Aware Beam Center Placement with Conditional Divergence Refinement for Multi Aperture FSO Receivers')
P[0].style='Title'
for r in P[0].runs:font_run(r,16,bold=True)
replace(P[7],
    'Asymmetric blockage and anisotropic pointing jitter can create deep power losses in multi-aperture free-space optical (FSO) receivers. We study risk-aware beam-center placement by maximizing the worst-state conditional mean of the lowest 5% of collected-power realizations, with divergence as an optional refinement. Exact Gaussian-to-aperture coupling is combined with correlated unit-mean lognormal scintillation. In the original three-pair scenario set, ten independent trials give a tail mean of 0.10074 versus 0.04598 for mean-power optimization, a 119.07% increase with 22.29% less mean collection. The 1.25% increment over center-only optimization is not statistically established. Expanding the design set to all six adjacent outer-aperture pairs exposes scenario dependence and motivates re-optimization. Five independent searches give similar validation objectives; the retained expanded-set design achieves 0.08741 (95% CI 0.08687–0.08795), versus 0.08069 (0.07997–0.08142) for center-only optimization at 100 μrad. Its paired increment is 8.32%, showing that divergence can matter under different blockage geometry. Tail levels of 1%, 5%, and 10% produce different footprint trade-offs. An illustrative equal-weight direct-detection link model provides SNR, BER, and outage comparisons, without equating power-tail gains with universal communication gains. Transmitter mapping, model checks, seed records, convergence histories, and timings define a reproducible simulation workflow; field availability remains outside its scope.')

before(P[13]._p,
    'Robustness is also studied in related optical, wireless, and network settings. Guo et al. [34] use an optical intelligent reflecting surface for visible-light positioning in indoor non-line-of-sight vehicular scenarios. Guo et al. [35] jointly allocate power and active/passive beamforming in an RF UAV–RIS communication system. Wei et al. [36] protect network connectivity against failures using robust optimization and random search. These works provide context for obstruction-aware optical design, coupled controls, and failure-aware objectives, respectively. Their propagation models, controls, and performance criteria differ from the present seven-aperture FSO footprint and worst-state lower-tail power objective; none is treated as a direct numerical FSO benchmark.')
replace(P[13],
    'This work studies a static footprint for one central and six peripheral circular receiver apertures. The risk criterion maximizes the minimum statewise lower-tail collected-power mean. Center placement is the principal nominal mechanism; divergence is a condition-dependent additional control. The original compact scenario set S3 contains no pair blockage and three asymmetric adjacent outer-pair losses. To test its limits, we evaluate all two-aperture losses and independently optimize an expanded set S6 containing no pair blockage and all six adjacent outer-pair losses. Neither set is assigned empirical occurrence probabilities. The reward-side lower-tail definition is related to conditional value-at-risk [22] and is stated directly to avoid loss/reward sign ambiguity.')
replace(P[14],
    'The study contributes a statewise lower-tail objective for risk-aware center placement; an exact Gaussian coupling model with anisotropic pointing, correlated scintillation, and explicit retained-power blockage factors; ablations that distinguish nominal center-placement gains from condition-dependent divergence gains; and reproducibility studies spanning optimizer seeds, SAA seeds, tail levels, expanded blockage states, and computational cost. An illustrative link model compares the resulting footprints against mean-power, average-SNR, and outage criteria adapted to this receiver. All results are simulations, with no claim of experimental validation or globally certified optimization.')

# Normalization and sampling definitions are adjacent to the model equations.
equation_before(P[34]._p,expr('E[',sub('h','i'),'] = exp(',frac([ss('σ','ln','2'),'(',sub('C','ii'),' − 1)'],'2'),') = 1,     ',sub('C','ii'),' = 1'),'6a')
replace(P[36],
    'Here Dij is aperture separation, σln = 0.38, and rc = 55 mm. Equation (6a) verifies unit mean for Cii = 1. The unregularized seven-aperture matrix in (7) has eigenvalues 0.65153, 0.70299, 0.80275 (twice), 1.16692 (twice), and 1.70613; its minimum eigenvalue is positive and its condition number is 2.619. Cholesky reconstruction has maximum absolute error 1.11 × 10⁻¹⁶. The source code adds 10⁻¹²I for numerical safety, making E[hi] differ from one by less than 8 × 10⁻¹⁴. In 200,000 independent normalization draws, aperture means range from 0.99961 to 1.00152, at most 1.73 standard errors from one. Cij is latent Gaussian correlation; irradiance correlation is [exp(σln²Cij) − 1]/[exp(σln²) − 1]. This phenomenological kernel is not a measured phase-coherence function or a Fried-parameter model [3,23,32].')
replace(P[37],
    'The original scenario set S3 = {none, (3,4), (1,2), (5,6)} uses three disjoint adjacent outer pairs, covering each peripheral aperture once. It is a compact stress set and is insufficient to establish robustness to other shadow orientations. In Fig. 1, the midpoint of pair (1,2) lies at azimuth 30°, along the major pointing axis; the omitted opposite pair (4,5), at 210°, lies on that same axis. Anisotropic jitter and a shifted beam therefore prevent a rotational-equivalence argument. Structural supports or moving obstacles could obscure adjacent apertures, but no measured structural-shadow map is available. Accordingly, S6 = {none, (1,2), (2,3), (3,4), (4,5), (5,6), (6,1)} includes every canonical adjacent outer-pair orientation and is optimized in Section 5.10. The unblocked label means no pair blockage; micro-obscuration still applies. Larger sectors, central-aperture loss, and arbitrary shadows require additional states; the retained S3 design is also stress-tested on all 21 pairs in Section 5.7.')
before(P[38]._p,
    'Every blockage factor denotes retained optical power, not fractional power lost or a field-amplitude factor. For aperture i in the selected pair, Bs,i is Uniform[0,0.12]; otherwise Bs,i = 1. For every aperture i and every Monte Carlo realization conditional on any state s, an independent event Ms,i is Bernoulli(0.01). Thus Pmicro = 0.01 = 1% per aperture per realization, rather than a probability of entering a receiver state. Given an event, Us,i is Uniform[0,0.4]. The total retained-power multiplier is:')
equation_before(P[38]._p,expr(sub('b','s,i'),' = ',sub('B','s,i'),'[1 − ',sub('M','s,i'),' + ',sub('M','s,i'),sub('U','s,i'),']'),'7a')
before(P[38]._p,
    'The variables Bs,i, Ms,i, and Us,i are independent across aperture and realization and independent of pointing and scintillation, conditional on s; only hi is spatially correlated. Pair membership is deterministic within a fixed state, but each blocked aperture draws its own multiplier. Pair blockage and micro-obscuration multiply when they coincide. The normalization check observes event frequencies 0.9765%–1.0320% across apertures, consistent with the prescribed 1%. These are controlled stressors without a calibrated temporal arrival process.')
replace(P[39],
    'Here bs,i is the retained-power multiplier in (7a), and ηs is dimensionless collected power before fixed terminal and atmospheric losses. Since hi describes local scintillation, ηs can exceed its no-fade geometric value; it is not a whole-wavefront energy-conservation model. Section 2.5 defines the assumed receiver summation and a conditional link-performance mapping. Table 1 summarizes notation.')

# Correct the notation and Table 2 roles consistently.
for row in T1.rows:
    symbol=''.join(row.cells[0]._tc.xpath('.//w:t/text() | .//m:t/text()'))
    if symbol=='bs,i':set_cell(row.cells[1],'State-dependent retained-power multiplier')
for row in T2.rows:
    key=row.cells[0].text
    vals=None
    if key=='Strong-block attenuation':vals=['Pair-block retained power','bblock','Uniform [0,0.12]','—','Independent retained-power multiplier per blocked aperture and realization']
    elif key=='Micro-obscuration probability':vals=['Micro-obscuration probability','Pmicro','0.01 = 1%','—','Per aperture per realization, conditional on each state; Bernoulli event']
    elif key=='Micro-obscuration attenuation':vals=['Micro-event retained power','bmicro','Uniform [0,0.4]','—','Independent retained-power multiplier, applied only on a micro-event']
    elif key=='Optimization state set':vals=['Scenario sets','S3; S6','4; 7 states','—','S3: none + 3 pairs; S6: none + all 6 adjacent outer pairs']
    elif key=='Tail probability':vals=['Tail probability','α','0.05; 0.01, 0.10','—','5% nominal; 1% and 10% sensitivity in Section 5.11']
    elif key=='DE population and generations':vals=['DE population and generations','—','33; ≤55 or ≤65','—','S3: at most 55; S6: at most 65 generations; center-only population 22']
    elif key=='SAA base seed':vals=['SAA base seeds','seedSAA','20273264; 20261024400','—','S3 and retained S6, respectively; state increments 1009']
    if vals:
        for c,t in zip(row.cells,vals):set_cell(c,t)
format_table(T2,[1.35,.70,1.10,.52,2.68])

# An explicit electrical combining assumption and a reproducible illustrative link.
anchor=P[43]._p
before(anchor,'2.5 Receiver combining and illustrative link performance','Heading 2')
before(anchor,
    'We assume intensity modulation/direct detection with seven photodetectors and an equal-weight sum of their electrical currents. The aperture optical powers add noncoherently; there is no optical phase combining, branch selection, or maximal-ratio weighting. All seven front ends remain active and contribute independent receiver noise even when their signal is blocked. Scintillation can remain correlated as specified in (7). Let Pavg be the equiprobable OOK average launch power, Tf a fixed power transmission, R the responsivity, and iN the per-branch input-current noise density over electrical equivalent noise bandwidth B. After removal of a constant background current, the OFF and ON mean currents are 0 and 2RPavgTfηs, with summed noise variance 7iN²B. Assuming signal-independent Gaussian noise and an adaptive midpoint decision threshold gives:')
equation_before(anchor,expr(ss('σ','Σ','2'),' = 7',ss('i','N','2'),'B,     κ = ',frac(['R',sub('P','avg'),sub('T','f')],sub('σ','Σ')),',     ',sub('γ','s'),' = ',sup(expr('(κ',sub('η','s'),')'),'2')),'8d')
equation_before(anchor,expr(sub('P','b'),'(',sub('η','s'),') = Q(κ',sub('η','s'),'),     ',sub('η','th'),' = ',frac([sup('Q','−1'),'(',sub('P','b,target'),')'],'κ')),'8e')
equation_before(anchor,expr(sub('P','out,s'),' = Pr(',sub('η','s'),' < ',sub('η','th'),'),     ',sub('P̄','b,s'),' = ',sub('E','s'),'[Q(κ',sub('η','s'),')]'),'8f')
before(anchor,
    'Here Q is the standard Gaussian upper-tail function; γs is the squared decision-distance-to-noise ratio, fixing the SNR convention. The illustrative operating point uses Pavg = −15 dBm = 31.62 μW, ON power 63.25 μW, Tf = 0.10 (10 dB fixed loss), R = 0.9 A/W, iN = 5 pA/√Hz, and B = 100 MHz. The assumed fixed loss represents atmospheric and terminal transmission; it is applied once, because geometric collection is already included in ηs. No additional free-space path-loss factor is inserted. These values give σΣ = 0.1323 μA, κ = 21.514, and ηth = 0.14364 for conditional BER target 10⁻³, equivalent to γth = 9.5495. The comparisons hold launch power and all losses/noise settings fixed across footprints.')
before(anchor,
    'Worst-state outage and average BER are maxs Pout,s and maxs P̄b,s, respectively. They depend on the whole distribution relative to the decision threshold, whereas L5 averages its lowest 5%; an L5 gain alone does not imply an equal SNR, BER, or outage gain. This thermal-noise example omits coding, background variability, bandwidth-dependent receiver design, and signal-dependent shot noise. At η = 1 and the nominal ON power, signal shot-noise variance is about 1% of the summed electronic-noise variance; the lognormal model is unbounded, so this is a reference-point check, not a uniform bound. Stronger signals or background require a fuller noise model. The example connects the footprint to defined communication metrics without claiming a calibrated deployment link budget.')

# Position center placement first, and explicitly state the optional joint extension.
replace(P[47],'For a selected tail fraction α, nominally 0.05, define the lower-tail conditional mean as')
replace(P[49],
    'This quantity averages the weakest α fraction of realizations within a state. We choose α = 0.05 as a design compromise between rare-fade emphasis and finite-sample precision, not as an empirical outage target. At Nopt = 5,000, its estimate averages 250 observations per state, compared with 50 at α = 0.01 and 500 at α = 0.10. Section 5.11 independently re-optimizes at all three levels. As a reward-side counterpart of CVaR [22], Lα,s has the same dimensionless power normalization as ηs and should not be interpreted as a channel availability probability.')
replace(P[50],
    'Risk-aware center placement optimizes (x₀,y₀) at a prescribed feasible divergence θref = 100 μrad. The optional three-variable extension also selects θ:')
equation_before(P[51]._p,expr(sub('max',expr(sub('x','0'),', ',sub('y','0'))),' ',sub('min',expr('s ∈ ',upright('S'))),' ',sub('L','α,s'),'(',sub('x','0'),', ',sub('y','0'),', ',sub('θ','ref'),')'),'10a')
eq11=equation_before(P[51]._p,expr(sub('max',expr(sub('x','0'),', ',sub('y','0'),', θ')),' ',sub('min',expr('s ∈ ',upright('S'))),' ',sub('L','α,s'),'(',sub('x','0'),', ',sub('y','0'),', θ)'),'11')
P[51]._p.getparent().remove(P[51]._p)
equation_before(P[53]._p,expr('−0.10 ',upright('m'),' ≤ ',sub('x','0'),', ',sub('y','0'),' ≤ 0.10 ',upright('m'),',     60 ',upright('μrad'),' ≤ θ ≤ 350 ',upright('μrad')),'12')
P[53]._p.getparent().remove(P[53]._p)
before(P[55]._p if 55 in P else P[56]._p,
    'Expanded S6 searches use the same 5,000 samples per state, bounds, mutation, recombination, tolerance, and one-worker immediate-update DE, with maxiter = 65 (33 population members for three variables, 22 for center-only). Five joint searches independently vary both SAA seed and optimizer seed; Section 5.6 separately isolates these sources for S3. Generation histories are saved before polishing; reported final objectives and nfev include polishing, while history callbacks are excluded from nfev. Independent evaluation uses ten new 20,000-sample/state trials with common random numbers across candidate designs within each trial. The first S6 search is predeclared as retained; independent validation never selects the best restart. The intervals summarize the mean of a fixed design’s empirical worst-state statistic across Monte Carlo trials, not uncertainty in the population-optimal controls or a global-optimality certificate.')
before(P[57]._p,
    'For S6 we also adapt two objective types from related FSO work: average-SNR beam-coordinate design [11] and outage-based beam-parameter design [13]. Under Section 2.5, the former maximizes the arithmetic state average of E[ηs²], proportional to average γs; the latter minimizes the arithmetic state average Pr(ηs < 0.14364) at the nominal power. All optimize the same three footprint controls and bounds with Nopt = 5,000 and the same fixed S6 SAA. We additionally re-optimize mean power over S6. Equal weights define comparison criteria, not measured state probabilities. These are objective-level adaptations to our model, not reproductions of the MIMO or UAV architectures, algorithms, or numerical results of [11,13]. Section 5.12 reports both the criteria and their resulting worst-state performance.')
replace(P[58],
    'Figure 2 summarizes the original S3 workflow. After defining the aperture geometry, pointing distribution, and correlated lognormal field, random realizations are cached within each state. A candidate center and optionally its divergence determine exact Gaussian aperture coupling; multiplication by retained-power and scintillation factors followed by noncoherent summation gives ηs. The S6 extension follows the same workflow with seven states rather than four, and the additional tail levels and link criteria change only the scoring function.')
replace(P[59],
    'For each state the lowest ceil(αN) power samples are averaged, and the minimum statewise value scores a robust candidate. Fixed-seed common random numbers compare candidates consistently. Independent trials then quantify Monte Carlo variability of retained designs. Generation histories, separate seed studies, SAA-budget checks, expanded blockage patterns, and link-performance curves test distinct aspects of the resulting footprint; termination of the heuristic does not imply mathematical global optimality.')

# Make the older extended-set discussion point to the actual new re-optimization.
replace(findp('These tests directly delimit'),
    'These tests delimit the scenario dependence of the original 119.07% gain. Protection against omitted pairs requires including them during design. Section 5.10 therefore re-optimizes on S6 rather than carrying the original S3 floor into the expanded setting. The all-21-pair checks in this subsection apply to the retained S3 designs; they do not establish all-pair optimality of the new S6 design. No empirical obstruction probabilities or arbitrary-shadow guarantees are inferred.')
replace(findp('The additional numerical studies in Sections'),
    'Sections 5.6–5.9 report additional tests of the original S3 designs separately from the historical nominal tables. Sections 5.10–5.13 add S6 re-optimization, tail-level sensitivity, link comparisons, and control quantization. No retained footprint is selected using independent validation outcomes. Exact scripts, seeds, settings, software versions, source snapshot, raw outputs, and timing records are supplied in FSO_Light_Green_Revision_Evidence.zip. The S6 random pools contain independently drawn realizations conditional on each state; common random numbers are shared across candidate designs within a state, not identical physical time series across states.')

# Correct the original scatter label using its unchanged numerical table.
original_metrics=pd.read_csv(ROOT/'repository/design_metrics.csv')
f,ax=plt.subplots(figsize=(6.7,3.7))
positions=[(-18,10),(7,-3),(7,-10),(8,6),(-112,9),(-8,-12)]
for j,(_,r) in enumerate(original_metrics.iterrows()):
    ax.scatter(r.equal_state_mean_eta,r.worst_state_lcvar5,s=38,zorder=3,color=plt.cm.tab10(j))
    ax.annotate(r['design'],(r.equal_state_mean_eta,r.worst_state_lcvar5),
        xytext=positions[j],textcoords='offset points',fontsize=8,ha='right' if j==5 else 'left')
ax.set(xlabel='Equal-state mean collected-power factor',ylabel='Worst-state 5% lower-tail mean',
    xlim=(.216,.352),ylim=(.040,.108))
ax.grid(alpha=.2);f.tight_layout();f.savefig(FIG/'Fig3_mean_tail_corrected.png');plt.close(f)
P[75].clear();P[75].alignment=WD_ALIGN_PARAGRAPH.CENTER
P[75].add_run().add_picture(str(FIG/'Fig3_mean_tail_corrected.png'),width=Inches(5.5))
replace(P[76],'Fig. 3 Mean–tail comparison under S3 for six footprints evaluated with 60,000 samples per state. Risk-aware center placement accounts for most nominal lower-tail improvement.')

# Added figures use numerical data only.
f,axs=plt.subplots(1,2,figsize=(6.7,2.85))
colors=['#234f70','#9b5b28','#46775b','#75528b','#8c3f4a']
for j,(_,r) in enumerate(seed6.iterrows()):
    h=pd.read_csv(PDATA/f'history_{r.run}.csv')
    axs[0].plot(h.generation,h.sample_score,lw=1.1,color=colors[j],label=f'Run {j+1}')
    nm='S6_joint_alpha05' if j==0 else r.run
    su=sm6(nm)
    axs[1].errorbar(j+1,su['mean'],yerr=[[su['mean']-su.ci95_low],[su.ci95_high-su['mean']]],fmt='o',color=colors[j],capsize=3,ms=4)
axs[0].set(xlabel='Generation before polishing',ylabel='SAA worst-state L₅',ylim=(.035,.095))
axs[0].legend(fontsize=7,ncol=2,loc='lower right')
axs[1].set(xlabel='Independent optimization run',ylabel='Evaluation mean worst-state L₅',xticks=range(1,6),ylim=(.084,.090))
for ax in axs:ax.grid(alpha=.2)
f.tight_layout();f.savefig(FIG/'Fig10_S6_convergence.png');plt.close(f)
f,axs=plt.subplots(1,2,figsize=(6.7,2.9))
for name,label,color,ls in [
    ('S3_legacy_mean','S3 mean power','#9a6528',':'),
    ('S3_legacy_joint','S3 joint tail','#888888','--'),
    ('S6_center_alpha05','S6 center tail','#46775b','-.'),
    ('S6_joint_alpha05','S6 joint tail','#234f70','-'),
    ('S6_outage','S6 outage criterion','#8b4968','--')]:
    g=curve6[curve6.design==name]
    axs[0].semilogy(g.average_tx_dBm,g.worst_outage,color=color,ls=ls,lw=1.4,label=label)
    axs[1].semilogy(g.average_tx_dBm,g.worst_mean_BER,color=color,ls=ls,lw=1.4,label=label)
axs[0].set(ylabel='Worst-state outage',ylim=(.002,1.3))
axs[1].set(ylabel='Worst-state mean BER',ylim=(5e-5,.5))
for ax in axs:
    ax.set(xlabel='Average launch power (dBm)',xticks=range(-23,-10,3));ax.grid(alpha=.22,which='both')
axs[1].legend(fontsize=7,loc='lower left');f.tight_layout();f.savefig(FIG/'Fig11_link_power.png');plt.close(f)

anchor=P[101]._p
before(anchor,'5.10 Reoptimization over all six adjacent pair losses','Heading 2')
before(anchor,
    'S6 includes one state with no pair blockage and six states with adjacent outer-pair blockage. Five independent joint optimizations use SAA base seeds 20,261,024,400; 20,261,124,400; 20,261,224,400; 20,261,324,400; and 20,261,424,400, paired with optimizer seeds 701–705. Each state seed increments by 1009. Table 8 reports every final parameter vector and training objective, followed by its mean independent-trial worst-state L5. The ten evaluation trials have 20,000 samples per state, base seeds 20,261,024,100 + 100,000(t − 1), and the same 1009 state increment. Run 1 is the retained S6 joint design by a predeclared rule, not because its test score is largest.')
rows=[]
for j,(_,r) in enumerate(seed6.iterrows()):
    nm='S6_joint_alpha05' if j==0 else r.run
    rows.append([str(j+1),f'{r.x0_m*1000:.2f}',f'{r.y0_m*1000:.2f}',f'{r.theta_urad:.2f}',f'{r.training_score:.5f}',ci6(nm)])
table_before(anchor,'Table 8 Independent S6 joint searches at α = 0.05 and their ten-trial evaluation intervals',
    ['Run','x₀ (mm)','y₀ (mm)','θ (μrad)','SAA L₅','Evaluation L₅ (95% CI)'],rows,[.40,.68,.68,.72,.78,2.29])
before(anchor,
    f'All five searches meet the SciPy termination criterion after {int(seed6.nit.min())}–{int(seed6.nit.max())} generations. Parameter means ± sample SD are x₀ = {seed6.x0_m.mean()*1000:.2f} ± {seed6.x0_m.std()*1000:.2f} mm, y₀ = {seed6.y0_m.mean()*1000:.2f} ± {seed6.y0_m.std()*1000:.2f} mm, and θ = {seed6.theta_urad.mean():.2f} ± {seed6.theta_urad.std():.2f} μrad. Training objectives range from {seed6.training_score.min():.5f} to {seed6.training_score.max():.5f}, while independent-trial mean objectives range from 0.08637 to 0.08763. These five-run summaries describe this finite search protocol, not a sampling distribution of the exact population optimizer. Figure 10 includes the full per-generation training histories and the independent-trial intervals. The S6 joint runs take {seed6.elapsed_s.min():.1f}–{seed6.elapsed_s.max():.1f} s and {int(seed6.nfev.min())}–{int(seed6.nfev.max())} objective calls; three separate one-worker optimization jobs shared the host, so these elapsed times are workload-specific. All eleven S6 searches, including the other criteria, terminate successfully; their individual logs are supplied.')
before(anchor,
    f'The retained S6 joint vector is (−29.76, 15.14) mm with θ = 124.76 μrad. Its worst-state L5 is {ci6("S6_joint_alpha05")}, compared with {ci6("S3_legacy_joint")} for the retained S3 joint vector evaluated on S6, a {100*(value6("S6_joint_alpha05")/value6("S3_legacy_joint")-1):.2f}% increase. The S6 mean-power design gives {ci6("S6_mean")}; the S6 joint tail design increases this metric by {100*(value6("S6_joint_alpha05")/value6("S6_mean")-1):.2f}% while reducing equal-state mean collection from {value6("S6_mean","mean_eta"):.5f} to {value6("S6_joint_alpha05","mean_eta"):.5f} ({100*(1-value6("S6_joint_alpha05","mean_eta")/value6("S6_mean","mean_eta")):.2f}%). These expanded-set figures replace no historical S3 values and are not a claim that one footprint dominates every state.')
before(anchor,
    f'The separately optimized S6 center-only vector is (−30.68, 16.18) mm at θ = 100 μrad, with L5 = {ci6("S6_center_alpha05")}. Common trial streams support a paired comparison: the joint-minus-center difference is {cmp6["mean_difference"]:.5f}, with 95% CI [{cmp6["ci95_low"]:.5f}, {cmp6["ci95_high"]:.5f}], or {cmp6["mean_relative_percent"]:.2f}% relative to the center-only mean. This interval excludes zero for these fixed retained designs under S6. It does not overturn the original S3 result, where the 1.25% increment was not established, nor does it prove a positive increment for all channel conditions. Divergence is therefore a condition-dependent refinement, whose value changes when omitted blockage orientations are protected.')
before(anchor,
    f'A fresh 60,000-realization randomized-sector test draws one of the six adjacent pairs uniformly for each realization. The retained S6 joint and center-only mixture L5 values are {mix6.loc["S6_joint_alpha05","L5"]:.5f} and {mix6.loc["S6_center_alpha05","L5"]:.5f}, respectively. The center-only mixture is larger even though its conditional worst-state floor is lower. A lower tail of a mixture is a different objective from the minimum of conditional lower tails, so this ordering is consistent with the design criterion. The test changes orientation per realization without introducing new pair geometries; it is not evidence for arbitrary shadow shapes or a measured occurrence distribution.')
image_before(anchor,FIG/'Fig10_S6_convergence.png',6.25,
    'Fig. 10 S6 joint optimization at α = 0.05. Left: generation-by-generation best SAA objective before final polishing for five independent SAA/optimizer seed pairs. Right: mean worst-state L5 and 95% Student-t intervals over ten independent 20,000-sample/state trials for each fixed optimized vector.')

before(anchor,'5.11 Sensitivity to the selected lower tail fraction','Heading 2')
before(anchor,
    'Independent S6 joint searches at α = 0.01, 0.05, and 0.10 use the same 5,000-sample/state SAA and optimizer seed 701; each score averages 50, 250, or 500 samples per state. Table 9 evaluates every resulting footprint at all three tail levels using the ten independent evaluation trials. Thus a footprint is not evaluated only at the quantile used to train it. The table reports mean empirical worst-state metrics; per-cell 95% intervals and statewise 60,000-sample checks are included in the evidence files.')
rows=[]
for a,nm,run in [(.01,'S6_joint_alpha01','S6_joint_alpha01'),(.05,'S6_joint_alpha05','S6_joint_seed1'),(.10,'S6_joint_alpha10','S6_joint_alpha10')]:
    r=opt6[opt6.run==run].iloc[0]
    rows.append([f'{a:.2f}',f'{r.x0_m*1000:.2f}',f'{r.y0_m*1000:.2f}',f'{r.theta_urad:.2f}',f'{value6(nm,"worst_L1"):.5f}',f'{value6(nm):.5f}',f'{value6(nm,"worst_L10"):.5f}'])
table_before(anchor,'Table 9 Cross evaluation of S6 footprints optimized for three tail fractions',
    ['Training α','x₀ (mm)','y₀ (mm)','θ (μrad)','Test L₁','Test L₅','Test L₁₀'],rows,[.75,.75,.75,.80,1.0,1.0,1.0])
before(anchor,
    'The 1% design selects a wider footprint (151.26 μrad) and raises the evaluated 1% tail mean to 0.06085, compared with 0.05622 for the 5% design, while reducing the evaluated 10% tail mean. The 10% design narrows the divergence to 112.69 μrad and trades away the deepest-tail protection. This observed change in footprint and cross-tail ranking shows that the choice of risk fraction matters; it does not establish one universally preferred α. The 1% optimization has only 50 training tail observations per state, so independent SAA replication or a larger optimization budget would be especially valuable before using that design. The main 5% choice remains a stated design preference, not a fitted physical reliability requirement.')

before(anchor,'5.12 Link performance and comparisons of optimization criteria','Heading 2')
before(anchor,
    'Table 10 evaluates the retained S3 and S6 footprints on the same seven S6 states under Section 2.5. Mean SNR is the arithmetic state average of E[γs] in linear units; outage and average BER are the maximum statewise values. Each entry is the mean over ten independent trials. The average-SNR and outage baselines adapt objective types from [11,13], respectively, while keeping our propagation model, footprint variables, power, noise, state weights, and optimizer budget fixed. The empirical outage objective is discontinuous and its DE solution is only a heuristic baseline, not a certified optimum. Full per-metric intervals are provided in the evidence files.')
rows=[]
for nm,label in [('S3_legacy_joint','S3 joint tail'),('S6_mean','S6 mean power'),('S6_snr','S6 average SNR [11]'),('S6_outage','S6 average outage [13]'),('S6_center_alpha05','S6 center tail'),('S6_joint_alpha05','S6 joint tail')]:
    rows.append([label,f'{value6(nm):.5f}',f'{value6(nm,"mean_SNR"):.2f}',f'{100*value6(nm,"worst_outage"):.2f}%',f'{value6(nm,"worst_mean_BER"):.3e}'])
table_before(anchor,'Table 10 Illustrative S6 link comparison at average launch power −15 dBm',
    ['Design criterion','Worst L₅','Mean SNR','Worst outage','Worst mean BER'],rows,[2.0,.85,.85,.95,1.40])
before(anchor,
    'The average-SNR and mean-power designs select θ = 60 μrad and nearly cancel the pointing bias, producing the largest average SNR but weak worst-state tails. The average-outage criterion selects θ = 82.74 μrad and gives 10.56% worst-state outage, versus 13.49% for the S6 joint tail design at the stated power. The center-only tail design gives 10.30% outage. Conversely, worst-state average BER is 0.00247 for the S6 joint tail design, versus 0.00330 for center-only and 0.00506 for the outage criterion. A lower conditional power-tail score is therefore not a universal communication-performance ranking. These comparisons expose rather than suppress the consequences of the selected objective.')
before(anchor,
    'Figure 11 varies average launch power from −23 to −11 dBm while retaining every footprint, the 10 dB fixed transmission loss, and the receiver noise parameters. Each curve uses a separate common 60,000-sample/state evaluation, not the ten-trial means of Table 10. The outage curves cross: for example, the S6 joint and center-only values are 13.52% and 10.42% at −15 dBm, but 2.37% and 2.93% at −13 dBm. Power-dependent ranking is expected because the BER threshold moves through different parts of the power distribution. The outage baseline was optimized only at −15 dBm and is not re-optimized along the sweep. This illustrates why the original 119.07% S3 power-tail gain must not be reported as a corresponding BER reduction or field-availability gain.')
image_before(anchor,FIG/'Fig11_link_power.png',6.0,
    'Fig. 11 S6 worst-state outage for conditional BER target 10⁻³ (left) and worst-state mean uncoded OOK BER (right). Each fixed footprint uses 60,000 samples per state at each launch power, with the receiver assumptions of Section 2.5.')

before(anchor,'5.13 Practical steering and divergence resolution','Heading 2')
before(anchor,
    'For small angles at 1 km, an offset command (x₀,y₀) corresponds to beam-direction angles approximately (x₀/z,y₀/z). The original retained command is (−39.08,11.95) μrad; the S6 command is (−29.76,15.14) μrad, referenced to nominal boresight before residual pointing bias. A 1 mm receiver-plane control step is a 1 μrad beam-direction step. These are optical steering angles, not a specified actuator or mirror-rotation calibration. The S6 divergence 124.76 μrad corresponds to M² ≈ 2.529 at wT = 10 mm under Section 2.1; a practical terminal must calibrate waist, focus, beam quality, and divergence together within its clear aperture. Transmitted power is held fixed rather than increased to compensate for a broader footprint.')
before(anchor,
    'As a numerical control-resolution check, the S6 center-only and joint vectors are rounded to 1 mm in each offset and to 1 μrad in divergence. The resulting commands are (−31,16) mm at 100 μrad and (−30,15) mm at 125 μrad. Across the same ten independent trial pools, mean worst-state L5 changes from 0.08069 to 0.08066 (−0.043%) and from 0.08741 to 0.08756 (+0.179%), respectively. Small positive rounding changes reflect finite SAA/search precision and Monte Carlo sampling; they do not imply that coarse control improves the population optimum. This local check does not validate tracking dynamics, aperture clipping during adjustment, steering range, jitter spectra, receiver saturation, or hardware power constraints. Those quantities require a specified terminal and measured calibration.')

replace(P[102],
    'The results are conditional on a design-level simulation. The lognormal field and 55 mm covariance scale are prescribed rather than measured, and the Rytov conversion is an interpretive equivalent. S6 covers all adjacent outer-pair orientations but not central loss, larger sectors, or arbitrary structural shadows; the all-21-pair tests of the historical S3 vector do not establish an S6 guarantee for those cases. Independent optimizer/SAA restarts and successful termination show numerical repeatability within the tested protocol, not global optimality. The optional divergence increment is unestablished in S3 and positive for the retained fixed vectors in S6, so it is condition-dependent. The link example assumes equal-weight current summation, adaptive midpoint OOK detection, fixed losses, and constant Gaussian noise; it omits coding, signal-dependent noise, and a measured system budget. The calibrated Gaussian envelope and centroid model omit wavefront distortion and joint tilt–spreading–scintillation statistics; separate perturbation checks are not wave-optics validation. Finally, every footprint is static, and the reported intervals describe Monte Carlo estimation for retained designs rather than parameter or deployment uncertainty.')
replace(P[103],
    'Further work should use measured obstruction maps and aperture covariance, independently replicate the extreme-tail optimization with larger budgets, and compare against phase-screen propagation and outdoor or hardware-in-the-loop measurements. A practical design should select its receiver combining, noise model, reliability threshold, and permitted blockage set before choosing between a tail, average-SNR, or outage objective. Adaptive steering and integrated or structured-light receiver controls [14,21,27–31] can be evaluated when terminal calibration and feedback are specified.')
replace(P[105],
    'Risk-aware center placement provides the main nominal robustness gain for the modeled seven-aperture FSO receiver. In S3, the retained joint footprint improves repeated-trial worst-state L5 from 0.04598 to 0.10074, a 119.07% increase with 22.29% less mean collection; its 1.25% increment over center-only design is not statistically established. Expanded blockage tests expose a lower floor for that footprint. Re-optimizing on all six adjacent outer-pair states yields L5 = 0.08741 (95% CI 0.08687–0.08795), versus 0.08069 for center-only optimization at 100 μrad, with paired increment CI [0.00638,0.00705]. Divergence can therefore provide a useful refinement when the protected geometry changes. Multiple optimization seeds, explicit normalization and positive-definiteness checks, convergence histories, and cross-evaluation at three tail fractions make the numerical scope reproducible. The illustrative link comparison shows that better deep-tail power can coexist with worse outage at a particular threshold, so no universal BER or availability gain follows from the tail objective. The framework supports scenario-conditioned footprint design with explicitly stated physical and receiver assumptions.')
replace(P[111],
    f'Code and Data Availability: The original simulation code and numerical CSV outputs are available in the FSO repository cited in [26]: {REPO}. The Code Availability URL and [26] refer to this same repository. Supplementary archive FSO_Light_Green_Revision_Evidence.zip contains the retrieved original-model snapshot, pink_experiments.py and earlier revision scripts, exact optimization and evaluation seeds/settings, all generation histories, raw CSV/JSON outputs, link assumptions, software versions, figures, source-verification notes, and checksums. These extension files are supplied with this revision; no claim is made that they have already been pushed to the public repository. All reported data are Monte Carlo simulations.')

for num,txt in enumerate([
    'Guo, Y.; Wang, F.; Li, R.; Shi, S.; Li, X.; da Costa, D. B.: Optical IRS Assisted-Visible Light Positioning in Indoor Non-LOS IoV Scenarios. IEEE Internet of Things Journal 12(14), 27686–27698 (2025). https://doi.org/10.1109/JIOT.2025.3563896',
    'Guo, X.; Liu, L.; Zhang, J.; Chen, S.: Joint Optimization of Power Allocation and Beamforming for UAV-RIS Communication System. IEEE Wireless Communications Letters 14(6), 1733–1737 (2025). https://doi.org/10.1109/LWC.2025.3554304',
    'Wei, W.; Sun, G.; Zhang, Q.: A random searching algorithm for efficiently solving the connectivity-oriented robust optimization problem on large-scale networked systems. Applied Soft Computing 173, 112924 (2025). https://doi.org/10.1016/j.asoc.2025.112924',
],34):
    p=make_p();P[139]._p.addprevious(p._p);reference(p,num,txt,all_highlight=True)

# Normalize captions and keep short numeric tables together.
for tab in D.tables:
    if len(tab.rows)<=10:
        for ri,row in enumerate(tab.rows):
            for cell in row.cells:
                for p in cell.paragraphs:p.paragraph_format.keep_with_next=(ri<len(tab.rows)-1)
P[57].style='Heading 2'
for p in D.paragraphs:
    if p.text.startswith(('Fig.','Table ')):
        p.style='Caption';p.paragraph_format.keep_with_next=p.text.startswith('Table ')
        p.paragraph_format.keep_together=True
        for r in p.runs:font_run(r,9.5)
    if p._p.xpath('.//w:drawing'):
        p.paragraph_format.keep_with_next=True;p.paragraph_format.space_after=Pt(2)
    if not p.text.strip() and not p._p.xpath('.//w:drawing | .//w:pict | .//m:oMath'):
        p.paragraph_format.space_after=Pt(0);p.paragraph_format.space_before=Pt(0);p.paragraph_format.line_spacing=Pt(2)
        for r in p.runs:r.font.size=Pt(2)
if not P[139].text.strip() and not P[139]._p.xpath('.//w:drawing | .//m:oMath'):
    P[139]._p.getparent().remove(P[139]._p)
D.core_properties.title=P[0].text
