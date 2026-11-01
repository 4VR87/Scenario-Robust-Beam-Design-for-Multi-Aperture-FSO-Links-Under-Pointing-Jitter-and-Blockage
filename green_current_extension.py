"""Four-comment revision: contribution, parameter evidence, states, validation."""
GDATA=ROOT/'green_revision_data'
gsensitivity=pd.read_csv(GDATA/'parameter_sensitivity_summary.csv')
gvalidation=pd.read_csv(GDATA/'independent_validation_summary.csv')
gpaired=pd.read_csv(GDATA/'paired_power_checks.csv')
gquadrature=pd.read_csv(GDATA/'quadrature_checks.csv')
gconfig=json.loads((GDATA/'independent_validation_config.json').read_text())
gsettings=json.loads((GDATA/'parameter_sensitivity_config.json').read_text())

def green_metric(frame,design,**filters):
    selected=frame[frame.design==design]
    for key,value in filters.items():selected=selected[selected[key]==value]
    if len(selected)!=1:raise RuntimeError((design,filters,len(selected)))
    return selected.iloc[0]
def green_ci(row):return f'{row["mean"]:.5f} ({row.ci95_low:.5f}–{row.ci95_high:.5f})'
def green_scientific(value):
    mantissa,exponent=f'{value:.2e}'.split('e')
    return mantissa+' × 10'+str(int(exponent)).translate(str.maketrans('-0123456789','⁻⁰¹²³⁴⁵⁶⁷⁸⁹'))

replace(P[0],'Scenario Robust Beam Center Design for Multi Aperture FSO Receivers with Conditional Divergence Refinement')
for r in P[0].runs:font_run(r,16,bold=True)
replace(P[7],
    'Multi-aperture free-space optical (FSO) receivers can suffer severe fades when pointing jitter and partial receiver loss interact. We formulate static beam-center placement as maximization of the minimum statewise conditional mean of the lowest 5% of collected-power realizations, with divergence as an optional refinement. The technical contribution is a receiver-availability-conditioned lower-tail design criterion, evaluated against mean-power, average-SNR, and outage objectives motivated by recent FSO beam optimization. Exact Gaussian aperture coupling is combined with anisotropic pointing and correlated unit-mean lognormal scintillation. For the original three blocked-pair scenarios, ten independent trials give a worst-state tail mean of 0.10074 versus 0.04598 for mean-power design, a 119.07% increase with a 22.29% mean-collection reduction. Expanded tests cover all 21 aperture pairs; re-optimization over all six adjacent outer pairs yields 0.08741 (95% CI 0.08687–0.08795). Independent direct aperture integration, spectral covariance factors, and scrambled-Sobol sampling reproduce the expanded-set joint-design tail mean as 0.08782, compared with 0.08791 for analytic coupling and conventional Monte Carlo at 65,536 samples per state. Literature establishes physical context for the fading family, while sensitivity tests disclose dependence on uncalibrated scintillation correlation and blockage assumptions. The novelty concerns scenario-conditioned risk allocation, not a new optimizer or propagation law. Numerical agreement supports implementation reliability; experimental calibration and field availability remain unestablished.')

# Explicitly compare the technical objective with recent FSO controls and evidence.
replace(P[12],
    'Beam-coordinate optimization, divergence tuning, and joint multi-beam control are established FSO techniques [11–13,25]. Table 11 compares their controls, objectives, and reported validation with the present formulation. Our contribution is the explicit receiver-state criterion minₛ L0.05,s: each failure state first receives its own lower-tail score, and the weakest score drives static footprint allocation. Neither exact Gaussian coupling, differential evolution, CVaR, nor joint center/divergence control is claimed as individually novel. Blockage-aware FSO analysis also predates this work [39]. The specific question here is how to distribute one Gaussian footprint across separated apertures when their conditional availability differs. Numerical comparisons adapt objective types to the same receiver model; they do not reproduce the hardware architectures of the cited papers.')
comparison_rows=[
    ['Cha et al. [11], 2024','Transmit beam coordinates','Received SNR; statistical-CSI robust pointing; numerical study','Robust pointing already exists; our additional control target is the weakest receiver-state tail.'],
    ['Park et al. [12], 2025','Divergence angle','Gain/link-budget optimization under elliptical pointing; analytical/numerical study','Anisotropy and divergence are established; our criterion conditions on aperture availability.'],
    ['Ling et al. [13], 2026','Multi-beam positions and divergence','Average received power and outage; simulation and UAV-assisted hardware','Joint beam control is established; objective-level comparisons here use one static seven-aperture model.'],
    ['Lim et al. [25], 2026','Phase-plate mode and relative angle','Power efficiency and turbulence resilience; 600 km uplink wave-optics simulation','Mode adaptation uses different optical controls; no cross-system performance ratio is asserted.'],
    ['Present study','Static center; optional divergence','Minimum conditional L₅ over prescribed receiver-loss states; analytic and independent numerical checks','Contribution is state-conditioned lower-tail allocation and quantified mean–tail trade-offs; no field validation.'],
]
table_before(P[13]._p,'Table 11 Positioning against representative recent FSO beam optimization studies',
    ['Work','Beam controls','Criterion and evidence','Relation to this contribution'],comparison_rows,[.98,1.10,2.03,2.24])
replace(P[14],
    'The specific contribution is a receiver-availability-conditioned minimax lower-tail formulation for static beam-center allocation, with divergence evaluated as a condition-dependent refinement. We compare its objective with mean-power, average-SNR, and outage criteria on the same receiver and power/noise assumptions; test additional failure patterns and re-optimize an expanded state set; quantify the effect of uncertain channel stress parameters; and independently verify the numerical channel and tail estimates using direct intensity integration and randomized quasi-Monte Carlo. The contribution concerns the formulation and its numerical evidence; the optimizer and physical propagation assumptions follow established methods.')

# Source evidence is separated from nominal, deliberately uncalibrated stress values.
anchor=P[38]._p
before(anchor,
    'Published measurements support the choice of a weak-scintillation model family without fixing this study’s numerical parameters. Nazari et al. [37] fitted normalized irradiance on a 220 m FSO link with a lognormal law. Liu et al. [38] studied an 865 m four-channel link and reported a scintillation-index reduction from 0.36 to 0.04 through multi-aperture reception, while comparing several fading distributions. These are external experimental contexts, not measurements of our seven-aperture channel. Their aperture averaging and coherence quantities do not identify the phenomenological rc in (7). Our σln = 0.38 gives scintillation index 0.15535 and the weak-Rytov equivalent in Section 2.4; it is a representative design assumption rather than a fitted site parameter.')
before(anchor,
    'Likewise, the Uniform retained-power distributions and 1% micro-event probability are deliberately selected severity controls. Existing blockage-aware FSO work [39] motivates including obstruction risk; it does not establish these distributions or their exact endpoints. Table 12 identifies the evidence status of each quantity. Section 5.14 varies every listed parameter separately using retained footprints, so the conclusions are conditional on the tested ranges rather than on an implied experimental calibration.')
table_before(anchor,'Table 12 Provenance and tested ranges of scintillation and blockage assumptions',
    ['Quantity','Nominal setting','Evidence status','Sensitivity settings'],[
      ['σln','0.38','Lognormal family: [32,37,38]; nominal value is assumed, SI = 0.15535.','0.15, 0.38, 0.65'],
      ['rc','55 mm','Phenomenological latent-correlation scale; no measured covariance fit.','20, 55, 120 mm'],
      ['bblock','Uniform [0,0.12]','Assumed retained-power stress; [39] supports obstruction modelling, not this law.','Upper endpoint 0, 0.12, 0.30'],
      ['Pmicro','0.01 = 1%','Assumed per-aperture, per-realization event rate; no arrival-time calibration.','0%, 1%, 5%'],
      ['bmicro','Uniform [0,0.4]','Assumed retained-power multiplier conditional on a micro-event.','Upper endpoint 0.1, 0.4, 0.8'],
    ],[.72,1.20,2.76,1.67])
before(anchor,
    'Terminal calibration would require simultaneous aperture-power and pointing records: normalize clear-path power for geometric coupling, estimate log-irradiance variances and cross-covariances, measure conditional transmission under representative masks, and estimate micro-event frequency per aperture exposure. The obstruction map should then determine the state set before re-optimization. The present parameter sweep is a sensitivity analysis, not a confidence region for those unmeasured quantities.')
font_run(P[99].add_run(' Section 5.15 extends this local coupling check to a separately implemented end-to-end numerical evaluator, including a different sampling scheme and covariance factorization.'),highlight=True)
before(P[68]._p,
    'The additional checks in Sections 5.14 and 5.15 use new random streams and never reselect the retained footprints. Parameter sensitivity has six independent 20,000-sample/state trials with common random numbers across cases within each trial. Independent numerical validation compares analytic coupling with pseudorandom Monte Carlo against a separate direct-integration/scrambled-Sobol implementation, at 4,096, 16,384, and 65,536 samples per state and eight independent replications per budget. Both evaluations use the same seven S6 states. Their intervals describe repeated numerical estimators for fixed designs; they do not certify the physical fading/blockage laws.')

# Numerical figures contain recorded results, with intervals, rather than illustrations.
fields=[('sigma_ln','σln',1),('corr_length_m','rc (mm)',1000),('block_max','Pair retained-power upper bound',1),
        ('micro_probability','Micro-event probability (%)',100),('micro_max','Micro retained-power upper bound',1)]
f,axs=plt.subplots(2,3,figsize=(6.7,4.35));axs=axs.ravel()
for ax,(field,label,scale) in zip(axs,fields):
    settings={c['case']:c for c in gsettings['scenarios'] if c['family'] in ['nominal',field]}
    for name,caption,color in [('S6_mean','S6 mean power','#99652c'),('S6_joint_alpha05','S6 joint tail','#2e674a')]:
        selected=gsensitivity[(gsensitivity.design==name)&gsensitivity.case.isin(settings)].copy()
        selected['x']=selected.case.map(lambda c:settings[c][field]*scale);selected=selected.sort_values('x')
        ax.errorbar(selected.x,selected['mean'],yerr=[selected['mean']-selected.ci95_low,selected.ci95_high-selected['mean']],
            marker='o',ms=3,lw=1,capsize=2,label=caption,color=color)
    ax.set(xlabel=label,ylim=(.025,.105));ax.tick_params(labelsize=7.5);ax.xaxis.label.set_size(8)
    ax.grid(alpha=.2)
for ax in [axs[0],axs[3]]:ax.set_ylabel('Worst-state L₅',fontsize=8.5)
handles,labels=axs[0].get_legend_handles_labels();axs[5].axis('off')
axs[5].legend(handles,labels,loc='upper left',frameon=False,fontsize=8)
axs[5].text(.02,.47,'Six trials per setting\n20,000 samples per state\n95% t intervals\nFixed footprints',transform=axs[5].transAxes,fontsize=8,va='top',linespacing=1.5)
f.tight_layout();f.savefig(FIG/'Fig12_parameter_sensitivity.png');plt.close(f)

f,axs=plt.subplots(1,2,figsize=(6.7,2.9))
for ax,name,title in zip(axs,['S6_joint_alpha05','S3_legacy_joint'],['S6 joint footprint','S3 joint footprint on S6']):
    for method,label,color,marker in [('Analytic_MC','Analytic coupling + MC','#235782','o'),('Independent_RQMC','2D integration + RQMC','#2e674a','s')]:
        s=gvalidation[(gvalidation.design==name)&(gvalidation.method==method)&(gvalidation.metric=='worst_L5')].sort_values('n')
        ax.errorbar(s.n,s['mean'],yerr=[s['mean']-s.ci95_low,s.ci95_high-s['mean']],marker=marker,ms=4,lw=1.1,capsize=3,label=label,color=color)
    ax.set_xscale('log',base=2);ax.set_xticks([4096,16384,65536],['4,096','16,384','65,536'])
    ax.set(xlabel='Samples per state',ylabel='Mean worst-state L₅',title=title);ax.tick_params(labelsize=8);ax.grid(alpha=.2)
axs[1].legend(fontsize=7,loc='lower right');f.tight_layout();f.savefig(FIG/'Fig13_independent_validation.png');plt.close(f)

anchor=P[101]._p
before(anchor,'5.14 Sensitivity to assumed channel and blockage parameters','Heading 2')
before(anchor,
    'The nominal S6 joint-tail, S6 mean-power, and historical S3 joint footprints are held fixed while changing one parameter at a time. Six independent trials use 20,000 samples in every S6 state, base seed 202610240900 + 1,000,000(t − 1), and state increment 1009. Within a trial/state, the same underlying pointing, Gaussian, and uniform draws are shared across parameter cases and footprints. This isolates changes in σln, latent correlation length, pair-block retained-power endpoint, micro-event probability, and micro-event retained-power endpoint. Each altered case keeps the other four parameters nominal. Table 13 compares the S6 mean-power and joint-tail designs; all three footprints and the raw trial values are retained in the supplementary evidence.')
rows=[];gains=[]
case_labels={
 'Nominal':'Nominal', 'sigma_ln=0.15':'σln = 0.15','sigma_ln=0.65':'σln = 0.65',
 'corr_length_m=0.02':'rc = 20 mm','corr_length_m=0.12':'rc = 120 mm',
 'block_max=0':'Block upper = 0','block_max=0.3':'Block upper = 0.30',
 'micro_probability=0':'Micro event = 0%','micro_probability=0.05':'Micro event = 5%',
 'micro_max=0.1':'Micro upper = 0.1','micro_max=0.8':'Micro upper = 0.8'}
for case in gsettings['scenarios']:
    name=case['case'];mean=green_metric(gsensitivity,'S6_mean',case=name);tail=green_metric(gsensitivity,'S6_joint_alpha05',case=name)
    gain=100*(tail['mean']/mean['mean']-1);gains.append(gain)
    rows.append([case_labels[name],green_ci(mean),green_ci(tail),f'{gain:.1f}%'])
table_before(anchor,'Table 13 Fixed-design parameter sensitivity over S6 with six-trial 95% intervals',
    ['Altered setting','Mean-power design L₅','Joint-tail design L₅','Relative gain'],rows,[1.25,2.03,2.03,1.04])
before(anchor,
    f'The joint-tail mean ranges from {gsensitivity[gsensitivity.design=="S6_joint_alpha05"]["mean"].min():.5f} to {gsensitivity[gsensitivity.design=="S6_joint_alpha05"]["mean"].max():.5f} across these settings, while its ratio-of-means advantage over the mean-power footprint ranges from {min(gains):.1f}% to {max(gains):.1f}%. Increasing σln from 0.15 to 0.65 reduces its tail mean from 0.09707 to 0.06934. Increasing rc from 20 to 120 mm reduces it from 0.09078 to 0.07874, consistent with less spatial-diversity benefit when aperture fluctuations are more correlated. Thus both the estimated gain and the absolute tail floor depend on the assumed scintillation strength and spatial correlation.')
before(anchor,
    'Complete loss of each selected pair (block upper endpoint zero) lowers the joint-tail metric to 0.08104, whereas a 0.30 endpoint gives 0.09631 because more power is retained. Raising micro-event probability from zero to 5% lowers it from 0.08907 to 0.08208. Changing only the micro-event endpoint from 0.1 to 0.8 changes it from 0.08695 to 0.08800. Figure 12 displays the intervals and trends. The nominal ranking persists in this finite set of one-at-a-time tests; simultaneous shifts, different marginal fading laws, time correlation, and site-specific obstruction distributions remain untested by this sweep. No tested interval is an empirical confidence bound or a deployment guarantee.')
image_before(anchor,FIG/'Fig12_parameter_sensitivity.png',6.25,
    'Fig. 12 Sensitivity of fixed S6 footprints to five assumed channel/blockage parameters. Only the horizontal-axis parameter changes in each panel. Points are six-trial mean worst-state L5 values with 95% Student-t intervals; each trial uses 20,000 samples per state. Connecting lines guide the eye.')

before(anchor,'5.15 Independent numerical implementation and validation','Heading 2')
before(anchor,
    'An independent evaluator is implemented in green_independent_solver.py without importing the original simulator, its random-number generator, Marcum-Q/noncentral-chi-square coupling, or partition-based tail function. It integrates the normalized two-dimensional Gaussian intensity directly over each circular aperture. For displacement d, use u = ρ²/a² to write:')
def integral_node(lower,upper,body):
    n=OxmlElement('m:nary');pr=OxmlElement('m:naryPr')
    for tag,val in [('chr','∫'),('limLoc','subSup'),('subHide','0'),('supHide','0')]:
        e=OxmlElement('m:'+tag);e.set(qn('m:val'),val);pr.append(e)
    n.append(pr);n.append(wrap('sub',expr(lower)));n.append(wrap('sup',expr(upper)));n.append(wrap('e',body));return n
exponent=frac(expr('−2[',sup('d','2'),' + ',sup('a','2'),'u − 2da',radical(expr('u')),' cos ψ]'),expr(sup('w','2')))
integrand=expr('exp(',exponent,') dψ du')
equation_before(anchor,expr(sub('g','Q'),'(d,w) = ',frac(expr(sup('a','2')),expr('π',sup('w','2'))),
    integral_node('0','1',expr(integral_node('0','2π',integrand)))),'13')
before(anchor,
    'The numerical rule uses 32 Gauss–Legendre nodes in u and 128 equally spaced azimuths. For each fixed footprint, a shape-preserving cubic interpolant tabulates the directly integrated coupling at 0.2 mm displacement steps over 0–550 mm; extrapolation is rejected. Independent scrambled-Sobol draws in 30 dimensions generate the pointing, scintillation, pair-retention, micro-event, and micro-retention variables. Gaussian covariance factors use symmetric eigendecomposition instead of the source Cholesky path. The independent tail routine fully sorts the samples, rather than invoking the source partition routine. The physical assumptions and retained parameter vectors are identical; the numerical routes differ.')
before(anchor,
    f'At 276 off-grid displacement points for each of four footprints, the largest direct-integration difference between 32 × 128 and 64 × 256 rules is {green_scientific(gquadrature.mid_high_max_abs.max())}. The refined rule differs from the analytic coupling by at most {green_scientific(gquadrature.high_analytic_max_abs.max())}; tabulation/interpolation contributes at most {green_scientific(gquadrature.interpolated_analytic_max_abs.max())} in coupling fraction on those checks. A paired 8,192-sample test in each of seven states gives maximum collected-power difference {green_scientific(gpaired.max_abs_power_error.max())} and maximum statewise L5 difference {green_scientific(float((gpaired.primary_L5-gpaired.independent_L5).abs().max()))}. Geometry and both covariance matrices agree entrywise between implementations. These checks address transcription, integration, interpolation, and tail-statistic consistency.')
before(anchor,
    'Fresh analytic/Monte Carlo and independent quadrature/randomized-quasi-Monte Carlo evaluations then use eight replications at each of 4,096, 16,384, and 65,536 samples per state. Monte Carlo base seeds are 202610240600 and Sobol base seeds 202610240700, with replication increments 1,000,000 and state increments 1009. Sobol budgets are powers of two and preserve the first points of the same scramble across budget levels. Each method shares its draws across retained designs. Table 14 gives the largest-budget estimates; Fig. 13 shows both budget sequences. Student-t intervals use seven degrees of freedom across independent trials or independent scrambles, respectively.')
rows=[];relative=[];independent_relative_percent={}
for name,label in [('S3_legacy_mean','S3 mean-power footprint'),('S3_legacy_joint','S3 joint-tail footprint'),('S6_center_alpha05','S6 center-only footprint'),('S6_joint_alpha05','S6 joint-tail footprint')]:
    a=green_metric(gvalidation,name,n=65536,method='Analytic_MC',metric='worst_L5')
    b=green_metric(gvalidation,name,n=65536,method='Independent_RQMC',metric='worst_L5')
    percent=100*(b['mean']/a['mean']-1);relative.append(abs(percent))
    independent_relative_percent[name]=percent
    rows.append([label,green_ci(a),green_ci(b),f'{percent:+.3f}%'])
table_before(anchor,'Table 14 Independent numerical comparison on S6 at 65,536 samples per state',
    ['Fixed footprint','Analytic coupling and MC L₅','Direct integration and RQMC L₅','Relative difference'],rows,[1.45,1.95,1.95,1.0])
before(anchor,
    f'The independent joint-design mean is 0.08782, compared with 0.08791 for analytic coupling and ordinary Monte Carlo; their relative difference is {independent_relative_percent["S6_joint_alpha05"]:.3f}%. Across all four footprints the largest absolute relative difference is {max(relative):.3f}%, and every independent mean lies inside the corresponding analytic-Monte-Carlo 95% interval at the largest budget. Both approaches preserve the footprint ordering at all three budgets. Agreement and smaller observed scramble dispersion support numerical reliability in this test; no general quasi-Monte-Carlo convergence rate is asserted for the non-smooth tail/minimum functional.')
before(anchor,
    'This is independent computational validation of the stated Gaussian/lognormal/blockage model. It cannot validate whether that model or its parameter values describe a real atmosphere, obstruction process, transmitter, or receiver. The external experiments in [13,37,38] are literature evidence and are not re-labelled as validation data for our terminal. No matching measured multi-aperture record or hardware data are used here. Wave-optics and calibrated field validation therefore remain a separate requirement for deployment claims.')
image_before(anchor,FIG/'Fig13_independent_validation.png',6.25,
    'Fig. 13 Independent numerical validation of two fixed footprints on all S6 states. Analytic coupling with conventional Monte Carlo is compared with direct two-dimensional aperture integration and scrambled-Sobol sampling. Error bars show 95% intervals over eight independent replications at each sample budget; no footprint is retuned.')

# Update headline scope and conclusions to include the new, bounded evidence.
replace(P[102],P[102].text+
    ' The independent implementation in Section 5.15 verifies numerical calculations under the same assumed laws; it does not replace experimental or alternative-propagation-model validation. Exact nominal blockage probabilities and correlation lengths remain uncalibrated, and the one-at-a-time parameter study does not cover simultaneous uncertainty or temporal blockage dynamics.')
replace(P[105],
    'The technical contribution is a receiver-availability-conditioned lower-tail criterion for static FSO beam-center allocation, with condition-dependent divergence refinement. Recent FSO studies already optimize beam positions, divergence, or adaptive modes; here each prescribed aperture-loss state has its own lower-tail score, and the weakest score determines the footprint. The original S3 comparison gives a 119.07% repeated-trial tail improvement with a 22.29% mean-collection reduction, while its 1.25% divergence increment is not statistically established. All-21-pair stress tests expose the limits of that scenario set. Re-optimization on S6 yields L5 = 0.08741 (95% CI 0.08687–0.08795), and separate parameter sweeps quantify dependence on scintillation correlation and assumed blockage severity. A new independent direct-integration/Sobol evaluator gives 0.08782, compared with 0.08791 for analytic coupling/Monte Carlo at 65,536 samples per state. Agreement supports implementation reliability without validating the physical assumptions. The model families have literature support, but the exact correlation and blockage parameters remain design stressors. These results support scenario-conditioned optical design and transparent objective comparisons, rather than experimental performance, universal communication gains, or field-availability predictions.')
replace(P[111],
    f'Code and Data Availability: The original model and numerical outputs are available in the FSO repository [26]: {REPO}. Supplementary archive FSO_Light_Green_Revision_Evidence.zip contains the original-model snapshot, retained design vectors, all added optimization/evaluation scripts and raw outputs, green_independent_solver.py, green_experiments.py, parameter-sensitivity trials, independent-validation trials, quadrature checks, exact settings and seeds, software versions, figures, and source-verification notes. The added studies are supplied with this revision and have not been pushed to the public repository. All new performance records are simulations; no measured terminal data are claimed.')

for num,txt in enumerate([
 'Nazari, Z.; Gholami, A.; Vali, Z.; Sedghi, M.; Ghassemlooy, Z.: Experimental investigation of scintillation effect on FSO channel. In: 2016 24th Iranian Conference on Electrical Engineering (ICEE), pp. 1629–1633. IEEE (2016). https://doi.org/10.1109/IranianCEE.2016.7585782',
 'Liu, S.; Tian, S.; Li, B.; Liu, Z.; Yao, H.: Statistical Study of Free-Space Optical Transmission Using Multi-Aperture Receivers Under Real-Measured Atmospheric Turbulence. Photonics 13(1), 63 (2026). https://doi.org/10.3390/photonics13010063',
 'Djordjevic, G. T.; Petkovic, M. I.; Spasic, M.; Antic, D. S.: Outage capacity of FSO link with pointing errors and link blockage. Optics Express 24(1), 219–230 (2016). https://doi.org/10.1364/OE.24.000219',
],37):
    p=D.add_paragraph();reference(p,num,txt,all_highlight=True)

# Number comparison and provenance tables in their actual manuscript order.
table_number_map={}
for p in D.paragraphs:
    match=re.match(r'^Table (\d+)\b',p.text)
    following=p._p.getnext()
    if match and following is not None and following.tag==qn('w:tbl'):
        old=int(match.group(1))
        if old in table_number_map:raise RuntimeError(f'Duplicate caption number {old}')
        table_number_map[old]=len(table_number_map)+1
math_before_renumber=len(D.element.xpath('.//m:oMath'))
for p in D.paragraphs:
    # Match across Word runs, but change only the numerical text spans. This
    # retains inline Office Math, hyperlinks, fields, and run formatting.
    nodes=p._p.xpath('.//w:t');starts=[];offset=0
    for node in nodes:
        starts.append(offset);offset+=len(node.text or '')
    combined=''.join(node.text or '' for node in nodes)
    for match in reversed(list(re.finditer(r'\bTable\s+(\d+)\b',combined))):
        old=match.group(1);new=str(table_number_map[int(old)])
        if old==new:continue
        lo,hi=match.span(1)
        indices=[i for i,node in enumerate(nodes) if starts[i]<hi and starts[i]+len(node.text or '')>lo]
        first,last=indices[0],indices[-1]
        if first==last:
            node=nodes[first];value=node.text or ''
            node.text=value[:lo-starts[first]]+new+value[hi-starts[first]:]
        else:
            nodes[first].text=(nodes[first].text or '')[:lo-starts[first]]+new
            nodes[last].text=(nodes[last].text or '')[hi-starts[last]:]
            for i in indices[1:-1]:nodes[i].text=''
        for i in indices:
            parent=nodes[i].getparent()
            if parent.tag==qn('w:r'):
                props=parent.find(qn('w:rPr'))
                if props is None:props=OxmlElement('w:rPr');parent.insert(0,props)
                marker=OxmlElement('w:highlight');marker.set(qn('w:val'),'yellow');props.append(marker)
assert len(D.element.xpath('.//m:oMath'))==math_before_renumber
for p in D.paragraphs:
    following=p._p.getnext()
    is_table_caption=p.text.startswith('Table ') and following is not None and following.tag==qn('w:tbl')
    if p.text.startswith('Fig.') or is_table_caption:
        p.style='Caption';p.paragraph_format.keep_together=True;p.paragraph_format.keep_with_next=p.text.startswith('Table ')
        for r in p.runs:font_run(r,9.5)
    elif p.text.startswith('Table ') and p.style.name=='Caption':
        p.style='Normal';p.paragraph_format.keep_together=False;p.paragraph_format.keep_with_next=False
        for r in p.runs:font_run(r,11)
    if p._p.xpath('.//w:drawing'):
        p.paragraph_format.keep_with_next=True;p.paragraph_format.space_after=Pt(2)
for tab in D.tables:
    if len(tab.rows)<=13:
        for ri,row in enumerate(tab.rows):
            for cell in row.cells:
                for p in cell.paragraphs:p.paragraph_format.keep_with_next=(ri<len(tab.rows)-1)
D.core_properties.title=P[0].text
(OUT/'table_number_map.json').write_text(json.dumps(table_number_map,indent=2))
