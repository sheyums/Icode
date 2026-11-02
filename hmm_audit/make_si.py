# Builds SupplementaryInformation_draft.docx: checks, controls and robustness tests with their outcomes (2026-11-02).
# Every number is copied from the ledgers (hmm_audit/RESULTS_LEDGER.md, RESULTS_LEDGER_planB.md, PORT_VALIDATION.md).
import os, docx
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
F = '/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit/figures/'
C8 = '/home/user/Icode/analysis_20r/figC8/figC8_wake_hazard.png'
d = docx.Document()
st = d.styles['Normal']; st.font.name = 'Arial'; st.font.size = Pt(10)
nfig = [0]; ntab = [0]
def H(t, l=1): d.add_heading(t, level=l)
def P(t, it=False):
    p = d.add_paragraph(); r = p.add_run(t); r.italic = it; return p
def FIG(path, cap, w=6.3):
    nfig[0] += 1; d.add_picture(path, width=Inches(w)); d.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
    p = d.add_paragraph(); r = p.add_run(f'Fig. S{nfig[0]}. '); r.bold = True; p.add_run(cap).font.size = Pt(9)
def TAB(cap, head, rows):
    ntab[0] += 1; p = d.add_paragraph(); r = p.add_run(f'Table S{ntab[0]}. '); r.bold = True; p.add_run(cap)
    t = d.add_table(rows=1, cols=len(head)); t.style = 'Light Grid Accent 1'
    for i, h in enumerate(head): t.rows[0].cells[i].text = h
    for row in rows:
        c = t.add_row().cells
        for i, v in enumerate(row): c[i].text = str(v)
    for row in t.rows:
        for c in row.cells:
            for pp in c.paragraphs:
                for rr in pp.runs: rr.font.size = Pt(8.5)
    d.add_paragraph()
d.add_heading('Supplementary Information (draft, 2 November 2026)', 0)
P('Checks, controls and robustness tests behind the main results, with their outcomes. Each analysis was specified in a '
  'written plan before it was run; post hoc analyses are labelled as such. Null results are reported alongside positive ones.', it=True)

H('S1. Three behavioural states, and the hidden Markov model')
P('Models with two, three and four hidden states were fitted to every fly. The two-state model is the classic sleep/wake '
  'description: one wake state, in which quiet and active minutes are noisy readings, and one sleep state. In the '
  'four-state model the data choose the extra state. Three states had the lowest BIC in 76–96% of flies of every genotype.')
FIG(F + 'figS_hmm_order_bic.png', 'Hidden-state order selected by BIC per fly, with n = minutes (left) and n = symbol changes (right). 187 flies.')
TAB('Order selection and the identity of the extra state in the four-state model.',
    ['Genotype', 'Best K = 2 / 3 / 4 (n = minutes)', 'n = changes', 'Median BIC(2) − BIC(3)', 'Median BIC(4) − BIC(3)', '4-state model: extra state'],
    [['CS', '0 / 96 / 4%', '0 / 88 / 12%', '+334', '+71', 'second active wake (17/25)'],
     ['per0 LD', '15 / 78 / 7%', '4 / 89 / 7%', '+273', '+81', 'second active wake (17/27)'],
     ['per0 DD', '8 / 76 / 16%', '4 / 56 / 40%', '+213', '+22', 'second sleep (16/25)'],
     ['fumin', '6 / 85 / 9%', '2 / 72 / 26%', '+221', '+60', 'second active wake (47/53)'],
     ['sleepless', '5 / 95 / 0%', '5 / 91 / 4%', '+187', '+77', 'second active wake (49/57)']])
P('This establishes three dynamically distinct states. It does not test whether activity within wake falls into two '
  'discrete levels: the quiet/active boundary is a threshold (θ = 0.5 of each fly\'s typical wake activity).')
FIG(F + 'fig1_hmm_emissions.png', 'Median emission matrices of the three-state hidden Markov model: each hidden state emits mostly one symbol (quiet 80–94%, active 91–96%, sleep 100%).')
P('Whether the three-state hidden model improves on the plain observed chain is a separate question: it does not in '
  '72–77% of flies by BIC, its likelihood gain (about 3 nats per 1000 min) is matched by one extra minute of memory, and '
  'allowing slow (6-hour) drift of the rates removes almost none of it.')

H('S2. Robustness of the quiet/active split')
P('Every direction of circulation that was significant at θ = 0.5 had the same sign at θ = 0.33 and 0.75; only '
  'near-zero, non-significant blocks changed sign, and magnitudes grew with θ. The reversal in sleepless by day held at '
  'all three thresholds. Recomputing the symbols on a 30-s grid preserved the directions (Fig. S3).')
FIG(F + 'fig17_resolution_30s_vs_1min.png', 'Circulation at 30-s against 1-min resolution.')

H('S3. Agreement between the observed chain and the hidden model')
FIG(F + 'figA1_direction_flux_SI_hmm.png', 'Observable current, cycle flux and entropy production per fly, with the hidden model\'s median (grey triangles). For sleepless, whose median is zero because many flies have no qualifying transitions, the mean is shown (diamonds).')
P('Direction: the hidden model\'s circulation has the same sign as the observed one in 81–100% of flies, including the '
  'reversal in sleepless by day (Fig. S5).')
FIG(F + 'fig_sss_reversal.png', 'Sleepless runs the cycle backwards in the light phase. Per fly, observed-chain flux (circle) joined to the hidden model\'s (triangle). Flies whose observed chain carries no cycle (|J| < 10⁻¹²) are plotted at 0 and excluded from the percentage and the test: light phase 13 of 50, dark 5 of 53.', w=5.2)
TAB('Onsager conductance and affinity per turn, observed and hidden (fly bootstrap, 1000 resamples).',
    ['Genotype', 'L observed [95% CI]', 'L hidden (path counts)', 'L hidden (model rates)', 'A_eff obs / hidden (nats/cycle)', 'Chain affinity A'],
    [['CS', '0.27 [0.23, 0.30]', '0.22 [0.17, 0.23]', '0.24 [0.19, 0.25]', '0.57 / 1.10', '0.47'],
     ['CS DD', '0.26 [0.22, 0.31]', '—', '—', '0.51 / —', '0.50'],
     ['per0 LD', '0.29 [0.26, 0.32]', '0.25 [0.23, 0.27]', '0.25 [0.23, 0.29]', '0.43 / 0.95', '0.42'],
     ['per0 DD', '0.35 [0.29, 0.40]', '0.33 [0.20, 0.36]', '0.34 [0.22, 0.37]', '0.64 / 0.97', '0.64'],
     ['fumin', '0.24 [0.19, 0.26]', '0.24 [0.17, 0.28]', '0.22 [0.17, 0.24]', '1.11 / 1.45', '1.09'],
     ['sleepless', '0.15 [0.12, 0.18]', '0.14 [0.12, 0.18]', '0.11 [0.08, 0.15]', '0.70 / 1.10', '0.79']])
P('L in cycles² nats⁻¹ h⁻¹. No hidden models were fitted to CS in constant darkness. Sleepless has the lowest conductance '
  'at both levels. The observed A_eff agrees with the affinity of the fitted chain, as single-cycle thermodynamics requires; '
  'the hidden values are about twice as large because observed dissipation is a lower bound. Genotype differences in L '
  'were tested by a fly bootstrap: sleepless < CS (p < 0.001, BH q 0.003); per0 DD > CS only nominally (p 0.026, q 0.065). '
  'A Wald test on the slope was not used: in 400 simulated pairs of genotypes with equal L it rejected in 22% of cases at '
  'nominal 5%, against 5.2% for the bootstrap difference.')
TAB('Energy landscape, observed chain against hidden model, per fly.',
    ['Genotype, phase', 'E_S − E_Q median obs / hidden', 'Spearman ρ', 'Routing prediction: same sign / ρ'],
    [['CS light', '−0.94 / −1.08', '0.95', '84% / 0.71'], ['CS dark', '−1.78 / −1.85', '0.83', '88% / 0.90'],
     ['per0 LD light', '−1.74 / −1.73', '0.92', '93% / 0.58'], ['per0 LD dark', '−1.58 / −1.61', '0.97', '93% / 0.79'],
     ['per0 DD', '−1.75 / −1.85', '0.88', '100% / 0.69'], ['fumin light', '0.64 / 0.74', '0.90', '91% / 0.78'],
     ['fumin dark', '0.97 / 1.07', '0.91', '86% / 0.68'], ['sleepless light', '1.19 / 1.13', '0.92', '70% / 0.56'],
     ['sleepless dark', '1.38 / 1.38', '0.94', '70% / 0.51']])
P('The energies are robust to the quiet/active boundary in every genotype. The branching-corrected routing prediction, which '
  'depends on the rare sleep-to-quiet transition, is robust except in sleepless, where it is not interpreted.')

H('S4. Memory of sleep length')
FIG(F + 'figS_C7_memory_with_hmm_check.png', 'Exit route against sleep-bout length, observed (a) and hidden (b); CS at night (c).')
P('In CS at night, sleep bouts that end through quiet wake are shorter than those ending in active wake (mean difference in '
  'log length −0.29; within-fly permutation p < 0.001, Holm < 0.001). The same holds at the hidden level (Holm 0.045, at the '
  'resolution limit of 400 simulations). per0 LD in the light phase is significant but marginal: Holm 0.04 at the '
  'pre-registered 2000 permutations, and 0.033 with 20000.')

H('S5. Wake-bout hazard and its controls')
FIG(C8, 'Wake-bout hazard at the four pre-set bin counts, per0 (top) and CS (bottom) in constant darkness; dashed: best two-exponential model, whose hazard can only fall.')
TAB('Rise test and controls.', ['Analysis', 'Outcome'],
    [['per0 DD rise (999 null datasets)', 'RiseNullP 0.001 / 0.009 / 0.059 / 0.086 at 12 / 16 / 20 / 24 bins: supported but bin-sensitive'],
     ['CS DD rise', 'RiseNullP ≥ 0.14 at all counts; no rise in either 12-h half'],
     ['Power to see a per0-sized rise in CS DD (12 bins)', '1.00 at n = 3309; 0.83 and 0.77 at the half sizes'],
     ['Renewal (burst–pause) null', 'p 0.015 / 0.015 / 0.050 / 0.050: inconclusive by the 3-of-4 rule; renewal makes less of a rise'],
     ['Time-reversible hidden-model null', 'not usable: best model (≤ 6 hidden states) gave 22.5% more wake bouts than observed (gate 15%)']])

H('S6. Beneath the observed cycle')
FIG(F + 'figA4_beneath_the_surface.png', 'Hidden layer against the observed chain (a, b), extra irreversibility from 3-minute memory (c), and the thermodynamic uncertainty relation (d).')
P('Extra irreversibility from memory beyond one minute qualified in per0 DD (+0.06 nats/h, also in sequences without sleep) '
  'and in fumin at night (+0.04 nats/h, only in sequences containing sleep). The thermodynamic uncertainty relation showed no '
  'dissipation beyond that of the observed chain (ratio 0.90 in per0 DD, 0.75 in CS DD; no fly above its own simulated null).')
P('Waiting times: dwell times in quiet and active wake depend on the direction of passage beyond the observed chain in per0 DD '
  '(+0.083 nats/h; Holm 0.0006; 28% of flies individually) and CS DD (+0.078; Holm 0.0002; 32%). The scoring is exactly '
  'time-symmetric (the reversed 1-s record gave the reversed symbols in all 275,040 minutes), so the sleep criterion cannot '
  'produce the effect from reversible activity. A stricter control with time-reversible models of the 1-s activity could not '
  'be built: such models gave 22.5% (per0 DD) and 40% (CS DD) more wake bouts than observed.')
FIG(F + 'fig_dwell_direction_CS_night.png', 'Where the direction dependence lies (descriptive): CS at night. Active-wake stays after waking last longer than the mirror-image stays that end in sleep (4.4 against 3.5 min, longer in 80% of flies); quiet-wake stays do not differ.')
FIG(F + 'fig_planB21_profiles_primary.png', 'Activity at the two ends of wake bouts, 1-s resolution (descriptive). Flies wind down over about 5 minutes before sleep and return to their usual activity within about 20 s of waking. The asymmetry was significant against per-bout and day-block sign-flip nulls (per0 DD Holm 0.014; CS DD 4 × 10⁻⁵). It confirms the observed cycle in raw counts and is not evidence of irreversibility beyond it.')

H('S7. Clock-locked modulation of the cycle')
P('The test compares the daily amplitude of the pooled cycle flux with a null that rotates each fly\'s record by an '
  'independent random offset. Calibration on 200 constant-rate datasets matched to per0 DD: 4.5% rejections at nominal 5%. '
  'Power at per0 DD\'s sample size against CS LD\'s amplitude: 0.80. Results: CS LD amplitude 0.077 (Holm 0.015), CS DD 0.063 '
  '(Holm 0.030), per0 LD 0.071 (Holm 0.059, not detected), per0 DD 0.006 (p 0.96). The first pre-specified summary, the range '
  'of the twelve bins, had no power (its band lay above the point estimate even in per0 DD) and was not used for inference.')

H('S8. Sleep deprivation: a null result')
P('26 CS males, 9 h of deprivation (ZT 15–24 of day 4) after three baseline days. Morning sleep (ZT 0–6) rose on the first '
  'recovery day (+42 min), but it also rose between baseline days before any deprivation (+28 min, p 0.007), so recovery could '
  'not be separated from this drift (excess over the drift +11.6 min, p 0.094). No readout of the cycle changed beyond the '
  'day-to-day variation of the same flies: observable current, flux and entropy production (paired tests against the placebo '
  'contrast p 0.62–0.96), the six transition probabilities (post hoc, all Holm ≥ 0.2), and wake activity (−0.017 after '
  'deprivation and −0.017 between baseline days). A repeat with undeprived controls recorded in parallel is needed.')

H('S9. Silencing of dFB and VNC neurons')
FIG(F + 'figC12_kir_silencing.png', 'Constitutive Kir2.1 silencing against the enhancerless-GAL4 control (19, 25 and 34 flies; one run, one incubator).')
P('Sleep fell by about 40% after VNC silencing (light 38.5 → 23.9, dark 42.6 → 25.3 min/h; Holm 2 × 10⁻⁸) and by 19% at night '
  'after dFB silencing (42.6 → 34.7; Holm 0.001). Silenced flies fell asleep from quiet wake less readily and woke more '
  'readily. No change in observable current, flux or entropy production was detected (12 tests, all Holm ≥ 0.66), in a '
  'background whose cycle is weak even in controls. Monitor check: differences between monitors within a genotype were small '
  '(largest 3.6 min/h, against 14.6 for the VNC effect). Excluding the first day changed nothing. The dependence of the '
  'exit route on sleep length was present in all three groups and unchanged by silencing (all Holm ≥ 0.32).')

H('S10. Statistical and numerical notes')
P('No-cycle fly–phases (|J| < 10⁻¹²) are counted, set to zero and never assigned a sign. Infinite values (a transition '
  'seen in one direction only) are ranked highest in medians. Edges seen in one direction only enter the entropy production '
  'through a rule-of-three bound. Wilcoxon signed-rank tests with ties use the normal approximation; the exact and '
  'approximate p differ by less than 0.02 in the one affected comparison (sleepless dark, p ≈ 0.66).')

H('S11. Independent reimplementation')
P('All analyses of the observed chain, the hidden model, the thermodynamic quantities, block irreversibility, the uncertainty '
  'relation, sleep-length memory and the wake-bout hazard were reimplemented independently in MATLAB and checked against the '
  'original computation: 21 of 21 reference tests and every reported group median of the observable analysis; per-fly '
  'hidden-model log-likelihoods to 2.5 × 10⁻¹⁰; all 160 deterministic numbers of the irreversibility and uncertainty-relation '
  'analyses exactly, and 180 Monte Carlo statistics within sampling error, with every verdict identical.')
d.save('SupplementaryInformation_draft.docx'); print('figures', nfig[0], 'tables', ntab[0])
