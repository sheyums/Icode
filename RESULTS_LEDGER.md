# Results ledger

**Append-only.** One entry per result, **including null and inconclusive
results**. Do not delete or rewrite an entry; add a correcting entry that
references it.

The point of recording nulls is that selection cannot be audited from the
successes alone. A ledger with only positive findings gives no way to know
whether a p-value of 0.03 was the first analysis run or the twentieth.

Rule 3 in `CLAUDE.md` governs this file. Each entry names the `PRESPEC.md`
entry it answers.

## Template

```
### YYYY-MM-DD  <short name>

Commit          Hash of the code that produced it (git log -1 --format=%H).
Prespec         Which PRESPEC.md entry this answers.
Time window     The narrowest block the number was computed from.
n               Sample size actually used, after truncation at xmin.
Result          The number, with its uncertainty. State the units.
Matched null    The same statistic under a sample-size-matched null.
Verdict         Support / against / inconclusive, per the pre-specified rule.
Not reportable  Anything that must NOT be quoted from this run, and why.
```

## Entries

<!-- Append below. Never edit an entry above. -->

### (none yet)

Earlier per0 DD results are in `FINDINGS.md`, written before this ledger
existed. They are not retrofitted here: back-filling a ledger from memory
would give it the appearance of a contemporaneous record without the
substance, which is worse than an honestly empty file.

### 2026-09-15  per0 DD wake: does any single fly carry the hazard rise?

```
Commit          35ee73a358c741f6a1b7209f8b7f35a58585d7d1 (Icode main; code unchanged since 8565082)
Script          wake_loo_rise.m, session scratch, not in the repo
Prespec         2026-09-15 "per0 DD wake: does any single fly carry the hazard rise?"
Time window     The whole analysed record of each fly (2.8-4.0 days, constant darkness). Every
                number was computed on 24 flies pooled, or on 3989 - n_i bouts after random removal.
n               3989 wake bouts from 25 flies; xmin 2 s, dt 1 s, 1-s bridged scoring. 3888 bouts
                after the one influential fly was removed.
Result          Full data RiseRatio 1.405, 1.388, 1.376, 1.430 at 12, 16, 20 and 24 bins.
                One fly, row 14 (n = 101), fell below the 2.5th percentile of its own random
                removals at all four counts; no other fly fell below at more than two of them.
                Removing fly 14 left RiseRatio 1.351, 1.326, 1.323, 1.370 against the null medians
                of that same run, 1.090, 1.139, 1.188, 1.240, so above the median at every count.
                RiseNullP there was 0.003, 0.025, 0.137 and 0.187, with NullDropped 0 throughout.
                Removing fly 10 (n = 239) instead RAISED the ratio to 1.684-1.936, so that fly
                damps the rise rather than causing it.
Matched null    For influence: 199 removals of the same number of bouts as the fly, drawn uniformly
                from all 3989 (RandomSeed 1000 + fly), on the same bins. For the p-values after
                removal: 999 draws from the two-exponential fit to the remaining bouts.
Verdict         Per the pre-specified rule: some flies matter, but the rise does not depend on them
                alone, and no single fly carries it. For the design question, a pooled wake hazard
                in the fuller datasets can be read as such, with per-fly influence checked rather
                than assumed.
Not reportable  The p-values after removing fly 14 are not a test of the rise in the full data, and
                they straddle 0.05 across bin counts (0.003 to 0.187), so no single one should be
                quoted. Nothing here says whether the rise lives in the flies or is produced by
                scoring bursty movement into wake bouts.
```

### 2026-09-15  per0 DD sleep: pooled K=3 against a per-fly null

```
Commit          35ee73a358c741f6a1b7209f8b7f35a58585d7d1 (Icode main; code unchanged since 8565082)
Script          sleep_pooled_null.m, session scratch, not in the repo
Prespec         2026-09-15 "per0 DD sleep: pooled K=3 against a per-fly null"
Time window     Per-fly fits: the whole analysed record of each fly (2.8-4.0 days, constant
                darkness). The tested statistic: all 25 flies pooled.
n               3987 sleep bouts from 25 flies; xmin 300 s, dt 1 s, 1-s bridged scoring. Per fly
                80 to 249 bouts, median 159.
Result          Observed logL -31626.8558 at K = 2 and -31618.4952 at K = 3, so LR = 16.72,
                reproducing FINDINGS. The observed K = 3 fit is the FINDINGS one: tau 54.6, 1040.5
                and 3141.2 s with observed weights 0.0496, 0.9298 and 0.0206. Fitted one fly at a
                time, only 10 of the 25 flies supported two timescales; the other 15 came back
                degenerate at K = 2 and entered the null with one. Pooling 25 such simulated flies
                gave LR well ABOVE the observed value: median 26.68, 95th percentile 48.14, maximum
                58.94 over B = 199, with no negative LR and no failed fit. p = 0.86. The slowest of
                the three fitted components held under 5 percent of bouts in 11 percent of the
                replicates.
Matched null    199 pooled datasets, each with the same 25 flies at the same bout counts, every fly
                drawn from its own fitted model. The staged rule stopped at 199 because p > 0.12.
Verdict         The pre-specified reading is supported: differences between flies are enough to
                produce the pooled third component, and more than enough. A pooled mixture fit on
                these data says nothing about how many timescales an individual fly has, so the
                fuller datasets are to be analysed fly by fly or with a hierarchical model.
Not reportable  Not evidence that no fly has a third timescale: the test asks what pooling can
                manufacture, not what individuals do. The null is if anything too generous, since it
                simulates from unshrunk per-fly fits whose estimation noise widens the spread
                between flies by an unknown amount; a hierarchical, shrunk null would be the sharper
                comparison. The bouts are bridged, so merging sits inside both the data and the
                per-fly fits.
```

### 2026-09-15  Correction to "per0 DD sleep: pooled K=3 against a per-fly null"

```
Commit          35ee73a358c741f6a1b7209f8b7f35a58585d7d1 (no new computation; wording only)
Corrects        2026-09-15 "per0 DD sleep: pooled K=3 against a per-fly null", recorded in 84a5187
What was wrong  That entry said a pooled mixture fit "says nothing about how many timescales an
                individual fly has". Too strong, as Sheyum pointed out.
Correct reading A pooled fit does constrain the population: K > 1 means the flies are not all one
                exponential at a common rate. What it cannot do is attribute the extra timescales,
                because differences between flies and structure within a fly produce the same pooled
                distribution. On these data the third component needs nothing beyond differences
                between flies (p = 0.86), so that component cannot be read as a third timescale in
                any fly. One direction does carry: mixing cannot make a hazard rise, so a
                non-monotone pooled hazard does imply non-monotone hazards in some flies. That is
                what the wake entry rests on.
Separate        The per-fly evidence is a different analysis from the pooled fit. Fitted one at a
evidence        time, 10 of the 25 flies supported two timescales. That count is the fitter
                degeneracy flag at about 160 bouts a fly, not a test, and at that bout count an
                unresolvable second timescale and a genuine rare component look alike, as the 219-fly
                LD work measured.
Unchanged       The numbers and the verdict on the pre-specified question stand: p = 0.86, and the
                unit of analysis for the fuller datasets is the fly or a hierarchical model.
```

### 2026-09-15  per0 DD wake: the hazard rise at four fixed bin counts

```
Commit          172315b (Icode main; analysis code unchanged since 8565082)
Script          wake_rise_fixed_bins.m, session scratch, not in the repo
Prespec         2026-09-15 "per0 DD wake: the hazard rise at four fixed bin counts"
Time window     The whole analysed record of each of the 25 flies (2.8-4.0 days, constant
                darkness), pooled. That pooled set is the narrowest block here.
n               3989 wake bouts; xmin 2 s, dt 1 s, 1-s bridged scoring.
Result          Null K = 2 fit: tau 8.4 and 661.9 s, observed weights 0.079 and 0.921,
                logL -29328.2460.
                bins 12: trough 119.8 s, peak 1094.0 s, RiseRatio 1.405, null median 1.083,
                         95th 1.210, RiseNullP 0.001
                bins 16: trough 109.3 s, peak  997.7 s, RiseRatio 1.388, null median 1.136,
                         95th 1.285, RiseNullP 0.010
                bins 20: trough 161.0 s, peak  944.3 s, RiseRatio 1.376, null median 1.186,
                         95th 1.390, RiseNullP 0.068
                bins 24: trough 144.1 s, peak  910.0 s, RiseRatio 1.430, null median 1.235,
                         95th 1.512, RiseNullP 0.101
                NullDropped 0 at every count; one bin dropped and one saturated at each, as
                expected. The values at 16 and 20 bins reproduce FINDINGS exactly.
                The observed ratio barely moves across bin counts (1.376 to 1.430). What moves is
                the NULL: its median climbs from 1.083 to 1.235 and its 95th percentile from 1.210
                to 1.512, because more bins put fewer bouts in each and a spurious peak-over-trough
                ratio grows.
Matched null    999 datasets of n = 3989 drawn from that two-exponential fit, each put through the
                whole estimator on its own bin edges.
Verdict         Supported but bin-sensitive, the wording fixed on 2026-09-14: a rise at 12 and 16
                bins, none at 20 and 24. Every count reported, not the smallest p.
Not reportable  Do not quote 0.001 on its own. The bin sensitivity here is a loss of power rather
                than a wobble in the data, since the observed ratio is stable while the null widens;
                that reading is an observation from these four counts, not a separate test. The test
                is against a falling-hazard null only: it says nothing about whether the rise is in
                the flies or is produced by scoring bursty movement into wake bouts.
```

### 2026-10-31  20r gate: per0 DD re-extraction FAILS (16 of 25 flies retained)

```
Commit          49ffe9e05e8019cb4f41962f2ba595559d14edd0 (Icode branch 20r-lclaude)
Script          analysis_20r/gate_20r.m; output analysis_20r/gate_20r_report.txt
Prespec         2026-10-31 "20r gate amendment: stored bout files not found" (replaces the 20r-A gate)
Code            Rockefeller Matlab scripts, git HEAD 001deca, current versions (SHA-256 prefixes):
                process_LocomToSleep 7aec365d5b79c859, sleeptrace 0deb702eed5cbfcc,
                ironout f22f91c262f8c370, segment_sleep_phases a627f698159bc97a,
                extractBoutsForModelCompetition 64f6f15985d04e5e. MATLAB R2026a Update 5.
                Settings: timeInterval 1, bridgeWakeSeconds 1, Teliminate 24, deadThresh 12,
                minDays 3, phase 'entire'.
Input           The 25 per0 DD *_from1stLD.txt files; SHA-256 per file in
                analysis_20r/per0DD_input_sha256.txt.
Time window     Each retained fly's whole record after the pipeline's trimming, pooled.
n               16 flies retained, not 25. 2717 wake bouts (xmin 2 s), 2715 sleep bouts (xmin 300 s);
                61.45 fly-days covered by bouts (FINDINGS: ~90).
Result          GATE FAIL on every item except bin centres:
                (1) flies 16 / 25, wake 2717 / 3989, sleep 2715 / 3987: FAIL.
                (2) 16-bin table: 15 kept bins, centres match to printed precision, but at-risk and
                    event counts do not (first bin 2717 at risk, 26 events vs 3989, 46): FAIL.
                (3) null K = 2: tau 10.31 and 642.55 s, q 0.0948 / 0.9052, logL -19828.8228 vs
                    -29328.2460: FAIL (different data).
                (4) RiseRatio 1.225 / 1.044 / 1.371 / 1.495 vs 1.405 / 1.388 / 1.376 / 1.430: FAIL.
                Cause: the 9 discarded flies (2_18_2018 #54, #58; 6_25_2018_2boardmb #4, #8, #9, #14,
                #15, #27, #28) have files of exactly 3.00 days (259201 s). Teliminate = 24 h leaves
                48 h < minDays 3 x 24 h, so PIPELINE:InsufficientData drops them. The other 16 files
                are 4-7 days long; FINDINGS reports per-fly records of 2.8-4.0 days and ~90
                fly-days for 25 flies, which these settings cannot produce from these files. The
                reference set was therefore made with different settings or differently trimmed
                input; which, is not established.
Saved           analysis_20r/per0DD_bouts_20r.mat, SHA-256
                ab9c327be079804525f6c2c10ddd79624c666916dcc425d8914157a1c4d3ce2f.
                NOT the reference set: the failed extraction, kept so the failure is checkable.
Matched null    Not applicable (reproduction check).
Verdict         Gate FAILS. Per the amendment: no CS DD number has been looked at, 20r-A/B/C not run.
                No alternative settings were tried: choosing settings until the counts match would be
                fitting the gate, and needs a new pre-specified entry.
Not reportable  Every number above except the cause. The 16-fly RiseRatios in particular (1.044 at 16
                bins, trough in the second-to-last bin) are a different, smaller dataset and say
                nothing about the per0 DD rise.
```

### 2026-10-31  20r gate, amendment 2 (pre-cut per0 DD files, Teliminate 0): FAILS, close

```
Commit          a3611dd00187abb9bbd89790817bb85799cb8de1 (Icode branch 20r-lclaude)
Script          analysis_20r/gate_20r.m (Teliminate 0); output analysis_20r/gate_20r_report_T0.txt
Prespec         2026-10-31 "20r gate amendment 2: the pre-cut per0 DD files, Teliminate 0"
Code            Rockefeller Matlab scripts, HEAD 001deca, current versions (hashes as in the first
                gate entry). MATLAB R2026a Update 5. timeInterval 1, bridgeWakeSeconds 1,
                Teliminate 0, deadThresh 12, minDays 3, phase 'entire'.
Input           per0DD activity files_1sec res.zip (SHA-256 778599db...), 25 cutdata files.
Time window     Each fly's whole cut record (16 x 4.000 d + 9 x 3.000 d), pooled.
n               25 flies. 4017 wake bouts (xmin 2 s), 4022 sleep bouts (xmin 300 s); 90.45
                fly-days covered by bouts.
Result          (1) flies 25: pass. Wake 4017 vs 3989 (+28), sleep 4022 vs 3987 (+35): FAIL.
                (2) 16-bin centres match; at-risk and event counts do not. Every bin holds the same
                    or slightly more bouts than FINDINGS (events +3, +1, +3, 0, +2, +2, 0, 0, 0, 0,
                    +2, +4, +9, +1, +1; at risk +28 in the first bin, +11 in the last): FAIL.
                (3) null K = 2: tau 8.27 and 663.45 s, q 0.0813 / 0.9188, logL -29523.4966 vs
                    -29328.2460: FAIL (a different n, so logL cannot match).
                (4) RiseRatio 1.416 / 1.402 / 1.386 / 1.435 vs 1.405 / 1.388 / 1.376 / 1.430: FAIL
                    at 3 decimals. Trough and peak bins identical to the 2026-09-15 entry at all four
                    counts (119.8/1094.0, 109.3/997.7, 161.0/944.3, 144.1/910.0 s).
                Reading: the right flies and the right records, with 28 extra wake and 35 extra sleep
                bouts spread over the whole duration range. 4017 wake equals the count cClaude
                reports for the Python S/Q/W definition on these files, so the current MATLAB
                pipeline and the Python one agree on wake bouts here; the reference 3989/3987 came
                from something else. Candidate, NOT tested: the pre-09-22 pipeline (04e0508^:
                ironout b54ade9, segment_sleep_phases f6c81f5, process_LocomToSleep 31d341d), the
                code that existed when the reference set was made on or before 2026-09-15.
Saved           analysis_20r/per0DD_bouts_20r_T0.mat, SHA-256 b3655e6bdaed7c18b6b6768e2aba1260bb8a964ba60e0d9ac6398cb473ad5eee. NOT the reference set.
Matched null    Not applicable (reproduction check).
Verdict         Gate FAILS. No CS DD number looked at; 20r-A/B/C not run. No other code version or
                setting tried.
Not reportable  The RiseRatios here are close to the recorded ones but are not a replication: the bout
                set differs by 63 bouts, and the gate is a reproduction check, not a test.
```

### 2026-10-31  20r gate amendment 3 (pre-09-22 code, C1 and C2): both FAIL; stopping rule reached

```
Commit          f5a50a49b03fcdf924a3cdb2c8e31381bf9da67b (Icode branch 20r-lclaude)
Script          analysis_20r/gate_20r.m; outputs analysis_20r/gate_20r_report_C1.txt, _C2.txt
Prespec         2026-10-31 "20r gate amendment 3: pre-09-22 pipeline, two pre-listed configurations,
                stopping rule"
Code            Rockefeller repo at 04e0508^ = ef1bce38626b71216c795d436cc6ad89be2f2af3, exported by
                git archive to a scratch folder, added with addpath -begin. SHA-256 prefixes of the
                exported files (LF line endings): process_LocomToSleep 9d9d64ab39da231b, ironout
                f1db073793d96bf6, segment_sleep_phases 5cd6b76ec3ab2840, sleeptrace fd84c3c4bcd0f001,
                extractBoutsForModelCompetition 64f6f15985d04e5e. requiredFilesAndProducts listed
                exactly these five; which -all for each (both runs): first entry the export,
                second the current Rockefeller copy, e.g.
                  ironout: <scratch>\rock_pre0922\ironout.m
                           C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts\ironout.m
                After each run which('ironout') resolved to the current Rockefeller copy again
                (path restored).
Input           The 25 per0 DD cutdata files (zip SHA-256 778599db...). Teliminate 0, timeInterval 1,
                bridgeWakeSeconds 1, deadThresh 12; C1 minDays 2, C2 minDays 2.8.
Time window     Each fly's whole cut record, pooled.
n               Both: 25 flies, 4017 wake / 4022 sleep bouts.
Result          C1: no match. 4017 / 4022 vs 3989 / 3987; 16-bin centres match, at-risk and events do
                not (same table as amendment 2). Items (3)-(4) not computed (amendment 3: only for a
                matching configuration).
                C2: no match, same counts.
                The C1 bouts are bit-identical to amendment 2's (current code), and C2's to C1's
                (isequaln on the full wake and sleep matrices). So the 2026-09-22 recoding does not
                change the bouts on these data; the "bit-identical" claim in 04e0508's message,
                previously unchecked, holds here.
Diagnostic      Per fly: wake 81-249 (recorded 81-248), sleep 81-249 (recorded 80-249). First bout of
                the record kept in 0 of 25 flies, last bout kept in 0 of 25, in both configurations:
                segment_sleep_phases drops every bout touching the first or last sample as censored
                (20 sleep and 14-ish wake per run, reported by the extractor as "censored by
                segmenter"). Per-fly counts in gate_20r_report_C1.txt.
Saved           per0DD_bouts_20r_C1.mat 922749f1ad8bf9dfc7fc98dcc09dfb02c666237058f727536e55577e0d75c754
                per0DD_bouts_20r_C2.mat 62f652a47a4f935745728989dedf93b909ac2cf85f10de9004dd1b40c832cb5f
                (identical contents to per0DD_bouts_20r_T0.mat; the files differ only in metadata).
Matched null    Not applicable (reproduction check).
Verdict         No configuration reproduces the reference set. STOPPING RULE reached: no further
                settings are tried. 20r moves to the replication route, which must be written as its
                own PRESPEC entry before anything runs. No CS DD number has been looked at.
Not reportable  Nothing here is a hazard result. The stored 3989 / 3987 set is not reproducible from
                the code and data available; its provenance is unresolved.
```

### 2026-10-31  Correction to "20r gate amendment 3 (pre-09-22 code, C1 and C2)"

```
Commit          0d007e8 (no new computation; read from per0DD_bouts_20r_C1.mat, meta.NumCensoredByPipeline)
Corrects        The Diagnostic line, which gave the segmenter's dropped wake bouts as "14-ish".
Correct         30 wake and 20 sleep bouts dropped as boundary-censored, over the 25 flies (C1; C2
                identical). Whole-trace counts, as the extractor reports them.
Unchanged       Everything else in that entry.
```

### 2026-10-31  20r: MATLAB scoring chain against the G1 golden bridging cases

```
Commit          ee138cfe59ebf0123496eec9045d4a1a653c6bec (Icode branch 20r-lclaude)
Script          analysis_20r/validate_G1_20r.m; output analysis_20r/G1_validation.txt
Prespec         None of its own: a code check required by the replication route (Sheyum/cClaude,
                2026-10-31) before that route's PRESPEC entry names the extractor.
Code            Current Rockefeller pipeline (HEAD 001deca): sleeptrace -> ironout(.,1,1) ->
                segment_sleep_phases(timeInterval=1), as process_LocomToSleep calls them.
Golden          hmm_audit/golden/G1_sleep_bouts.json (SHA-256 prefix 49ecca22e49f6643), from sqw_ref.sleep_bouts
                (Python); 0-based starts shifted +1.
Result          8 cases x {bridged, unbridged} = 16 comparisons. 14 AGREE, 2 DIFFER.
                Every bridging case agrees: 301-1-300 bridged to 602 s; 301-1-100, 300-2-300 and
                299-1-400 not bridged; the chain 300-1-300-1-300 bridged to 902 s; exactly 300 s is a
                bout; 299 s is not.
                The 2 differences are one case, at_start_no_left_flank (400 s inactive from the first
                sample, 1 active, 400 inactive, 5 active):
                  Python bridged    [1 801]          MATLAB kept []
                  Python unbridged  [1 400] [402 400]  MATLAB kept [402 400]
                Before segment_sleep_phases' censoring, MATLAB's scored wave holds exactly Python's
                bouts ([1 801] and [1 400] [402 400]). So the bridging agrees there too; the
                difference is EDGE HANDLING ONLY: MATLAB drops a bout touching the first or last
                sample as censored, Python keeps it.
                Consistent with the per0 DD counts: Python 4042 sleep vs MATLAB 4022 is 20, exactly
                the 20 sleep bouts MATLAB drops as censored there (C1 correction entry); wake 4017 in
                both.
Matched null    Not applicable (code check).
Verdict         Neither side changed. The lab convention for bouts at the record edges is Sheyum's to
                designate; the 20r extractor will be the designated one, stated in the replication
                route's PRESPEC entry.
Not reportable  Not a test of either convention's correctness. Neither convention reproduces the lost
                3989 / 3987 set.
```

### 2026-10-31  20r replication: per0 DD wake hazard rise at the four fixed bin counts

```
Commit          d9334dd63768bb576308f5aa611ab1583215a29f (Icode branch 20r-lclaude)
Script          analysis_20r/extract_20r.m, analysis_20r/rise_test_20r.m; log
                analysis_20r/replication_per0DD.txt
Prespec         2026-10-31 "20r amendment: stored per0 DD bouts lost (replication route)", step (2):
                the 2026-09-15 "per0 DD wake: the hazard rise at four fixed bin counts", unchanged.
Extractor       Current Rockefeller pipeline (HEAD 001deca), edge bouts dropped (Sheyum's convention),
                Teliminate 0, on the 25 cutdata files (zip 778599db...).
Bouts saved     analysis_20r/per0DD_bouts_20r_replication.mat, SHA-256
                b14ba34047b04d028154da4ae492959fca1528eeac2c3f0714b84e71b0591b30
                (wake and sleep per fly, NaN-padded, seconds, with file list and trace lengths).
                This is now the per0 DD reference set for 20r.
Time window     Each fly's whole cut record (16 x 4.000 d + 9 x 3.000 d), 25 flies pooled. The pooled
                set is the narrowest block.
n               4017 wake bouts (xmin 2 s, dt 1 s); 30 wake and 20 sleep bouts dropped at record edges.
                The 2026-09-15 set had 3989.
Result          Null K = 2: tau 8.27 and 663.45 s, observed weights 0.081 / 0.919, logL -29523.4966
                (2026-09-15: 8.4 and 661.9 s, 0.079 / 0.921).
                bins 12: trough 119.8 s, peak 1094.0 s, RiseRatio 1.416, null median 1.089, 95th 1.208,
                         RiseNullP 0.001   (2026-09-15: 1.405, 1.083, 1.210, 0.001)
                bins 16: trough 109.3 s, peak  997.7 s, RiseRatio 1.402, null median 1.136, 95th 1.283,
                         RiseNullP 0.009   (2026-09-15: 1.388, 1.136, 1.285, 0.010)
                bins 20: trough 161.0 s, peak  944.3 s, RiseRatio 1.385, null median 1.185, 95th 1.401,
                         RiseNullP 0.059   (2026-09-15: 1.376, 1.186, 1.390, 0.068)
                bins 24: trough 144.1 s, peak  910.0 s, RiseRatio 1.435, null median 1.231, 95th 1.503,
                         RiseNullP 0.086   (2026-09-15: 1.430, 1.235, 1.512, 0.101)
                NullDropped 0 at every count; one bin dropped and one saturated at each.
                Trough and peak bins identical to 2026-09-15 at all four counts.
Matched null    999 datasets of n = 4017 from the fitted two-exponential model, each through the whole
                estimator on its own bin edges.
Verdict         Supported but bin-sensitive (rise at 12 and 16 bins, not at 20 and 24): the 2026-09-15
                verdict REPLICATES on the re-extracted bouts, count by count.
Not reportable  Do not quote 0.001 alone; every count is the result. Not an independent replication: the
                same flies and records, re-scored, differing from the lost set by 28 wake bouts. Says
                nothing about whether the rise is in the flies or made by scoring (20r-B).
```

### 2026-10-31  20r-A: CS DD wake hazard rise (blind), influence, sleep control, model library

```
Commit          0197c70 for run_20rA.m / influence_20r.m (code unchanged through 7eae6a4); extractor
                and rise test as in the per0 DD replication entry. Log analysis_20r/run_20rA_log.txt,
                results analysis_20r/run_20rA_results.mat.
Prespec         2026-10-31 20r-A, with the extractor of "20r amendment: stored per0 DD bouts lost",
                step (3). Run after the per0 DD replication, as that entry requires; blind until then.
Data            CS DD, 25 flies: CS_DD_activfiles.zip (SHA-256 33495ce4...), each file cut to its first
                345600 samples (hashes of raw and cut files in analysis_20r/CSDD_input_sha256.txt).
                25 of 25 flies retained, 100.00 fly-days of record. Bouts saved to
                analysis_20r/CSDD_bouts_20r.mat, SHA-256
                80402183169d4c5a92dbda2a06b2ba3c7f78f536bd6480fa857ecc33aed2dfa5.
Time window     Each fly's whole cut record (4.000 d), 25 flies pooled; the pooled set is the
                narrowest block.
n               3309 wake bouts (xmin 2 s), 3314 sleep bouts (xmin 300 s), dt 1 s; 30 wake and 20
                sleep bouts dropped at record edges.

(1) Rise test, wake. Null K = 2: tau 9.83 and 1356.0 s, q 0.089 / 0.911, logL -26525.3468.
                bins 12: trough 3255.9 s, peak  7628.8 s, RiseRatio 0.327, null median 1.099, 95th
                         1.247, RiseNullP 1.000
                bins 16: trough 5543.1 s, peak 10497.3 s, RiseRatio 0.482, null median 1.159, 95th
                         1.342, RiseNullP 1.000
                bins 20: trough 7628.4 s, peak 12714.0 s, RiseRatio 0.862, null median 1.218, 95th
                         1.448, RiseNullP 0.994
                bins 24: trough 9437.8 s, peak 14446.6 s, RiseRatio 1.453, null median 1.276, 95th
                         1.575, RiseNullP 0.140
                NullDropped 0 at every count; one bin dropped and one saturated at each.
                Verdict: NOT SUPPORTED (RiseNullP >= 0.05 at all four counts).
                The trough lies in the far tail (3256-9438 s) at every count: the CS DD wake hazard
                falls through the 100-1000 s range where per0 DD's rises (per0 troughs 109-161 s,
                peaks 910-1094 s). RiseRatio < 1 at 12, 16 and 20 bins: the only kept bins after the
                trough are lower than it. BoutHazard's header says RiseRatio is ">= 1 for ANY curve";
                that does not hold when the trough is the second-to-last kept bin. The verdict does
                not depend on it (p ~ 1), but the header claim is wrong; flagged, not changed.
(2) Per-fly influence (2026-09-15 design). No fly influential: fly 15 below its 2.5th percentile at
                20 bins only, fly 21 at 24 bins only, every other fly at none. Verdict by the rule:
                no single fly carries the (absent) rise. With no rise to carry, this is vacuous.
(3) Sleep control. Null K = 2: tau 154.4 and 1199.0 s, q 0.146 / 0.854. RiseRatio 0.731 / 1.091 /
                0.577 / 1.927, null median 1.200 / 1.326 / 1.430 / 1.568, RiseNullP 0.982 / 0.893 /
                0.997 / 0.210, NullDropped 0. Verdict: not supported, as in per0 DD sleep.
(4) CompareBoutModels(x, 2, SamplingInterval=1, MaxComponents=5), default library.
                R.Skipped: empty. Table (AICc): hyperexp K=4 (k 7, logL -26283.56, w 0.794); K=5
                (dAICc 3.83); hyper_erlang and hyperexp K=3 tied (dAICc 5.76, logL -26288.45 both);
                weibull_mix 35.1; exp_weibull 265; weibull 267; gamma = powerlaw_cutoff 301 (logL
                identical, nesting check holds); K=2 476; K=1 = erlang (logL identical) 1158;
                powerlaw 6141; chisquared not interpretable (not unit-invariant). No row Degenerate.
                R.Selected / SelectionPath: hyperexp K=4, the first row tested, G = 32.80 on 40 bins,
                p in [0.428, 0.748]: accepted. Rows below it are untested, not accepted.
                Guards: Degenerate 0, smallest branch q x n = 114 bouts, AtRateBound false,
                WeightsUntruncatedReliable true. q = 0.061 / 0.035 / 0.842 / 0.062 at tau 5.7 /
                45.9 / 1017 / 6093 s.
                hyper_erlang: the fitted shapes are [1 1 1], so it IS hyperexp K=3 (logL equal to 4
                decimals). ShapeSweep logL by series-branch shape m = 1..6: -26288.448, -26288.751,
                -26289.682, -26291.002, -26292.497, -26294.060: falls from m = 1. No series branch is
                supported in CS DD wake bouts (per0 DD: +15.8 nats from m = 1 to 2, FINDINGS).
Matched null    (1) and (3): 999 datasets of the observed n from that set's fitted two-exponential
                model, each through the whole estimator on its own bin edges. (2): 199 random removals
                of n_i bouts per fly (RandomSeed 1000 + fly).
Verdict         20r-A: the per0 DD wake-hazard rise is NOT seen in CS DD at any of the four bin counts,
                and the hyper-Erlang series branch collapses to an exponential. Per 20r-B's rule, the
                CS DD renewal control is NOT run (20r-A did not give "supported" or "supported but
                bin-sensitive"). Per 20r-C's fixed choices, CS DD 20r-C is NOT run (hyper_erlang not
                selected, and identified only with no series branch).
Not reportable  Not "CS DD has four timescales": AICc picks K=4 but BIC prefers K=3 (dBIC 6.43), the
                order LRT was not triggered, and a pooled mixture cannot attribute components to
                flies (2026-09-15 correction). Not "the rise needs a broken clock": this is one CS
                cohort against one per0 cohort, different flies and runs, and the per0 rise is itself
                bin-sensitive. The influence result is vacuous here and must not be cited as
                robustness of anything.
```

### 2026-10-31  20r-C (EXPLORATORY): per0 DD series-branch share per fly vs 20p excess irreversibility

```
Commit          61228262d5b1c4fcd984c5a9c86a289dd4fa9538 (analysis_20r/run_20rC.m; log
                analysis_20r/run_20rC_per0DD_log.txt; results analysis_20r/run_20rC_per0DD.mat)
Prespec         2026-10-31 20r-C, with step (4) of "20r amendment: stored per0 DD bouts lost". CS DD
                NOT run: 20r-A's hyper_erlang collapsed to shapes [1 1 1] (no series branch).
Time window     Each fly's whole cut record; the fit pools the 25 flies.
n               25 flies (all matched in hmm_audit/e_perfly_20p.csv by file name, all with a Q/W
                reference); 4017 wake bouts in the pooled fit.
Result          Pooled fit, FitHyperErlangMLE Shapes [1 1 5], discrete, RandomSeed 1: logL -29503.5874,
                GuardOK 1; rates 1.514e-3, 0.1064, 6.88e-3 /s; observed q 0.803 / 0.081 / 0.116
                (FINDINGS on the lost set: 1.515e-3, 0.1053, 6.949e-3; q 0.807 / 0.079 / 0.114).
                Against the K = 2 two-exponential on the same bouts (replication entry, -29523.4966):
                +19.9 nats for two more parameters.
                Per-fly mean responsibility of the 5-stage branch: median 0.111, range 0.056-0.168.
                Spearman with E_nats_per_h:          rho 0.338, p 0.099, n 25
                Spearman with E_wakeonly_nats_per_h: rho 0.139, p 0.505, n 25
Matched null    Not applicable (exploratory).
Verdict         None, by design. Exploratory and hypothesis-generating only.
Not reportable  Neither correlation is evidence of anything: post hoc, n = 25, noisy per-fly E, two
                tests run, neither below 0.05. Must not be cited unless replicated by a pre-specified
                test on new flies. The responsibilities come from ONE pooled fit applied to every fly,
                so they measure where a fly's bouts fall on the pooled mixture, not a per-fly fit.
                The +19.9 nats is a comparison with K = 2, not the matched-k comparison of FINDINGS.
```

### 2026-10-31  20r-A follow-up (b): CS DD wake hazard rise within the two 12-h record-clock halves

```
Commit          043e13b7deff16001c0fc67bce88ba2e9f00f8fc (analysis_20r/run_20rA_followups.m 'halves',
                rise_test_20r.m with NullReplicates 999); log analysis_20r/followup_halves_log.txt,
                results analysis_20r/followup_halves_CSDD.mat
Prespec         2026-10-31 "20r-A follow-ups" (b), and its amendment (power at the half sizes),
                committed in 9389016 before this output was read.
Time window     Wake bouts STARTING in record hours 0-12 (resp. 12-24) of each of the 4 days, 25 flies
                pooled; each half is the narrowest block. Not clipped at the half boundary.
n               1536 (hours 0-12) and 1773 (hours 12-24); xmin 2 s, dt 1 s.
Result          Hours 0-12. Null K = 2: tau 9.19 and 1949.7 s, q 0.073 / 0.927.
                bins 12: trough 3255.9 s, RiseRatio 0.326, null median 1.176, 95th 1.474, p 1.000
                bins 16: trough 5543.1 s, RiseRatio 0.487, null median 1.270, 95th 1.703, p 1.000
                bins 20: trough 7628.4 s, RiseRatio 0.859, null median 1.386, 95th 1.907, p 0.998
                bins 24: trough 9437.8 s, RiseRatio 1.613, null median 1.501, 95th 2.189, p 0.345
                Hours 12-24. Null K = 2: tau 7.87 and 811.7 s, q 0.089 / 0.911.
                bins 12: trough 1845.4 s, RiseRatio 0.671, null median 1.142, 95th 1.357, p 1.000
                bins 16: trough 3026.6 s, RiseRatio 0.674, null median 1.226, 95th 1.506, p 1.000
                bins 20: trough 4073.1 s, RiseRatio 0.448, null median 1.304, 95th 1.664, p 1.000
                bins 24: trough 3341.9 s, RiseRatio 0.797, null median 1.383, 95th 1.822, p 0.993
                NullDropped 0 throughout. Troughs again in the far tail (1845-9438 s), never in the
                100-1000 s range where per0 DD rises. The slow wake time constant differs between
                halves (1950 vs 812 s): the halves are not the same law, as the free-running rhythm
                predicts.
Matched null    Per half, 999 datasets of that half's n from the half's own two-exponential fit, each
                through the whole estimator on its own bin edges.
Verdict         Not supported in either half (p >= 0.345 at every count). By the pre-specified reading:
                no rise within a half either, so masking by the rhythm does not explain the pooled CS
                DD null. Two halves tested, no correction; correction could only weaken a rise, and
                there is none. How much this negative is worth depends on the per0 power at n = 1536
                and 1773 (amendment), still running; to be read together with that entry.
Not reportable  Not "CS DD has no rise" on its own: until the half-size power is in, low power is not
                excluded. The censoring gap applies (bouts assigned by start, not clipped at the
                boundary; long bouts carry other-half time, direction of the bias not established).
                The halves approximate subjective day/night only to ~1-2 h; tau not estimated.
```

### 2026-10-31  20r-B: per0 DD wake hazard rise against a burst-pause renewal null

```
Commit          7eae6a4fbe63997dfb6af88699012923b4b837e8 (analysis_20r/run_20rB.m; extractor
                extract_20r.m); log analysis_20r/run_20rB_per0DD_log.txt, results
                analysis_20r/run_20rB_per0DD.mat
Prespec         2026-10-31 20r-B, with the extractor of "20r amendment: stored per0 DD bouts lost".
Data            per0 DD reference set (per0DD_bouts_20r_replication.mat, b14ba340...) and its 25 cut
                count files.
Time window     Each fly's whole cut record for the null construction; pooled over the 25 flies for
                the statistic. The pooled set is the narrowest block.
n               Observed 4017 wake bouts (xmin 2 s, dt 1 s). Synthetic pooled wake bouts per
                replicate: median 4616 (4499-4763), so all 199 were subsampled to 4017.
Pre-run check   Fly 1 (12_29_2017 #10), one synthetic trace (rng 2000): active-run quantiles (10/25/
                50/75/90/99%) 1 1 1 1 2 5 s real and synthetic; inactive 1 1 2 5 20 1066 s real vs
                1 1 2 5 22 1135 s synthetic; 10381 vs 9802 runs. Passed (read before the run).
Result          bins 12: observed RiseRatio 1.416; renewal median 1.027, 95th 1.145; p_ren 0.015
                bins 16: observed RiseRatio 1.402; renewal median 1.047, 95th 1.197; p_ren 0.015
                bins 20: observed RiseRatio 1.385; renewal median 1.088, 95th 1.358; p_ren 0.050
                bins 24: observed RiseRatio 1.435; renewal median 1.114, 95th 1.418; p_ren 0.050
                (p_ren = (1 + #{RiseRatio* >= observed}) / 200: 2, 2, 9 and 9 of 199 replicates
                reached the observed ratio.)
                Descriptive, as pre-specified: the renewal null's median RiseRatio (1.027, 1.047,
                1.088, 1.114) is BELOW the two-exponential null medians of 2026-09-15 (1.083, 1.136,
                1.186, 1.235) and of the replication (1.089, 1.136, 1.185, 1.231) at every count.
                Renewal plus scoring does not itself manufacture a rise; it makes less of one than a
                falling-hazard mixture. Renewal also yields ~15% MORE wake bouts than observed from
                the same run-length marginals.
Matched null    199 pooled datasets, each the 25 flies' own alternating renewals of their active and
                inactive runs (same record lengths, same run-length marginals, no correlation between
                successive runs), scored by the identical pipeline, subsampled to n = 4017;
                RandomSeed 2000 + b.
Verdict         INCONCLUSIVE (bin-sensitive), per the rule: p_ren < 0.05 at 2 counts (12, 16), >= 0.05
                at 2 (20, 24). The two counts that fail sit exactly at 0.050 (9 of 199), the boundary
                of the rule. CS DD 20r-B not run (20r-A not supported).
Not reportable  Not "the rise is beyond renewal" and not "renewal explains it": the pre-specified rule
                gives neither. Do not quote the 12/16-bin p_ren alone. The 0.050s are at the
                resolution limit of B = 199 (steps of 0.005); a larger B could move them either way,
                and choosing B after seeing this would be rule 4's failure. The ~15% bout-count excess
                means subsampling discards synthetic bouts, so the null's bin occupancy differs from
                the data's; its effect on the ratio is not measured.
```

### 2026-10-31  20r-A follow-up (a), amendment: power to see the per0 DD rise at n = 1536

```
Commit          7bd0e99f92554a6a71eeda4604c0eb0236fcef49 (analysis_20r/run_20rA_followups.m 'power',
                rise_test_20r.m); log analysis_20r/followup_power_n1536_log.txt (a snapshot of the
                running half-size log at the end of this set), results
                analysis_20r/followup_power_per0DD_n1536.mat
Prespec         "Amendment to 20r-A follow-ups: power also at each half's n" (9389016).
Data            per0 DD reference wake bouts (per0DD_bouts_20r_replication.mat, b14ba340...).
Time window     Each subsample, pooled over the flies it holds.
n               1536 per subsample (the CS DD record-hours 0-12 count); 199 subsamples, median 9 whole
                flies each plus part of one more; rng(3000 + b).
Result          Power (fraction with RiseNullP < 0.05), bins 12 / 16 / 20 / 24:
                  0.774 / 0.598 / 0.432 / 0.377
                Median RiseRatio 1.537 / 1.601 / 1.664 / 1.792; median RiseNullP 0.010 / 0.025 /
                0.065 / 0.080. Verdict over the four counts: "supported" 0.302; "supported" or
                "supported but bin-sensitive" 0.809. NullDropped > 0 in none.
Matched null    Each subsample's own two-exponential fit, 199 null datasets of n = 1536, each through
                the whole estimator on its own bins (NullReplicates 199, Sheyum 2026-10-31).
Verdict         Descriptive (no rule). At n = 1536 a per0-like rise is detected at 12 bins in 77% of
                subsamples and gives a "supported" or "bin-sensitive" verdict in 81%. CS DD hours 0-12
                (n = 1536) gave p = 1.000 / 1.000 / 0.998 / 0.345: not a near miss.
Not reportable  The power is for THESE 25 per0 flies, resampled with overlap, not for a new cohort, and
                assumes CS DD's rise, if any, would be per0-sized. 199 null replicates per subsample:
                p in steps of 0.005.
```
