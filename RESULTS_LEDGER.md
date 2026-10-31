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
