# Rockefeller bout extractor, as used in the 20r runs

Exported 2026-11-01 for PRESPEC 20u test 1 (Sheyum approved sharing). The five files are **unchanged, byte-for-byte**
copies of `C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts` at that folder's git HEAD **001deca**
(current versions, the ones the 20r replication route and 20r-A / 20r-B ran). `.gitattributes` here stores them
without line-ending conversion (they are CRLF files), so a checkout on any platform hashes the same.

| file | SHA-256 | 20r gate entry prefix |
|---|---|---|
| `process_LocomToSleep.m` | `7aec365d5b79c859943ec4590f28eb9bc25e788d5a4ed211f661f1e8ac5570a5` | 7aec365d5b79c859 |
| `sleeptrace.m` | `0deb702eed5cbfcc35c03ec6caf918b4e903435f749e01e032600c8d4d5d5371` | 0deb702eed5cbfcc |
| `ironout.m` | `f22f91c262f8c37015067c7daed4b2e5336ac70c66122f06063d10241fcdd590` | f22f91c262f8c370 |
| `segment_sleep_phases.m` | `a627f698159bc97ad51d9901836ea57e1aa9916b08dcca8b6d596a2c76ff92a5` | a627f698159bc97a |
| `extractBoutsForModelCompetition.m` | `64f6f15985d04e5e5f4850c20dfe2f9e54bc9508d52e7a495aa1fb0adb4141dc` | 64f6f15985d04e5e |

**Helpers:** none. `matlab.codetools.requiredFilesAndProducts` on the five files (R2026a Update 5) lists only these
five files and the product MATLAB (no toolbox).

**Settings used in 20r** (the call site is `analysis_20r/extract_20r.m`; PRESPEC "20r amendment: stored per0 DD bouts
lost (replication route)"):

| setting | value |
|---|---|
| timeInterval | 1 (s) |
| bridgeWakeSeconds | 1 (`ironout` bridges wake runs of 1 s, so wake xmin = 2 s) |
| Teliminate | 0 (h; the cut files already start at subjective day, first 24 h of DD removed before export) |
| deadThresh | 12 |
| minDays | 3 |
| phase | 'entire' |

Edge-touching bouts are dropped as censored (Sheyum's convention, `segment_sleep_phases`); sleep xmin = 300 s.
Reference output: `analysis_20r/per0DD_bouts_20r_replication.mat` (SHA-256 b14ba340..., 25 flies, 4017 wake bouts
>= 2 s, 4022 sleep bouts >= 300 s); input hashes per file in `analysis_20r/per0DD_input_sha256.txt`.
