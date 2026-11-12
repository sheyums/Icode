"""Parse a raw multibeam week file into a compact .npz (int8 positions, float64 time).
Exploration loader (E0); format facts it relies on are asserted, not assumed."""
import sys, numpy as np
src, dst = sys.argv[1], sys.argv[2]
rows, times, rowtype = [], [], []
with open(src, 'r', newline='') as f:
    hdr = [next(f) for _ in range(3)]
    for line in f:
        p = line.rstrip('\r\n').split('\t')
        assert len(p) == 65, (len(p), line[:80])
        t = p[64]
        rows.append(np.array(p[:64], dtype=np.int8))
        times.append(float(t) if t.strip() else np.nan)
X = np.vstack(rows); T = np.array(times)
np.savez_compressed(dst, X=X, T=T, header=np.array(hdr))
print(hdr); print(X.shape)
