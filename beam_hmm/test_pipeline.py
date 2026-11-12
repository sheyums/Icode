"""Unit checks for pipeline.py (synthetic, known answers)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pipeline as pl
# spike filter
x = np.array([5, 5, 12, 5, 5, 6, 5, 9, 9]); y, k = pl.spike_filter(x)
assert k == 1 and list(y) == [5, 5, 5, 5, 5, 6, 5, 9, 9], y
# bridging: 301 inactive + 1 active + 300 inactive -> one 602-s bout; 301+1+100 -> not bridged
a = np.r_[np.ones(10), np.zeros(301), 1, np.zeros(300), np.ones(10)].astype(np.int8)
b = pl.sleep_bouts(a); assert b == [(10, 602, False)], b
a = np.r_[np.ones(10), np.zeros(301), 1, np.zeros(100), np.ones(10)].astype(np.int8)
b = pl.sleep_bouts(a); assert b == [(10, 301, False)], b
# excluded seconds break and censor a bout
a = np.r_[np.ones(10), np.zeros(400), -np.ones(5), np.zeros(400), np.ones(10)].astype(np.int8)
b = pl.sleep_bouts(a); assert b == [(10, 400, True), (415, 400, True)], b
# activity_1s: samples at 0.4 s; change at t=1.2 -> second 1 active
t = np.arange(10) * 0.4; x = np.array([3, 3, 3, 4, 4, 4, 4, 4, 4, 4]); v = np.ones(10, bool)
assert list(pl.activity_1s(x, t, v)) == [0, 1, 0, 0], pl.activity_1s(x, t, v)
# edge flicker detection: A-B-A with |B-A| = 1 counts both changes
t = np.arange(300) * 0.4; x = np.full(300, 7); x[50] = 8; x[200] = 12; v = np.ones(300, bool)
D = pl.minute_descriptors(x, t, v)
assert D['flick'][0] == 2 and D['changes'][1] == 2 and D['flickonly_s'][0] == 1, D
print('pipeline unit checks: all passed')
