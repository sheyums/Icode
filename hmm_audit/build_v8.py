# V8: V7 with all tracked changes accepted, sections reordered into a logical structure, group headings added,
# equations renumbered in the new order (four minor ones left unnumbered), in-text references updated.
import docx, re, copy, sys
SRC = '/tmp/claude-0/acc3.docx'; OUT = 'MaterialsMethod_V8_reorganized.docx'
D = docx.Document(SRC); body = D.element.body
W = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'; M = '{http://schemas.openxmlformats.org/officeDocument/2006/math}'
def ptext(el): return ''.join(t.text or '' for t in el.iter(W+'t'))
def is_head(el):
    if el.tag != W+'p': return False
    runs = [r for r in el.iter(W+'r') if ''.join(t.text or '' for t in r.iter(W+'t')).strip()]
    if not runs: return False
    bi = lambda r: (r.find(W+'rPr') is not None and r.find(W+'rPr').find(W+'b') is not None and r.find(W+'rPr').find(W+'i') is not None)
    return all(bi(r) for r in runs) and len(ptext(el)) < 130
els = [e for e in body if e.tag != W+'sectPr']; sect = body.find(W+'sectPr')
# ---- 1. equation renumbering (document order of V7 -> new numbers); None = unnumbered
MAP = {20:1, 21:2, 22:3, 23:4, 24:5, 25:6, 5:7, 6:8, 7:9, 8:10, 9:11, 10:None, 11:12, 12:13, 13:None,
       14:14, 15:15, 16:16, 17:17, 18:18, 19:19, 1:20, 2:21, 3:None, 4:None}
eq_done = set()
for el in els:
    t = ptext(el).strip()
    is_eq = (el.tag == W+'tbl') or any(True for _ in el.iter(M+'oMath')) or len(t) < 12
    if not is_eq: continue
    for node in el.iter(W+'t'):
        def rep(m):
            n = int(m.group(1)); eq_done.add(n); new = MAP.get(n, n)
            return '' if new is None else f'({new})'
        if node.text and re.fullmatch(r'\s*(\(\d+\)\s*)+', node.text):
            node.text = re.sub(r'\((\d+)\)', rep, node.text)
# references: explicit forms found in V7 (checked by hand), replaced within single text nodes
REF = [('Eqs. 5–25', 'Eqs. 1–19'), ('Eqs. 5–19', 'Eqs. 7–19'), ('Eqs. (7)-(9)', 'Eqs. (9)–(11)'), ('Eqs. 11, 6', 'Eqs. 12, 8'),
       ('Eqs. 14–15', 'Eqs. 14–15'), ('Eq. 14, 15', 'Eqs. 14, 15'), ('Eqs. 14, 15', 'Eqs. 14, 15'), ('Eqs. 18, 19', 'Eqs. 18, 19'),
       ('Eqs. 20, 22', 'Eqs. 1, 3'), ('Eqs. 5 ', 'Eqs. 7 ')]
SINGLE = {20:1, 21:2, 22:3, 23:4, 24:5, 25:6, 5:7, 6:8, 7:9, 9:11, 11:12, 12:13, 16:16, 14:14}
nref = 0; missed = []
for el in els:
    for node in el.iter(W+'t'):
        if not node.text: continue
        s = node.text
        for a, b in REF:
            if a in s: s = s.replace(a, b); nref += 1
        def rs(m):
            global nref
            n = int(m.group(2)); nref += 1
            return f"{m.group(1)}{SINGLE[n]}" if n in SINGLE else m.group(0)
        s = re.sub(r'(Eq\. \(?)(\d+)(?=[\)\.,;: ])', rs, s) if 'Eqs.' not in s else re.sub(r'(Eq\. \(?)(\d+)(?=[\)\.,;: ])', rs, s)
        node.text = s
# ---- 2. sections
secs = {}; order0 = []; cur = None
for el in els:
    if is_head(el): cur = ptext(el); secs[cur] = [el]; order0.append(cur)
    else: secs.setdefault(cur, []).append(el)
GROUPS = [
 ('Animals, recordings and genetic manipulations', ['Fly husbandry and genetic backcrossing', 'Neurogenetics', 'Activity and sleep measurements', 'Sleep deprivator']),
 ('Behavioural states and their models', ['Hidden Markov model of sleep-wake dynamics', 'Observed S/Q/W chain', 'Hidden Markov model: robustness check and model-based estimate', 'Validation of the inference']),
 ('Direction of circulation', ['Model-free cycle current', 'Round trips around sleep bouts', 'Changes in the direction of circulation']),
 ('Thermodynamic quantities', ['Thermodynamic quantities: what is computed, and from what', 'Net cycle flux', 'Cycle affinity', 'Entropy production rate', 'Onsager linear response and thermodynamic conductance', 'Energy landscape quantification', 'Routing fraction and routing deviation']),
 ('Beyond the observed chain', ['Kullback–Leibler divergence and time irreversibility', 'Observable and hidden levels compared']),
 ('Circadian modulation', ['Circadian modulation of the cycle']),
 ('Memory and bout durations', ['Dependence of routing on bout duration', 'Fitting probability distributions to sleep and wake bouts', 'Wake-bout hazard and its controls']),
 ('Perturbations', ['Sleep deprivation: within-fly contrasts with a day-to-day control', 'Silencing of dorsal fan-shaped body and ventral nerve cord neurons']),
 ('Statistics and software', ['Group comparisons of derived quantities', 'Cliff’s delta effect size', 'Software']),
 ('Sections from the earlier version not used in this study (to be deleted)', ['Overlap between hidden states and observable sleep/wake', 'Identification of the circadian-encoding hidden state via synthetic lesion', 'Multi-exponential fitting using log-transformed histogram method', 'Detection of routing subgroups within per0'])]
used = [h for _, hs in GROUPS for h in hs]
missing = [h for h in order0 if h not in used]; extra = [h for h in used if h not in secs]
assert not missing and not extra, (missing, extra)
RENAME = {'Hidden Markov model of sleep-wake dynamics': 'Sleep, quiet wake and active wake'}
template = secs['Fly husbandry and genetic backcrossing'][0]
def group_head(text):
    p = copy.deepcopy(template)
    rs = [r for r in p.iter(W+'r')]
    for r in rs[1:]: r.getparent().remove(r)
    r = rs[0]; rp = r.find(W+'rPr')
    if rp is not None and rp.find(W+'i') is not None: rp.remove(rp.find(W+'i'))
    if rp is not None and rp.find(W+'iCs') is not None: rp.remove(rp.find(W+'iCs'))
    if rp is not None and rp.find(W+'caps') is None: rp.append(rp.makeelement(W+'caps', {}))
    ts = list(r.iter(W+'t')); ts[0].text = text
    for t in ts[1:]: t.text = ''
    return p
for el in els: body.remove(el)
DROP = 'Sections from the earlier version not used in this study (to be deleted)'   # deleted at Sheyum's request, 2026-11-02
dropped_ids = set()
for g, hs in GROUPS:
    if g == DROP:
        for h in hs:
            for el in secs[h]:
                for cr in el.iter(W+'commentReference'): dropped_ids.add(cr.get(W+'id'))
GROUPS = [(g, hs) for g, hs in GROUPS if g != DROP]
for gi, (g, hs) in enumerate(GROUPS, 1):
    body.insert(len(body) - 1, group_head(f'{gi}. {g}'))
    for h in hs:
        for el in secs[h]:
            if el is secs[h][0] and h in RENAME:
                ts = list(el.iter(W+'t')); ts[0].text = RENAME[h]
                for t in ts[1:]: t.text = ''
            body.insert(len(body) - 1, el)
# remove the dropped sections' comments from the comments part
for rel in D.part.rels.values():
    if rel.reltype.endswith('/comments'):
        cp = rel.target_part; from lxml import etree; root = etree.fromstring(cp.blob)
        for c in list(root):
            if c.get(W+'id') in dropped_ids: root.remove(c)
        cp._blob = etree.tostring(root, xml_declaration=True, encoding='UTF-8', standalone=True)
D.save(OUT)
print('dropped comments:', sorted(dropped_ids))
print('equation numbers seen:', sorted(eq_done)); print('reference replacements:', nref); print('sections:', len(used))
