# Builds MaterialsMethod_V7_Claude_tracked.docx from the user's V6: tracked insertions/deletions (author "Claude")
# for the S/Q/W HMM analysis, plus comments on sections that rest on the binary per-window model.
# Rerunnable: always starts from mm_v6.docx.
import re, os, subprocess, shutil, json, sys
from xml.sax.saxutils import escape
HERE = os.path.dirname(os.path.abspath(__file__))
SK = '/root/.claude/skills/synced/093a3a8c-05d9-4e67-ab6f-64868e3645dc_9e28ec93-5c7e-480f-b37c-a97e2b130e27/docx/scripts'
AUTH = 'Claude'; DATE = '2026-10-24T08:00:00Z'
TXT = json.load(open(os.path.join(HERE, 'methods_text.json')))
U = os.path.join(HERE, 'v7'); shutil.rmtree(U, ignore_errors=True); os.makedirs(U)
subprocess.run(['unzip', '-q', os.path.join(HERE, 'mm_v6.docx'), '-d', U], check=True)
subprocess.run(['find', U, '-type', 'l', '-delete'], check=True)
subprocess.run(['python3', os.path.join(SK, 'merge_runs.py'), U], check=True, capture_output=True)
DOC = os.path.join(U, 'word', 'document.xml')
rid = [6000]
def nid():
    rid[0] += 1; return str(rid[0])
def attrs(): return f'w:id="{nid()}" w:author="{AUTH}" w:date="{DATE}"'
FONT = '<w:rFonts w:ascii="Arial" w:eastAsiaTheme="majorEastAsia" w:hAnsi="Arial" w:cs="Arial"/>'
def run(text, italic=False, bold=False, sub=False, sup=False):
    rp = FONT + ('<w:b/><w:bCs/>' if bold else '') + ('<w:i/><w:iCs/>' if italic else '') + \
         ('<w:vertAlign w:val="subscript"/>' if sub else '') + ('<w:vertAlign w:val="superscript"/>' if sup else '')
    return f'<w:r><w:rPr>{rp}</w:rPr><w:t xml:space="preserve">{escape(text)}</w:t></w:r>'
def parse_markup(s):
    """Light markup: *italic*, _{sub}, ^{sup}. Returns runs."""
    out = []; i = 0; buf = ''; ital = False
    while i < len(s):
        if s[i] == '*':
            if buf: out.append(run(buf, ital)); buf = ''
            ital = not ital; i += 1
        elif s.startswith('_{', i) or s.startswith('^{', i):
            if buf: out.append(run(buf, ital)); buf = ''
            j = s.index('}', i); out.append(run(s[i+2:j], ital, sub=s[i] == '_', sup=s[i] == '^')); i = j + 1
        else:
            buf += s[i]; i += 1
    if buf: out.append(run(buf, ital))
    return ''.join(out)
def ins_para(text, heading=False):
    rp = FONT + ('<w:b/><w:bCs/><w:i/><w:iCs/>' if heading else '')
    ppr = f'<w:pPr><w:spacing w:before="240" w:after="80"/><w:jc w:val="both"/><w:rPr><w:ins {attrs()}/>{rp}</w:rPr></w:pPr>'
    body = run(text, italic=True, bold=True) if heading else parse_markup(text)
    return f'<w:p>{ppr}<w:ins {attrs()}>{body}</w:ins></w:p>'
def del_para(p):
    """Track-delete a whole paragraph: every run, every math zone, and the paragraph mark."""
    def del_run(m):
        r = m.group(0).replace('<w:t>', '<w:delText>').replace('<w:t ', '<w:delText ').replace('</w:t>', '</w:delText>')
        return f'<w:del {attrs()}>{r}</w:del>'
    # math zones: wrap inner content
    parts = re.split(r'(<m:oMath>.*?</m:oMath>)', p, flags=re.S); out = []
    for part in parts:
        if part.startswith('<m:oMath>'):
            out.append(f'<m:oMath><w:del {attrs()}>{part[len("<m:oMath>"):-len("</m:oMath>")]}</w:del></m:oMath>')
        else:
            out.append(re.sub(r'<w:r>(?:(?!</w:r>).)*</w:r>|<w:r [^>]*>(?:(?!</w:r>).)*</w:r>', del_run, part, flags=re.S))
    p = ''.join(out)
    # paragraph mark
    if '<w:pPr>' in p and re.search(r'<w:pPr>.*?</w:pPr>', p, re.S):
        ppr = re.search(r'<w:pPr>.*?</w:pPr>', p, re.S).group(0)
        if '<w:rPr>' in ppr: nppr = ppr.replace('<w:rPr>', f'<w:rPr><w:del {attrs()}/>', 1)
        else: nppr = ppr.replace('</w:pPr>', f'<w:rPr><w:del {attrs()}/></w:rPr></w:pPr>')
        p = p.replace(ppr, nppr, 1)
    return p
x = open(DOC, encoding='utf8').read()
b0 = x.index('<w:body>')
paras = [(m.start() + b0, m.end() + b0) for m in re.finditer(r'<w:p[ >].*?</w:p>', x[b0:], re.S)]
def ptext(i): return ''.join(re.findall(r'<w:t[^>]*>([^<]*)</w:t>', x[paras[i][0]:paras[i][1]]))
# anchors by text, so a re-run on a slightly edited V6 fails loudly instead of editing the wrong place
def find(prefix):
    hits = [i for i in range(len(paras)) if ptext(i).startswith(prefix)]
    assert len(hits) == 1, (prefix, hits); return hits[0]
ops = []   # (para index, 'before'|'after'|'delete', xml)
for e in TXT['edits']:
    i = find(e['anchor'])
    if e['op'] == 'delete': ops.append((i, 'delete', None))
    else:
        xml = ''.join(ins_para(t['text'], t.get('heading', False)) for t in e['paras'])
        ops.append((i, e['op'], xml))
for i, op, xml in sorted(ops, key=lambda o: (o[0], {'before': 0, 'delete': 1, 'after': 2}[o[1]]), reverse=True):
    s, e = paras[i]
    if op == 'delete': x = x[:s] + del_para(x[s:e]) + x[e:]
    elif op == 'after': x = x[:e] + xml + x[e:]
    else: x = x[:s] + xml + x[s:]
open(DOC, 'w', encoding='utf8').write(x)
# comments: anchored on the whole heading/paragraph text of the target paragraph
for c in TXT['comments']:
    r = subprocess.run(['python3', os.path.join(SK, 'comment.py'), U, c['text'], '--author', AUTH, '--initials', 'C'],
                       capture_output=True, text=True)
    if r.returncode:
        r = subprocess.run(['python3', os.path.join(SK, 'comment.py'), U, c['text']], capture_output=True, text=True)
    cid = re.search(r'w:id="(\d+)"', r.stdout).group(1)
    x = open(DOC, encoding='utf8').read(); b0 = x.index('<w:body>')
    paras = [(m.start() + b0, m.end() + b0) for m in re.finditer(r'<w:p[ >].*?</w:p>', x[b0:], re.S)]
    i = find(c['anchor']); s, e = paras[i]; p = x[s:e]
    pe = p.index('</w:pPr>') + len('</w:pPr>') if '</w:pPr>' in p else p.index('>') + 1
    p = (p[:pe] + f'<w:commentRangeStart w:id="{cid}"/>' + p[pe:-len('</w:p>')] +
         f'<w:commentRangeEnd w:id="{cid}"/><w:r><w:rPr><w:rStyle w:val="CommentReference"/></w:rPr><w:commentReference w:id="{cid}"/></w:r></w:p>')
    x = x[:s] + p + x[e:]; open(DOC, 'w', encoding='utf8').write(x)
OUT = os.path.join(HERE, 'MaterialsMethod_V7_Claude_tracked.docx')
if os.path.exists(OUT): os.remove(OUT)
subprocess.run(f'cd {U} && zip -qXr {OUT} .', shell=True, check=True)
print('written', OUT, 'revision ids used', rid[0] - 6000)
