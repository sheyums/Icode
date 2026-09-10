#!/usr/bin/env python3
"""Generate FitHyperexponentialMLE_oct.m, an Octave-runnable twin of
FitHyperexponentialMLE.m.

Octave (through 8.x) cannot parse a MATLAB "arguments" validation block.
This script replaces that block -- and only that block -- with an
equivalent defaults + name-value parsing shim, leaving the body byte for
byte identical, so test_FitHyperexponentialMLE.m can be run without a
MATLAB licence:

    python3 tools/make_octave_twin.py
    octave-cli test_FitHyperexponentialMLE.m

The twin is a build artifact. Regenerate it after every edit to
FitHyperexponentialMLE.m; never edit it directly.
"""
import re, sys, os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
src=open(os.path.join(ROOT, 'FitHyperexponentialMLE.m')).read()
m=re.search(r'^arguments\n(.*?)^end\n', src, re.S|re.M)
body=m.group(1)
lines=[l.strip() for l in body.strip().split('\n') if l.strip().startswith('options.')]
defs=[]
for l in lines:
    name=re.match(r'options\.(\w+)',l).group(1)
    default=l.rsplit('=',1)[1].strip() if '=' in l.rsplit('options.',1)[-1] else '[]'
    defs.append(f"options.{name} = {default};")
shim="% AUTO-GENERATED Octave twin: 'arguments' block replaced by varargin parsing.\n"
shim+="\n".join(defs)+"\n"
shim+="""for ii = 1:2:numel(varargin)
    options.(varargin{ii}) = varargin{ii+1};
end
validateattributes(eventseries,{'numeric'},{'real'});
validateattributes(xmin,{'numeric'},{'scalar','positive'});
"""
out=src[:m.start()]+shim+src[m.end():]
out=out.replace('function H = FitHyperexponentialMLE(eventseries, xmin, options)',
                'function H = FitHyperexponentialMLE_oct(eventseries, xmin, varargin)',1)
open(os.path.join(ROOT, 'FitHyperexponentialMLE_oct.m'), 'w').write(out)
print("generated; options parsed:", len(defs))
