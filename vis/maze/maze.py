from __future__ import print_function
import re
from subprocess import check_output

m_pat = re.compile(r'\s*(\d+)\s+(\d+)\s+(\d+)')

def mummer_matches(fn1, fn2, l, m='mumreference', verbose=False):
    cmd = 'mummer -{} -l {} -b -c {} {} 2>/dev/null'.format(m, l, fn1, fn2)
    if verbose:
        print('#', cmd)
    mummer_out = check_output(cmd, shell=True)
    fwd = []
    rev = []
    for l in mummer_out.splitlines():
        if l.startswith('>'):
            is_rev = True if 'Reverse' in l else False
        else:
            s1, s2, l = [int(x) for x in m_pat.match(l).group(1,2,3)]
            if is_rev:
                rev.append([s1, s1+l-1, s2, s2-l+1])
            else:
                fwd.append([s1, s1+l-1, s2, s2+l-1])
    return {'fwd': fwd, 'rev': rev}
