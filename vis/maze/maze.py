from __future__ import print_function
import re, os, glob, sys
from tempfile import NamedTemporaryFile
from subprocess import check_output, call, Popen, PIPE


################################################################################
#
# CONFIG
#
# Specify the paths to external tools, if not in PATH
#
tool = {}
tool['mummer']     = 'mummer'
tool['lastdb']     = 'lastdb' # last >= 584 is needed. Todo(meiers): write requirements
tool['lastal']     = 'lastal'
tool['last-split'] = 'last-split'
#
# END CONFIG
#
################################################################################


m_pat = re.compile(r'\s*(\d+)\s+(\d+)\s+(\d+)')

def mummer_matches(fn1, fn2, l, m='mumreference', verbose=False):
    cmd = '{} -{} -l {} -b -c {} {} 2>/dev/null'.format(tool['mummer'], m, l, fn1, fn2)
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


################################################################################
# maze detail 

def LASTsplit_matches(ref, seq):
    '''Given two DNA sequences, compute last-split alignment
    and return matches'''
    matches = []
    with NamedTemporaryFile(delete=False, suffix='.fa') as f_ref, \
         NamedTemporaryFile(delete=False, suffix='.fa') as f_read:
        fn_ref = f_ref.name
        fn_read = f_read.name
        print('>tmp_ref', file=f_ref)
        print(ref, file=f_ref)
        print('>tmp_read', file=f_read)
        print(seq, file=f_read)
    try: 
        cmd = '{} {} {} 2>/dev/null'.format(tool['lastdb'], fn_ref, fn_ref)
        call(cmd, shell=True)
        p1 = Popen([tool['lastal'], '-a', '3', fn_ref, fn_read], stdout=PIPE)
        p2 = Popen(tool['last-split'], stdin=p1.stdout, stdout=PIPE)
        last_split = p2.communicate()[0]
    except Exception as e:
        print ("FATAL: LAST crashed.", e, sep='\n', file=sys.stderr)
        last_split=''
    finally: 
        #for x in glob.glob(str(fn_ref + '*')):
        #    os.remove(x)
        #os.remove(fn_read)
        pass
    matches = []
    it = iter(last_split.splitlines())
    while True:
        try:
            ln = it.next()
        except StopIteration, e:
            break
        if ln.startswith("s"):
            ln2 = it.next()
            match=dict()
            # reference entry
            s = ln.split()
            match['d1'] = int(s[2])
            match['d2'] = int(s[2])+int(s[3])
            match['ali1'] = s[6]
            # query entry
            s = ln2.split()
            match['q1'] = int(s[2])
            match['q2'] = int(s[2])+int(s[3])
            match['ali2'] = s[6]
            match['strand'] = s[4]
            match['qlen'] = int(s[5])
            # mismatches
            match['sim'] = len(match['ali1']) - sum([a==b for (a,b) in 
                            zip(match['ali1'].upper(), match['ali2'].upper())])
            match['record'] = _format_last_alignment(match)
            matches.append(match)

    return matches

# Todo(meiers): Do we really need this?
def _format_last_alignment(m, w=60):
    a1 = m['ali1']
    a2 = m['ali2']
    am = ''.join([ (':' if c.upper()==C.upper() else ' ') for c,C in zip(a1,a2) ])
    l1 = []
    l2 = []
    lm = []
    while a1:
        l1.append(a1[:w])
        l2.append(a2[:w])
        lm.append(am[:w])
        a1 = a1[w:]
        a2 = a2[w:]
        am = am[w:]
    final = ''        
    for i in range(len(l1)):
        final += 'ref'.format(m['strand']).ljust(8) + l1[i] + '\n'
        final += ''.ljust(8) + lm[i] + '\n'
        final += 'read{}'.format(m['strand']).ljust(8) + l2[i] + '\n'
        final += '\n'
    return (final)

def _transform_coords(match):
    '''Transform coordinates from LAST type (0-based, minus strand are 
        coordinates in the reverse complement sequence) to traditional
        coordinates (1-based, end is last based covered)'''
    m = match.copy()
    m['d1'] = match['d1'] + 1
    if m['strand'] == '+':
        m['q1'] = match['q1'] + 1
    else:
        m['q1'] = match['qlen'] - match['q2'] + 1
        m['q2'] = match['qlen'] - match['q1']
    return m